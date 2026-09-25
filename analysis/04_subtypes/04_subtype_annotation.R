#!/usr/bin/env Rscript
################################################################################
# 04_subtypes_consolidated_publication.R  [Seurat v5 / publication-ready]
#
# Pipeline:
#   PART 1: 7-celltype scoring with 5 methods (ALL on same cell set; no resample bias)
#   PART 3: Cerebrovascular-cell subtypes (PC/SMC/Art/Cap/Ven) - if vascular cells selected
#   PART 4: Astrocyte subtypes (Homeostatic/Intermediate/Reactive)
#   PART 5: Microglia subtypes (Homeostatic/DAM/Inflammatory)
#   PART 6: Bulk deconvolution back-validation (ssGSEA + stats + survival)
#
# Key fixes vs initial:
#   - All 5 methods run on SAME sub-sampled cells (true comparability)
#   - sc variable name conflict eliminated (renamed to scoring score)
#   - AverageExpression column-name (gN vs N) handled robustly
#   - Bulk deconvolution: ssGSEA + Wilcoxon + FDR + boxplot with p-values
#   - Per-figure source_data CSV (publication output)
#   - Unified publication theme + npg/glasbey palette
#   - PDF + PNG (300dpi) + TIFF (600dpi)
################################################################################

set.seed(123456)
options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 200 * 1024^3)
paths_to_add <- "/path/to/R_libraries"
for (p in paths_to_add) {
  if (dir.exists(p) && !(p %in% .libPaths())) .libPaths(c(p, .libPaths()))
}

suppressPackageStartupMessages({
  library(Seurat); library(qs)
  library(ggplot2); library(patchwork); library(ggpubr)
  library(dplyr); library(tidyr); library(tibble)
  library(AUCell); library(GSVA); library(GSEABase)
  library(harmony); library(pheatmap); library(lme4)
  library(RColorBrewer); library(scales); library(viridis)
  library(singscore); library(speckle)
})

has_ucell    <- requireNamespace("UCell",     quietly = TRUE)
has_singsc   <- requireNamespace("singscore", quietly = TRUE)
has_ggsci    <- requireNamespace("ggsci",     quietly = TRUE)
has_pals     <- requireNamespace("pals",      quietly = TRUE)
has_ggrastr  <- requireNamespace("ggrastr",   quietly = TRUE)
has_ComplexH <- requireNamespace("ComplexHeatmap", quietly = TRUE)

# ============================================================================
# USER PARAMETERS
# ============================================================================
TEST_MODE <- TRUE
if (TEST_MODE) {
  IN_QS    <- "/path/to/project/results/02_annotation/rds/scRNA_annotated.qs"
  BULK_TAB <- "/path/to/project/results/03_bulk/tables"
  BULK_DAT <- "/path/to/public_data/bulk"
  OUT      <- "/path/to/project/results/04_subtypes"
} else {
  IN_QS    <- "/path/to/project/results/02_annotation/rds/scRNA_annotated.qs"
  BULK_TAB <- "/path/to/project/results/03_bulk/tables"
  BULK_DAT <- "/path/to/public_data/bulk"
  OUT      <- "/path/to/project/results/04_subtypes"
}

# Sub-sampling for heavy methods (ssGSEA/AUCell). Set to NA to disable.

OUT_FIG <- file.path(OUT, "figures_main")
OUT_ED  <- file.path(OUT, "figures_extended_data")
OUT_TAB <- file.path(OUT, "tables")
OUT_RDS <- file.path(OUT, "rds")
OUT_LOG <- file.path(OUT, "logs")
OUT_SRC <- file.path(OUT, "source_data")
for (d in c(OUT_FIG, OUT_ED, OUT_TAB, OUT_RDS, OUT_LOG, OUT_SRC))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

LOG_FILE <- file.path(OUT_LOG, "node04_run.log")
log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", ts, msg)
  cat(line, "\n"); cat(line, "\n", file = LOG_FILE, append = TRUE)
}

# ============================================================================
# publication THEME / PALETTE / HELPERS
# ============================================================================
publication_theme <- function(base = 7) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      axis.text   = element_text(color = "black", size = base),
      axis.title  = element_text(color = "black", size = base + 0.5),
      axis.line   = element_line(color = "black", linewidth = 0.35),
      axis.ticks  = element_line(color = "black", linewidth = 0.35),
      legend.text = element_text(size = base - 0.5),
      legend.title= element_text(size = base),
      legend.key.size = unit(0.32, "cm"),
      strip.background = element_blank(),
      strip.text  = element_text(size = base, face = "bold"),
      plot.title  = element_text(size = base + 1, face = "bold", hjust = 0.5),
      panel.grid  = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
}

publication_palette <- function(n) {
  if (n <= 10 && has_ggsci) ggsci::pal_npg("nrc")(n)
  else if (has_pals) as.character(pals::glasbey(n))
  else colorRampPalette(brewer.pal(8, "Set1"))(n)
}

# 7-celltype palette (consistent with node02)
ct_pal <- c(
  "Oligodendrocytes" = "#0072B2", "Excitatory" = "#999999",
  "Inhibitory" = "#56B4E9",       "Astrocytes" = "#CC79A7",
  "Microglia"  = "#D55E00",       "OPCs" = "#009E73",
  "Cerebrovascular cells" = "#F0E442",
  "Endothelial"= "#F0E442"
)
group_palette <- c("CN" = "#4DBBD5", "AD" = "#E64B35",
                   "NC" = "#4DBBD5", "UC" = "#E64B35",
                   "Ctrl" = "#4DBBD5",
                   "Control" = "#4DBBD5", "Disease" = "#E64B35",
                   "ND" = "#4DBBD5")

standardize_group <- function(x) {
  v <- tolower(trimws(as.character(x)))
  v[v %in% c("cn", "nc", "nd", "ctrl", "control", "normal", "healthy")] <- "CN"
  v[v %in% c("ad", "uc", "disease", "case", "alz", "alzheimers")] <- "AD"
  v
}

save_fig <- function(p, name, dir, w = 6, h = 4) {
  pdf_f  <- file.path(dir, paste0(name, ".pdf"))
  png_f  <- file.path(dir, paste0(name, ".png"))
  svg_f  <- file.path(dir, paste0(name, ".svg"))
  tiff_f <- file.path(dir, paste0(name, ".tiff"))
  ggsave(pdf_f, p, width = w, height = h, units = "in", device = cairo_pdf)
  ggsave(png_f, p, width = w, height = h, units = "in", dpi = 300)
  tryCatch(ggsave(svg_f, p, width = w, height = h, units = "in"),
           error = function(e) log_msg(sprintf("  SVG skipped (%s)", name)))
  tryCatch(ggsave(tiff_f, p, width = w, height = h, units = "in", dpi = 600,
                  compression = "lzw", device = "tiff"),
           error = function(e) log_msg(sprintf("  TIFF skipped (%s)", name)))
  log_msg(sprintf("  Saved: %s", name))
}
safe_write_csv <- function(x, f)
  tryCatch(write.csv(x, f, row.names = FALSE),
           error = function(e) log_msg(sprintf("  CSV save failed: %s", f)))
save_source <- function(df, fig_name)
  safe_write_csv(df, file.path(OUT_SRC, paste0(fig_name, "_source_data.csv")))

# Robust column lookup for AverageExpression (handles "g0" vs "0")
ave_col <- function(mat, key) {
  key <- as.character(key)
  cn <- colnames(mat)
  # Try exact match
  if (key %in% cn) return(key)
  # Try g-prefix (Seurat v5 prepends "g" to numeric keys)
  gk <- paste0("g", key)
  if (gk %in% cn) return(gk)
  # Try numeric match (if colnames are numbers)
  if (key %in% as.character(as.numeric(cn))) return(key)
  # Debug: print available columns
  cat("  [ave_col] key=", key, " cols=", paste(head(cn, 10)), "...\n")
  return(NA_character_)
}

prepare_v5 <- function(obj, a = "RNA") {
  DefaultAssay(obj) <- a
  if (inherits(obj[[a]], "Assay5")) obj <- JoinLayers(obj, assay = a)
  obj
}

harmony_emb <- function(obj, gbv = "orig.ident", du = 1:15, th = 2) {
  pca_emb <- Embeddings(obj, "pca")[, du]
  harm_emb <- harmony::RunHarmony(pca_emb, obj@meta.data,
                                  vars_use = gbv, theta = th,
                                  verbose = FALSE)
  colnames(harm_emb) <- paste0("harmony_", 1:ncol(harm_emb))
  rownames(harm_emb) <- colnames(obj)
  obj[["harmony"]] <- CreateDimReducObject(harm_emb, key = "harmony_", assay = "RNA")
  obj
}

geom_pt <- function(...) {
  if (has_ggrastr) ggrastr::geom_point_rast(..., raster.dpi = 300)
  else geom_point(...)
}

# ============================================================================
# PART 1: LOAD + 5-METHOD SCORING (all on FULL seu_all)
# Methods 1-4: AddModuleScore / AUCell / UCell / singscore on ALL 184k cells
# Method 5: ssGSEA via pseudo-bulk (celltype x sample), then back-fill to cells
# ============================================================================
log_msg("================ PART 1: 7-celltype scoring (5 methods) ================")
stopifnot(file.exists(IN_QS))
seu_all <- qread(IN_QS); seu_all <- prepare_v5(seu_all)
log_msg(sprintf("Loaded: %d cells x %d features", ncol(seu_all), nrow(seu_all)))

ct_order <- intersect(
  c("Oligodendrocytes","Excitatory","Inhibitory","Astrocytes",
    "Microglia","OPCs","Cerebrovascular cells","Endothelial"),
  unique(as.character(seu_all$celltype))
)
log_msg(sprintf("Cell types found: %s", paste(ct_order, collapse = ", ")))

if ("Group" %in% colnames(seu_all[[]])) {
  seu_all$Group <- standardize_group(seu_all$Group)
} else if ("group" %in% colnames(seu_all[[]])) {
  seu_all$Group <- standardize_group(seu_all$group)
}

# Bulk up-genes
up_file <- file.path(BULK_TAB, "bulk_up_geneset_AD_vs_CN.txt")
if (!file.exists(up_file)) up_file <- file.path(BULK_TAB, "bulk_up_geneset.txt")
stopifnot(file.exists(up_file))
up_genes <- readLines(up_file)
up_genes_use <- intersect(up_genes, rownames(seu_all))
log_msg(sprintf("Bulk up-genes: %d total, %d available in scRNA", length(up_genes), length(up_genes_use)))
if (length(up_genes_use) < 10) stop("Too few bulk up-genes overlap with scRNA features.")

log_msg("  Method 1: AddModuleScore on ALL cells ...")
seu_all <- AddModuleScore(seu_all, features = list(up_genes_use),
                          name = "AD_AMS_", assay = "RNA", seed = 1)
score_ams <- setNames(seu_all$AD_AMS_1, colnames(seu_all))

log_msg("  Method 2: AUCell on ALL cells ...")
expr_all <- GetAssayData(seu_all, assay = "RNA", layer = "data")
ranks_all <- AUCell_buildRankings(expr_all, nCores = 1, plotStats = FALSE, verbose = FALSE)
auc_all <- AUCell_calcAUC(list(AD_up = up_genes_use), ranks_all,
                          aucMaxRank = ceiling(nrow(ranks_all) * 0.05), verbose = FALSE)
score_aucell <- setNames(as.numeric(SummarizedExperiment::assay(auc_all)[1, ]), colnames(seu_all))
rm(ranks_all, auc_all); gc()

# ---- Method 3: UCell (chunked to avoid segfault on large matrix) ----
if (has_ucell) {
  log_msg("  Method 3: UCell on ALL cells (chunked 30k/batch)...")
  suppressPackageStartupMessages(library(UCell))
  chunk_size <- 30000
  cell_chunks <- split(colnames(seu_all), ceiling(seq_along(colnames(seu_all)) / chunk_size))
  score_ucell_list <- list()
  for (i in seq_along(cell_chunks)) {
    cells <- cell_chunks[[i]]
    log_msg(sprintf("    UCell chunk %d/%d (%d cells)...", i, length(cell_chunks), length(cells)))
    seu_chunk <- seu_all[, cells]
    res <- tryCatch({
      seu_chunk <- UCell::AddModuleScore_UCell(seu_chunk, features = list(AD_up = up_genes_use),
                                               name = "", assay = "RNA")
      list(scores = setNames(seu_chunk$AD_up, cells), ok = TRUE)
    }, error = function(e) {
      log_msg(sprintf("    UCell chunk %d FAILED: %s", i, e$message))
      list(scores = setNames(rep(NA_real_, length(cells)), cells), ok = FALSE)
    })
    score_ucell_list[[i]] <- res$scores
    rm(seu_chunk, res); gc()
  }
  score_ucell <- unlist(score_ucell_list)
  score_ucell <- score_ucell[colnames(seu_all)]
  rm(score_ucell_list, cell_chunks); gc()
  n_ok <- sum(!is.na(score_ucell))
  log_msg(sprintf("  UCell done: %d/%d cells scored", n_ok, length(score_ucell)))
} else {
  log_msg("  UCell not installed - skipped.")
  score_ucell <- setNames(rep(NA_real_, ncol(seu_all)), colnames(seu_all))
}

if (has_singsc) {
  log_msg("  singscore will be computed in post-processing step.")
  score_ss <- setNames(rep(NA_real_, ncol(seu_all)), colnames(seu_all))
} else {
  log_msg("  singscore not installed - skipped.")
  score_ss <- setNames(rep(NA_real_, ncol(seu_all)), colnames(seu_all))
}
rm(expr_all); gc()

# ---- Method 5: ssGSEA via pseudo-bulk ----
log_msg("  Method 5: ssGSEA via pseudo-bulk (celltype x sample) ...")
if (!"orig.ident" %in% colnames(seu_all[[]])) seu_all$orig.ident <- "Sample1"
pb_meta <- seu_all[[]] %>% dplyr::select(orig.ident, celltype) %>%
  mutate(pb_id = paste0(celltype, "__", orig.ident))
pb_ids <- unique(pb_meta$pb_id); pb_ids <- pb_ids[!is.na(pb_ids)]
log_msg(sprintf("  Building %d pseudo-bulk profiles ...", length(pb_ids)))
expr_l <- GetAssayData(seu_all, assay = "RNA", layer = "data")
pb_mat <- do.call(cbind, lapply(pb_ids, function(pid) {
  cells_in <- colnames(seu_all)[pb_meta$pb_id == pid]
  if (length(cells_in) == 1) return(as.matrix(expr_l[, cells_in, drop = FALSE]))
  Matrix::rowMeans(expr_l[, cells_in, drop = FALSE])
}))
colnames(pb_mat) <- pb_ids; pb_mat <- as.matrix(pb_mat)
gp_pb <- ssgseaParam(pb_mat, list(AD_up = up_genes_use), minSize = 10, maxSize = 500, normalize = TRUE)
ssgsea_pb <- gsva(gp_pb, verbose = FALSE)
score_gsva <- setNames(as.numeric(ssgsea_pb[1, pb_meta$pb_id]), colnames(seu_all))
rm(pb_mat, ssgsea_pb, gp_pb); gc()

score_methods <- c("AddModuleScore", "AUCell", "UCell", "singscore", "ssGSEA")
score_df <- data.frame(
  cell = colnames(seu_all),
  celltype = as.character(seu_all$celltype),
  AddModuleScore = score_ams[colnames(seu_all)],
  AUCell = score_aucell[colnames(seu_all)],
  UCell = score_ucell[colnames(seu_all)],
  singscore = score_ss[colnames(seu_all)],
  ssGSEA = score_gsva[colnames(seu_all)],
  stringsAsFactors = FALSE
) %>% dplyr::filter(!is.na(celltype) & celltype %in% ct_order)
score_df$celltype <- factor(score_df$celltype, levels = ct_order)

score_df_norm <- score_df
for (m in score_methods) {
  x <- score_df_norm[[m]]
  if (all(is.na(x))) next
  rng <- range(x, na.rm = TRUE)
  score_df_norm[[m]] <- if (diff(rng) == 0) 0 else (x - rng[1]) / diff(rng)
}
score_df_norm$Combined <- rowMeans(score_df_norm[, score_methods], na.rm = TRUE)
score_df$Combined <- score_df_norm$Combined

meta_add <- score_df_norm[, c("cell", score_methods, "Combined")]
rownames(meta_add) <- meta_add$cell; meta_add$cell <- NULL
seu_all <- AddMetaData(seu_all, meta_add)

# ---- Compute singscore in PART 1 (chunked, in-memory, not post-hoc) ----
log_msg("  Computing singscore (25000-cell chunks)...")
sc_u <- intersect(up_genes_use, rownames(seu_all))
sc_chunks <- split(colnames(seu_all), ceiling(seq_along(colnames(seu_all)) / 25000))
sc_list <- list()
for (i in seq_along(sc_chunks)) {
  e <- as.matrix(GetAssayData(seu_all[, sc_chunks[[i]]], assay = "RNA", layer = "data"))
  r <- rankGenes(e)
  ss <- simpleScore(r, upSet = sc_u, knownDirection = TRUE)
  sc_list[[i]] <- setNames(ss$TotalScore, sc_chunks[[i]])
  log_msg(sprintf("    singscore chunk %d/%d done", i, length(sc_chunks)))
  rm(e, r, ss); gc()
}
sc_all <- unlist(sc_list)
log_msg(sprintf("  singscore: %d cells, NA=%d", length(sc_all), sum(is.na(sc_all))))

# Update seu_all + score_df_norm with real singscore, re-normalize
seu_all$singscore <- sc_all[colnames(seu_all)]
score_df_norm$singscore <- sc_all[score_df_norm$cell]
for (m in score_methods) {
  x <- score_df_norm[[m]]; if (all(is.na(x))) next
  rng <- range(x, na.rm = TRUE)
  score_df_norm[[m]] <- if (diff(rng) == 0) 0 else (x - rng[1]) / diff(rng)
}
score_df_norm$Combined <- rowMeans(score_df_norm[, score_methods], na.rm = TRUE)
seu_all$Combined <- score_df_norm$Combined
rm(sc_all, sc_list, sc_chunks, sc_u); gc()

qs::qsave(seu_all, file.path(OUT_RDS, "scoring_intermediate_all.qs"), preset = "fast")
safe_write_csv(score_df_norm, file.path(OUT_TAB, "scoring_per_cell.csv"))
log_msg(sprintf("PART 1 done (singscore NA=%d)", sum(is.na(score_df_norm$singscore))))

singscore_csv <- file.path(OUT_TAB, "scoring_per_cell.csv")
if (file.exists(singscore_csv)) {
  score_df_norm <- read.csv(singscore_csv)
  score_df_norm$celltype <- factor(score_df_norm$celltype, levels = intersect(ct_order, unique(score_df_norm$celltype)))
  if ("Combined" %in% colnames(score_df_norm)) {
    seu_all$Combined <- score_df_norm$Combined[match(colnames(seu_all), score_df_norm$cell)]
  }
  log_msg(sprintf("Score reloaded: %d cells, singscore NA=%d", nrow(score_df_norm), sum(is.na(score_df_norm$singscore))))
}

log_msg("================ PART 2: Statistical ranking ================")

# Per-method per-celltype mean
score_summary <- score_df_norm %>%
  group_by(celltype) %>%
  summarise(across(all_of(c(score_methods, "Combined")), ~mean(.x, na.rm = TRUE)),
            n = dplyr::n(), .groups = "drop")
safe_write_csv(score_summary, file.path(OUT_TAB, "scoring_summary_by_celltype.csv"))

# Inter-method consensus: rank each method, then average rank per celltype
ranks_mat <- sapply(score_methods, function(m) {
  v <- score_summary[[m]]
  if (all(is.na(v))) return(rep(NA, length(v)))
  rank(-v, ties.method = "average")
})
rownames(ranks_mat) <- score_summary$celltype
mean_rank <- rowMeans(ranks_mat, na.rm = TRUE)
rank_df <- data.frame(
  celltype = score_summary$celltype,
  mean_rank = mean_rank,
  combined_score = score_summary$Combined,
  stringsAsFactors = FALSE
) %>% arrange(mean_rank)
safe_write_csv(rank_df, file.path(OUT_TAB, "celltype_consensus_rank.csv"))

# Statistical test: each celltype vs all others (Wilcoxon on Combined)
wilcox_df <- do.call(rbind, lapply(ct_order, function(ct) {
  in_grp  <- score_df_norm$Combined[score_df_norm$celltype == ct]
  out_grp <- score_df_norm$Combined[score_df_norm$celltype != ct]
  if (length(in_grp) < 5 || length(out_grp) < 5)
    return(data.frame(celltype = ct, median_in = NA, median_out = NA,
                      delta = NA, p = NA, log10p = NA))
  w <- wilcox.test(in_grp, out_grp, alternative = "greater")
  data.frame(
    celltype   = ct,
    median_in  = median(in_grp, na.rm = TRUE),
    median_out = median(out_grp, na.rm = TRUE),
    delta      = median(in_grp, na.rm = TRUE) - median(out_grp, na.rm = TRUE),
    p          = w$p.value,
    log10p     = -log10(pmax(w$p.value, 1e-300))
  )
}))
wilcox_df$padj <- p.adjust(wilcox_df$p, method = "BH")
wilcox_df <- wilcox_df %>% arrange(desc(delta))
safe_write_csv(wilcox_df, file.path(OUT_TAB, "scoring_wilcox_per_celltype.csv"))

# Selection rule: padj < 0.05 AND delta > 0 AND in top by mean_rank
sig_ct <- wilcox_df %>% filter(padj < 0.05 & delta > 0) %>% pull(celltype)
top_ct <- rank_df %>% filter(mean_rank <= median(mean_rank)) %>% pull(celltype)
selected_for_subtype <- intersect(sig_ct, top_ct)
# Always-include from biology priors (NVU paper)
nvu_priors <- intersect(c("Cerebrovascular cells","Endothelial","Astrocytes","Microglia"), ct_order)
selected_for_subtype <- union(selected_for_subtype, nvu_priors)
log_msg(sprintf("Cell types selected for subtype analysis: %s",
                paste(selected_for_subtype, collapse = ", ")))

selection_log <- data.frame(
  celltype = ct_order,
  mean_rank = rank_df$mean_rank[match(ct_order, rank_df$celltype)],
  combined  = rank_df$combined_score[match(ct_order, rank_df$celltype)],
  delta     = wilcox_df$delta[match(ct_order, wilcox_df$celltype)],
  padj      = wilcox_df$padj[match(ct_order, wilcox_df$celltype)],
  selected  = ct_order %in% selected_for_subtype
)
safe_write_csv(selection_log, file.path(OUT_TAB, "celltype_selection_log.csv"))

# ----- Figures for scoring -----
ct_cols_used <- ct_pal[intersect(names(ct_pal), ct_order)]

# Fig 1a: Combined score violin + box + stats
p_vln <- ggplot(score_df_norm, aes(celltype, Combined, fill = celltype)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.85, linewidth = 0.3) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white",
               alpha = 0.9, linewidth = 0.3) +
  scale_fill_manual(values = ct_cols_used) +
  labs(x = NULL, y = "AD signature score (Combined)",
       title = "AD signature across cell types") +
  publication_theme(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
save_fig(p_vln, "Fig1a_AD_signature_violin", OUT_FIG, w = 5, h = 4)
save_source(score_df_norm[, c("cell","celltype","Combined")], "Fig1a_AD_signature_violin")

# Fig 1b: Method heatmap
sh <- score_summary %>%
  dplyr::select(celltype, all_of(c(score_methods, "Combined"))) %>%
  pivot_longer(-celltype, names_to = "method", values_to = "score")
sh$method <- factor(sh$method, levels = rev(c(score_methods, "Combined")))

p_hm <- ggplot(sh, aes(celltype, method, fill = score)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", score)), size = 2.2) +
  scale_fill_gradientn(colors = c("#3B4CC0", "white", "#B40426"),
                       name = "Score") +
  labs(x = NULL, y = NULL, title = "Five-method scoring matrix") +
  publication_theme(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_hm, "Fig1b_method_heatmap", OUT_FIG, w = 5.5, h = 3.5)
save_source(sh, "Fig1b_method_heatmap")

# Fig 1c: Wilcoxon delta + significance
wilcox_df$celltype <- factor(wilcox_df$celltype, levels = wilcox_df$celltype[order(wilcox_df$delta)])
p_w <- ggplot(wilcox_df, aes(celltype, delta, fill = celltype)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = ifelse(padj < 0.001, "***",
                               ifelse(padj < 0.01, "**",
                                      ifelse(padj < 0.05, "*", "ns")))),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = ct_cols_used) +
  labs(x = NULL, y = "Delta score (vs others)",
       title = "Statistical separation per cell type") +
  publication_theme(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
save_fig(p_w, "Fig1c_wilcoxon_delta", OUT_FIG, w = 5, h = 4)
save_source(wilcox_df, "Fig1c_wilcoxon_delta")

# Fig 1 panel
fig1 <- (p_vln | p_hm | p_w) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9))
save_fig(fig1, "Fig1_AD_scoring_overview", OUT_FIG, w = 14, h = 4.2)

# UMAP colored by Combined score
if ("umap" %in% names(seu_all@reductions)) {
  tryCatch({
    umap_df <- as.data.frame(Embeddings(seu_all, "umap"))
    colnames(umap_df) <- c("UMAP_1", "UMAP_2")
    umap_df$Combined <- seu_all$Combined
    umap_df$celltype <- as.character(seu_all$celltype)
    p_um <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = Combined)) +
      geom_point(size = 0.3, alpha = 0.85) +
      scale_color_viridis_c(option = "magma", name = "AD score") +
      labs(title = "AD signature on UMAP") +
      publication_theme(7) + theme(aspect.ratio = 1)
    save_fig(p_um, "Fig1d_AD_score_UMAP", OUT_FIG, w = 5, h = 4.5)
    save_source(umap_df, "Fig1d_AD_score_UMAP")
  }, error = function(e) log_msg(sprintf("Fig1d skipped: %s", e$message)))
}
gc()

# ED Fig 1: Each method's per-cell distribution
ed1_data <- score_df_norm %>%
  pivot_longer(all_of(score_methods), names_to = "method", values_to = "value")
p_ed1 <- ggplot(ed1_data, aes(celltype, value, fill = celltype)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.85, linewidth = 0.25) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
               alpha = 0.9, linewidth = 0.25) +
  scale_fill_manual(values = ct_cols_used) +
  facet_wrap(~ method, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Score (min-max scaled)",
       title = "Per-method scoring distributions") +
  publication_theme(6.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
        legend.position = "none")
save_fig(p_ed1, "ED_Fig1_per_method_violin", OUT_ED, w = 9, h = 6)
save_source(ed1_data, "ED_Fig1_per_method_violin")

# ED Fig 2: Method correlation
cor_mat <- cor(score_df_norm[, score_methods], use = "pairwise.complete.obs", method = "spearman")
cor_long <- cor_mat %>% as.data.frame() %>% rownames_to_column("m1") %>%
  pivot_longer(-m1, names_to = "m2", values_to = "rho")
p_ed2 <- ggplot(cor_long, aes(m1, m2, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.5) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
  labs(x = NULL, y = NULL, title = "Method-method concordance") +
  publication_theme(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_ed2, "ED_Fig2_method_correlation", OUT_ED, w = 5, h = 4.2)
save_source(cor_long, "ED_Fig2_method_correlation")

log_msg("PART 2 done.")

# Save checkpoint for worker subprocesses
qs::qsave(seu_all, file.path(OUT_RDS, "seu_all_checkpoint.qs"), preset = "fast")
rm(seu_all); gc()
log_msg("seu_all checkpoint saved. Memory freed.")
log_msg("DEBUG: entering subprocess dispatch...")

# ============================================================================
# PART 3-5: SUBTYPE ANALYSIS (per selected celltype)
# ============================================================================
# Common subclustering wrapper
subcluster_one <- function(obj, ct_name, res = 0.5, dims_pca = 1:30, dims_h = 1:15) {
  obj <- prepare_v5(obj)
  obj <- NormalizeData(obj, verbose = FALSE) %>%
    FindVariableFeatures(nfeatures = 2000, verbose = FALSE)
  hvg <- setdiff(VariableFeatures(obj),
                 grep("^MT-|^RPS|^RPL|^MALAT1|^FOS|^JUN|^HSP",
                      rownames(obj), value = TRUE))
  VariableFeatures(obj) <- hvg
  npcs <- min(30, max(5, ncol(obj) - 1))
  obj <- ScaleData(obj, verbose = FALSE) %>% RunPCA(npcs = npcs, verbose = FALSE)
  if ("orig.ident" %in% colnames(obj[[]]) && length(unique(obj$orig.ident)) > 1) {
    dh <- intersect(seq_len(npcs), dims_h)
    obj <- harmony_emb(obj, "orig.ident", dh, 2)
    red <- "harmony"
  } else {
    red <- "pca"
  }
  use_dims <- intersect(seq_len(npcs), 1:10)
  obj <- FindNeighbors(obj, reduction = red, dims = use_dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = res, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = red, dims = use_dims, verbose = FALSE)
  obj
}

# Generic subtype assignment by marker mean (with margin reporting)
assign_subtype <- function(obj, marker_list, cluster_col = NULL) {
  if (is.null(cluster_col)) {
    snn_cols <- grep("_snn_res", colnames(obj[[]]), value = TRUE)
    cluster_col <- tail(snn_cols, 1)
  }
  Idents(obj) <- obj[[cluster_col]][[1]]
  avg <- AverageExpression(obj, assays = "RNA", layer = "data", verbose = FALSE)$RNA
  res <- list()
  for (cl in levels(Idents(obj))) {
    cn <- ave_col(avg, cl)
    if (is.na(cn)) { res[[cl]] <- list(top = "Unassigned", margin = 0, scores = rep(0, length(marker_list))); next }
    sc <- sapply(marker_list, function(g) {
      gx <- intersect(g, rownames(avg))
      if (length(gx) == 0) return(NA_real_)
      mean(avg[gx, cn], na.rm = TRUE)
    })
    sc[is.na(sc)] <- 0
    o <- order(sc, decreasing = TRUE)
    res[[cl]] <- list(top = names(sc)[o[1]],
                      margin = sc[o[1]] - sc[o[2]],
                      scores = sc)
  }
  res
}

run_subtype_analysis <- function(seu_in, ct_name, marker_list, sub_col, sub_pal,
                                 part_label) {
  log_msg(sprintf("================ %s: %s subtypes ================", part_label, ct_name))
  obj <- subset(seu_in, celltype == ct_name)
  if (ncol(obj) < 30) {
    log_msg(sprintf("  Too few cells (%d) - skip", ncol(obj))); return(NULL)
  }
  obj <- subcluster_one(obj, ct_name)
  
  # Filter markers to genes present
  marker_present <- lapply(marker_list, function(g) intersect(g, rownames(obj)))
  marker_present <- marker_present[sapply(marker_present, length) >= 1]
  for (nm in names(marker_present))
    log_msg(sprintf("  %s markers available: %d (%s)", nm,
                    length(marker_present[[nm]]),
                    paste(marker_present[[nm]], collapse = ",")))
  
  ass <- assign_subtype(obj, marker_present)
  
  if (ct_name == "Astrocytes") {
    for (cl in names(ass)) {
      s <- ass[[cl]]$scores
      h <- s["Homeostatic"]; if(is.na(h)||h<=0) h <- 0.01
      r <- s["Reactive"]
      ratio <- r / h
      if (ratio > 0.5) ass[[cl]]$top <- "Reactive"
      else if (ratio > 0.25) ass[[cl]]$top <- "Intermediate"
      else ass[[cl]]$top <- "Homeostatic"
      s2 <- sort(s, decreasing=TRUE)
      ass[[cl]]$margin <- s2[1] - s2[2]
    }
    log_msg("  Using ratio-based assignment for Astrocytes (H/I/R)")
  }
  
  scores_mat <- do.call(rbind, lapply(ass, function(x) x$scores))
  rownames(scores_mat) <- names(ass)
  margin_vec <- sapply(ass, function(x) x$margin)
  top_vec    <- sapply(ass, function(x) x$top)
  
  # Save cluster-subtype log
  cl_log <- data.frame(
    cluster = names(ass),
    assigned = top_vec,
    margin   = margin_vec,
    confidence = cut(margin_vec, breaks = c(-Inf, 0.1, 0.3, Inf),
                     labels = c("Low","Medium","High"))
  )
  cl_log <- cbind(cl_log, scores_mat)
  safe_write_csv(cl_log, file.path(OUT_TAB,
                                   sprintf("subtype_assignment_%s.csv", ct_name)))
  
  # Map to cells
  cluster_col <- grep("_snn_res", colnames(obj[[]]), value = TRUE)
  cluster_col <- tail(cluster_col, 1)
  cl_ids <- as.character(obj[[cluster_col]][[1]])
  obj[[sub_col]] <- unname(top_vec[cl_ids])
  obj[[sub_col]][is.na(obj[[sub_col]])] <- "Unassigned"
  
  used <- intersect(names(sub_pal), unique(obj[[sub_col]][[1]]))
  pal_use <- sub_pal[used]
  if (!"Unassigned" %in% names(pal_use) && "Unassigned" %in% obj[[sub_col]][[1]])
    pal_use <- c(pal_use, "Unassigned" = "grey75")
  
  # ---- UMAP ----
  umap_df <- as.data.frame(Embeddings(obj, "umap"))
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  umap_df$subtype <- obj[[sub_col]][[1]]
  if ("Group" %in% colnames(obj[[]])) umap_df$group <- as.character(obj$Group)
  
  p_um <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = subtype)) +
    geom_point(size = 0.5, alpha = 0.85) +
    scale_color_manual(values = pal_use) +
    guides(color = guide_legend(override.aes = list(size = 2.5))) +
    labs(title = sprintf("%s subtypes", ct_name), x = "UMAP 1", y = "UMAP 2", color = NULL) +
    publication_theme(7) + theme(aspect.ratio = 1)
  save_fig(p_um, sprintf("Fig_%s_UMAP", ct_name), OUT_FIG, w = 5.2, h = 4.5)
  save_source(umap_df, sprintf("Fig_%s_UMAP", ct_name))
  
  present_subtypes <- unique(obj[[sub_col]][[1]])
  marker_for_dot <- marker_present[intersect(names(marker_present), present_subtypes)]
  marker_for_dot <- rev(marker_for_dot)
  if (length(marker_for_dot) > 0 && length(present_subtypes) > 1) {
    Idents(obj) <- obj[[sub_col]][[1]]
    obj <- SetIdent(obj, value = factor(obj[[sub_col]][[1]], levels = intersect(names(pal_use), present_subtypes)))
    p_dot <- DotPlot(obj, features = marker_for_dot, dot.scale = 4) +
      scale_color_gradientn(colors = c("grey92", "#FDAE61", "#D73027", "#67001F")) +
      scale_size(range = c(0.5, 5)) +
      labs(x = NULL, y = NULL, title = sprintf("%s markers", ct_name)) +
      publication_theme(6.5) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
    total_features <- length(unique(unlist(marker_for_dot)))
    save_fig(p_dot, sprintf("Fig_%s_DotPlot", ct_name), OUT_FIG,
             w = max(5, 0.25 * total_features + 2.5), h = 3.5)
    save_source(p_dot$data, sprintf("Fig_%s_DotPlot", ct_name))
    
    # ---- Heatmap (z-score) ----
    mat_avg <- as.matrix(AverageExpression(obj, assays = "RNA", layer = "data",
                                           features = pg, group.by = sub_col,
                                           verbose = FALSE)$RNA)
    if (nrow(mat_avg) >= 2 && ncol(mat_avg) >= 2) {
      ms <- t(scale(t(mat_avg)))
      ms[ms > 2] <- 2; ms[ms < -2] <- -2
      ms[is.na(ms)] <- 0
      
      pal_heat <- pal_use
      names(pal_heat) <- gsub("_", "-", names(pal_heat))
      ann_col <- data.frame(Subtype = colnames(ms), row.names = colnames(ms))
      pheatmap(ms, cluster_rows = FALSE, cluster_cols = FALSE,
               color = colorRampPalette(c("#053061","#2166AC","#F7F7F7","#B2182B","#67001F"))(100),
               annotation_col = ann_col,
               annotation_colors = list(Subtype = pal_heat),
               annotation_names_col = FALSE, border_color = NA,
               cellwidth = 40, cellheight = 14,
               filename = file.path(OUT_FIG, sprintf("Fig_%s_Heatmap.pdf", ct_name)),
               width = 6, height = max(4, 0.22 * length(pg) + 1.5))
      log_msg(sprintf("  Heatmap saved: Fig_%s_Heatmap.pdf", ct_name))
      save_source(as.data.frame(ms) %>% rownames_to_column("gene"),
                  sprintf("Fig_%s_Heatmap", ct_name))
    }
    
    # ---- Composition by group ----
    if ("Group" %in% colnames(obj[[]])) {
      comp <- as.data.frame(obj[[]]) %>%
        dplyr::count(Group, !!sym(sub_col)) %>%
        group_by(Group) %>%
        mutate(prop = n / sum(n)) %>% ungroup()
      colnames(comp)[2] <- "subtype"
      p_comp <- ggplot(comp, aes(Group, prop, fill = subtype)) +
        geom_col(color = "white", linewidth = 0.3) +
        scale_fill_manual(values = pal_use) +
        scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
        labs(x = NULL, y = "Proportion", fill = NULL,
             title = sprintf("%s composition by group", ct_name)) +
        publication_theme(7)
      save_fig(p_comp, sprintf("Fig_%s_Composition", ct_name), OUT_FIG, w = 4, h = 4)
      save_source(comp, sprintf("Fig_%s_Composition", ct_name))
      
      # ---- CN vs AD per subtype (Propeller FDR: pseudobulk + limma eBayes) ----
      comp_by_cell <- as.data.frame(obj[[]])[, c("orig.ident", "Group", sub_col), drop=FALSE]
      colnames(comp_by_cell)[3] <- "subtype"
      comp_by_cell$Group <- factor(standardize_group(comp_by_cell$Group), levels = c("CN", "AD"))
      sample_prop <- comp_by_cell %>%
        dplyr::count(orig.ident, Group, subtype) %>%
        group_by(orig.ident) %>%
        mutate(prop = n / sum(n)) %>% ungroup()
      
      # Propeller: pseudobulk limma test following a recent single-nucleus application
      prop_res <- tryCatch({
        suppressMessages(propeller(clusters = comp_by_cell$subtype,
                                   sample = comp_by_cell$orig.ident,
                                   group = comp_by_cell$Group,
                                   transform = "logit", robust = TRUE, trend = FALSE))
      }, error = function(e) NULL)
      
      # Build annotation: Delta% + Propeller FDR
      subtypes_present <- unique(sample_prop$subtype)
      anno_list <- lapply(subtypes_present, function(st) {
        dat <- sample_prop %>% filter(subtype == st)
        cn_mean <- mean(dat$prop[dat$Group == "CN"]) * 100
        ad_mean <- mean(dat$prop[dat$Group == "AD"]) * 100
        delta <- ad_mean - cn_mean
        if (!is.null(prop_res) && st %in% rownames(prop_res)) {
          fdr <- prop_res[st, "FDR"]
          label <- sprintf("%+.1f%%\nFDR=%.2f", delta, fdr)
        } else {
          label <- sprintf("%+.1f%%", delta)
        }
        data.frame(subtype = st, mean_CN = cn_mean, mean_AD = ad_mean,
                   delta = delta, label = label, stringsAsFactors = FALSE)
      })
      sig_df <- do.call(rbind, anno_list)
      
      p_sig <- ggplot(sample_prop, aes(Group, prop, fill = Group)) +
        geom_boxplot(width = 0.5, alpha = 0.7, outlier.size = 0.5) +
        geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
        geom_text(data = sig_df, aes(label = label, x = 1.5, y = Inf, vjust = 2),
                  inherit.aes = FALSE, size = 2.8, hjust = 0.5) +
        facet_wrap(~subtype, nrow = 1, scales = "free_y") +
        scale_fill_manual(values = group_palette) +
        labs(x = NULL, y = "Proportion (per sample)",
             title = paste0(ct_name, ": AD vs CN (Propeller FDR")) +
        publication_theme(6) +
        theme(strip.text = element_text(size = 5, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 5))
      save_fig(p_sig, sprintf("Fig_%s_ADvsCN", ct_name), OUT_FIG, w = 2*length(subtypes_present)+1, h = 3.5)
    
    # ---- DEGs per subtype ----
    Idents(obj) <- obj[[sub_col]][[1]]
    if (length(unique(Idents(obj))) >= 2) {
      degs <- tryCatch(
        FindAllMarkers(obj, only.pos = TRUE, logfc.threshold = 0.25,
                       min.pct = 0.1, verbose = FALSE),
        error = function(e) NULL
      )
      if (!is.null(degs) && nrow(degs) > 0) {
        safe_write_csv(degs, file.path(OUT_TAB,
                                       sprintf("subtype_DEGs_%s.csv", ct_name)))
      }
    }
  }
  
  qs::qsave(obj, file.path(OUT_RDS, sprintf("scRNA_%s_annotated.qs",
                                            tolower(substr(ct_name,1,5)))),
            preset = "fast")
  log_msg(sprintf("  %s subtype distribution:", ct_name))
  print(table(obj[[sub_col]][[1]]))
  obj
}
} # <- close run_subtype_analysis function body

# Marker definitions
ec_markers <- list(
  Pericyte  = c("PDGFRB","RGS5","NOTCH3"),
  SMC       = c("ACTA2","TAGLN","CNN1"),
  Arterial  = c("GJA5","EFNB2","BMX"),
  Capillary = c("MFSD2A","RGCC","CA4"),
  Venous    = c("NR2F2","ACKR1")
)
ec_pal <- c("Pericyte"="#8B4513","SMC"="#FF7F00","Arterial"="#E41A1C",
            "Capillary"="#4DAF4A","Venous"="#377EB8","Unassigned"="grey75")

astro_markers <- list(
  Homeostatic  = c("ALDH1L1","SLC1A2","GLUL"),
  Intermediate = c("VIM","CLU"),
  Reactive = c("GFAP","CD44")
)
astro_pal <- c("Homeostatic"="#66C2A5","Intermediate"="#FFD92F",
               "Reactive"="#FC8D62","Unassigned"="grey75")

micro_markers <- list(
  Homeostatic  = c("P2RY12","TMEM119","CX3CR1"),
  DAM          = c("APOE","LPL","CST7","SPP1"),
  Inflammatory = c("IL1B","CCL2","TNF")
)
micro_pal <- c("Homeostatic"="#66C2A5","DAM"="#FC8D62",
               "Inflammatory"="#E78AC3","Unassigned"="grey75")

# Run subtype analyses via isolated subprocess (100% memory release on exit)
sub_objs <- list()
SUBTYPE_WORKER <- file.path(Sys.getenv("AD_NVU_REPO_ROOT", "."), "analysis", "04_subtypes", "04a_subtype_worker.R")
RSCRIPT_BIN <- Sys.which("Rscript")
CKPT_QS <- file.path(OUT_RDS, "seu_all_checkpoint.qs")

run_subtype_subprocess <- function(ct_name) {
  log_msg(sprintf(">>> Launching subprocess for %s ...", ct_name))
  safe_ct <- gsub("[^A-Za-z0-9]+", "_", ct_name)
  wlog <- file.path(OUT_LOG, sprintf("worker_%s.log", safe_ct))
  # Write a .bat file to set env vars and launch worker (system2() env is unreliable on Windows)
  bat <- file.path(OUT_LOG, sprintf("launch_%s.bat", safe_ct))
  bat_lines <- c(
    sprintf('set "WORKER_CELLTYPE=%s"', ct_name),
    sprintf('set "WORKER_CKPT_QS=%s"', CKPT_QS),
    sprintf('set "WORKER_OUT=%s"', OUT),
    sprintf('set "WORKER_OUT_FIG=%s"', OUT_FIG),
    sprintf('set "WORKER_OUT_ED=%s"', OUT_ED),
    sprintf('set "WORKER_OUT_TAB=%s"', OUT_TAB),
    sprintf('set "WORKER_OUT_RDS=%s"', OUT_RDS),
    sprintf('set "WORKER_OUT_SRC=%s"', OUT_SRC),
    sprintf('set "WORKER_LOG_FILE=%s"', wlog),
    sprintf('"%s" --vanilla "%s" > "%s" 2>&1', RSCRIPT_BIN, SUBTYPE_WORKER, wlog)
  )
  writeLines(bat_lines, bat)
  exit_code <- shell(shQuote(bat), wait = TRUE)
  unlink(bat)
  log_msg(sprintf("  subprocess exit code: %s", as.character(exit_code)))
  if (is.null(exit_code) || length(exit_code) == 0 || exit_code == 0) {
    rds_path <- file.path(OUT_RDS, sprintf("subtype_%s.qs", ct_name))
    if (file.exists(rds_path)) {
      res <- qs::qread(rds_path)
      log_msg(sprintf("<<< %s subprocess OK", ct_name))
      return(res)
    }
  }
  log_msg(sprintf("!!! %s subprocess FAILED (exit=%s), see %s", ct_name, as.character(exit_code), wlog))
  return(NULL)
}

vascular_ct_name <- if ("Cerebrovascular cells" %in% ct_order) "Cerebrovascular cells" else "Endothelial"
if (vascular_ct_name %in% selected_for_subtype) {
  sub_objs[[vascular_ct_name]] <- run_subtype_subprocess(vascular_ct_name)
}
if ("Astrocytes" %in% selected_for_subtype) {
  sub_objs[["Astrocytes"]] <- run_subtype_subprocess("Astrocytes")
}
if ("Microglia" %in% selected_for_subtype) {
  sub_objs[["Microglia"]] <- run_subtype_subprocess("Microglia")
}

# ============================================================================
# PART 6: BULK DECONVOLUTION (ssGSEA + Wilcoxon + boxplot with p-values)
# ============================================================================
log_msg("================ PART 6: Bulk deconvolution back-validation ================")

bulk_eset_file <- file.path(BULK_DAT, "GSE132903_eSet.Rdata")
bulk_anno_file <- file.path(BULK_DAT, "GPL10558_bioc.rda")

if (!file.exists(bulk_eset_file) || !file.exists(bulk_anno_file)) {
  log_msg("  Bulk reference files missing - PART 6 skipped.")
} else {
  load(bulk_eset_file)
  if (exists("gset")) gse <- gset[[1]] else stop("Object 'gset' not found in eSet file.")
  pdata <- pData(gse)
  diag <- gsub("diagnosis: ", "", pdata$characteristics_ch1.3)
  group_list <- factor(ifelse(diag == "AD", "AD", "CN"), levels = c("CN", "AD"))
  expr_bulk <- exprs(gse)
  
  load(bulk_anno_file)
  anno_gpl <- GPL10558_bioc; colnames(anno_gpl) <- c("probe_id", "symbol")
  common_p <- intersect(rownames(expr_bulk), anno_gpl$probe_id)
  expr_bulk <- expr_bulk[common_p, , drop = FALSE]
  anno_gpl <- anno_gpl[match(common_p, anno_gpl$probe_id), , drop = FALSE]
  expr_bulk <- rowsum(expr_bulk, anno_gpl$symbol, na.rm = TRUE)
  expr_bulk <- expr_bulk[!rownames(expr_bulk) %in% c("", NA), , drop = FALSE]
  log_msg(sprintf("  Bulk: %d genes x %d samples (AD=%d, CN=%d)",
                  nrow(expr_bulk), ncol(expr_bulk),
                  sum(group_list == "AD"), sum(group_list == "CN")))
  
  # Build marker gene sets per subtype from each celltype subobject
  all_gs <- list()
  for (ct_name in names(sub_objs)) {
    obj <- sub_objs[[ct_name]]; if (is.null(obj)) next
    sub_col <- switch(ct_name,
                      "Cerebrovascular cells" = "vascular_subtype",
                      "Endothelial" = "vascular_subtype",
                      "Astrocytes"  = "astro_subtype",
                      "Microglia"   = "micro_subtype")
    if (!sub_col %in% colnames(obj[[]])) next
    Idents(obj) <- obj[[sub_col]][[1]]
    if (length(unique(Idents(obj))) < 2) next
    deg_file <- file.path(OUT_TAB, sprintf("subtype_DEGs_%s.csv", ct_name))
    if (file.exists(deg_file)) {
      degs <- read.csv(deg_file)
      # take top 50 per subtype, present in bulk
      gs <- split(degs, degs$cluster)
      gs <- lapply(gs, function(d) {
        d <- d[order(d$p_val_adj), ]
        head(intersect(d$gene, rownames(expr_bulk)), 50)
      })
      gs <- gs[sapply(gs, length) >= 5]
      if (length(gs) > 0)
        names(gs) <- paste0(ct_name, "_", names(gs))
      all_gs <- c(all_gs, gs)
    }
  }
  
  if (length(all_gs) >= 1) {
    log_msg(sprintf("  Subtype gene sets for ssGSEA: %d", length(all_gs)))
    
    gs_collection <- GeneSetCollection(mapply(function(g, nm)
      GeneSet(g, geneIdType = SymbolIdentifier(), setName = nm),
      all_gs, names(all_gs)))
    
    p_ssg <- ssgseaParam(as.matrix(expr_bulk), gs_collection,
                         minSize = 5, maxSize = 500, normalize = TRUE)
    bulk_ssgsea <- gsva(p_ssg, verbose = FALSE)
    
    safe_write_csv(
      cbind(geneset = rownames(bulk_ssgsea), as.data.frame(bulk_ssgsea)),
      file.path(OUT_TAB, "bulk_ssGSEA_subtype_scores.csv"))
    
    # Wilcoxon AD vs CN per geneset
    stat_df <- do.call(rbind, lapply(rownames(bulk_ssgsea), function(g) {
      v <- bulk_ssgsea[g, ]
      ad <- v[group_list == "AD"]; cn <- v[group_list == "CN"]
      if (length(ad) < 3 || length(cn) < 3)
        return(data.frame(geneset = g, mean_AD = NA, mean_CN = NA,
                          delta = NA, p = NA))
      w <- wilcox.test(ad, cn)
      data.frame(geneset = g,
                 mean_AD = mean(ad, na.rm = TRUE),
                 mean_CN = mean(cn, na.rm = TRUE),
                 delta = mean(ad, na.rm = TRUE) - mean(cn, na.rm = TRUE),
                 p = w$p.value)
    }))
    stat_df$padj <- p.adjust(stat_df$p, method = "BH")
    stat_df <- stat_df %>% arrange(padj)
    safe_write_csv(stat_df, file.path(OUT_TAB, "bulk_subtype_AD_vs_CN_stats.csv"))
    
    # Boxplot per geneset
    box_df <- as.data.frame(bulk_ssgsea) %>%
      rownames_to_column("geneset") %>%
      pivot_longer(-geneset, names_to = "sample", values_to = "score") %>%
      mutate(group = group_list[match(sample, colnames(bulk_ssgsea))])
    
    box_df$geneset <- factor(box_df$geneset,
                             levels = stat_df$geneset[order(stat_df$padj)])
    
    p_box <- ggplot(box_df, aes(group, score, fill = group)) +
      geom_boxplot(width = 0.6, outlier.size = 0.6, linewidth = 0.3) +
      geom_jitter(width = 0.15, size = 0.5, alpha = 0.6) +
      stat_compare_means(method = "wilcox.test",
                         label = "p.signif", size = 3,
                         comparisons = list(c("CN","AD"))) +
      facet_wrap(~ geneset, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = group_palette) +
      labs(x = NULL, y = "ssGSEA score",
           title = "Subtype gene-set enrichment in bulk (AD vs CN)") +
      publication_theme(7) + theme(legend.position = "none")
    n_sets <- nrow(stat_df)
    save_fig(p_box, "Fig_Bulk_Deconv_Boxplot", OUT_FIG,
             w = 9, h = 2.6 * ceiling(n_sets / 3))
    save_source(box_df, "Fig_Bulk_Deconv_Boxplot")
    
    # Heatmap of scores
    if (has_ComplexH) {
      suppressPackageStartupMessages(library(ComplexHeatmap))
      mat_z <- t(scale(t(bulk_ssgsea)))
      mat_z[mat_z > 2] <- 2; mat_z[mat_z < -2] <- -2
      ord <- order(group_list)
      mat_z <- mat_z[, ord]; gl_ord <- group_list[ord]
      
      ha <- HeatmapAnnotation(
        Group = gl_ord,
        col = list(Group = c("CN" = unname(group_palette["CN"]),
                             "AD" = unname(group_palette["AD"]))),
        annotation_name_side = "left",
        annotation_legend_param = list(Group = list(title = "Group")))
      ht <- Heatmap(mat_z,
                    name = "z-score",
                    top_annotation = ha,
                    cluster_columns = FALSE, cluster_rows = TRUE,
                    show_column_names = FALSE,
                    col = colorRampPalette(c("#3B4CC0","white","#B40426"))(100),
                    row_names_gp = gpar(fontsize = 7),
                    column_title = "Bulk ssGSEA z-scores")
      pdf(file.path(OUT_FIG, "Fig_Bulk_Deconv_Heatmap.pdf"), width = 8, height = 5)
      draw(ht); dev.off()
      png(file.path(OUT_FIG, "Fig_Bulk_Deconv_Heatmap.png"),
          width = 8, height = 5, units = "in", res = 300)
      draw(ht); dev.off()
      save_source(as.data.frame(mat_z) %>% rownames_to_column("geneset"),
                  "Fig_Bulk_Deconv_Heatmap")
    }
    
    # Volcano-like summary
    stat_df$signif <- ifelse(stat_df$padj < 0.05 & stat_df$delta > 0, "Up_AD",
                             ifelse(stat_df$padj < 0.05 & stat_df$delta < 0, "Down_AD","NS"))
    p_vol <- ggplot(stat_df, aes(delta, -log10(padj), color = signif)) +
      geom_point(size = 2.5, alpha = 0.85) +
      scale_color_manual(values = c("Up_AD" = "#E64B35",
                                    "Down_AD" = "#4DBBD5",
                                    "NS" = "grey70")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
      labs(x = "Delta ssGSEA (AD - CN)", y = expression(-log[10](FDR)),
           title = "Subtype-level AD enrichment in bulk", color = NULL) +
      publication_theme(7)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p_vol <- p_vol +
        ggrepel::geom_text_repel(
          data = subset(stat_df, padj < 0.05),
          aes(label = geneset), size = 2.2, max.overlaps = 20,
          show.legend = FALSE)
    }
    save_fig(p_vol, "Fig_Bulk_Deconv_Volcano", OUT_FIG, w = 5.5, h = 4)
    save_source(stat_df, "Fig_Bulk_Deconv_Volcano")
  } else {
    log_msg("  No subtype gene sets available - skipping ssGSEA.")
  }
}

# ============================================================================
# OUTPUT MANIFEST
# ============================================================================
manifest <- data.frame(
  file = c(
    "figures_main/Fig1_AD_scoring_overview.*",
    "figures_main/Fig1a_AD_signature_violin.*",
    "figures_main/Fig1b_method_heatmap.*",
    "figures_main/Fig1c_wilcoxon_delta.*",
    "figures_main/Fig1d_AD_score_UMAP.*",
    "figures_main/Fig_<CellType>_UMAP.*",
    "figures_main/Fig_<CellType>_DotPlot.*",
    "figures_main/Fig_<CellType>_Heatmap.pdf",
    "figures_main/Fig_<CellType>_Composition.*",
    "figures_main/Fig_Bulk_Deconv_Boxplot.*",
    "figures_main/Fig_Bulk_Deconv_Heatmap.*",
    "figures_main/Fig_Bulk_Deconv_Volcano.*",
    "figures_extended_data/ED_Fig1_per_method_violin.*",
    "figures_extended_data/ED_Fig2_method_correlation.*",
    "tables/scoring_per_cell.csv",
    "tables/scoring_summary_by_celltype.csv",
    "tables/scoring_wilcox_per_celltype.csv",
    "tables/celltype_consensus_rank.csv",
    "tables/celltype_selection_log.csv",
    "tables/subtype_assignment_<CellType>.csv",
    "tables/subtype_DEGs_<CellType>.csv",
    "tables/bulk_ssGSEA_subtype_scores.csv",
    "tables/bulk_subtype_AD_vs_CN_stats.csv",
    "source_data/*_source_data.csv",
    "rds/scoring_intermediate_sub.qs",
    "rds/scRNA_<CellType>_annotated.qs"
  ),
  description = c(
    "Main combined scoring overview (a-c)",
    "Combined AD-score violin (Fig 1a)",
    "5-method x celltype heatmap (Fig 1b)",
    "Wilcoxon delta with significance stars (Fig 1c)",
    "AD score on UMAP (Fig 1d)",
    "Subtype UMAP per cell type",
    "Subtype marker dot plot",
    "Subtype z-score heatmap",
    "Subtype composition by group",
    "Bulk subtype boxplot AD vs CN",
    "Bulk ssGSEA z-score heatmap",
    "Bulk-level volcano of subtype gene sets",
    "Per-method violin (ED Fig 1)",
    "Method-method Spearman heatmap (ED Fig 2)",
    "Per-cell 5-method scores (normalized)",
    "Mean per celltype x method",
    "Wilcoxon test per celltype (vs others)",
    "Inter-method consensus rank",
    "Selection rationale for subtype analysis",
    "Per-cluster subtype assignment + margin",
    "Per-subtype DEGs",
    "Bulk ssGSEA matrix",
    "Bulk Wilcoxon stats per subtype set",
    "Per-figure source data (publication output)",
    "Sub-sampled scored object",
    "Annotated subtype Seurat objects"
  ),
  stringsAsFactors = FALSE
)
safe_write_csv(manifest, file.path(OUT_TAB, "output_manifest.csv"))

capture.output(sessionInfo(), file = file.path(OUT_LOG, "sessionInfo_node04.txt"))

# ===== FINAL STEP: Compute singscore + regenerate dependent figures =====
# Reload seu_all (was freed before PART 3 for memory safety)
if (!exists("seu_all")) {
  ckpt_qs <- file.path(OUT_RDS, "seu_all_checkpoint.qs")
  if (file.exists(ckpt_qs)) {
    seu_all <- qs::qread(ckpt_qs)
    log_msg("seu_all reloaded for FINAL STEP.")
  }
}
singscore_csv <- file.path(OUT_TAB, "scoring_per_cell.csv")
need_singscore <- TRUE
if (file.exists(singscore_csv)) {
  csv_check <- read.csv(singscore_csv, nrow = 100)
  if (sum(is.na(csv_check$singscore)) == 0) need_singscore <- FALSE
}
if (need_singscore) {
  log_msg("Computing singscore from loaded seu_all (no qread)...")
  suppressPackageStartupMessages(library(singscore))
  up_final <- file.path(BULK_TAB, "bulk_up_geneset_AD_vs_CN.txt")
  if (!file.exists(up_final)) up_final <- file.path(BULK_TAB, "bulk_up_geneset.txt")
  up <- readLines(up_final)
  u <- intersect(up, rownames(seu_all))
  chunks <- split(colnames(seu_all), ceiling(seq_along(colnames(seu_all)) / 25000))
  sc_list <- list()
  for (i in seq_along(chunks)) {
    e <- as.matrix(GetAssayData(seu_all[, chunks[[i]]], assay = "RNA", layer = "data"))
    r <- rankGenes(e); ss <- simpleScore(r, upSet = u, knownDirection = TRUE)
    sc_list[[i]] <- setNames(ss$TotalScore, chunks[[i]])
    log_msg(sprintf("    singscore chunk %d/%d", i, length(chunks)))
    rm(e, r, ss); gc()
  }
  sc <- unlist(sc_list)
  log_msg(sprintf("  singscore done: NA=%d/%d", sum(is.na(sc)), length(sc)))
  csv <- read.csv(singscore_csv, stringsAsFactors = FALSE)
  csv$singscore <- sc[csv$cell]
  for (m in c("AddModuleScore","AUCell","UCell","singscore","ssGSEA")) {
    x <- csv[[m]]; if (all(is.na(x))) next
    rng <- range(x, na.rm = TRUE)
    csv[[m]] <- if (diff(rng) == 0) 0 else (x - rng[1]) / diff(rng)
  }
  csv$Combined <- rowMeans(csv[, c("AddModuleScore","AUCell","UCell","singscore","ssGSEA")], na.rm = TRUE)
  seu_all$singscore <- sc[colnames(seu_all)]
  seu_all$Combined <- csv$Combined[match(colnames(seu_all), csv$cell)]
  write.csv(csv, singscore_csv, row.names = FALSE)
  log_msg(paste("  singscore CSV saved: NA=", sum(is.na(csv$singscore)), "/", nrow(csv)))
  rm(sc, sc_list, chunks); gc()
  
  # Regenerate singscore-dependent figures
  log_msg("Regenerating scoring figures with singscore...")
  score_df_norm <- read.csv(singscore_csv)
  score_df_norm$celltype <- factor(score_df_norm$celltype, levels = intersect(ct_order, unique(score_df_norm$celltype)))
  
  # Fig1a
  p1a <- ggplot(score_df_norm, aes(celltype, Combined, fill = celltype)) +
    geom_violin(trim = T, scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.5) +
    scale_fill_manual(values = ct_cols_used) + theme_classic(base_size = 11) +
    labs(x = "", y = "AD Signature (Combined)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), legend.position = "none")
  save_fig(p1a, "Fig1a_AD_signature_violin", OUT_FIG, 6, 4.5)
  
  # Fig1b
  sm <- score_df_norm %>% group_by(celltype) %>%
    summarise(across(c(AddModuleScore, AUCell, UCell, singscore, ssGSEA, Combined), mean, .names = "{.col}"), .groups = "drop")
  sh <- sm %>% select(celltype, AddModuleScore, AUCell, UCell, singscore, ssGSEA, Combined) %>%
    pivot_longer(-celltype, names_to = "method", values_to = "score")
  sh$method <- factor(sh$method, rev(c("AddModuleScore","AUCell","UCell","singscore","ssGSEA","Combined")))
  p1b <- ggplot(sh, aes(celltype, method, fill = score)) + geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = round(score, 3)), size = 2.8) +
    scale_fill_gradientn(colors = c("#3B4CC0","white","#B40426")) +
    theme_classic(base_size = 11) + labs(x = "", y = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), panel.grid = element_blank())
  save_fig(p1b, "Fig1b_method_heatmap", OUT_FIG, 5.5, 3.5)
  
  # Fig1c
  pvals <- sapply(levels(score_df_norm$celltype), function(ct) {
    x <- score_df_norm$Combined[score_df_norm$celltype == ct]
    y <- score_df_norm$Combined[score_df_norm$celltype != ct]
    wilcox.test(x, y)$p.value
  })
  padj <- p.adjust(pvals, "BH")
  dd <- data.frame(celltype = names(padj), delta = sapply(names(padj), function(ct) {
    mean(score_df_norm$Combined[score_df_norm$celltype == ct]) - mean(score_df_norm$Combined)
  }), sig = cut(padj, c(-Inf, 0.001, 0.01, 0.05, Inf), c("***","**","*","ns")), stringsAsFactors = FALSE)
  p1c <- ggplot(dd, aes(reorder(celltype, delta), delta, fill = delta)) + geom_col(width = 0.7) + coord_flip() +
    scale_fill_gradient2(low = "#4DBBD5", mid = "grey95", high = "#E64B35") +
    geom_text(aes(label = paste0(sig, " p=", format(padj, scientific = T, digits = 2))), hjust = -0.02, size = 3) +
    expand_limits(y = max(dd$delta) * 1.35) + geom_hline(yintercept = 0) + theme_classic(base_size = 11) +
    labs(x = "", y = "\U0394 Combined (vs rest)") + theme(legend.position = "none")
  save_fig(p1c, "Fig1c_wilcoxon_delta", OUT_FIG, 6.5, 4)
  
  # Fig1_overview (patchwork)
  fig1 <- wrap_plots(p1a, p1b, p1c, ncol = 3, widths = c(1.8, 1.8, 1.2)) +
    plot_annotation(title = "AD Up-Regulated Gene Signature by Cell Type",
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10)))
  save_fig(fig1, "Fig1_AD_scoring_overview", OUT_FIG, 10.5, 4)
  
  # ED_Fig1
  ed1 <- score_df_norm %>% pivot_longer(c(AddModuleScore, AUCell, UCell, singscore, ssGSEA), names_to = "method", values_to = "value")
  p_ed1 <- ggplot(ed1, aes(celltype, value, fill = celltype)) +
    geom_violin(trim = T, scale = "width", alpha = 0.85, linewidth = 0.25) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.9, linewidth = 0.25) +
    scale_fill_manual(values = ct_cols_used) + facet_wrap(~method, scales = "free_y", ncol = 3) +
    labs(x = "", y = "Score (min-max scaled)", title = "Per-method scoring") +
    theme_classic(base_size = 6.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5), legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 9))
  ggsave(file.path(OUT_ED, "ED_Fig1_per_method_violin.pdf"), p_ed1, width = 9, height = 6, dpi = 300)
  
  # ED_Fig2
  cor_mat <- cor(score_df_norm[, c("AddModuleScore","AUCell","UCell","singscore","ssGSEA")], use = "pairwise", method = "spearman")
  cor_long <- as.data.frame(as.table(cor_mat)); names(cor_long) <- c("M1","M2","r")
  p_ed2 <- ggplot(cor_long, aes(M1, M2, fill = r)) + geom_tile(color = "white") +
    geom_text(aes(label = round(r, 2)), size = 3.5) +
    scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", limits = c(-1, 1)) +
    theme_classic(base_size = 10) + labs(x = "", y = "", title = "Method Spearman correlation") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
  ggsave(file.path(OUT_ED, "ED_Fig2_method_correlation.pdf"), p_ed2, width = 5.5, height = 4.5, dpi = 300)
  
  log_msg("Scoring figures regenerated with singscore.")
}

log_msg("================ ALL DONE ================")
log_msg(sprintf("Selected for subtype analysis: %s",
                paste(selected_for_subtype, collapse = ", ")))
log_msg("Outputs: figures_main/, figures_extended_data/, tables/, source_data/, rds/, logs/")




