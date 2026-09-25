set.seed(123456)  # reproducibility
options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 200 * 1024^3)
paths_to_add <- "/path/to/R_libraries"
for (p in paths_to_add) {
  if (dir.exists(p) && !(p %in% .libPaths())) .libPaths(c(p, .libPaths()))
}

# ---- TEST MODE (2026-05-20) ----
TEST_MODE <- TRUE
if (TEST_MODE) {
  INPUT_RDS <- "/path/to/project/results/01_atlas/rds"
  OUT      <- "/path/to/project/results/06_communication/outputs"
  RAW_DIR <- "/path/to/public_data"
  cat("[TEST MODE] Output to results/06_communication/outputs\n")
} else {
  INPUT_RDS <- "/path/to/project/results/01_atlas/rds"
  OUT      <- "/path/to/project/results/06_communication"
  RAW_DIR <- "/path/to/public_data"
}

# ==================================================================
# 06_communication_consolidated.R
# ==================================================================

OUT       <- if (TEST_MODE) "/path/to/project/results/06_communication/outputs" else "/path/to/project/results/06_communication"
OUT_FIGS  <- file.path(OUT, "main_figures")
OUT_TAB   <- file.path(OUT, "tables")
OUT_RDS   <- file.path(OUT, "rds")
OUT_METH  <- file.path(OUT, "methods")
OUT_RES   <- file.path(OUT, "results")
OUT_ATT   <- file.path(OUT, "attachments")
OUT_SUP   <- file.path(OUT, "supp_figures")

dir.create(OUT,      showWarnings=FALSE, recursive=TRUE)
dir.create(OUT_FIGS, showWarnings=FALSE, recursive=TRUE)
dir.create(OUT_TAB,  showWarnings=FALSE, recursive=TRUE)
dir.create(OUT_RDS,  showWarnings=FALSE, recursive=TRUE)
dir.create(OUT_METH, showWarnings=FALSE, recursive=TRUE)
dir.create(OUT_RES,  showWarnings=FALSE, recursive=TRUE)
dir.create(OUT_ATT,  showWarnings=FALSE, recursive=TRUE)
dir.create(OUT_SUP,  showWarnings=FALSE, recursive=TRUE)

setwd(OUT)

INPUT_RDS <- "/path/to/project/results/02_annotation/rds"
# ========== Part1: LIANA+GO+Reactome ==========
# File: scripts_final/node06_GO_Reactome_LIANA_SCTour_AD_NVU.R

################################################################################
# node06_FINAL_GO_Reactome_LIANA_CapEC_SCTour.R
#
# Project: AD-NVU endothelial communication analysis
#
# Purpose:
#   1. GO / Reactome enrichment for cerebrovascular / astrocyte / microglia subtypes
#   2. CN / AD LIANA communication landscape
#   3. Differential cerebrovascular-to-NVU communication: AD vs CN
#   4. Complete capillary endothelial functional module scoring
#   5. Capillary endothelial dysfunction / activation composite score
#   6. Capillary ActHigh / ActLow / ActMid LIANA
#   7. Curated publication-style figures
#   8. SCTour export for Endo / Astro / Micro / CCC / GO results
#
# Key design:
#   - Seurat v5 compatible
#   - LIANA 0.1.11 compatible as much as possible
#   - Convert Seurat to SingleCellExperiment before LIANA
#   - No hard re-clustering of CapEC
#   - Interpret CapEC as disease-state continuum, not transdifferentiation
#   - Hard-coded PROJECT_DIR
#   - No getwd()
#   - No cells.1 / cells.2 in FindMarkers
#   - No useDingbats argument
#   - No scale_size(na.value = ...)
#
################################################################################

################################################################################
# Module 0: Environment setup
################################################################################

PROJECT_DIR <- "/path/to/R_libraries"

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(enrichplot)
  library(patchwork)
})

if (!requireNamespace("liana", quietly = TRUE)) {
  stop("Package 'liana' is required. Please install LIANA before running this script.")
}

# Use sequential plan for stability on Windows
suppressPackageStartupMessages(library(future))
plan("sequential")

# Capture errors properly
options(error = function() {
  cat("\n!!! FATAL ERROR at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Last error message:", geterrmessage(), "\n")
  sink(NULL)
  q(save = "no", status = 1, runLast = FALSE)
})

# Fix namespace conflicts: BiocGenerics::select shadows dplyr::select
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate
arrange <- dplyr::arrange

set.seed(20260516)

CONFIG <- list(
  project_dir = PROJECT_DIR,

  full_rds  = "/path/to/project/results/02_annotation/rds/scRNA_annotated.qs",
  endo_rds  = "/path/to/project/results/04_subtypes/rds/scRNA_cereb_annotated.qs",
  astro_rds = "/path/to/project/results/04_subtypes/rds/subtype_Astrocytes.qs",
  micro_rds = "/path/to/project/results/04_subtypes/rds/subtype_Microglia.qs",
  out_tables = file.path(OUT, "tables"),
  out_plots  = file.path(OUT, "main_figures"),
  out_supp   = file.path(OUT, "supplementary_figures"),
  out_rds    = file.path(OUT, "rds"),
  out_sctour = file.path(OUT, "sctour"),
  out_logs   = file.path(OUT, "logs"),

  assay = "RNA",

  group_col = "Group",
  group_levels = c("CN", "AD"),

  major_col = "celltype",
  endo_sub_col = "vascular_subtype",
  astro_sub_col = "astro_subtype",
  micro_sub_col = "micro_subtype",

  min_cells_marker = 10,
  min_pct_marker = 0.10,
  logfc_marker = 0.25,
  padj_marker = 0.05,

  min_cells_liana = 5,
  expr_prop_liana = 0.05,
  max_cells_per_label = 200000,
  target_total = 500000,

  capec_high_quantile = 0.67,
  capec_low_quantile  = 0.33,

  top_n_liana_dotplot = 25,
  top_n_lr_barplot = 30,
  top_n_state_lr = 40,
  top_n_curated_lr = 25,

  run_sctour_export = TRUE,

  ## ----- Checkpoint config -----
  ckpt_dir  = file.path(OUT, "rds/checkpoints"),
  use_checkpoint = TRUE         # TRUE = resume from checkpoints if they exist
)

dir.create(CONFIG$out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG$out_plots,  recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG$out_supp,   recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG$out_rds,    recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG$out_sctour, recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG$out_logs,   recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG$ckpt_dir,   recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(CONFIG$out_logs, "node06_FINAL_run_log.txt")
sink(LOG_FILE, split = TRUE)

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf(...))
  cat("\n")
}

msg("Parallel plan: sequential (multisession was unstable on Windows)")

check_file <- function(path, label = path) {
  if (!file.exists(path)) {
    stop("Missing file: ", label, "\nPath: ", path)
  }
}

theme_publication <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text = element_text(face = "bold"),
      legend.key = element_blank(),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1)
    )
}

p_value_label <- function(p_value) {
  if (is.na(p_value)) return("n.s.")
  if (p_value < 1e-10) return("p<1e-10")
  paste0("p=", format.pval(p_value, digits = 2, eps = 1e-10))
}

wilcox_group_annotation <- function(df, value_col, group_col = "Group", facet_col = NULL, y_expand = 0.08) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  if (!group_col %in% colnames(df) || !value_col %in% colnames(df)) return(NULL)
  if (!all(c("CN", "AD") %in% unique(as.character(df[[group_col]])))) return(NULL)

  if (is.null(facet_col)) {
    tmp <- df %>%
      select(group = all_of(group_col), value = all_of(value_col)) %>%
      filter(!is.na(group), !is.na(value))
    if (nrow(tmp) == 0) return(NULL)
    pval <- tryCatch(wilcox.test(value ~ group, data = tmp)$p.value, error = function(e) NA_real_)
    ymax <- max(tmp$value, na.rm = TRUE)
    ymin <- min(tmp$value, na.rm = TRUE)
    if (!is.finite(ymax) || !is.finite(ymin)) return(NULL)
    data.frame(
      x = 1.5,
      y = ymax + (ymax - ymin) * y_expand,
      label = p_value_label(pval),
      stringsAsFactors = FALSE
    )
  } else {
    out <- df %>%
      select(facet = all_of(facet_col), group = all_of(group_col), value = all_of(value_col)) %>%
      filter(!is.na(facet), !is.na(group), !is.na(value)) %>%
      group_by(facet) %>%
      summarise(
        p_value = tryCatch(wilcox.test(value ~ group)$p.value, error = function(e) NA_real_),
        ymax = max(value, na.rm = TRUE),
        ymin = min(value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        x = 1,
        y = ymax + (ymax - ymin) * y_expand,
        label = vapply(p_value, p_value_label, character(1))
      )
    out
  }
}

safe_add_stat_layer <- function(plot, stat_df, facet_value_col = NULL) {
  if (is.null(stat_df) || nrow(stat_df) == 0) return(plot)
  if (is.null(facet_value_col)) {
    plot + annotate("text", x = 1.5, y = stat_df$y[[1]], label = stat_df$label[[1]], size = 2.8)
  } else {
    plot + geom_text(
      data = stat_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 2.6,
      vjust = 0
    )
  }
}

safe_ggsave_all <- function(prefix, plot, width = 7, height = 5, dpi = 300) {
  ggsave(paste0(prefix, ".png"), plot, width = width, height = height, dpi = dpi)
  ggsave(paste0(prefix, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(prefix, ".svg"), plot, width = width, height = height)
}


sanitize_label <- function(x) {
  x <- as.character(x)
  x <- ifelse(is.na(x) | x == "", "Unknown", x)
  x <- gsub("[[:space:]/\\-]+", "_", x)
  x <- gsub("[^A-Za-z0-9_]", "", x)
  x
}

harmonize_group_labels <- function(obj) {
  if (CONFIG$group_col %in% colnames(obj@meta.data)) {
    g <- as.character(obj@meta.data[[CONFIG$group_col]])
    g[g == "NC"] <- "CN"
    g[g == "UC"] <- "AD"
    obj@meta.data[[CONFIG$group_col]] <- factor(g, levels = CONFIG$group_levels)
  }
  if (CONFIG$major_col %in% colnames(obj@meta.data)) {
    ct <- as.character(obj@meta.data[[CONFIG$major_col]])
    ct[ct %in% c("Endothelial", "Endothelial_cells", "EC")] <- "Cerebrovascular cells"
    obj@meta.data[[CONFIG$major_col]] <- ct
  }
  obj
}

merge_subtype_to_full <- function(full_obj, subset_obj, subtype_col) {
  if (!subtype_col %in% colnames(subset_obj@meta.data)) {
    stop("Missing subtype column in subset object: ", subtype_col)
  }
  vals <- as.character(subset_obj@meta.data[[subtype_col]])
  names(vals) <- rownames(subset_obj@meta.data)
  full_obj@meta.data[[subtype_col]] <- NA_character_
  common <- intersect(rownames(full_obj@meta.data), names(vals))
  full_obj@meta.data[common, subtype_col] <- vals[common]
  full_obj
}

detect_sample_col <- function(obj) {
  candidates <- c("sample_id", "Sample", "sample", "orig.ident", "donor", "subject", "patient")
  hit <- candidates[candidates %in% colnames(obj@meta.data)]
  if (length(hit) == 0) {
    obj$sample_id_auto <- obj$orig.ident
    return("sample_id_auto")
  }
  hit[1]
}

get_data_matrix <- function(obj, assay = "RNA") {
  DefaultAssay(obj) <- assay
  mat <- tryCatch(
    GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) NULL
  )
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    mat <- GetAssayData(obj, assay = assay, layer = "counts")
  }
  mat
}

zscore_vec <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - m) / s
}

# ============================================================
# Checkpoint system: allows resume after interruption
# ============================================================

checkpoint_file <- function(module_name) {
  file.path(CONFIG$ckpt_dir, paste0("ckpt_", module_name, ".qs"))
}

checkpoint_exists <- function(module_name) {
  file.exists(checkpoint_file(module_name))
}

checkpoint_save <- function(module_name, ..., env = parent.frame()) {
  if (!CONFIG$use_checkpoint) return(invisible())

  obj_names <- as.character(substitute(...()))
  if (length(obj_names) == 0 || is.null(obj_names) || (length(obj_names) == 1 && obj_names[1] == "")) {
    writeLines("done", file.path(CONFIG$ckpt_dir, paste0("flag_", module_name, ".txt")))
    msg("Checkpoint SAVED (flag): %s", module_name)
    return(invisible())
  }

  to_save <- list()
  for (nm in obj_names) {
    to_save[[nm]] <- tryCatch(get(nm, envir = env), error = function(e) NULL)
  }
  to_save$.module <- module_name
  to_save$.timestamp <- Sys.time()

  qsave(to_save, checkpoint_file(module_name))
  msg("Checkpoint SAVED: %s (%s objects)", module_name, length(setdiff(names(to_save), c(".module", ".timestamp"))))
}

checkpoint_load <- function(module_name, env = parent.frame()) {
  f <- checkpoint_file(module_name)
  if (!file.exists(f)) return(FALSE)

  msg("Checkpoint LOADED: %s", module_name)
  data <- qread(f)
  for (nm in setdiff(names(data), c(".module", ".timestamp"))) {
    assign(nm, data[[nm]], envir = env)
  }
  msg("  -> restored objects: %s", paste(setdiff(names(data), c(".module", ".timestamp")), collapse = ", "))
  TRUE
}

checkpoint_flag_done <- function(module_name) {
  f <- file.path(CONFIG$ckpt_dir, paste0("flag_", module_name, ".txt"))
  writeLines("done", f)
  msg("Checkpoint FLAG set: %s", module_name)
}

################################################################################
# Module 1: Load objects and metadata harmonization
################################################################################

msg("=====================================================================")
msg("Module 1: Loading Seurat objects")
msg("=====================================================================")

check_file(CONFIG$full_rds,  "full Seurat object")
check_file(CONFIG$endo_rds,  "endo Seurat object")
check_file(CONFIG$astro_rds, "astro Seurat object")
check_file(CONFIG$micro_rds, "micro Seurat object")

sc_full  <- qread(CONFIG$full_rds)
sc_endo  <- qread(CONFIG$endo_rds)
sc_astro <- qread(CONFIG$astro_rds)
sc_micro <- qread(CONFIG$micro_rds)

sc_full  <- harmonize_group_labels(sc_full)
sc_endo  <- harmonize_group_labels(sc_endo)
sc_astro <- harmonize_group_labels(sc_astro)
sc_micro <- harmonize_group_labels(sc_micro)

sc_full <- merge_subtype_to_full(sc_full, sc_endo,  CONFIG$endo_sub_col)
sc_full <- merge_subtype_to_full(sc_full, sc_astro, CONFIG$astro_sub_col)
sc_full <- merge_subtype_to_full(sc_full, sc_micro, CONFIG$micro_sub_col)

DefaultAssay(sc_full)  <- CONFIG$assay
DefaultAssay(sc_endo)  <- CONFIG$assay
DefaultAssay(sc_astro) <- CONFIG$assay
DefaultAssay(sc_micro) <- CONFIG$assay

sample_col <- detect_sample_col(sc_full)

################################################################################
# ???1???????? - ???????????
################################################################################
verify_data_integrity <- function(sc_full, sc_endo, sc_astro, sc_micro) {
  cat("\n===== [QC] node06 input integrity =====\n")
  errors <- character()
  if (CONFIG$endo_sub_col %in% colnames(sc_endo@meta.data)) {
    actual_endo <- sort(table(sc_endo@meta.data[[CONFIG$endo_sub_col]], useNA='ifany'), decreasing=TRUE)
    cat("[QC] sc_endo vascular_subtype:\n")
    print(actual_endo)
    total_endo <- sum(actual_endo)
    cat(sprintf("[QC] sc_endo cells: %d\n", total_endo))
  } else {
    msg("[ERROR] sc_endo missing %s column", CONFIG$endo_sub_col)
    errors <- c(errors, paste0("sc_endo missing ", CONFIG$endo_sub_col, " column"))
  }
  overlap <- intersect(colnames(sc_full), colnames(sc_endo))
  cat(sprintf("[QC] sc_full / sc_endo overlap: %d (expected: %d)\n", length(overlap), ncol(sc_endo)))
  if (length(overlap) != ncol(sc_endo)) {
    msg("[ERROR] sc_full ? sc_endo ???????!")
    errors <- c(errors, sprintf("sc_full/sc_endo overlap: expected %d, got %d", ncol(sc_endo), length(overlap)))
  }
  if (length(errors) > 0) msg("[WARN] Data integrity check failed: %s", paste(errors, collapse="; "))
  TRUE
}

tryCatch({ verify_data_integrity(sc_full, sc_endo, sc_astro, sc_micro) }, error = function(e) { msg("[FATAL] ??????: %s", e$message); stop(e) })
################################################################################
required_full_cols <- c(CONFIG$group_col, CONFIG$major_col)
missing_full_cols <- setdiff(required_full_cols, colnames(sc_full@meta.data))
if (length(missing_full_cols) > 0) {
  stop("Missing required columns in sc_full metadata: ", paste(missing_full_cols, collapse = ", "))
}

sc_full[[CONFIG$group_col]] <- factor(
  sc_full[[CONFIG$group_col]][, 1],
  levels = CONFIG$group_levels
)

if (CONFIG$group_col %in% colnames(sc_endo@meta.data)) {
  sc_endo[[CONFIG$group_col]] <- factor(sc_endo[[CONFIG$group_col]][, 1], levels = CONFIG$group_levels)
}
if (CONFIG$group_col %in% colnames(sc_astro@meta.data)) {
  sc_astro[[CONFIG$group_col]] <- factor(sc_astro[[CONFIG$group_col]][, 1], levels = CONFIG$group_levels)
}
if (CONFIG$group_col %in% colnames(sc_micro@meta.data)) {
  sc_micro[[CONFIG$group_col]] <- factor(sc_micro[[CONFIG$group_col]][, 1], levels = CONFIG$group_levels)
}

msg("Full cells:  %d", ncol(sc_full))
msg("Cerebrovascular cells: %d", ncol(sc_endo))
msg("Astro cells: %d", ncol(sc_astro))
msg("Micro cells: %d", ncol(sc_micro))
msg("Sample column used: %s", sample_col)

################################################################################
# Module 1.1: Construct LIANA labels
################################################################################

msg("Constructing liana_label")

meta <- sc_full@meta.data

meta$liana_label <- sanitize_label(meta[[CONFIG$major_col]])

if (CONFIG$endo_sub_col %in% colnames(meta)) {
  is_endo <- !is.na(meta[[CONFIG$endo_sub_col]]) &
    meta[[CONFIG$endo_sub_col]] != ""
  meta$liana_label[is_endo] <- paste0("Cerebrovascular_", sanitize_label(meta[[CONFIG$endo_sub_col]][is_endo]))
}

if (CONFIG$astro_sub_col %in% colnames(meta)) {
  is_astro <- !is.na(meta[[CONFIG$astro_sub_col]]) &
    meta[[CONFIG$astro_sub_col]] != ""
  meta$liana_label[is_astro] <- paste0("Astro_", sanitize_label(meta[[CONFIG$astro_sub_col]][is_astro]))
}

if (CONFIG$micro_sub_col %in% colnames(meta)) {
  is_micro <- !is.na(meta[[CONFIG$micro_sub_col]]) &
    meta[[CONFIG$micro_sub_col]] != ""
  meta$liana_label[is_micro] <- paste0("Micro_", sanitize_label(meta[[CONFIG$micro_sub_col]][is_micro]))
}

sc_full$liana_label <- meta$liana_label

write.csv(
  as.data.frame(table(sc_full$liana_label, sc_full[[CONFIG$group_col]][, 1])),
  file.path(CONFIG$out_tables, "node06_LIANA_label_cellcounts_before_downsampling.csv"),
  row.names = FALSE
)

################################################################################
# Module 2: Marker genes and GO / Reactome enrichment
################################################################################

msg("=====================================================================")
msg("Module 2: Markers and GO / Reactome enrichment")
msg("=====================================================================")

symbol_to_entrez <- function(genes) {
  genes <- unique(genes)
  genes <- genes[!is.na(genes) & genes != ""]
  if (length(genes) == 0) return(character(0))

  conv <- suppressMessages(
    bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )

  unique(conv$ENTREZID)
}

run_enrich_one <- function(gene_symbols, subtype, family) {
  entrez <- symbol_to_entrez(gene_symbols)

  if (length(entrez) < 5) {
    return(list(go = NULL, reactome = NULL))
  }

  ego <- tryCatch(
    enrichGO(
      gene          = entrez,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.20,
      readable      = TRUE
    ),
    error = function(e) NULL
  )

  ere <- tryCatch(
    enrichPathway(
      gene          = entrez,
      organism      = "human",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.20,
      readable      = TRUE
    ),
    error = function(e) NULL
  )

  go_df <- NULL
  re_df <- NULL

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    go_df <- as.data.frame(ego) %>%
      mutate(subtype = subtype, family = family, database = "GO_BP")
  }

  if (!is.null(ere) && nrow(as.data.frame(ere)) > 0) {
    re_df <- as.data.frame(ere) %>%
      mutate(subtype = subtype, family = family, database = "Reactome")
  }

  list(go = go_df, reactome = re_df)
}

run_subtype_enrichment <- function(obj, subtype_col, family) {
  if (!subtype_col %in% colnames(obj@meta.data)) {
    warning("Subtype column not found: ", subtype_col)
    return(list(go = NULL, reactome = NULL, markers = NULL))
  }

  obj <- subset(obj, cells = colnames(obj)[!is.na(obj[[subtype_col]][, 1]) & obj[[subtype_col]][, 1] != ""])
  if (ncol(obj) == 0) {
    return(list(go = NULL, reactome = NULL, markers = NULL))
  }

  Idents(obj) <- obj[[subtype_col]][, 1]

  tab <- table(Idents(obj))
  keep_levels <- names(tab)[tab >= CONFIG$min_cells_marker]
  obj <- subset(obj, idents = keep_levels)

  if (length(unique(Idents(obj))) < 2) {
    warning("Less than two valid subtypes for ", family)
    return(list(go = NULL, reactome = NULL, markers = NULL))
  }

  # Seurat V5: JoinLayers before FindAllMarkers
  obj <- JoinLayers(obj)

  markers <- tryCatch(
    FindAllMarkers(
      obj,
      assay = CONFIG$assay,
      only.pos = TRUE,
      min.pct = CONFIG$min_pct_marker,
      logfc.threshold = CONFIG$logfc_marker,
      test.use = "wilcox"
    ),
    error = function(e) {
      warning("FindAllMarkers failed for ", family, ": ", e$message)
      NULL
    }
  )

  if (is.null(markers) || nrow(markers) == 0) {
    return(list(go = NULL, reactome = NULL, markers = NULL))
  }

  markers <- markers %>%
    mutate(family = family, subtype = as.character(cluster)) %>%
    filter(p_val_adj <= CONFIG$padj_marker)

  enrich_list <- markers %>%
    group_by(subtype) %>%
    summarise(genes = list(unique(gene)), .groups = "drop") %>%
    mutate(enrich = map2(genes, subtype, ~ run_enrich_one(.x, .y, family)))

  go_df <- bind_rows(map(enrich_list$enrich, "go"))
  re_df <- bind_rows(map(enrich_list$enrich, "reactome"))

  list(go = go_df, reactome = re_df, markers = markers)
}

# Module 2 checkpoint skip
module2_done <- file.exists(file.path(CONFIG$out_tables, "node06_GO_BP_enrichment_all.csv")) &&
  file.exists(file.path(CONFIG$out_tables, "node06_Reactome_enrichment_all.csv")) &&
  file.exists(file.path(CONFIG$out_tables, "node06_subtype_positive_markers_all.csv"))
if (module2_done) {
  msg("  Module 2 enrichment files found - skipping FindAllMarkers/GO/Reactome")
  markers_all <- read.csv(file.path(CONFIG$out_tables, "node06_subtype_positive_markers_all.csv"), stringsAsFactors = FALSE)
  go_all <- read.csv(file.path(CONFIG$out_tables, "node06_GO_BP_enrichment_all.csv"), stringsAsFactors = FALSE)
  reactome_all <- read.csv(file.path(CONFIG$out_tables, "node06_Reactome_enrichment_all.csv"), stringsAsFactors = FALSE)
} else {
  res_endo  <- run_subtype_enrichment(sc_endo,  CONFIG$endo_sub_col,  "Cerebrovascular")
  res_astro <- run_subtype_enrichment(sc_astro, CONFIG$astro_sub_col, "Astro")
  res_micro <- run_subtype_enrichment(sc_micro, CONFIG$micro_sub_col, "Micro")
  go_all <- bind_rows(res_endo$go, res_astro$go, res_micro$go)
  reactome_all <- bind_rows(res_endo$reactome, res_astro$reactome, res_micro$reactome)
  markers_all <- bind_rows(res_endo$markers, res_astro$markers, res_micro$markers)
}

write.csv(go_all, file.path(CONFIG$out_tables, "node06_GO_BP_enrichment_all.csv"), row.names = FALSE)
write.csv(reactome_all, file.path(CONFIG$out_tables, "node06_Reactome_enrichment_all.csv"), row.names = FALSE)
write.csv(markers_all, file.path(CONFIG$out_tables, "node06_subtype_positive_markers_all.csv"), row.names = FALSE)

qsave(go_all, file.path(CONFIG$out_rds, "node06_GO_BP_enrichment_all.qs"))
qsave(reactome_all, file.path(CONFIG$out_rds, "node06_Reactome_enrichment_all.qs"))
qsave(markers_all, file.path(CONFIG$out_rds, "node06_subtype_positive_markers_all.qs"))
################################################################################
# Module 3: GO / Reactome dotplots
################################################################################

msg("=====================================================================")
msg("Module 3: GO / Reactome dotplots")
msg("=====================================================================")

make_enrich_dotplot <- function(df, family_name, title_text) {
  if (is.null(df) || nrow(df) == 0) {
    return(
      ggplot() +
        theme_void() +
        labs(title = paste0(title_text, " - no enrichment"))
    )
  }

  tmp <- df %>%
    filter(family == family_name) %>%
    group_by(subtype) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 8) %>%
    ungroup()

  if (nrow(tmp) == 0) {
    return(
      ggplot() +
        theme_void() +
        labs(title = paste0(title_text, " - no enrichment"))
    )
  }

  tmp <- tmp %>%
    mutate(
      Description = str_trunc(Description, 60),
      neglog10FDR = -log10(p.adjust + 1e-300)
    )

  ggplot(tmp, aes(x = subtype, y = reorder(Description, neglog10FDR))) +
    geom_point(aes(size = Count, color = neglog10FDR), alpha = 0.85) +
    scale_color_gradient(low = "#3C5488", high = "#E64B35") +
    theme_publication(base_size = 8) +
    labs(
      title = title_text,
      x = NULL,
      y = NULL,
      color = expression(-log[10]("FDR")),
      size = "Count"
    )
}

for (fam in c("Cerebrovascular", "Astro", "Micro")) {
  p_go <- make_enrich_dotplot(go_all, fam, paste0(fam, ": GO BP enrichment"))
  p_re <- make_enrich_dotplot(reactome_all, fam, paste0(fam, ": Reactome enrichment"))

  safe_ggsave_all(
    file.path(CONFIG$out_plots, paste0("Fig_node06_GO_DotPlot_", fam)),
    p_go,
    width = 7.2,
    height = 6.0
  )

  safe_ggsave_all(
    file.path(CONFIG$out_plots, paste0("Fig_node06_Reactome_DotPlot_", fam)),
    p_re,
    width = 7.2,
    height = 6.0
  )
}

msg("=====================================================================")
msg("Module 4: LIANA helper functions")
msg("=====================================================================")

downsample_by_label <- function(obj, label_col, max_cells_per_label = 2000, target_total = 30000, seed = 20260516, protect_labels = c("Cerebrovascular cells", "CapEC_ActHigh", "CapEC_ActMid", "CapEC_ActLow", "Cerebrovascular_Capillary", "Cerebrovascular_Pericyte", "Cerebrovascular_SMC", "Cerebrovascular_Arterial", "Cerebrovascular_Venous", "Capillary", "Pericyte", "SMC", "Arterial", "Venous")) {
  set.seed(seed)
  meta <- obj@meta.data
  labels <- meta[[label_col]]
  cells_by_label <- split(colnames(obj), labels)
  sampled <- unlist(lapply(names(cells_by_label), function(label) {
    cells <- cells_by_label[[label]]
    if (label %in% protect_labels) return(cells)
    if (length(cells) > max_cells_per_label) sample(cells, max_cells_per_label) else cells
  }))
  if (length(sampled) > target_total) {
    sampled <- sample(sampled, target_total)
  }
  subset(obj, cells = sampled)
}

seurat_to_sce_for_liana <- function(obj, label_col, assay = "RNA") {
  DefaultAssay(obj) <- assay
  
  counts <- tryCatch(
    GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) NULL
  )

  data <- tryCatch(
    GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) NULL
  )

  if (is.null(counts)) {
    counts <- data
  }

  if (is.null(data)) {
    data <- counts
  }

  sce <- SingleCellExperiment(
    assays = list(
      counts = counts,
      logcounts = data
    ),
    colData = obj@meta.data
  )

  colData(sce)[[label_col]] <- as.character(obj@meta.data[[label_col]])

  sce
}

normalize_liana_result <- function(df) {
  # LIANA 0.1.11 returns list with $aggregate element
  if (is.list(df) && !is.data.frame(df)) {
    if ("aggregate" %in% names(df)) {
      df <- as.data.frame(df$aggregate)
    } else if ("liana_scores" %in% names(df)) {
      df <- as.data.frame(df$liana_scores)
    } else {
      # Try first data.frame element
      for (elem in df) {
        if (is.data.frame(elem) && nrow(elem) > 0) {
          df <- as.data.frame(elem)
          break
        }
      }
    }
  } else {
    df <- as.data.frame(df)
  }

  if (nrow(df) == 0) stop("Empty LIANA result")

  # Standardize column names: LIANA uses ligand.complex/receptor.complex
  # Create ligand/receptor aliases for downstream compatibility
  if ("ligand.complex" %in% colnames(df)) df$ligand <- df$ligand.complex
  if ("receptor.complex" %in% colnames(df)) df$receptor <- df$receptor.complex
  if (!"source" %in% colnames(df) && "source_std" %in% colnames(df)) df$source <- df$source_std
  if (!"target" %in% colnames(df) && "target_std" %in% colnames(df)) df$target <- df$target_std

  required <- c("source", "target", "ligand", "receptor")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("LIANA result missing columns: ", paste(missing, collapse = ", "))
  }

  # Compute communication strength from aggregate_rank (primary) or magnitude_std (fallback)
  if ("aggregate_rank" %in% colnames(df)) {
    df$communication_strength <- -log10(pmax(df$aggregate_rank, 1e-12))
  } else if ("magnitude_std" %in% colnames(df)) {
    df$communication_strength <- df$magnitude_std
  } else if ("magnitude_rank" %in% colnames(df)) {
    df$communication_strength <- -log10(pmax(df$magnitude_rank, 1e-12))
  } else {
    df$communication_strength <- 1
  }

  df
}

run_liana_core <- function(obj, label_col) {
  sce <- seurat_to_sce_for_liana(obj, label_col = label_col, assay = CONFIG$assay)

  # LIANA: explicitly specify methods, then aggregate manually (v0.1.11 compatible)
  res <- tryCatch({
      raw_res <- liana::liana_wrap(
        sce,
        method = c("natmi", "connectome", "logfc", "sca", "cellphonedb"),
        idents_col = label_col,
        assay.type = "logcounts",
        expr_prop = CONFIG$expr_prop_liana,
        verbose = FALSE
      )

      # Check if raw_res has an aggregate element
      if (is.list(raw_res) && "aggregate" %in% names(raw_res)) {
        as.data.frame(raw_res$aggregate)
      } else if (is.list(raw_res) && !is.data.frame(raw_res)) {
        # Manually aggregate if aggregate not automatically returned
        agg <- tryCatch(liana::liana_aggregate(raw_res), error = function(e) NULL)
        if (!is.null(agg)) {
          msg("LIANA: manual aggregation successful")
          as.data.frame(agg)
        } else {
          # Fallback: use cellphonedb results as primary
          msg("LIANA: using cellphonedb as primary (aggregation unavailable)")
          df <- as.data.frame(raw_res$cellphonedb)
          if (nrow(df) > 0) df else NULL
        }
      } else {
        as.data.frame(raw_res)
      }
    }, error = function(e1) {
      msg("LIANA failed: %s", e1$message)
      NULL
    }
  )

  if (is.null(res)) stop("LIANA returned NULL after all attempts")

  msg("LIANA result: %d rows", nrow(res))

  # Standardize: create ligand/receptor from .complex columns
  if ("ligand.complex" %in% colnames(res) && !"ligand" %in% colnames(res)) res$ligand <- res$ligand.complex
  if ("receptor.complex" %in% colnames(res) && !"receptor" %in% colnames(res)) res$receptor <- res$receptor.complex

  # Compute communication_strength
  if ("aggregate_rank" %in% colnames(res)) {
    res$communication_strength <- -log10(pmax(res$aggregate_rank, 1e-12))
  } else if ("mean_rank" %in% colnames(res)) {
    res$communication_strength <- -log10(pmax(res$mean_rank, 1e-12))
  } else if ("lr.mean" %in% colnames(res)) {
    # cellphonedb-style
    res$communication_strength <- res$lr.mean
  } else if ("magnitude_rank" %in% colnames(res)) {
    res$communication_strength <- -log10(pmax(res$magnitude_rank, 1e-12))
  } else {
    msg("Warning: no aggregate_rank/mean_rank/lr.mean in LIANA output")
    res$communication_strength <- 1
  }

  res
}

run_liana_one_group <- function(obj, group_value, label_col = "liana_label", prefix = "node06") {
  msg("Running LIANA for group: %s; label_col: %s", group_value, label_col)

  cells_use <- rownames(obj@meta.data)[obj@meta.data[[CONFIG$group_col]] == group_value]
  obj_g <- subset(obj, cells = cells_use)

  counts_df <- obj_g@meta.data %>%
    count(.data[[label_col]], name = "n_cells") %>%
    rename(label = 1) %>%
    arrange(desc(n_cells))

  valid_labels <- counts_df %>%
    filter(n_cells >= CONFIG$min_cells_liana) %>%
    pull(label)

  obj_g <- subset(obj_g, cells = rownames(obj_g@meta.data)[obj_g@meta.data[[label_col]] %in% valid_labels])

  if (length(unique(obj_g@meta.data[[label_col]])) < 2) {
    warning("Less than two valid LIANA labels for group: ", group_value)
    return(list(raw_df = NULL, agg_df = NULL, endo_nvu_df = NULL, counts = counts_df))
  }

  obj_g <- downsample_by_label(
    obj_g,
    label_col = label_col,
    max_cells_per_label = CONFIG$max_cells_per_label,
    target_total = CONFIG$target_total,
    seed = 20260516 + ifelse(group_value == "AD", 1, 0)
  )

  counts_after <- obj_g@meta.data %>%
    count(.data[[label_col]], name = "n_cells") %>%
    rename(label = 1) %>%
    mutate(group = group_value)

  raw_df <- run_liana_core(obj_g, label_col = label_col)
  raw_df$group <- group_value

  agg_df <- raw_df %>%
    group_by(group, source, target) %>%
    summarise(
      n_lr = n(),
      mean_strength = mean(communication_strength, na.rm = TRUE),
      median_strength = median(communication_strength, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_strength))

  endo_nvu_df <- raw_df %>%
    filter(str_detect(source, "^Cerebrovascular_|^CapEC_")) %>%
    arrange(desc(communication_strength))

  write.csv(
    raw_df,
    file.path(CONFIG$out_tables, paste0(prefix, "_LIANA_raw_", group_value, ".csv")),
    row.names = FALSE
  )
  write.csv(
    agg_df,
    file.path(CONFIG$out_tables, paste0(prefix, "_LIANA_aggregate_all_", group_value, ".csv")),
    row.names = FALSE
  )
  write.csv(
    endo_nvu_df,
    file.path(CONFIG$out_tables, paste0(prefix, "_LIANA_Cerebrovascular_to_NVU_", group_value, ".csv")),
    row.names = FALSE
  )
  write.csv(
    counts_after,
    file.path(CONFIG$out_tables, paste0(prefix, "_LIANA_cellcounts_after_downsampling_", group_value, ".csv")),
    row.names = FALSE
  )

  list(raw_df = raw_df, agg_df = agg_df, endo_nvu_df = endo_nvu_df, counts = counts_after)
}

################################################################################
# Module 5: LIANA CN / AD full NVU (SLOW ~30-60 min)
################################################################################

msg("=====================================================================")
msg("Module 5: LIANA CN / AD full NVU")
msg("=====================================================================")

CKPT_LIANA_CN <- file.path(CONFIG$out_rds, "node06_LIANA_CN.qs")
CKPT_LIANA_AD <- file.path(CONFIG$out_rds, "node06_LIANA_AD.qs")

if (file.exists(CKPT_LIANA_CN) && file.exists(CKPT_LIANA_AD)) {
  msg("LIANA checkpoints found - loading cached results (to resume fresh, delete %s/*.qs)", CONFIG$out_rds)
  liana_CN <- qread(CKPT_LIANA_CN)
  liana_AD <- qread(CKPT_LIANA_AD)
} else {
  sc_full_typed <- sc_full
  cat(sprintf("[LIANA] sc_full_typed: %d cells\n", ncol(sc_full_typed)))

  liana_CN <- run_liana_one_group(sc_full_typed, "CN", label_col = "liana_label", prefix = "node06")
  liana_AD <- run_liana_one_group(sc_full_typed, "AD", label_col = "liana_label", prefix = "node06")
  qsave(liana_CN, CKPT_LIANA_CN)
  qsave(liana_AD, CKPT_LIANA_AD)
}

cellcounts_all <- bind_rows(liana_CN$counts, liana_AD$counts)

write.csv(
  cellcounts_all,
  file.path(CONFIG$out_tables, "node06_CellCounts_NVU_after_downsampling.csv"),
  row.names = FALSE
)

qsave(liana_CN, file.path(CONFIG$out_rds, "node06_LIANA_CN.qs"))
qsave(liana_AD, file.path(CONFIG$out_rds, "node06_LIANA_AD.qs"))

################################################################################
# Module 6: Differential CCC AD vs CN
################################################################################

msg("=====================================================================")
msg("Module 6: Differential cerebrovascular-to-NVU communication AD vs CN")
msg("=====================================================================")

diff_ccc <- NULL

if (!is.null(liana_CN$endo_nvu_df) && !is.null(liana_AD$endo_nvu_df)) {
  cn <- liana_CN$endo_nvu_df %>%
    select(source, target, ligand, receptor, communication_strength_CN = communication_strength)

  ad <- liana_AD$endo_nvu_df %>%
    select(source, target, ligand, receptor, communication_strength_AD = communication_strength)

  diff_ccc <- full_join(
    ad,
    cn,
    by = c("source", "target", "ligand", "receptor")
  ) %>%
    mutate(
      communication_strength_AD = replace_na(communication_strength_AD, 0),
      communication_strength_CN = replace_na(communication_strength_CN, 0),
      delta_AD_vs_CN = communication_strength_AD - communication_strength_CN,
      direction = case_when(
        delta_AD_vs_CN > 0 ~ "AD_up",
        delta_AD_vs_CN < 0 ~ "AD_down",
        TRUE ~ "unchanged"
      )
    ) %>%
    arrange(desc(abs(delta_AD_vs_CN)))

  write.csv(
    diff_ccc,
    file.path(CONFIG$out_tables, "node06_Differential_CCC_Cerebrovascular_to_NVU_AD_vs_CN.csv"),
    row.names = FALSE
  )
}

msg("=====================================================================")
msg("Module 7: Complete CapEC functional module scoring")
msg("=====================================================================")

if (!CONFIG$endo_sub_col %in% colnames(sc_endo@meta.data)) {
  stop(CONFIG$endo_sub_col, " column not found in sc_endo")
}

cap_cells <- rownames(sc_endo@meta.data)[
  sc_endo@meta.data[[CONFIG$endo_sub_col]] %in% c("Capillary", "CapEC", "Capillary_EC")
]

if (length(cap_cells) < 10) {
  stop("Too few capillary cells detected. Please check ", CONFIG$endo_sub_col, " values.")
}

sc_cap <- subset(sc_endo, cells = cap_cells)
DefaultAssay(sc_cap) <- CONFIG$assay

capec_gene_sets <- list(
  CapEC_BBB_maintenance = c(
    "CLDN5", "OCLN", "TJP1", "TJP2", "JAM2", "JAM3",
    "PECAM1", "CDH5", "ESAM", "SLC2A1", "MFSD2A",
    "ABCB1", "ABCG2", "LRP1", "KDR", "FLT1", "SOX17"
  ),
  CapEC_Endothelial_activation = c(
    "ICAM1", "VCAM1", "SELE", "SELP", "EGR1", "FOS",
    "JUN", "RELA", "NFKBIA", "TNFAIP3", "EDN1",
    "ANGPT2", "PLAU", "PLAUR", "SERPINE1"
  ),
  CapEC_Inflammatory = c(
    "CCL2", "CCL3", "CCL4", "CCL5", "CXCL1", "CXCL2",
    "CXCL3", "CXCL8", "CXCL10", "IL6", "PTGS2",
    "NFKBIA", "TNFAIP3", "IRAK2"
  ),
  CapEC_IFN_response = c(
    "ISG15", "IFI6", "IFI27", "IFIT1", "IFIT2", "IFIT3",
    "MX1", "MX2", "OAS1", "OAS2", "OAS3", "STAT1",
    "IRF7", "CXCL10"
  ),
  CapEC_ECM_remodeling = c(
    "COL4A1", "COL4A2", "COL18A1", "FN1", "LAMA4",
    "LAMB1", "LAMC1", "MMP2", "MMP14", "TIMP1",
    "TIMP2", "SPARC", "SERPINE1"
  ),
  CapEC_Permeability_angiogenic = c(
    "VEGFA", "KDR", "FLT1", "ANGPT2", "TEK", "PLVAP",
    "VWF", "ESM1", "APLN", "APLNR", "DLL4", "JAG1",
    "ADM", "RAMP2"
  ),
  CapEC_Transport = c(
    "SLC2A1", "MFSD2A", "ABCB1", "ABCG2", "SLC7A5",
    "SLC16A1", "SLC38A5", "TFRC", "LRP1", "INSR",
    "ATP1A1", "ATP1B1"
  )
)

score_module <- function(obj, genes, score_name, assay = "RNA") {
  mat <- get_data_matrix(obj, assay = assay)
  genes_present <- intersect(genes, rownames(mat))

  if (length(genes_present) == 0) {
    warning("No genes present for module: ", score_name)
    return(rep(NA_real_, ncol(obj)))
  }

  score <- Matrix::colMeans(mat[genes_present, , drop = FALSE])
  as.numeric(score)
}

for (nm in names(capec_gene_sets)) {
  raw_col <- paste0(nm, "_score")
  z_col <- paste0(nm, "_score_z")

  sc_cap[[raw_col]] <- score_module(sc_cap, capec_gene_sets[[nm]], raw_col, CONFIG$assay)
  sc_cap[[z_col]] <- zscore_vec(sc_cap[[raw_col]][, 1])
}

sc_cap$CapEC_BBB_loss_score_z <- -1 * rowMeans(
  sc_cap@meta.data[, c(
    "CapEC_BBB_maintenance_score_z",
    "CapEC_Transport_score_z"
  ), drop = FALSE],
  na.rm = TRUE
)

dys_cols <- c(
  "CapEC_Endothelial_activation_score_z",
  "CapEC_Inflammatory_score_z",
  "CapEC_IFN_response_score_z",
  "CapEC_ECM_remodeling_score_z",
  "CapEC_Permeability_angiogenic_score_z",
  "CapEC_BBB_loss_score_z"
)

sc_cap$CapEC_dysfunction_composite_score <- rowMeans(
  sc_cap@meta.data[, dys_cols, drop = FALSE],
  na.rm = TRUE
)

sc_cap$CapEC_dysfunction_composite_score_z <- zscore_vec(
  sc_cap$CapEC_dysfunction_composite_score
)

q_high <- quantile(
  sc_cap$CapEC_dysfunction_composite_score_z,
  probs = CONFIG$capec_high_quantile,
  na.rm = TRUE
)

q_low <- quantile(
  sc_cap$CapEC_dysfunction_composite_score_z,
  probs = CONFIG$capec_low_quantile,
  na.rm = TRUE
)

sc_cap$CapEC_activation_state <- case_when(
  sc_cap$CapEC_dysfunction_composite_score_z >= q_high ~ "CapEC_ActHigh",
  sc_cap$CapEC_dysfunction_composite_score_z <= q_low  ~ "CapEC_ActLow",
  TRUE ~ "CapEC_ActMid"
)

sc_cap$CapEC_activation_state <- factor(
  sc_cap$CapEC_activation_state,
  levels = c("CapEC_ActLow", "CapEC_ActMid", "CapEC_ActHigh")
)

capec_cell_scores <- sc_cap@meta.data %>%
  mutate(cell = rownames(sc_cap@meta.data))

score_cols_raw <- paste0(names(capec_gene_sets), "_score")
score_cols_z <- paste0(names(capec_gene_sets), "_score_z")

write.csv(
  capec_cell_scores,
  file.path(CONFIG$out_tables, "node06_CapEC_ModuleScore_CellLevel.csv"),
  row.names = FALSE
)

qsave(sc_cap, file.path(CONFIG$out_rds, "node06_CapEC_with_ModuleScores.qs"))

################################################################################
# Module 7.1: CapEC sample-level summary and statistics
################################################################################

msg("CapEC sample-level statistics")

if (!sample_col %in% colnames(capec_cell_scores)) {
  capec_cell_scores$sample_id_auto <- sc_cap$orig.ident
  sample_col_cap <- "sample_id_auto"
} else {
  sample_col_cap <- sample_col
}

capec_sample_scores <- capec_cell_scores %>%
  mutate(
    Group = .data[[CONFIG$group_col]],
    sample_id = .data[[sample_col_cap]]
  ) %>%
  group_by(Group, sample_id) %>%
  summarise(
    n_capec = n(),
    across(
      all_of(c(
        score_cols_z,
        "CapEC_BBB_loss_score_z",
        "CapEC_dysfunction_composite_score_z"
      )),
      ~ mean(.x, na.rm = TRUE)
    ),
    frac_ActHigh = mean(CapEC_activation_state == "CapEC_ActHigh", na.rm = TRUE),
    frac_ActLow  = mean(CapEC_activation_state == "CapEC_ActLow", na.rm = TRUE),
    frac_ActMid  = mean(CapEC_activation_state == "CapEC_ActMid", na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  capec_sample_scores,
  file.path(CONFIG$out_tables, "node06_CapEC_ModuleScore_SampleLevel.csv"),
  row.names = FALSE
)

test_one_score <- function(df, score_col) {
  tmp <- df %>%
    select(Group, value = all_of(score_col)) %>%
    filter(!is.na(value), !is.na(Group))

  groups <- unique(as.character(tmp$Group))
  if (!all(CONFIG$group_levels %in% groups)) {
    return(data.frame(
      score = score_col,
      n_CN = sum(tmp$Group == "CN"),
      n_AD = sum(tmp$Group == "AD"),
      mean_CN = NA_real_,
      mean_AD = NA_real_,
      delta_AD_vs_CN = NA_real_,
      p_value = NA_real_,
      test = NA_character_
    ))
  }

  n_cn <- sum(tmp$Group == "CN")
  n_ad <- sum(tmp$Group == "AD")

  mean_cn <- mean(tmp$value[tmp$Group == "CN"], na.rm = TRUE)
  mean_ad <- mean(tmp$value[tmp$Group == "AD"], na.rm = TRUE)

  pval <- tryCatch(
    wilcox.test(value ~ Group, data = tmp)$p.value,
    error = function(e) NA_real_
  )

  data.frame(
    score = score_col,
    n_CN = n_cn,
    n_AD = n_ad,
    mean_CN = mean_cn,
    mean_AD = mean_ad,
    delta_AD_vs_CN = mean_ad - mean_cn,
    p_value = pval,
    test = "Wilcoxon rank-sum"
  )
}

score_cols_for_stats <- c(
  "CapEC_BBB_maintenance_score_z",
  "CapEC_Endothelial_activation_score_z",
  "CapEC_Inflammatory_score_z",
  "CapEC_IFN_response_score_z",
  "CapEC_ECM_remodeling_score_z",
  "CapEC_Permeability_angiogenic_score_z",
  "CapEC_Transport_score_z",
  "CapEC_BBB_loss_score_z",
  "CapEC_dysfunction_composite_score_z",
  "frac_ActHigh",
  "frac_ActLow",
  "frac_ActMid"
)

score_cols_for_stats <- intersect(score_cols_for_stats, colnames(capec_sample_scores))

capec_score_stats <- map_dfr(
  score_cols_for_stats,
  ~ test_one_score(capec_sample_scores, .x)
) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj_BH)

write.csv(
  capec_score_stats,
  file.path(CONFIG$out_tables, "node06_CapEC_ModuleScore_SampleLevel_Stats.csv"),
  row.names = FALSE
)

################################################################################
# Module 7.2: CapEC module score figures
################################################################################

msg("CapEC module score figures")

plot_score_labels <- c(
  CapEC_BBB_maintenance_score_z = "BBB maintenance",
  CapEC_Endothelial_activation_score_z = "Endothelial activation",
  CapEC_Inflammatory_score_z = "Inflammatory",
  CapEC_IFN_response_score_z = "IFN response",
  CapEC_ECM_remodeling_score_z = "ECM remodeling",
  CapEC_Permeability_angiogenic_score_z = "Permeability / angiogenic",
  CapEC_Transport_score_z = "Transport",
  CapEC_BBB_loss_score_z = "BBB loss",
  CapEC_dysfunction_composite_score_z = "Dysfunction composite"
)

sample_plot_cols <- intersect(names(plot_score_labels), colnames(capec_sample_scores))

sample_long <- capec_sample_scores %>%
  select(Group, sample_id, all_of(sample_plot_cols)) %>%
  pivot_longer(
    cols = all_of(sample_plot_cols),
    names_to = "score_name",
    values_to = "score_value"
  ) %>%
  mutate(
    score_label = factor(plot_score_labels[score_name], levels = plot_score_labels[sample_plot_cols])
  )

p_capec_sample <- ggplot(sample_long, aes(x = Group, y = score_value, fill = Group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.25, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.85) +
  facet_wrap(~ score_label, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("CN" = "#4DBBD5", "AD" = "#E64B35")) +
  theme_publication(base_size = 8) +
  labs(
    title = "CapEC functional module scores at sample level",
    x = NULL,
    y = "Mean z-score"
  )

sample_stat_df <- sample_long %>%
  group_by(score_label) %>%
  summarise(
    p_value = tryCatch(wilcox.test(score_value ~ Group)$p.value, error = function(e) NA_real_),
    ymax = max(score_value, na.rm = TRUE),
    ymin = min(score_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = vapply(p_value, p_value_label, character(1)),
    y = ymax + (ymax - ymin) * 0.08,
    x = 1.5
  )

p_capec_sample <- p_capec_sample +
  geom_text(
    data = sample_stat_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.4,
    vjust = 0
  )


safe_ggsave_all(
  file.path(CONFIG$out_plots, "Fig_node06_CapEC_ModuleScores_SampleLevel"),
  p_capec_sample,
  width = 8.2,
  height = 6.4
)

key_violin_scores <- c(
  "CapEC_BBB_maintenance_score_z",
  "CapEC_Endothelial_activation_score_z",
  "CapEC_ECM_remodeling_score_z",
  "CapEC_dysfunction_composite_score_z"
)

key_violin_scores <- intersect(key_violin_scores, colnames(capec_cell_scores))

cell_long <- capec_cell_scores %>%
  mutate(Group = .data[[CONFIG$group_col]]) %>%
  select(cell, Group, all_of(key_violin_scores)) %>%
  pivot_longer(
    cols = all_of(key_violin_scores),
    names_to = "score_name",
    values_to = "score_value"
  ) %>%
  mutate(
    score_label = factor(plot_score_labels[score_name], levels = plot_score_labels[key_violin_scores])
  )

p_capec_violin <- ggplot(cell_long, aes(x = Group, y = score_value, fill = Group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.6, linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.2) +
  facet_wrap(~ score_label, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("CN" = "#4DBBD5", "AD" = "#E64B35")) +
  theme_publication(base_size = 8) +
  labs(
    title = "CapEC functional module scores at cell level",
    x = NULL,
    y = "Module score z"
  )

cell_stat_df <- cell_long %>%
  group_by(score_label) %>%
  summarise(
    p_value = tryCatch(wilcox.test(score_value ~ Group)$p.value, error = function(e) NA_real_),
    ymax = max(score_value, na.rm = TRUE),
    ymin = min(score_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = vapply(p_value, p_value_label, character(1)),
    y = ymax + (ymax - ymin) * 0.05,
    x = 1.5
  )

p_capec_violin <- p_capec_violin +
  geom_text(
    data = cell_stat_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.2,
    vjust = 0
  )

safe_ggsave_all(
  file.path(CONFIG$out_plots, "Fig_node06_CapEC_ModuleScores_CellLevel_Violin"),
  p_capec_violin,
  width = 8.5,
  height = 3.2
)

capec_state_comp <- capec_cell_scores %>%
  mutate(
    Group = .data[[CONFIG$group_col]],
    sample_id = .data[[sample_col_cap]]
  ) %>%
  filter(!is.na(CapEC_activation_state)) %>%
  count(Group, sample_id, CapEC_activation_state, name = "n") %>%
  group_by(Group, sample_id) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

write.csv(
  capec_state_comp,
  file.path(CONFIG$out_tables, "node06_CapEC_ActivationState_Composition.csv"),
  row.names = FALSE
)

p_capec_state <- ggplot(capec_state_comp, aes(x = Group, y = frac, fill = CapEC_activation_state)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.25, alpha = 0.75) +
  geom_jitter(aes(color = CapEC_activation_state), width = 0.12, size = 1.1, alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~ CapEC_activation_state, nrow = 1) +
  scale_fill_manual(values = c(
    "CapEC_ActLow" = "#4DBBD5",
    "CapEC_ActMid" = "#999999",
    "CapEC_ActHigh" = "#E64B35"
  )) +
  scale_color_manual(values = c(
    "CapEC_ActLow" = "#4DBBD5",
    "CapEC_ActMid" = "#999999",
    "CapEC_ActHigh" = "#E64B35"
  )) +
  theme_publication(base_size = 8) +
  labs(
    title = "CapEC activation-state composition",
    x = NULL,
    y = "Fraction per sample",
    fill = NULL
  )

state_stat_df <- capec_state_comp %>%
  group_by(CapEC_activation_state) %>%
  summarise(
    p_value = tryCatch(wilcox.test(frac ~ Group)$p.value, error = function(e) NA_real_),
    ymax = max(frac, na.rm = TRUE),
    ymin = min(frac, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = vapply(p_value, p_value_label, character(1)),
    y = ymax + (ymax - ymin) * 0.10,
    x = 1.5
  )

p_capec_state <- p_capec_state +
  geom_text(
    data = state_stat_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.4,
    vjust = 0
  )

safe_ggsave_all(
  file.path(CONFIG$out_plots, "Fig_node06_CapEC_ActivationState_Composition"),
  p_capec_state,
  width = 7.5,
  height = 3.0
)

################################################################################
# Module 7.3: Pericyte functional module scoring
# Pericyte shows strongest proportion shift (+9% UC)
################################################################################

msg("Module 7.3: Pericyte functional module scoring")
peri_cells <- rownames(sc_endo@meta.data)[
  sc_endo@meta.data[[CONFIG$endo_sub_col]] %in% c("Pericyte")]
if (length(peri_cells) >= 10) {
  sc_peri <- subset(sc_endo, cells = peri_cells)
  DefaultAssay(sc_peri) <- CONFIG$assay

  pericyte_gene_sets <- list(
    Pericyte_Contractile = c("ACTA2","TAGLN","CNN1","MYH11","MYLK","ACTG2","DES","CALD1","TPM1","TPM2"),
    Pericyte_PDGFRB_signaling = c("PDGFRB","RGS5","NOTCH3","ANPEP","CSPG4","COL4A1","COL4A2","LAMA4","FN1"),
    Pericyte_BBB_support = c("CLDN5","OCLN","TJP1","CDH5","PECAM1","MFSD2A","SLC2A1","ABCB1","KDR","FLT1","ANGPT1"),
    Pericyte_ECM = c("COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","FN1","LAMA4","LAMB1","SPARC","TIMP1","TIMP2","MMP2","MMP14")
  )

  for (nm in names(pericyte_gene_sets)) {
    scores_nm <- score_module(sc_peri, pericyte_gene_sets[[nm]], nm)
    sc_peri <- AddMetaData(sc_peri, scores_nm, col.name = paste0(nm, "_score_z"))
  }
  
  peri_score_names <- paste0(names(pericyte_gene_sets), "_score_z")
  peri_meta <- sc_peri@meta.data
  peri_meta$Group <- factor(peri_meta$Group, levels = c("CN", "AD"))

  peri_sample <- peri_meta %>% group_by(orig.ident, Group) %>%
    summarise(across(all_of(peri_score_names), mean, na.rm = TRUE), .groups = "drop")

  peri_stats <- do.call(rbind, lapply(peri_score_names, function(sn) {
    w <- wilcox.test(peri_sample[[sn]][peri_sample$Group == "AD"],
                     peri_sample[[sn]][peri_sample$Group == "CN"])
    data.frame(Module = sn, CN_mean = mean(peri_sample[[sn]][peri_sample$Group == "CN"]),
               AD_mean = mean(peri_sample[[sn]][peri_sample$Group == "AD"]), p = w$p.value)
  }))
  peri_stats$p_adj <- p.adjust(peri_stats$p, method = "BH")
  peri_stats$Module <- gsub("_score_z", "", peri_stats$Module)
  write.csv(peri_stats, file.path(CONFIG$out_tables, "node06_Pericyte_ModuleScore_Stats.csv"), row.names = FALSE)
  msg(sprintf("  Pericyte module scores saved. nCells=%d", ncol(sc_peri)))

  peri_long <- peri_meta %>% select(all_of(c("Group", peri_score_names))) %>%
    pivot_longer(cols = all_of(peri_score_names), names_to = "Module", values_to = "Score")
  peri_long$Module <- gsub("_score_z", "", peri_long$Module)

  p_peri <- ggplot(peri_long, aes(Group, Score, fill = Group)) +
    geom_violin(trim = TRUE, alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.15, alpha = 0.5, outlier.shape = NA) +
    facet_wrap(~ Module, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c(CN = "#4DBBD5", AD = "#E64B35")) +
    labs(x = NULL, y = "Module Score (z)", title = "Pericyte functional modules: CN vs AD") +
    theme_bw(base_size = 8) + theme(legend.position = "none")
  safe_ggsave_all(file.path(CONFIG$out_plots, "Fig_node06_Pericyte_ModuleScore"), p_peri, width = 10, height = 4)
  msg("  Pericyte module figure saved.")
} else {
  msg("  Too few Pericyte cells - skipping module scoring.")
}

################################################################################
# Module 7.4: Arterial/SMC vascular dysfunction module scoring
# strongest MiloR signal (N29/N36 Arterial nhood loss in Disease)
################################################################################

msg("Module 7.4: Arterial/SMC vascular contractile module scoring")
art_smc_cells <- rownames(sc_endo@meta.data)[
  sc_endo@meta.data[[CONFIG$endo_sub_col]] %in% c("Arterial", "SMC")]
if (length(art_smc_cells) >= 10) {
  sc_asmc <- subset(sc_endo, cells = art_smc_cells)
  DefaultAssay(sc_asmc) <- CONFIG$assay

  vsmc_gene_sets <- list(
    VSMC_Contractile = c("ACTA2","TAGLN","CNN1","MYH11","MYLK","ACTG2","CALD1","TPM1","TPM2","MYOCD","SRF"),
    VSMC_Arterial_specification = c("GJA5","GJA4","EFNB2","BMX","HEY1","HEY2","DLL4","JAG1","NOTCH1","NOTCH4","EPHB4","NRP1"),
    VSMC_Vascular_tone = c("ACE","AGT","AGTR1","NOS3","EDN1","EDNRA","PTGIS","PTGIR","KCNMA1","CACNA1C","ATP2A2"),
    VSMC_ECM = c("COL1A1","COL1A2","COL3A1","COL4A1","ELN","FBN1","FN1","SPARC","LOX","LOXL1","TIMP1","MMP2")
  )

  for (nm in names(vsmc_gene_sets)) {
    scores_nm <- score_module(sc_asmc, vsmc_gene_sets[[nm]], nm)
    sc_asmc <- AddMetaData(sc_asmc, scores_nm, col.name = paste0(nm, "_score_z"))
  }
  
  vsmc_score_names <- paste0(names(vsmc_gene_sets), "_score_z")
  vsmc_meta <- sc_asmc@meta.data
  vsmc_meta$Group <- factor(vsmc_meta$Group, levels = c("CN", "AD"))

  vsmc_sample <- vsmc_meta %>% group_by(orig.ident, Group) %>%
    summarise(across(all_of(vsmc_score_names), mean, na.rm = TRUE), .groups = "drop")

  vsmc_stats <- do.call(rbind, lapply(vsmc_score_names, function(sn) {
    w <- wilcox.test(vsmc_sample[[sn]][vsmc_sample$Group == "AD"],
                     vsmc_sample[[sn]][vsmc_sample$Group == "CN"])
    data.frame(Module = sn, CN_mean = mean(vsmc_sample[[sn]][vsmc_sample$Group == "CN"]),
               AD_mean = mean(vsmc_sample[[sn]][vsmc_sample$Group == "AD"]), p = w$p.value)
  }))
  vsmc_stats$p_adj <- p.adjust(vsmc_stats$p, method = "BH")
  vsmc_stats$Module <- gsub("_score_z", "", vsmc_stats$Module)
  write.csv(vsmc_stats, file.path(CONFIG$out_tables, "node06_VSMC_ModuleScore_Stats.csv"), row.names = FALSE)
  msg(sprintf("  VSMC module scores saved. nCells=%d", ncol(sc_asmc)))

  vsmc_long <- vsmc_meta %>% select(all_of(c("Group", vsmc_score_names))) %>%
    pivot_longer(cols = all_of(vsmc_score_names), names_to = "Module", values_to = "Score")
  vsmc_long$Module <- gsub("_score_z", "", vsmc_long$Module)

  p_vsmc <- ggplot(vsmc_long, aes(Group, Score, fill = Group)) +
    geom_violin(trim = TRUE, alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.15, alpha = 0.5, outlier.shape = NA) +
    facet_wrap(~ Module, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c(CN = "#4DBBD5", AD = "#E64B35")) +
    labs(x = NULL, y = "Module Score (z)", title = "Arterial/SMC vascular modules: CN vs AD") +
    theme_bw(base_size = 8) + theme(legend.position = "none")
  safe_ggsave_all(file.path(CONFIG$out_plots, "Fig_node06_VSMC_ModuleScore"), p_vsmc, width = 10, height = 4)
  msg("  VSMC module figure saved.")
} else {
  msg("  Too few Arterial/SMC cells - skipping module scoring.")
}

################################################################################
# Module 8: Save CapEC scores and state labels (lightweight, no OOM)
################################################################################

msg("=====================================================================")
msg("Module 8: Save CapEC scores and state labels (lightweight)")
msg("=====================================================================")

transfer_cols <- c(
  score_cols_raw,
  score_cols_z,
  "CapEC_BBB_loss_score_z",
  "CapEC_dysfunction_composite_score",
  "CapEC_dysfunction_composite_score_z",
  "CapEC_activation_state"
)
transfer_cols <- intersect(transfer_cols, colnames(sc_cap@meta.data))

capec_state_lookup <- sc_cap@meta.data[["CapEC_activation_state"]]
names(capec_state_lookup) <- rownames(sc_cap@meta.data)
capec_state_lookup <- capec_state_lookup[!is.na(capec_state_lookup) & capec_state_lookup != ""]
cat(sprintf("[Module 8] CapEC_activation_state lookup: %d cells\n", length(capec_state_lookup)))

common_endo_cells <- intersect(colnames(sc_endo), colnames(sc_cap))
for (cc in transfer_cols) {
  sc_endo@meta.data[common_endo_cells, cc] <- sc_cap@meta.data[common_endo_cells, cc]
}
sc_endo$liana_label_state <- NA_character_
is_cap_endo <- colnames(sc_endo) %in% names(capec_state_lookup)
sc_endo$liana_label_state[is_cap_endo] <- as.character(capec_state_lookup[colnames(sc_endo)[is_cap_endo]])

qsave(sc_endo, file.path(CONFIG$out_rds, "node06_scRNA_endo_with_CapEC_scores.qs"))
cat("[Module 8] sc_endo with CapEC saved (1,576 cells)\n")

qsave(capec_state_lookup, file.path(CONFIG$out_rds, "node06_CapEC_activation_state_lookup.qs"))
cat("[Module 8] CapEC state lookup saved\n")

state_df <- data.frame(
  cell = names(capec_state_lookup),
  state = as.character(capec_state_lookup),
  stringsAsFactors = FALSE
)
write.csv(state_df, file.path(CONFIG$out_tables, "node06_CapEC_activation_state_lookup.csv"), row.names = FALSE)

rm(sc_full); gc()
cat("processing...\n")

################################################################################
# Module 9: CapEC activation state LIANA (SLOW ~15-30 min)
################################################################################

msg("=====================================================================")
msg("Module 9: CapEC activation-state LIANA")
msg("=====================================================================")

CKPT_STATE_CN <- file.path(CONFIG$out_rds, "node06_State_LIANA_CN_subtyped.qs")
CKPT_STATE_AD <- file.path(CONFIG$out_rds, "node06_State_LIANA_AD_subtyped.qs")

if (file.exists(CKPT_STATE_CN) && file.exists(CKPT_STATE_AD)) {
  msg("CapEC state LIANA checkpoints found - loading cached results")
  liana_state_CN <- qread(CKPT_STATE_CN)
  liana_state_AD <- qread(CKPT_STATE_AD)
} else {
  cat("processing...\n")
  sc_full_state <- qread("/path/to/project/results/02_annotation/rds/scRNA_annotated.qs")
  cat(sprintf("[Module 9] sc_full_endofiltered: %d cells loaded\n", ncol(sc_full_state)))

  sc_full_state <- harmonize_group_labels(sc_full_state)
  sc_full_state <- merge_subtype_to_full(sc_full_state, sc_endo,  CONFIG$endo_sub_col)
  sc_full_state <- merge_subtype_to_full(sc_full_state, sc_astro, CONFIG$astro_sub_col)
  sc_full_state <- merge_subtype_to_full(sc_full_state, sc_micro, CONFIG$micro_sub_col)
  meta_s <- sc_full_state@meta.data
  meta_s$liana_label <- sanitize_label(meta_s[[CONFIG$major_col]])
  if (CONFIG$endo_sub_col %in% colnames(meta_s)) {
    is_endo <- !is.na(meta_s[[CONFIG$endo_sub_col]]) & meta_s[[CONFIG$endo_sub_col]] != ""
    meta_s$liana_label[is_endo] <- paste0("Cerebrovascular_", sanitize_label(meta_s[[CONFIG$endo_sub_col]][is_endo]))
  }
  if (CONFIG$astro_sub_col %in% colnames(meta_s)) {
    is_astro <- !is.na(meta_s[[CONFIG$astro_sub_col]]) & meta_s[[CONFIG$astro_sub_col]] != ""
    meta_s$liana_label[is_astro] <- paste0("Astro_", sanitize_label(meta_s[[CONFIG$astro_sub_col]][is_astro]))
  }
  if (CONFIG$micro_sub_col %in% colnames(meta_s)) {
    is_micro <- !is.na(meta_s[[CONFIG$micro_sub_col]]) & meta_s[[CONFIG$micro_sub_col]] != ""
    meta_s$liana_label[is_micro] <- paste0("Micro_", sanitize_label(meta_s[[CONFIG$micro_sub_col]][is_micro]))
  }

  capec_lookup <- qread(file.path(CONFIG$out_rds, "node06_CapEC_activation_state_lookup.qs"))
  meta_s$liana_label_state <- meta_s$liana_label
  has_state <- rownames(meta_s) %in% names(capec_lookup)
  meta_s$liana_label_state[has_state] <- as.character(capec_lookup[rownames(meta_s)[has_state]])

  sc_full_state@meta.data[["liana_label_state"]] <- meta_s$liana_label_state
  cat("processing...\n")
  rm(meta_s); gc()

  cat("processing...\n")
  print(sort(table(sc_full_state@meta.data[["liana_label_state"]], useNA="ifany"), decreasing=TRUE)[1:10])

  cat("processing...\n")
  liana_state_CN <- run_liana_one_group(sc_full_state, "CN", label_col = "liana_label_state", prefix = "node06_State")
  cat("processing...\n")
  liana_state_AD <- run_liana_one_group(sc_full_state, "AD", label_col = "liana_label_state", prefix = "node06_State")

  qsave(liana_state_CN, CKPT_STATE_CN)
  qsave(liana_state_AD, CKPT_STATE_AD)
  cat("processing...\n")

  rm(sc_full_state); gc()
  cat("processing...\n")
}

qsave(liana_state_CN, file.path(CONFIG$out_rds, "node06_State_LIANA_CN.qs"))
qsave(liana_state_AD, file.path(CONFIG$out_rds, "node06_State_LIANA_AD.qs"))
qsave(liana_state_CN, file.path(CONFIG$out_rds, "node06_State_LIANA_CN_subtyped.qs"))
qsave(liana_state_AD, file.path(CONFIG$out_rds, "node06_State_LIANA_AD_subtyped.qs"))

TARGET_FOCUS <- c(
  "Micro_DAM",
  "Micro_Inflammatory",
  "Astro_Reactive",
  "Astro_Intermediate",
  "Excitatory",
  "Inhibitory",
  "Oligodendrocytes",
  "OPCs"
)

extract_capec_state_lr <- function(df, group_value) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  tmp <- df %>%
    filter(source %in% c("CapEC_ActHigh", "CapEC_ActLow", "CapEC_ActMid")) %>%
    mutate(group = group_value)

  if (nrow(tmp) == 0) return(NULL)

  tmp
}

capec_state_lr_all <- bind_rows(
  extract_capec_state_lr(liana_state_CN$raw_df, "CN"),
  extract_capec_state_lr(liana_state_AD$raw_df, "AD")
)

write.csv(
  capec_state_lr_all,
  file.path(CONFIG$out_tables, "node06_CapEC_State_LIANA_all.csv"),
  row.names = FALSE
)

capec_state_lr_focus <- capec_state_lr_all %>%
  filter(target %in% TARGET_FOCUS | str_detect(target, "Astro|Micro|Excitatory|Inhibitory|OPC|Oligo"))

write.csv(
  capec_state_lr_focus,
  file.path(CONFIG$out_tables, "node06_CapEC_State_LIANA_TargetFocus.csv"),
  row.names = FALSE
)

################################################################################
# Module 9.1: CapEC ActHigh vs ActLow differential LR
################################################################################

msg("CapEC ActHigh vs ActLow LIANA comparison")

compare_high_low_one_group <- function(df, group_value) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  high <- df %>%
    filter(source == "CapEC_ActHigh") %>%
    select(target, ligand, receptor, strength_high = communication_strength)

  low <- df %>%
    filter(source == "CapEC_ActLow") %>%
    select(target, ligand, receptor, strength_low = communication_strength)

  out <- full_join(high, low, by = c("target", "ligand", "receptor")) %>%
    mutate(
      strength_high = replace_na(strength_high, 0),
      strength_low = replace_na(strength_low, 0),
      delta_high_vs_low = strength_high - strength_low,
      group = group_value,
      direction = case_when(
        delta_high_vs_low > 0 ~ "ActHigh_up",
        delta_high_vs_low < 0 ~ "ActHigh_down",
        TRUE ~ "unchanged"
      )
    ) %>%
    arrange(desc(abs(delta_high_vs_low)))

  out
}

state_diff_CN <- compare_high_low_one_group(liana_state_CN$raw_df, "CN")
state_diff_AD <- compare_high_low_one_group(liana_state_AD$raw_df, "AD")

state_diff_all <- bind_rows(state_diff_CN, state_diff_AD)

write.csv(
  state_diff_all,
  file.path(CONFIG$out_tables, "node06_CapEC_ActHigh_vs_ActLow_LIANA.csv"),
  row.names = FALSE
)

################################################################################
# Module 9.2: State LIANA figures
################################################################################

p_state_heat_uc <- NULL
p_curated <- NULL

if (!is.null(capec_state_lr_focus) && nrow(capec_state_lr_focus) > 0) {
  heat_df <- capec_state_lr_focus %>%
    filter(group == "AD") %>%
    group_by(source, target) %>%
    summarise(mean_strength = mean(communication_strength, na.rm = TRUE), .groups = "drop")

  if (nrow(heat_df) > 0) {
    p_state_heat_uc <- ggplot(heat_df, aes(x = target, y = source, fill = mean_strength)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_gradient(low = "grey95", high = "#E64B35") +
      theme_publication(base_size = 8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "AD capillary state-to-NVU communication",
        x = NULL,
        y = NULL,
        fill = "Strength"
      )

    safe_ggsave_all(
      file.path(CONFIG$out_plots, "Fig_node06_CapEC_State_LIANA_AD_Heatmap"),
      p_state_heat_uc,
      width = 7.5,
      height = 3.2
    )
  }
}

curated_lr_keywords <- c(
  "VEGFA", "KDR", "FLT1",
  "ANGPT2", "TEK",
  "ICAM1", "VCAM1",
  "CCL2", "CCR2",
  "CXCL10", "CXCR3",
  "TGFB1", "TGFBR",
  "JAG1", "NOTCH",
  "SEMA", "PLX",
  "COL4A1", "ITGA", "ITGB",
  "FN1"
)

if (!is.null(state_diff_AD) && nrow(state_diff_AD) > 0) {
  curated_df <- state_diff_AD %>%
    mutate(lr_pair = paste(ligand, receptor, sep = " - ")) %>%
    filter(
      str_detect(ligand, paste(curated_lr_keywords, collapse = "|")) |
        str_detect(receptor, paste(curated_lr_keywords, collapse = "|"))
    ) %>%
    arrange(desc(delta_high_vs_low)) %>%
    slice_head(n = CONFIG$top_n_curated_lr)

  write.csv(
    curated_df,
    file.path(CONFIG$out_tables, "node06_CapEC_ActHigh_Curated_LR_Candidates_AD.csv"),
    row.names = FALSE
  )

  if (nrow(curated_df) > 0) {
    curated_df <- curated_df %>%
      mutate(lr_pair = factor(lr_pair, levels = rev(unique(lr_pair))))

    p_curated <- ggplot(curated_df, aes(x = delta_high_vs_low, y = lr_pair, fill = target)) +
      geom_col(width = 0.75) +
      theme_publication(base_size = 8) +
      labs(
        title = "Activated capillary communication candidates in AD",
        x = "ActHigh - ActLow communication strength",
        y = NULL,
        fill = "Target"
      )

    safe_ggsave_all(
      file.path(CONFIG$out_plots, "Fig_node06_CapEC_ActHigh_Curated_LR_Candidates_AD"),
      p_curated,
      width = 7.2,
      height = 5.5
    )
  }
}

################################################################################
# Module 10: Microglia immune contamination QC
################################################################################

msg("=====================================================================")
msg("Module 10: Microglia immune contamination QC")
msg("=====================================================================")

immune_qc_genes <- c(
  "PTPRC", "LST1", "TYROBP", "AIF1", "CSF1R",
  "CD3D", "CD3E", "CD79A", "MS4A1", "NKG7",
  "LYZ", "S100A8", "S100A9", "FCGR3A"
)

micro_mat <- get_data_matrix(sc_micro, CONFIG$assay)
immune_present <- intersect(immune_qc_genes, rownames(micro_mat))

if (length(immune_present) > 0) {
  qc_score <- Matrix::colMeans(micro_mat[immune_present, , drop = FALSE])

  qc_df <- sc_micro@meta.data %>%
    mutate(
      cell = rownames(sc_micro@meta.data),
      immune_contamination_score = as.numeric(qc_score)
    )

  write.csv(
    qc_df,
    file.path(CONFIG$out_tables, "node06_Micro_ImmuneContamination_QC.csv"),
    row.names = FALSE
  )
}

################################################################################
# Module 11: SCTour export
################################################################################

msg("=====================================================================")
msg("Module 11: SCTour export")
msg("=====================================================================")

export_sctour_object <- function(obj, name, subtype_col = NULL) {
  out_dir <- file.path(CONFIG$out_sctour, name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  DefaultAssay(obj) <- CONFIG$assay

  counts <- tryCatch(
    GetAssayData(obj, assay = CONFIG$assay, layer = "counts"),
    error = function(e) NULL
  )

  data <- tryCatch(
    GetAssayData(obj, assay = CONFIG$assay, layer = "data"),
    error = function(e) NULL
  )

  if (is.null(counts)) counts <- data
  if (is.null(data)) data <- counts

  Matrix::writeMM(counts, file.path(out_dir, paste0(name, "_counts.mtx")))
  write.table(
    rownames(counts),
    file.path(out_dir, paste0(name, "_genes.tsv")),
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  write.table(
    colnames(counts),
    file.path(out_dir, paste0(name, "_barcodes.tsv")),
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  meta_out <- obj@meta.data %>%
    mutate(cell = rownames(obj@meta.data))

  if (!is.null(subtype_col) && subtype_col %in% colnames(meta_out)) {
    meta_out$sctour_state <- meta_out[[subtype_col]]
  }

  write.csv(
    meta_out,
    file.path(out_dir, paste0(name, "_metadata.csv")),
    row.names = FALSE
  )

  if ("umap" %in% names(obj@reductions)) {
    emb <- Embeddings(obj, "umap") %>%
      as.data.frame() %>%
      mutate(cell = rownames(.))
    write.csv(
      emb,
      file.path(out_dir, paste0(name, "_umap.csv")),
      row.names = FALSE
    )
  }

  invisible(TRUE)
}

if (CONFIG$run_sctour_export) {
  export_sctour_object(sc_endo,  "Cerebrovascular",  subtype_col = CONFIG$endo_sub_col)
  export_sctour_object(sc_astro, "Astro", subtype_col = CONFIG$astro_sub_col)
  export_sctour_object(sc_micro, "Micro", subtype_col = CONFIG$micro_sub_col)

  if (!is.null(diff_ccc)) {
    file.copy(
      file.path(CONFIG$out_tables, "node06_Differential_CCC_Cerebrovascular_to_NVU_AD_vs_CN.csv"),
      file.path(CONFIG$out_sctour, "Differential_CCC.csv"),
      overwrite = TRUE
    )
  }

  if (file.exists(file.path(CONFIG$out_tables, "node06_GO_BP_enrichment_all.csv"))) {
    file.copy(
      file.path(CONFIG$out_tables, "node06_GO_BP_enrichment_all.csv"),
      file.path(CONFIG$out_sctour, "GO_DEG_directional.csv"),
      overwrite = TRUE
    )
  }

  if (file.exists(file.path(CONFIG$out_tables, "node06_CapEC_ModuleScore_CellLevel.csv"))) {
    file.copy(
      file.path(CONFIG$out_tables, "node06_CapEC_ModuleScore_CellLevel.csv"),
      file.path(CONFIG$out_sctour, "CapEC_ModuleScore_CellLevel.csv"),
      overwrite = TRUE
    )
  }

  if (file.exists(file.path(CONFIG$out_tables, "node06_CapEC_ActHigh_vs_ActLow_LIANA.csv"))) {
    file.copy(
      file.path(CONFIG$out_tables, "node06_CapEC_ActHigh_vs_ActLow_LIANA.csv"),
      file.path(CONFIG$out_sctour, "CapEC_ActHigh_vs_ActLow_LIANA.csv"),
      overwrite = TRUE
    )
  }
}

################################################################################
# Module 12: Assemble publication-style summary figure
################################################################################

msg("=====================================================================")
msg("Module 12: Assemble summary figure")
msg("=====================================================================")

summary_plots <- list()

summary_plots[["A"]] <- p_capec_sample +
  labs(title = "A  CapEC functional programs")

summary_plots[["B"]] <- p_capec_state +
  labs(title = "B  CapEC activation-state composition")

if (!is.null(p_state_heat_uc)) {
  summary_plots[["C"]] <- p_state_heat_uc +
    labs(title = "C  CapEC state-to-NVU communication")
}

if (!is.null(p_curated)) {
  summary_plots[["D"]] <- p_curated +
    labs(title = "D  Activated CapEC communication candidates")
}

if (length(summary_plots) >= 3) {
  p_summary <- wrap_plots(
    summary_plots,
    ncol = 1,
    heights = c(1.3, 0.8, 0.9, 1.3)[seq_along(summary_plots)]
  )

  safe_ggsave_all(
    file.path(CONFIG$out_plots, "Fig_node06_Main_CapEC_State_LIANA_Summary"),
    p_summary,
    width = 9.0,
    height = 14.0
  )
}

msg("=====================================================================")
msg("Module 13: Final summary")
msg("=====================================================================")

output_summary <- data.frame(
  category = c(
    "tables",
    "plots",
    "supplementary",
    "rds_node06",
    "sctour",
    "logs"
  ),
  path = c(
    CONFIG$out_tables,
    CONFIG$out_plots,
    CONFIG$out_supp,
    CONFIG$out_rds,
    CONFIG$out_sctour,
    CONFIG$out_logs
  )
)

write.csv(
  output_summary,
  file.path(CONFIG$out_tables, "node06_FINAL_output_summary.csv"),
  row.names = FALSE
)

msg("Output summary:")
print(output_summary)

################################################################################
# Module 14: Session info
################################################################################

msg("=====================================================================")
msg("Module 14: Session info")
msg("=====================================================================")

sink(file.path(CONFIG$out_logs, "node06_FINAL_sessionInfo.txt"))
cat("node06 FINAL GO / Reactome / LIANA / CapEC / SCTour analysis\n")
cat("Date: ", as.character(Sys.time()), "\n\n")
print(sessionInfo())
sink()

msg("Analysis finished successfully.")

# Set final checkpoint flag for quality-control review step
checkpoint_flag_done("node06_FULL_COMPLETE")

sink()

source(file.path(Sys.getenv("AD_NVU_REPO_ROOT", "."), "analysis", "06_communication", "06a_communication_figure_helper.R"))

# ========== Part2: publication-grade figures ==========
# File: scripts_final/node06_Figures_publicationGrade_v2.R
if (FALSE) {


suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(viridis)
  library(pheatmap)
  library(RColorBrewer)
  library(cowplot)
  library(tibble)
})

PROJECT_DIR <- "/path/to/R_libraries"
OUT_DIR <- file.path(OUT, "main_figures_v2")
OUT_TABLES <- file.path(OUT, "tables")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(20260517)

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), sprintf(...), "\n")

group_colors <- c("NC" = "#0072B2", "UC" = "#D55E00")

wilcox_pval <- function(x, y) {
  if (length(x) < 3 || length(y) < 3) return(NA)
  wilcox.test(x, y, exact = FALSE)$p.value
}
format_pval <- function(p) {
  if (is.na(p) || p > 0.05) return("n.s.")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("n.s.")
}
format_pval_full <- function(p) {
  if (is.na(p) || p > 0.05) return(sprintf("P=%.2f (n.s.)", p))
  if (p < 0.001) return("P<0.001 ***")
  sprintf("P=%.3f %s", p, format_pval(p))
}

clean_label <- function(x) {
  complex_patterns <- c(
    "ACVR2A_BMPR1B", "ACVR2B_BMPR1B", "BMPR1A_BMPR2", "BMPR1B_BMPR2",
    "ITGA1_ITGB1", "ITGA5_ITGB1", "ITGA11_ITGB1", "ITGAV_ITGB1",
    "DLL4_NOTCH4", "HLA.DQA1_LAG3", "FLT1_KDR", "FLT1_FLT4",
    "VEGFA_KDR", "VEGFB_KDR", "VEGFC_KDR", "VEGFC_FLT4",
    "COL4A1_COL4A2", "COL4A2_COL4A1"
  )
  for (pat in complex_patterns) {
    replacement <- gsub("_", "-", pat)
    x <- gsub(pat, replacement, x, fixed = FALSE)
  }

  x <- gsub("_", " ", x)

  x <- gsub("reactive DAA", "Reactive DAA", x, ignore.case = TRUE)
  x <- gsub("homeostatic", "Homeostatic", x, ignore.case = TRUE)
  x <- gsub("intermediate", "Intermediate", x, ignore.case = TRUE)
  x <- gsub("inflammatory", "Inflammatory", x, ignore.case = TRUE)
  x <- gsub("capillary", "Capillary", x, ignore.case = TRUE)
  x <- gsub("arterial", "Arterial", x, ignore.case = TRUE)
  x <- gsub("venous", "Venous", x, ignore.case = TRUE)
  x <- gsub("\\bNC\\b", "NC", x)
  x <- gsub("\\bUC\\b", "UC", x)
  x <- gsub("\\.\\.", ".", x)
  trimws(x)
}

trunc_label <- function(x, max_len = 65) {
  ifelse(nchar(x) > max_len, paste0(substr(x, 1, max_len), "..."), x)
}

# ============================================================
# ============================================================
msg("========== 1. GO DotPlot (fix labels) ==========")

go_enrich <- read.csv(file.path(OUT_TABLES, "node06_GO_BP_enrichment_all.csv"))

fix_dotplot <- function(df, title_prefix, out_prefix, database = "GO") {
  df <- df %>% filter(p.adjust < 0.05)
  if (nrow(df) == 0) { msg("  skip: no sig"); return() }

  df_top <- df %>% group_by(subtype) %>% arrange(p.adjust) %>% slice_head(n = 3) %>% ungroup()
  df_top$subtype_clean <- clean_label(df_top$subtype)
  df_top$Description_clean <- df_top$Description
  df_top <- df_top %>% arrange(Description_clean)
  df_top$Description_clean <- factor(df_top$Description_clean, levels = rev(unique(df_top$Description_clean)))
  df_top$neglogFDR <- -log10(df_top$p.adjust + 1e-10)
  n_sub <- length(unique(df_top$subtype_clean))

  p <- ggplot(df_top, aes(x = subtype_clean, y = Description_clean, size = Count, color = neglogFDR)) +
    geom_point(alpha = 0.9) +
    scale_radius(range = c(2, 7), name = "Gene Count") +
    scale_color_viridis(option = "plasma", name = bquote(-log[10](FDR)),
                        breaks = c(2, 5, 10, 15), labels = c("2", "5", "10", "15")) +
    facet_wrap(~subtype_clean, scales = "free_x", ncol = min(3, n_sub)) +
    labs(title = sprintf("%s: %s enrichment", title_prefix, database),
         x = "", y = paste0(database, " Biological Process")) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(8, 8, 8, 80)
    )

  fig_h <- max(6, nrow(df_top) * 0.4)
  fig_w <- min(12, 4 * n_sub)
  ggsave(file.path(OUT_DIR, sprintf("Fig_node06_%s_DotPlot_%s.png", database, out_prefix)),
         p, width = fig_w, height = fig_h, dpi = 300)
  ggsave(file.path(OUT_DIR, sprintf("Fig_node06_%s_DotPlot_%s.pdf", database, out_prefix)),
         p, width = fig_w, height = fig_h)
  msg(sprintf("  %s DotPlot (%s) done [%d pathways, h=%.1f, w=%.1f]", database, out_prefix, nrow(df_top), fig_h, fig_w))
}

tryCatch({
  fix_dotplot(go_enrich, "Endo", "Endo", "GO")
  fix_dotplot(go_enrich, "Astro", "Astro", "GO")
  fix_dotplot(go_enrich, "Micro", "Micro", "GO")
}, error = function(e) msg("GO error: %s", e$message))

# ============================================================
# ============================================================
msg("========== 2. CapEC ModuleScores (exact p values) ==========")

sample_level <- read.csv(file.path(OUT_TABLES, "node06_CapEC_ModuleScore_SampleLevel.csv"))

modules <- c(
  "CapEC_BBB_maintenance_score_z" = "BBB Maintenance",
  "CapEC_Endothelial_activation_score_z" = "Endothelial Activation",
  "CapEC_IFN_response_score_z" = "IFN Response",
  "CapEC_ECM_remodeling_score_z" = "ECM Remodeling",
  "CapEC_Inflammatory_score_z" = "Inflammatory",
  "CapEC_Permeability_angiogenic_score_z" = "Permeability"
)

stat_results <- data.frame(Module=character(), Pval=numeric(), Pval_full=character(),
                           NC_mean=numeric(), UC_mean=numeric(), NC_n=integer(), UC_n=integer())
p_list <- list(); i <- 1

for (mod_col in names(modules)) {
  if (!(mod_col %in% colnames(sample_level))) next
  df <- sample_level[, c("sample_id", "Group", mod_col)]
  colnames(df) <- c("sample_id", "Group", "score")
  df <- na.omit(df)
  nc_vals <- df$score[df$Group == "NC"]
  uc_vals <- df$score[df$Group == "UC"]
  p <- wilcox_pval(nc_vals, uc_vals)
  pfull <- format_pval_full(p)
  stat_results <- rbind(stat_results, data.frame(
    Module=modules[mod_col], Pval=p, Pval_full=pfull,
    NC_mean=mean(nc_vals), UC_mean=mean(uc_vals),
    NC_n=length(nc_vals), UC_n=length(uc_vals)))

  y_max <- max(df$score) * 1.15
  y_bracket <- max(df$score) * 1.05
  pstar <- format_pval(p)

  p1 <- ggplot(df, aes(x = Group, y = score, fill = Group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.18, size = 2.5, alpha = 0.7, color = "black") +
    scale_fill_manual(values = group_colors) +
    annotate("segment", x = 1, xend = 2, y = y_bracket, yend = y_bracket, color = "black", linewidth = 0.7) +
    annotate("text", x = 1.5, y = y_bracket * 1.04,
             label = pfull, size = 3.2, fontface = "bold") +
    labs(title = modules[mod_col], y = "Mean z-score",
         subtitle = sprintf("NC (n=%d) vs UC (n=%d)", length(nc_vals), length(uc_vals))) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 8, color = "gray40"),
          legend.position = "none",
          axis.text = element_text(size = 9),
          plot.margin = margin(5, 5, 5, 5)) +
    ylim(NA, y_max)

  p_list[[i]] <- p1; i <- i + 1
}

write.csv(stat_results, file.path(OUT_TABLES, "node06_CapEC_ModuleScore_Statistics_v2.csv"), row.names=FALSE)

combined <- plot_grid(plotlist = p_list, ncol = 3, nrow = 2,
  labels = LETTERS[1:length(p_list)], label_size = 11, rel_heights = rep(1, length(p_list))) +
  labs(title = "CapEC Functional Module Scores: NC vs UC (Sample Level)") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))

ggsave(file.path(OUT_DIR, "Fig_node06_CapEC_ModuleScores_SampleLevel_stat_v3.png"),
       combined, width = 13, height = 8, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig_node06_CapEC_ModuleScores_SampleLevel_stat_v3.pdf"),
       combined, width = 13, height = 8)
msg("  SampleLevel done with exact P values")

# ============================================================
# ============================================================
msg("========== 3. LIANA Heatmap (fix gray: use raw score) ==========")

liana_uc <- read.csv(file.path(OUT_TABLES, "node06_State_LIANA_Endo_to_NVU_UC.csv"))
liana_uc$source <- clean_label(liana_uc$source)
liana_uc$target <- clean_label(liana_uc$target)

liana_use <- liana_uc %>%
  arrange(aggregate_rank) %>%
  slice_head(n = 35) %>%
  mutate(lr_pair = clean_label(paste0(ligand, "-", receptor)),
         st = paste0(source, " -> ", target)) %>%
  filter(!is.na(communication_strength), is.finite(communication_strength))

mat <- liana_use %>% select(st, lr_pair, communication_strength) %>% distinct() %>%
  pivot_wider(names_from = lr_pair, values_from = communication_strength, values_fn = mean) %>%
  column_to_rownames("st") %>% as.matrix()
mat[is.na(mat)] <- 0

mat_log <- log1p(mat)
mat_z <- t(scale(t(mat_log)))
mat_z[is.na(mat_z)] <- 0

pdf(file.path(OUT_DIR, "Fig_node06_CapEC_State_LIANA_UC_Heatmap_v4.pdf"), width = 13, height = 8)
pheatmap(mat_z,
  cluster_rows = TRUE, cluster_cols = TRUE,
  clustering_method = "ward.D2",
  color = hcl.colors(100, "Purple-Orange", rev = TRUE),
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize_row = 8, fontsize_col = 8,
  legend_breaks = c(-2,-1,0,1,2), legend_labels = c("-2","-1","0","1","2"),
  main = sprintf("CapEC Activation State LIANA: Top %d LR pairs (UC group)\nZ-score normalized log(1+comm. strength)", ncol(mat_z))
)
dev.off()

png(file.path(OUT_DIR, "Fig_node06_CapEC_State_LIANA_UC_Heatmap_v4.png"),
    width = 1300, height = 800, res = 150)
pheatmap(mat_z,
  cluster_rows = TRUE, cluster_cols = TRUE,
  clustering_method = "ward.D2",
  color = hcl.colors(100, "Purple-Orange", rev = TRUE),
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize_row = 8, fontsize_col = 8,
  legend_breaks = c(-2,-1,0,1,2), legend_labels = c("-2","-1","0","1","2"),
  main = sprintf("CapEC Activation State LIANA: Top %d LR pairs (UC group)\nZ-score normalized log(1+comm. strength)", ncol(mat_z))
)
dev.off()
msg("  LIANA Heatmap v4 done (z-score on log1p, v1 logic restored)")

# ============================================================
# ============================================================
msg("========== 4. Differential CCC (clean labels + color) ==========")

diff <- read.csv(file.path(OUT_TABLES, "node06_Differential_CCC_Endo_to_NVU_UC_vs_NC.csv"))
if (nrow(diff) > 0) {
  diff$source <- clean_label(diff$source)
  diff$target <- clean_label(diff$target)
  diff$LR <- clean_label(paste0(diff$ligand, "-", diff$receptor))
  
  source_colors <- c(
    "Capillary" = "#E69F00", "Arterial" = "#56B4E9",
    "Venous" = "#009E73", "EC" = "#F0E442"
  )
  
  diff_top <- diff %>% arrange(desc(abs(delta_UC_vs_NC))) %>% slice_head(n = 20) %>%
    mutate(direction = ifelse(delta_UC_vs_NC > 0, "UC_enhanced", "NC_enhanced"))
  
  p_diff <- ggplot(diff_top,
                  aes(x = reorder(LR, delta_UC_vs_NC),
                      y = delta_UC_vs_NC, fill = source)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = source_colors, name = "Source Cell") +
    labs(title = sprintf("Top %d Differential LR pairs: UC vs NC (EC to NVU)", length(unique(diff_top$LR))),
         subtitle = "Bar color = source cell type | Position = UC (right) vs NC (left)",
         x = "Ligand-Receptor Pair", y = bquote(Delta~"(UC - NC)")) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          legend.position = "right", legend.text = element_text(size = 9),
          plot.margin = margin(8, 8, 8, 100)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")
  
  ggsave(file.path(OUT_DIR, "Fig_node06_Differential_CCC_Top20_v4.png"),
         p_diff, width = 11, height = 7, dpi = 300)
  ggsave(file.path(OUT_DIR, "Fig_node06_Differential_CCC_Top20_v4.pdf"),
         p_diff, width = 11, height = 7)
  msg("  Differential CCC v4 done (source cell type coloring, v1 logic restored)")
} else {
  msg("  No differential LR pairs found - skipping Differential CCC plot.")
}

# ============================================================
# ============================================================
msg("========== 5. ActHigh LR (clean labels) ==========")

acthigh <- read.csv(file.path(OUT_TABLES, "node06_CapEC_ActHigh_Curated_LR_Candidates_UC.csv"))
acthigh$target_clean <- clean_label(acthigh$target)
acthigh$lr_pair_clean <- clean_label(acthigh$lr_pair)

target_colors <- c(
  "Astro Reactive DAA" = "#E69F00",
  "Micro DAM" = "#56B4E9",
  "Micro Homeostatic" = "#009E73",
  "Excitatory" = "#F0E442",
  "Inhibitory" = "#CC79A7",
  "Astro Intermediate" = "#D55E00"
)

acthigh_top <- acthigh %>% arrange(desc(abs(delta_high_vs_low))) %>% slice_head(n = 15)

p_act <- ggplot(acthigh_top,
                aes(x = reorder(lr_pair_clean, delta_high_vs_low),
                    y = delta_high_vs_low, fill = target_clean)) +
  geom_bar(stat = "identity", width = 0.65) +
  coord_flip() +
  scale_fill_manual(values = target_colors, name = "Target Cell Type") +
  labs(title = "Top LR pairs: ActHigh vs ActLow Capillary EC (UC)",
       subtitle = "Positive: stronger in ActHigh | Target cell types colored",
       x = "Ligand-Receptor Pair", y = "Delta communication strength") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10),
        legend.position = "right",
        legend.text = element_text(size = 8),
        plot.margin = margin(8, 8, 8, 100)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")

ggsave(file.path(OUT_DIR, "Fig_node06_CapEC_ActHigh_LR_v4.png"),
       p_act, width = 10, height = 7, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig_node06_CapEC_ActHigh_LR_v4.pdf"),
       p_act, width = 10, height = 7)
msg("  ActHigh LR v3 done")

# ============================================================
# ============================================================
n_files <- length(list.files(OUT_DIR, "\\.(png|pdf)$"))
msg("========== ALL DONE ==========")
msg("Output: %s", OUT_DIR)
msg("Total files: %d", n_files)
cat("\n=== Module Score Exact P Values ===\n")
print(stat_results, digits = 4)


# Part3 SKIPPED - redundant with Part1 analysis
}
if (FALSE) {
# ========== Part3: Integrated stats ==========
# File: scripts_final/node06_integrated.R

################################################################################
#
#
#   Rscript node06_integrated.R
#
#
################################################################################

# ################################################################################
# ################################################################################

PROJECT_DIR <- "/path/to/R_libraries"

cat("==================================================\n")
cat("processing...\n")
cat("==================================================\n\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(ggpubr)
})

cat("processing...\n")


capec_ms_file <- "rds/node06/node06_CapEC_with_ModuleScores.qs"
stopifnot("CapEC score file is missing" = file.exists(capec_ms_file))

capec_ms <- qread(capec_ms_file)
cat("processing...\n", nrow(capec_ms), "processing...\n")
cat("processing...\n")
print(table(capec_ms$Group, useNA = "ifany"))

stopifnot("CapEC cell count is insufficient" = nrow(capec_ms) > 1000)
cat("processing...\n")

# ################################################################################
# ################################################################################
cat("processing...\n")

state_comp_file <- "tables/node06_CapEC_ActivationState_Composition.csv"
if (file.exists(state_comp_file)) {
  capec_state <- read.csv(state_comp_file)
  cat("processing...\n", state_comp_file, "\n")
} else {
  stop("Run the communication workflow before generating the state-composition table")
}

run_wilcoxon <- function(df, state_col, group_col, val_col, state_val) {
  d <- df[df[[state_col]] == state_val, ]
  nc_vals <- d[[val_col]][d[[group_col]] == "NC"]
  uc_vals <- d[[val_col]][d[[group_col]] == "UC"]
  if (length(nc_vals) < 3 || length(uc_vals) < 3) return(NA_real_)
  p <- wilcox.test(nc_vals, uc_vals, exact = FALSE)$p.value
  wilcox.test(nc_vals, uc_vals, exact = FALSE)$p.value  # raw p-value
}

get_sig <- function(p) {
  if (is.na(p) || length(p) > 1 || p >= 0.05) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  return("*")
}

states <- c("CapEC_ActLow", "CapEC_ActMid", "CapEC_ActHigh")
p_vals_state <- sapply(states, function(st)
  run_wilcoxon(capec_state, "CapEC_activation_state", "Group", "frac", st))

cat("processing...\n")
print(p_vals_state)

theme_n <- theme_bw(base_size = 8) + theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(linewidth = 0.4, colour = "black"),
  axis.ticks = element_line(linewidth = 0.4),
  plot.title = element_text(size = 8, face = "bold", hjust = 0)
)

state_plots <- lapply(seq_along(states), function(i) {
  st <- states[i]
  d <- capec_state[capec_state$CapEC_activation_state == st, ]
  pv <- as.numeric(p_vals_state[[st]])
  sg <- get_sig(pv)
  ggplot(d, aes(x = Group, y = frac, fill = Group, color = Group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.25, alpha = 0.75) +
    geom_jitter(width = 0.12, size = 1.1, alpha = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = c("NC" = "#4DBBD5", "UC" = "#E64B35")) +
    scale_color_manual(values = c("NC" = "#4DBBD5", "UC" = "#E64B35")) +
    theme_n +
    labs(title = gsub("CapEC_", "", st), x = NULL,
         y = if(i == 1) "Fraction per sample") +
    ylim(0, 1.25) +
    annotate("text", x = 1.5, y = 1.12, label = sg,
             size = 5, fontface = "bold", color = "black")
})

p_state <- wrap_plots(state_plots, nrow = 1) +
  plot_annotation(
    title = "B  CapEC activation-state (Wilcoxon, * p<0.05, ** p<0.01, *** p<0.001)",
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0))
  )

ggsave(file.path(OUT, "main_figures", "Fig_node06_CapEC_ActivationState_Composition_integrated.png"),
       p_state, width = 7.5, height = 3.0, dpi = 300, units = "cm")
cat("processing...\n")

# ################################################################################
# ################################################################################
cat("processing...\n")

cell_score_file <- "tables/node06_CapEC_ModuleScore_CellLevel.csv"
capec_cell <- read.csv(cell_score_file)
cat("processing...\n", nrow(capec_cell), "processing...\n")

score_cols_z <- c(
  "CapEC_BBB_maintenance_score_z",
  "CapEC_Endothelial_activation_score_z",
  "CapEC_Inflammatory_score_z",
  "CapEC_IFN_response_score_z",
  "CapEC_ECM_remodeling_score_z",
  "CapEC_Permeability_angiogenic_score_z",
  "CapEC_Transport_score_z"
)

score_labels <- c(
  "BBB maintenance",
  "Endothelial activation",
  "Inflammatory",
  "IFN response",
  "ECM remodeling",
  "Permeability/angiogenic",
  "Transport"
)
names(score_labels) <- score_cols_z

score_long <- reshape2::melt(
  capec_cell[, c("Group", score_cols_z)],
  id.vars = "Group",
  measure.vars = score_cols_z,
  variable.name = "score_name",
  value.name = "score_value"
)
score_long$score_label <- score_labels[as.character(score_long$score_name)]

score_pvals <- sapply(score_cols_z, function(sn) {
  d <- score_long[score_long$score_name == sn, ]
  nc <- d$score_value[d$Group == "NC"]
  uc <- d$score_value[d$Group == "UC"]
  if (length(nc) < 3 || length(uc) < 3) return(NA_real_)
  p <- wilcox.test(nc, uc, exact = FALSE)$p.value
  wilcox.test(nc, uc, exact = FALSE)$p.value  # raw p-value
})
cat("processing...\n")
score_padj <- p.adjust(score_pvals, method = "BH")
print(score_pvals)

score_plots <- lapply(score_cols_z, function(sn) {
  d <- score_long[score_long$score_name == sn, ]
  pv <- as.numeric(score_pvals[sn])
  sg <- get_sig(pv)
  max_y <- max(d$score_value, na.rm = TRUE)
  ggplot(d, aes(x = Group, y = score_value, fill = Group)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.6, linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.2) +
    scale_fill_manual(values = c("NC" = "#4DBBD5", "UC" = "#E64B35")) +
    theme_n +
    labs(title = score_labels[sn], x = NULL, y = NULL) +
    ylim(range(d$score_value) * c(1, 1.18)) +
    annotate("text", x = 1.5, y = max_y * 1.1, label = sg,
             size = 4, fontface = "bold", color = "black")
})

p_violin <- wrap_plots(score_plots, nrow = 2, ncol = 4) +
  plot_annotation(
    title = "CapEC module scores (Wilcoxon, * p<0.05, ** p<0.01, *** p<0.001)",
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0))
  )

ggsave(file.path(OUT, "main_figures", "Fig_node06_CapEC_ModuleScores_integrated.png"),
       p_violin, width = 9, height = 5.5, dpi = 300, units = "cm")
cat("processing...\n")

# ################################################################################
# ################################################################################
cat("processing...\n")

go_df <- read.csv(file.path(OUT, "tables", "node06_GO_BP_enrichment_all.csv"))
react_df <- read.csv(file.path(OUT, "tables", "node06_Reactome_enrichment_all.csv"))

noise_pat <- paste0(c(
  "\\bkidney\\b", "\\breal\\b", "\\burogenital\\b",
  "\\bbone development\\b", "\\bskeletal\\b", "\\bcartilage\\b",
  "\\bmuscle contraction\\b", "\\bmuscle organ\\b",
  "\\bcell cycle\\b", "\\bDNA replication\\b",
  "\\bchromosome segregation\\b", "\\bmitotic\\b",
  "\\bcell division\\b", "\\bspindle\\b",
  "\\bB cell\\b", "\\bT cell\\b",
  "\\bxenobiotic\\b", "\\bdrug catabolic\\b",
  "\\bcobalamin\\b", "\\bretinol\\b"
), collapse = "|")

is_noise_vec <- function(x) {
  x <- as.character(x)
  idx <- is.na(x) | nchar(x) == 0
  res <- logical(length(x))
  res[!idx] <- grepl(noise_pat, x[!idx], ignore.case = TRUE)
  res
}

celltype_groups <- list(
  "EC" = c("Capillary", "Venous", "Arterial"),
  "Astrocyte"   = c("Homeostatic", "Intermediate", "Reactive_DAA"),
  "Microglia"   = c("DAM", "Homeostatic", "Inflammatory")
)

make_dotplot <- function(df, ct_name, subtypes, top_n = 8, color_low, color_high) {
  d <- df[df$subtype %in% subtypes & !is_noise_vec(df$Description) & df$p.adjust < 0.05, ]
  if (nrow(d) == 0) { cat("  [WARN]", ct_name, "processing...\n"); return(NULL) }
  d <- d[order(d$p.adjust), ]
  d <- do.call(rbind, lapply(split(d, d$subtype), function(x) head(x, top_n)))
  if (nrow(d) == 0) return(NULL)
  d$score <- -log10(d$p.adjust + 1e-12)
  d$Description <- substring(d$Description, 1, 60)
  ggplot(d, aes(x = subtype, y = reorder(Description, score),
                size = Count, color = score)) +
    geom_point() +
    scale_color_gradient(low = color_low, high = color_high) +
    scale_size(range = c(2, 5)) +
    theme_n +
    labs(title = paste0(ct_name, " GO (BP)"), x = NULL, y = NULL,
         color = "-log10(FDR)", size = "Gene count") +
    theme(axis.text.y = element_text(size = 6),
          plot.margin = margin(5, 5, 5, 120))
}

for (ct_name in names(celltype_groups)) {
  p_go <- make_dotplot(go_df, ct_name, celltype_groups[[ct_name]],
                        color_low = "#4DBBD5", color_high = "#E64B35")
  if (!is.null(p_go)) {
    ggsave(file.path(OUT, "main_figures", paste0("Fig_node06_GO_", ct_name,
                  "_DotPlot_integrated.png")),
           p_go, width = 7.5, height = 6, dpi = 300, units = "cm")
  }
}

react_pat <- paste0(c(
  "\\bkidney\\b", "\\breal\\b", "\\bmuscle\\b",
  "\\bcell cycle\\b", "\\bDNA replication\\b", "\\bchromosome\\b",
  "\\bxenobiotic\\b", "\\bdrug metabolic\\b"
), collapse = "|")

is_react_noise_vec <- function(x) {
  x <- as.character(x)
  idx <- is.na(x) | nchar(x) == 0
  res <- logical(length(x))
  res[!idx] <- grepl(react_pat, x[!idx], ignore.case = TRUE)
  res
}

make_react_dotplot <- function(df, ct_name, subtypes, top_n = 8) {
  d <- df[df$subtype %in% subtypes & !is_react_noise_vec(df$Description) & df$p.adjust < 0.05, ]
  if (nrow(d) == 0) return(NULL)
  d <- d[order(d$p.adjust), ]
  d <- do.call(rbind, lapply(split(d, d$subtype), function(x) head(x, top_n)))
  if (nrow(d) == 0) return(NULL)
  d$score <- -log10(d$p.adjust + 1e-12)
  d$Description <- substring(d$Description, 1, 60)
  ggplot(d, aes(x = subtype, y = reorder(Description, score),
                size = Count, color = score)) +
    geom_point() +
    scale_color_gradient(low = "#4DBBD5", high = "#3C5488") +
    scale_size(range = c(2, 5)) +
    theme_n +
    labs(title = paste0(ct_name, " Reactome"), x = NULL, y = NULL,
         color = "-log10(FDR)", size = "Gene count") +
    theme(axis.text.y = element_text(size = 6),
          plot.margin = margin(5, 5, 5, 120))
}

for (ct_name in names(celltype_groups)) {
  p_r <- make_react_dotplot(react_df, ct_name, celltype_groups[[ct_name]])
  if (!is.null(p_r)) {
    ggsave(file.path(OUT, "main_figures", paste0("Fig_node06_Reactome_", ct_name,
                  "_DotPlot_integrated.png")),
           p_r, width = 7.5, height = 6, dpi = 300, units = "cm")
  }
}

cat("processing...\n")

# ################################################################################
# ################################################################################
cat("processing...\n")

cat("==================================================\n")
cat("processing...\n")
cat("==================================================\n")
cat("processing...\n")
cat("processing...\n")
for (st in states) {
  cat("  ", st, ": p =", round(p_vals_state[[st]], 4),
      "->", get_sig(p_vals_state[[st]]), "\n")
}
cat("processing...\n")
for (nm in names(score_pvals)) {
  cat("  ", score_labels[nm], ": p =", round(score_pvals[nm], 4),
      "->", get_sig(score_pvals[nm]), "\n")
}
cat("processing...\n")
cat("processing...\n")
sig_modules <- names(which(score_pvals < 0.05))
if (length(sig_modules) > 0) {
  cat("processing...\n")
  for (m in sig_modules) cat("     -", score_labels[m], "(p=", round(score_pvals[m], 4), ")\n")
} else {
  cat("processing...\n")
}
cat("processing...\n")
cat("==================================================\n")
}
cat("Part3 skipped (redundant with Part1)\n")
