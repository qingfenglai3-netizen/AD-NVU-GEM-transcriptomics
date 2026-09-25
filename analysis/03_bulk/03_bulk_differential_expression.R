#!/usr/bin/env Rscript
################################################################################
# node03 bulk differential expression, publication-ready AD-NVU version
#
# Purpose
#   Analyze public MTG bulk RNA-seq data (GSE132903) as the bulk transcriptomic
#   anchor for downstream AD-NVU candidate prioritization.
#
# Main decisions
#   - Diagnosis labels are standardized as CN and AD for all manuscript outputs.
#   - Differential expression uses limma with AD versus CN.
#   - The main node figure emphasizes sample separation, DE effect size, and
#     top-gene expression patterns. Gene lists are retained as supplementary
#     tables rather than early-node main tables.
#
# Reproducibility
#   R 4.4.1-compatible; fixed seed; source data exported for each figure.
################################################################################

set.seed(123456)
options(stringsAsFactors = FALSE)

paths_to_add <- "/path/to/R_libraries"
for (p in paths_to_add) {
  if (dir.exists(p) && !(p %in% .libPaths())) .libPaths(c(p, .libPaths()))
}

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggrepel)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
  library(grid)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
TEST_MODE <- TRUE
if (TEST_MODE) {
  RAW_DIR <- "/path/to/public_data"
  OUT <- "/path/to/project/results/03_bulk"
} else {
  RAW_DIR <- "/path/to/public_data"
  OUT <- "/path/to/project/results/03_bulk"
}

IN_ESET <- "/path/to/public_data/bulk/GSE132903_eSet.Rdata"
IN_ANNOT <- "/path/to/public_data/bulk/GPL10558_bioc.rda"

OUT_MAIN <- file.path(OUT, "main_figures")
OUT_SUPP <- file.path(OUT, "supplementary_figures")
OUT_EXT <- file.path(OUT, "extended_data_figures")
OUT_TAB <- file.path(OUT, "tables")
OUT_RDS <- file.path(OUT, "rds")
OUT_LOG <- file.path(OUT, "logs")
OUT_SRC <- file.path(OUT, "source_data")
OUT_METHODS <- file.path(OUT, "methods")
node01_MARKERS <- "/path/to/project/results/01_atlas/outputs/supplementary_tables/Supplementary_Table_7_node01_all_cluster_markers.csv"
node02_ANNOT <- "/path/to/project/results/02_annotation/outputs/supplementary_tables/Supplementary_Table_2_node02_cluster_annotation_initial.csv"

OUT_FINAL <- file.path(OUT, "outputs")
OUT_F_MAIN <- file.path(OUT_FINAL, "main_figures")
OUT_F_SUPP <- file.path(OUT_FINAL, "supplementary_figures")
OUT_F_EXT <- file.path(OUT_FINAL, "extended_data_figures")
OUT_T_MAIN <- file.path(OUT_FINAL, "main_tables")
OUT_T_SUPP <- file.path(OUT_FINAL, "supplementary_tables")
OUT_METH <- file.path(OUT_FINAL, "methods_results")
OUT_EDIT <- file.path(OUT_FINAL, "editable_panels")

for (d in c(OUT_MAIN, OUT_SUPP, OUT_EXT, OUT_TAB, OUT_RDS, OUT_LOG, OUT_SRC,
            OUT_METHODS, OUT_F_MAIN, OUT_F_SUPP, OUT_F_EXT, OUT_T_MAIN,
            OUT_T_SUPP, OUT_METH, OUT_EDIT)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(OUT_LOG, "03_bulk_run.log")
log_msg <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n")
  tryCatch(cat(line, "\n", file = LOG_FILE, append = TRUE), error = function(e) NULL)
}

# ------------------------------------------------------------------------------
# Plot helpers
# ------------------------------------------------------------------------------
publication_theme <- function(base = 7) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      axis.text = element_text(color = "black", size = base),
      axis.title = element_text(color = "black", size = base + 0.5),
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      legend.text = element_text(size = base - 0.5),
      legend.title = element_text(size = base),
      legend.key.size = unit(0.32, "cm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base, face = "bold"),
      plot.title = element_text(size = base + 1, face = "bold", hjust = 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
}

group_palette <- c("CN" = "#4DBBD5", "AD" = "#E64B35")
de_palette <- c("Up in AD" = "#D55E00", "Down in AD" = "#0072B2", "NS" = "grey82")

save_fig <- function(p, name, dir, w = 6, h = 4) {
  files <- c(pdf = file.path(dir, paste0(name, ".pdf")),
             png = file.path(dir, paste0(name, ".png")),
             svg = file.path(dir, paste0(name, ".svg")),
             tiff = file.path(dir, paste0(name, ".tiff")))
  tryCatch(ggsave(files["pdf"], p, width = w, height = h, units = "in", device = cairo_pdf),
           error = function(e) log_msg(sprintf("PDF save failed for %s: %s", name, e$message)))
  tryCatch(ggsave(files["png"], p, width = w, height = h, units = "in", dpi = 300),
           error = function(e) log_msg(sprintf("PNG save failed for %s: %s", name, e$message)))
  tryCatch(ggsave(files["svg"], p, width = w, height = h, units = "in"),
           error = function(e) log_msg(sprintf("SVG save failed for %s: %s", name, e$message)))
  tryCatch(ggsave(files["tiff"], p, width = w, height = h, units = "in", dpi = 600,
                  compression = "lzw", device = "tiff"),
           error = function(e) log_msg(sprintf("TIFF save failed for %s: %s", name, e$message)))
  log_msg(sprintf("Saved %s.{pdf,png,svg,tiff}", name))
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE)
}

save_source <- function(x, name) {
  safe_write_csv(x, file.path(OUT_SRC, paste0(name, "_source_data.csv")))
}

copy_multi_format <- function(src_base, dst_base, src_dir, dst_dir) {
  for (ext in c("pdf", "png", "svg", "tiff")) {
    src <- file.path(src_dir, paste0(src_base, ".", ext))
    dst <- file.path(dst_dir, paste0(dst_base, ".", ext))
    if (file.exists(src)) file.copy(src, dst, overwrite = TRUE)
  }
}

copy_editable <- function(src_base, dst_base, src_dir, dst_dir) {
  for (ext in c("pdf", "svg")) {
    src <- file.path(src_dir, paste0(src_base, ".", ext))
    dst <- file.path(dst_dir, paste0(dst_base, ".", ext))
    if (file.exists(src)) file.copy(src, dst, overwrite = TRUE)
  }
}

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
log_msg("===== node03 bulk differential expression =====")
stopifnot(file.exists(IN_ESET), file.exists(IN_ANNOT))
load(IN_ESET)  # gset
gse <- gset[[1]]
pdata <- pData(gse)
expr_probe <- exprs(gse)

diag_raw <- gsub("diagnosis: ", "", pdata$characteristics_ch1.3)
group <- ifelse(diag_raw == "AD", "AD", ifelse(diag_raw == "ND", "CN", NA))
if (anyNA(group)) stop("Unknown diagnosis labels detected.")
group <- factor(group, levels = c("CN", "AD"))

sample_meta <- pdata %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(diagnosis_raw = diag_raw, group = as.character(group))

log_msg(sprintf("Samples: CN=%d | AD=%d", sum(group == "CN"), sum(group == "AD")))
log_msg(sprintf("Probe expression matrix: %d probes x %d samples", nrow(expr_probe), ncol(expr_probe)))

load(IN_ANNOT)  # GPL10558_bioc
anno <- GPL10558_bioc
colnames(anno) <- c("probe_id", "symbol")
anno$symbol <- trimws(as.character(anno$symbol))
common_probe <- intersect(rownames(expr_probe), anno$probe_id)
expr_probe <- expr_probe[common_probe, , drop = FALSE]
anno <- anno[match(common_probe, anno$probe_id), , drop = FALSE]
keep <- !is.na(anno$symbol) & anno$symbol != ""
expr_probe <- expr_probe[keep, , drop = FALSE]
anno <- anno[keep, , drop = FALSE]
log_msg(sprintf("Annotated probes: %d", nrow(expr_probe)))

# Probe-level limma first, then keep the strongest absolute logFC probe per gene.
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "AD_vs_CN")
fit <- eBayes(lmFit(expr_probe, design))
all_probe <- topTable(fit, coef = "AD_vs_CN", adjust.method = "fdr", number = Inf)
all_probe$probe_id <- rownames(all_probe)
all_probe <- left_join(all_probe, anno, by = "probe_id")

all_diff <- all_probe %>%
  filter(!is.na(symbol), symbol != "") %>%
  group_by(symbol) %>%
  slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(adj.P.Val, desc(abs(logFC)))

selected_probes <- all_diff$probe_id
expr_gene <- expr_probe[selected_probes, , drop = FALSE]
rownames(expr_gene) <- all_diff$symbol

TH_FC <- 0.58
TH_FDR <- 0.01
all_diff <- all_diff %>%
  mutate(
    neg_log10_fdr = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
    category = case_when(
      logFC > TH_FC & adj.P.Val < TH_FDR ~ "Up in AD",
      logFC < -TH_FC & adj.P.Val < TH_FDR ~ "Down in AD",
      TRUE ~ "NS"
    )
  )
up_genes <- all_diff %>% filter(category == "Up in AD")
down_genes <- all_diff %>% filter(category == "Down in AD")
log_msg(sprintf("DE genes: Up in AD=%d | Down in AD=%d", nrow(up_genes), nrow(down_genes)))

safe_write_csv(sample_meta, file.path(OUT_TAB, "bulk_GSE132903_sample_metadata.csv"))
safe_write_csv(all_diff, file.path(OUT_TAB, "bulk_GSE132903_DE_all_genes.csv"))
safe_write_csv(up_genes, file.path(OUT_TAB, "bulk_GSE132903_DE_up_in_AD.csv"))
safe_write_csv(down_genes, file.path(OUT_TAB, "bulk_GSE132903_DE_down_in_AD.csv"))
writeLines(up_genes$symbol, file.path(OUT_TAB, "bulk_up_geneset_AD_vs_CN.txt"))
save(all_diff, up_genes, down_genes, sample_meta, expr_gene, group,
     file = file.path(OUT_RDS, "bulk_GSE132903_node03_objects.Rdata"))

# ------------------------------------------------------------------------------
# Figure data
# ------------------------------------------------------------------------------
expr_gene_z <- t(scale(t(expr_gene)))
expr_gene_z[is.na(expr_gene_z)] <- 0
pca <- prcomp(t(expr_gene_z), center = FALSE, scale. = FALSE)
pca_var <- (pca$sdev^2) / sum(pca$sdev^2)
pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  group = as.character(group)
)

label_genes <- all_diff %>%
  filter(category != "NS") %>%
  arrange(adj.P.Val, desc(abs(logFC))) %>%
  group_by(category) %>%
  slice_head(n = 8) %>%
  ungroup()

top_heat_genes <- all_diff %>%
  filter(category != "NS") %>%
  arrange(adj.P.Val, desc(abs(logFC))) %>%
  slice_head(n = 40) %>%
  pull(symbol)
top_heat_genes <- intersect(top_heat_genes, rownames(expr_gene_z))
heat_df <- as.data.frame(expr_gene_z[top_heat_genes, , drop = FALSE]) %>%
  rownames_to_column("symbol") %>%
  pivot_longer(-symbol, names_to = "sample_id", values_to = "z") %>%
  left_join(sample_meta[, c("sample_id", "group")], by = "sample_id")
heat_df$symbol <- factor(heat_df$symbol, levels = rev(top_heat_genes))
sample_order <- pca_df %>% arrange(group, PC1) %>% pull(sample_id)
heat_df$sample_id <- factor(heat_df$sample_id, levels = sample_order)

summary_df <- data.frame(
  metric = c("Samples_CN", "Samples_AD", "Genes_tested", "Up_in_AD", "Down_in_AD",
             "log2FC_cutoff", "FDR_cutoff"),
  value = c(sum(group == "CN"), sum(group == "AD"), nrow(all_diff),
            nrow(up_genes), nrow(down_genes), TH_FC, TH_FDR)
)
safe_write_csv(summary_df, file.path(OUT_TAB, "bulk_GSE132903_DE_summary.csv"))

# ------------------------------------------------------------------------------
# Main figures
# ------------------------------------------------------------------------------
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 1.9, alpha = 0.9) +
  stat_ellipse(linewidth = 0.35, level = 0.68, show.legend = FALSE) +
  scale_color_manual(values = group_palette) +
  labs(
    title = "Bulk sample structure",
    x = sprintf("PC1 (%.1f%%)", 100 * pca_var[1]),
    y = sprintf("PC2 (%.1f%%)", 100 * pca_var[2]),
    color = NULL
  ) +
  publication_theme(7) +
  theme(aspect.ratio = 1)
save_fig(p_pca, "Fig3A_bulk_PCA", OUT_MAIN, w = 4.2, h = 3.8)
save_source(pca_df, "Fig3A_bulk_PCA")

p_volcano <- ggplot(all_diff, aes(logFC, neg_log10_fdr, color = category)) +
  geom_point(size = 0.65, alpha = 0.55, stroke = 0) +
  geom_vline(xintercept = c(-TH_FC, TH_FC), linetype = "dashed",
             color = "grey40", linewidth = 0.28) +
  geom_hline(yintercept = -log10(TH_FDR), linetype = "dashed",
             color = "grey40", linewidth = 0.28) +
  geom_text_repel(
    data = label_genes,
    aes(label = symbol),
    size = 2,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.25,
    segment.size = 0.18,
    show.legend = FALSE
  ) +
  scale_color_manual(values = de_palette, breaks = c("Up in AD", "Down in AD", "NS")) +
  labs(
    title = "AD versus CN differential expression",
    x = "log2 fold change (AD/CN)",
    y = "-log10(FDR)",
    color = NULL
  ) +
  publication_theme(7)
save_fig(p_volcano, "Fig3B_bulk_volcano", OUT_MAIN, w = 5.2, h = 4.2)
save_source(all_diff, "Fig3B_bulk_volcano")

p_heat <- ggplot(heat_df, aes(sample_id, symbol, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                       midpoint = 0, limits = c(-2.5, 2.5), oob = squish,
                       name = "z-score") +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  labs(title = "Top bulk DE genes", x = NULL, y = NULL) +
  publication_theme(5.5) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 4.5),
    strip.text = element_text(size = 6, face = "bold"),
    legend.position = "right"
  )
save_fig(p_heat, "Fig3C_bulk_top_DE_heatmap", OUT_MAIN, w = 6.2, h = 5.8)
save_source(heat_df, "Fig3C_bulk_top_DE_heatmap")

bar_df <- data.frame(
  category = factor(c("Up in AD", "Down in AD"), levels = c("Up in AD", "Down in AD")),
  n = c(nrow(up_genes), nrow(down_genes))
)
p_bar <- ggplot(bar_df, aes(category, n, fill = category)) +
  geom_col(width = 0.62, color = "black", linewidth = 0.25) +
  geom_text(aes(label = n), vjust = -0.35, size = 2.5) +
  scale_fill_manual(values = de_palette) +
  labs(title = "Significant DE genes", x = NULL, y = "Genes, n") +
  publication_theme(7) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
save_fig(p_bar, "Fig3D_bulk_DE_gene_counts", OUT_MAIN, w = 3.2, h = 3.4)
save_source(bar_df, "Fig3D_bulk_DE_gene_counts")

fig3_final <- (wrap_elements(full = p_pca) | wrap_elements(full = p_volcano)) /
  (wrap_elements(full = p_heat) | wrap_elements(full = p_bar)) +
  plot_layout(widths = c(1.05, 1), heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11))
save_fig(fig3_final, "Figure_3_node03_bulk_DE_FINAL", OUT_MAIN, w = 13, h = 10)

# ------------------------------------------------------------------------------
# Supplementary / extended figures
# ------------------------------------------------------------------------------
p_density <- as.data.frame(expr_gene) %>%
  rownames_to_column("symbol") %>%
  pivot_longer(-symbol, names_to = "sample_id", values_to = "expr") %>%
  left_join(sample_meta[, c("sample_id", "group")], by = "sample_id") %>%
  ggplot(aes(expr, group = sample_id, color = group)) +
  geom_density(linewidth = 0.25, alpha = 0.25) +
  scale_color_manual(values = group_palette) +
  labs(title = "Bulk expression density by sample", x = "Expression", y = "Density", color = NULL) +
  publication_theme(7)
save_fig(p_density, "FigS9_bulk_expression_density", OUT_SUPP, w = 5.5, h = 4)

p_pval <- ggplot(all_diff, aes(P.Value)) +
  geom_histogram(bins = 50, fill = "grey70", color = "white", linewidth = 0.1) +
  labs(title = "Nominal P-value distribution", x = "P value", y = "Genes, n") +
  publication_theme(7)
save_fig(p_pval, "FigS10_bulk_pvalue_distribution", OUT_SUPP, w = 4.8, h = 3.6)
save_source(all_diff[, c("symbol", "P.Value", "adj.P.Val", "logFC")], "FigS10_bulk_pvalue_distribution")

# ------------------------------------------------------------------------------
# Bulk-to-single-cell bridge: are AD bulk DE genes enriched in annotated cell types?
# ------------------------------------------------------------------------------
overlap_results <- NULL
if (file.exists(node01_MARKERS) && file.exists(node02_ANNOT)) {
  log_msg("Computing bulk DE overlap with node02 single-cell marker programs")
  marker_tbl <- read.csv(node01_MARKERS, check.names = FALSE)
  annot_tbl <- read.csv(node02_ANNOT, check.names = FALSE)
  marker_tbl$cluster_key <- paste0("g", as.character(marker_tbl$cluster))
  marker_tbl <- marker_tbl %>%
    left_join(annot_tbl[, c("cluster", "initial_celltype", "confidence")],
              by = c("cluster_key" = "cluster")) %>%
    filter(!is.na(initial_celltype), !is.na(gene), gene != "",
           avg_log2FC > 0.25, p_val_adj < 0.05) %>%
    group_by(initial_celltype, gene) %>%
    summarise(best_log2FC = max(avg_log2FC, na.rm = TRUE),
              best_adj_p = min(p_val_adj, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(initial_celltype, best_adj_p, desc(best_log2FC)) %>%
    group_by(initial_celltype) %>%
    slice_head(n = 150) %>%
    ungroup()

  universe <- unique(all_diff$symbol)
  de_sets <- list("Up in AD" = unique(up_genes$symbol),
                  "Down in AD" = unique(down_genes$symbol))
  overlap_results <- lapply(names(de_sets), function(direction) {
    de_gene_set <- intersect(de_sets[[direction]], universe)
    lapply(sort(unique(marker_tbl$initial_celltype)), function(ct) {
      marker_set <- intersect(unique(marker_tbl$gene[marker_tbl$initial_celltype == ct]), universe)
      a <- length(intersect(de_gene_set, marker_set))
      b <- length(setdiff(de_gene_set, marker_set))
      c <- length(setdiff(marker_set, de_gene_set))
      d <- length(setdiff(universe, union(de_gene_set, marker_set)))
      ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
      data.frame(
        direction = direction,
        celltype = ct,
        overlap_genes = a,
        de_genes = length(de_gene_set),
        marker_genes = length(marker_set),
        universe_genes = length(universe),
        odds_ratio = unname(ft$estimate),
        p_value = ft$p.value,
        overlap_gene_symbols = paste(sort(intersect(de_gene_set, marker_set)), collapse = ";"),
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()
  }) %>% bind_rows() %>%
    mutate(
      fdr = p.adjust(p_value, method = "BH"),
      neg_log10_fdr = -log10(pmax(fdr, .Machine$double.xmin)),
      direction = factor(direction, levels = c("Up in AD", "Down in AD")),
      celltype = factor(celltype, levels = rev(c(
        "Excitatory", "Inhibitory", "Oligodendrocytes", "OPCs",
        "Astrocytes", "Microglia", "Cerebrovascular cells"
      )))
    )

  safe_write_csv(overlap_results,
                 file.path(OUT_TAB, "bulk_scRNA_celltype_marker_overlap_enrichment.csv"))

  p_overlap <- ggplot(overlap_results,
                      aes(direction, celltype, size = overlap_genes,
                          color = neg_log10_fdr)) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(range = c(1.5, 7), breaks = c(5, 10, 20, 40),
                          name = "Overlap genes") +
    scale_color_gradient(low = "grey80", high = "#B40426",
                         name = "-log10(FDR)") +
    labs(
      title = "Bulk DE genes overlap single-cell marker programs",
      x = NULL,
      y = NULL
    ) +
    publication_theme(7) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      legend.position = "right"
    )
  save_fig(p_overlap, "FigS11_bulk_scRNA_celltype_marker_overlap", OUT_SUPP, w = 5.4, h = 4.4)
  save_source(overlap_results, "FigS11_bulk_scRNA_celltype_marker_overlap")
} else {
  log_msg("Skipped bulk-to-single-cell marker overlap: node01 marker or annotation table not found")
}

# ------------------------------------------------------------------------------
# Final package
# ------------------------------------------------------------------------------
copy_multi_format("Figure_3_node03_bulk_DE_FINAL", "Figure_3_node03_bulk_DE_FINAL", OUT_MAIN, OUT_F_MAIN)
copy_multi_format("Fig3A_bulk_PCA", "Figure_3_node03_panelA_bulk_PCA", OUT_MAIN, OUT_F_MAIN)
copy_multi_format("Fig3B_bulk_volcano", "Figure_3_node03_panelB_bulk_volcano", OUT_MAIN, OUT_F_MAIN)
copy_multi_format("Fig3C_bulk_top_DE_heatmap", "Figure_3_node03_panelC_bulk_top_DE_heatmap", OUT_MAIN, OUT_F_MAIN)
copy_multi_format("Fig3D_bulk_DE_gene_counts", "Figure_3_node03_panelD_bulk_DE_gene_counts", OUT_MAIN, OUT_F_MAIN)

copy_multi_format("FigS9_bulk_expression_density", "Figure_S9_03_bulk_expression_density", OUT_SUPP, OUT_F_SUPP)
copy_multi_format("FigS10_bulk_pvalue_distribution", "Figure_S10_03_bulk_pvalue_distribution", OUT_SUPP, OUT_F_SUPP)
copy_multi_format("FigS11_bulk_scRNA_celltype_marker_overlap", "Figure_S11_node03_bulk_scRNA_celltype_marker_overlap", OUT_SUPP, OUT_F_SUPP)

for (nm in c("Figure_3_node03_bulk_DE_FINAL", "Fig3A_bulk_PCA", "Fig3B_bulk_volcano",
             "Fig3C_bulk_top_DE_heatmap", "Fig3D_bulk_DE_gene_counts",
             "FigS9_bulk_expression_density", "FigS10_bulk_pvalue_distribution",
             "FigS11_bulk_scRNA_celltype_marker_overlap")) {
  src_dir <- ifelse(grepl("^FigS", nm), OUT_SUPP, OUT_MAIN)
  copy_editable(nm, paste0("node03_", nm, "_edit"), src_dir, OUT_EDIT)
}

file.copy(file.path(OUT_TAB, "bulk_GSE132903_DE_all_genes.csv"),
          file.path(OUT_T_SUPP, "Supplementary_Table_6_node03_bulk_DE_all_genes.csv"),
          overwrite = TRUE)
file.copy(file.path(OUT_TAB, "bulk_GSE132903_DE_up_in_AD.csv"),
          file.path(OUT_T_SUPP, "Supplementary_Table_7_node03_bulk_up_in_AD.csv"),
          overwrite = TRUE)
file.copy(file.path(OUT_TAB, "bulk_GSE132903_DE_down_in_AD.csv"),
          file.path(OUT_T_SUPP, "Supplementary_Table_8_node03_bulk_down_in_AD.csv"),
          overwrite = TRUE)
file.copy(file.path(OUT_TAB, "bulk_GSE132903_sample_metadata.csv"),
          file.path(OUT_T_SUPP, "Supplementary_Table_9_node03_bulk_sample_metadata.csv"),
          overwrite = TRUE)
file.copy(file.path(OUT_TAB, "bulk_GSE132903_DE_summary.csv"),
          file.path(OUT_T_SUPP, "Supplementary_Table_10_node03_bulk_DE_summary.csv"),
          overwrite = TRUE)
if (file.exists(file.path(OUT_TAB, "bulk_scRNA_celltype_marker_overlap_enrichment.csv"))) {
  file.copy(file.path(OUT_TAB, "bulk_scRNA_celltype_marker_overlap_enrichment.csv"),
            file.path(OUT_T_SUPP, "Supplementary_Table_11_node03_bulk_scRNA_marker_overlap_enrichment.csv"),
            overwrite = TRUE)
}

manifest <- data.frame(
  item = c("Figure 3", "Figure S9", "Figure S10", "Figure S11", "Supplementary Tables 6-11"),
  description = c(
    "Bulk MTG DE overview: PCA, volcano plot, top DE heatmap, and DE gene counts.",
    "Per-sample expression density QC.",
    "Nominal P-value distribution for bulk DE model diagnostics.",
    "Overlap enrichment between bulk AD DE genes and node01 single-cell cell-type marker programs.",
    "Bulk DE all genes, up/down genes, sample metadata, summary metrics, and bulk-to-single-cell marker overlap enrichment."
  )
)
safe_write_csv(manifest, file.path(OUT_FINAL, "03_bulk_output_manifest.csv"))

method_lines <- c(
  "# node03 bulk differential expression",
  "",
  sprintf("Samples: CN=%d, AD=%d.", sum(group == "CN"), sum(group == "AD")),
  sprintf("Genes tested after probe-to-gene selection: %d.", nrow(all_diff)),
  sprintf("Significance threshold: |log2FC| > %.2f and FDR < %.2g.", TH_FC, TH_FDR),
  sprintf("Significant genes: %d up in AD and %d down in AD.", nrow(up_genes), nrow(down_genes)),
  "Bulk-to-single-cell bridge: significant bulk DE genes were tested for overlap with node02 cell-type marker programs using one-sided Fisher exact tests and Benjamini-Hochberg correction.",
  "",
  "Main figure recommendation: Figure 3; no main table is recommended at this early node.",
  "The full DE table, significant gene sets, and marker-overlap statistics are retained as supplementary tables."
)
writeLines(method_lines, file.path(OUT_METH, "node03_bulk_methods_results_summary.md"))

capture.output(sessionInfo(), file = file.path(OUT_LOG, "sessionInfo_03_bulk.txt"))
log_msg("===== node03 bulk differential expression complete =====")
