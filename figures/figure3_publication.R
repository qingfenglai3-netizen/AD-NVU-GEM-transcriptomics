################################################################################
# Rebuild MANUSCRIPT Figure 3 and related supplementary figures
# Claim: the bulk AD-up transcriptional signature can be localized to specific
# major cell compartments in the integrated snRNA-seq atlas.
################################################################################

options(stringsAsFactors = FALSE)
set.seed(20260531)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
})

ROOT <- "/path/to/project/results"
NODE03 <- file.path(ROOT, "03_bulk")
NODE04 <- file.path(ROOT, "04_subtypes")
OUT <- file.path(ROOT, "publication_outputs")
DIR_MAIN <- file.path(OUT, "main_figures")
DIR_SUPP <- file.path(OUT, "supplementary_figures")
DIR_TABLE <- file.path(OUT, "supplementary_tables")
DIR_LOG <- file.path(OUT, "decision_logs")
DIR_SRC <- file.path(OUT, "source_data")
for (d in c(DIR_MAIN, DIR_SUPP, DIR_TABLE, DIR_LOG, DIR_SRC)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cm_to_in <- function(cm) cm / 2.54

theme_publication <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial", face = "bold", colour = "black"),
      axis.text = element_text(size = max(base_size - 0.2, 6.0), face = "bold", colour = "black"),
      axis.title = element_text(size = base_size, face = "bold", colour = "black"),
      axis.title.y = element_text(size = base_size, face = "bold", margin = margin(r = 1), vjust = 0.5),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      legend.text = element_text(size = max(base_size - 0.5, 6.0), face = "bold", colour = "black"),
      legend.title = element_text(size = max(base_size - 0.3, 6.0), face = "bold", colour = "black"),
      legend.key.size = unit(3.0, "mm"),
      strip.text = element_text(size = base_size, face = "bold", colour = "black"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.35),
      plot.title = element_text(size = base_size + 0.2, face = "bold", hjust = 0, colour = "black"),
      panel.grid = element_blank(),
      plot.margin = margin(3, 3, 3, 3, unit = "mm")
    )
}

save_publication <- function(p, file, width_cm = 17.0, height_cm = 17.8, dpi = 600) {
  w <- cm_to_in(width_cm)
  h <- cm_to_in(height_cm)
  ggsave(paste0(file, ".pdf"), p, width = w, height = h, device = cairo_pdf, bg = "white")
  ggsave(paste0(file, ".svg"), p, width = w, height = h, device = "svg", bg = "white")
  ggsave(paste0(file, ".png"), p, width = w, height = h, dpi = 600, bg = "white")
  ggsave(paste0(file, ".tiff"), p, width = w, height = h, dpi = dpi, bg = "white",
         compression = "lzw", device = "tiff")
}

panel_label <- function(p, lab) {
  label <- ggplot() +
    theme_void(base_family = "Arial") +
    annotate("text", x = 0, y = 0.5, label = lab, hjust = 0, vjust = 0.5,
             family = "Arial", fontface = "bold", size = 12 / .pt) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme(plot.margin = margin(0, 0, 0, 0))
  wrap_plots(label, wrap_elements(full = p), ncol = 1, heights = c(0.075, 1))
}

read_csv <- function(path) read.csv(path, check.names = FALSE)

ct_levels <- c("Oligodendrocytes", "Excitatory", "Inhibitory", "Astrocytes",
               "Microglia", "OPCs", "Cerebrovascular cells")
ct_cols <- c(
  "Oligodendrocytes" = "#0072B2",
  "Excitatory" = "#999999",
  "Inhibitory" = "#56B4E9",
  "Astrocytes" = "#CC79A7",
  "Microglia" = "#D55E00",
  "OPCs" = "#009E73",
  "Cerebrovascular cells" = "#B79F00"
)
group_cols <- c("CN" = "#4DBBD5", "AD" = "#E64B35")
de_cols <- c("Up in AD" = "#D55E00", "Down in AD" = "#0072B2", "NS" = "#D0D0D0")

copy_source <- function(from, name) {
  to <- file.path(DIR_SRC, name)
  file.copy(from, to, overwrite = TRUE)
  to
}

volcano <- read_csv(file.path(NODE03, "source_data", "Fig3B_bulk_volcano_source_data.csv"))
overlap <- read_csv(file.path(NODE03, "source_data", "FigS11_bulk_scRNA_celltype_marker_overlap_source_data.csv"))
score_violin <- read_csv(file.path(NODE04, "source_data", "Fig1a_AD_signature_violin_source_data.csv"))
score_heat <- read_csv(file.path(NODE04, "source_data", "Fig1b_method_heatmap_source_data.csv"))
score_delta <- read_csv(file.path(NODE04, "source_data", "Fig1c_wilcoxon_delta_source_data.csv"))
score_umap <- read_csv(file.path(NODE04, "source_data", "Fig1d_AD_score_UMAP_source_data.csv"))

pca <- read_csv(file.path(NODE03, "source_data", "Fig3A_bulk_PCA_source_data.csv"))
bulk_heat <- read_csv(file.path(NODE03, "source_data", "Fig3C_bulk_top_DE_heatmap_source_data.csv"))
bulk_counts <- read_csv(file.path(NODE03, "source_data", "Fig3D_bulk_DE_gene_counts_source_data.csv"))
pvals <- read_csv(file.path(NODE03, "source_data", "FigS10_bulk_pvalue_distribution_source_data.csv"))
method_corr <- read_csv(file.path(NODE04, "source_data", "ED_Fig2_method_correlation_source_data.csv"))

copy_source(file.path(NODE03, "source_data", "Fig3B_bulk_volcano_source_data.csv"),
            "Figure_3A_bulk_volcano_source_data.csv")
copy_source(file.path(NODE03, "source_data", "FigS11_bulk_scRNA_celltype_marker_overlap_source_data.csv"),
            "Figure_3B_bulk_celltype_marker_overlap_source_data.csv")
copy_source(file.path(NODE04, "source_data", "Fig1a_AD_signature_violin_source_data.csv"),
            "Figure_3D_AD_signature_celltype_scores_source_data.csv")
copy_source(file.path(NODE04, "source_data", "Fig1c_wilcoxon_delta_source_data.csv"),
            "Figure_3E_AD_signature_celltype_delta_source_data.csv")
copy_source(file.path(NODE04, "source_data", "Fig1d_AD_score_UMAP_source_data.csv"),
            "Figure_3F_AD_signature_UMAP_source_data.csv")
copy_source(file.path(NODE04, "source_data", "Fig1b_method_heatmap_source_data.csv"),
            "Figure_3C_method_score_matrix_source_data.csv")
copy_source(file.path(NODE04, "source_data", "ED_Fig2_method_correlation_source_data.csv"),
            "Supplementary_Figure_5B_method_correlation_source_data.csv")

volcano$category2 <- ifelse(volcano$category %in% c("Up in AD", "Down in AD"), volcano$category, "NS")
volcano$category2 <- factor(volcano$category2, levels = c("Down in AD", "NS", "Up in AD"))
volcano$neg_log10_fdr <- pmin(volcano$neg_log10_fdr, 16)
volcano$label <- ""
for (cat in c("Up in AD", "Down in AD")) {
  idx <- which(volcano$category2 == cat)
  idx <- idx[order(volcano$`adj.P.Val`[idx], -abs(volcano$logFC[idx]))]
  keep <- head(idx, 4)
  volcano$label[keep] <- volcano$symbol[keep]
}

pA <- ggplot(volcano, aes(logFC, neg_log10_fdr)) +
  geom_point(aes(colour = category2), size = 0.35, alpha = 0.75) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", linewidth = 0.35) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", linewidth = 0.35) +
  ggrepel::geom_text_repel(
    data = subset(volcano, label != ""),
    aes(label = label),
    size = 2.1,
    fontface = "bold",
    min.segment.length = 0,
    max.overlaps = 40,
    box.padding = 0.22,
    segment.size = 0.25
  ) +
  scale_colour_manual(values = de_cols, labels = c("Down in AD", "NS", "Up in AD"), name = NULL, drop = FALSE) +
  labs(x = "log2 fold change (AD/CN)", y = "-log10(FDR)") +
  guides(colour = "none") +
  theme_publication(7.2) +
  theme(legend.position = "none")

overlap$celltype <- factor(overlap$celltype, levels = rev(ct_levels))
overlap$direction <- factor(overlap$direction, levels = c("Up in AD", "Down in AD"))
overlap$neg_log10_fdr_cap <- pmin(overlap$neg_log10_fdr, 7.2)

pB <- ggplot(overlap, aes(direction, celltype)) +
  geom_point(aes(size = overlap_genes, fill = neg_log10_fdr_cap),
             shape = 21, colour = "black", stroke = 0.18, alpha = 0.92) +
  scale_y_discrete(labels = function(x) sub("Cerebrovascular cells", "Cerebrovascular\ncells", x, fixed = TRUE)) +
  scale_fill_gradient(low = "#E5E5E5", high = "#B2182B", name = "-log10(FDR)") +
  scale_size_area(max_size = 6.4, name = "Overlap genes", breaks = c(2, 6, 10, 12)) +
  labs(x = NULL, y = NULL) +
  guides(
    fill = guide_colorbar(title.position = "top", barwidth = unit(16, "mm"), barheight = unit(2.2, "mm")),
    size = guide_legend(title.position = "top", nrow = 1, override.aes = list(fill = "white"))
  ) +
  theme_publication(7.2) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.key.size = unit(2.8, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

score_violin$celltype <- factor(score_violin$celltype, levels = ct_levels)
score_heat$celltype <- factor(score_heat$celltype, levels = ct_levels)
score_heat$method <- factor(score_heat$method,
                            levels = c("AddModuleScore", "AUCell", "UCell", "singscore", "ssGSEA", "Combined"))
pC <- ggplot(score_heat, aes(celltype, method, fill = score)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", score)), size = 2.0, fontface = "bold") +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B2182B",
                       midpoint = 0.50, name = "Mean\nscore") +
  labs(x = NULL, y = NULL) +
  theme_publication(7.0) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

pD <- ggplot(score_violin, aes(celltype, Combined, fill = celltype)) +
  geom_violin(width = 0.86, linewidth = 0.25, colour = "black", scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.25, fill = "white") +
  scale_fill_manual(values = ct_cols, guide = "none") +
  labs(x = NULL, y = "AD-up signature score") +
  theme_publication(7.2) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

score_delta$celltype <- factor(score_delta$celltype, levels = score_delta$celltype[order(score_delta$delta)])
pE <- ggplot(score_delta, aes(delta, celltype, fill = celltype)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  scale_fill_manual(values = ct_cols, guide = "none") +
  labs(x = "Median difference vs other cells", y = NULL) +
  theme_publication(7.2)

score_umap$celltype <- factor(score_umap$celltype, levels = ct_levels)
score_umap_plot <- score_umap %>%
  group_by(celltype) %>%
  group_modify(~ .x[sample(seq_len(nrow(.x)), min(nrow(.x), 16000)), , drop = FALSE]) %>%
  ungroup()
pF <- ggplot(score_umap_plot, aes(UMAP_1, UMAP_2, colour = Combined)) +
  geom_point(size = 0.035, alpha = 0.75) +
  scale_colour_gradientn(
    colours = c("#2D004B", "#542788", "#F1A340", "#FFF176"),
    limits = quantile(score_umap$Combined, c(0.01, 0.99), na.rm = TRUE),
    oob = squish,
    name = "AD-up\nscore"
  ) +
  coord_fixed(ratio = 0.68) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(colour = guide_colorbar(title.position = "top", barwidth = unit(18, "mm"), barheight = unit(2.4, "mm"))) +
  theme_publication(7.2) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

main_fig <- ((panel_label(pA, "A") | panel_label(pB, "B")) /
             (panel_label(pC, "C") | panel_label(pF, "D"))) +
  plot_layout(heights = c(1, 1.05), guides = "keep") &
  theme(plot.margin = margin(3, 3, 3, 3, unit = "mm"))

save_publication(main_fig, file.path(DIR_MAIN, "Figure_3_bulk_AD_signature_celltype_mapping"),
         width_cm = 17.0, height_cm = 14.2, dpi = 600)

## Supplementary Figure 4: bulk differential-expression evidence and diagnostics.
pca$group <- factor(pca$group, levels = c("CN", "AD"))
pS4A <- ggplot(pca, aes(PC1, PC2, colour = group)) +
  stat_ellipse(linewidth = 0.35, alpha = 0.8, show.legend = FALSE) +
  geom_point(size = 1.45, alpha = 0.85) +
  scale_colour_manual(values = group_cols, name = NULL) +
  labs(x = "PC1", y = "PC2") +
  theme_publication(7.2) +
  theme(legend.position = "right")

bulk_heat$group <- factor(bulk_heat$group, levels = c("AD", "CN"))
bulk_heat$sample_id <- factor(bulk_heat$sample_id, levels = unique(bulk_heat$sample_id[order(bulk_heat$group)]))
bulk_heat$symbol <- factor(bulk_heat$symbol, levels = rev(unique(bulk_heat$symbol)))
pS4B <- ggplot(bulk_heat, aes(sample_id, symbol, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                       midpoint = 0, limits = c(-2.5, 2.5), oob = squish, name = "z-score") +
  labs(x = NULL, y = NULL) +
  theme_publication(6.0) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6.0, lineheight = 0.95, margin = margin(r = 1.5)),
    legend.position = "right"
  )

bulk_counts$category <- factor(bulk_counts$category, levels = c("Up in AD", "Down in AD"))
pS4C <- ggplot(bulk_counts, aes(category, n, fill = category)) +
  geom_col(width = 0.58, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = n), vjust = -0.35, size = 2.5, fontface = "bold") +
  scale_fill_manual(values = c("Up in AD" = "#D55E00", "Down in AD" = "#0072B2"), guide = "none") +
  labs(x = NULL, y = "Genes, n") +
  theme_publication(7.2) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

pS4D <- ggplot(pvals, aes(P.Value)) +
  geom_histogram(bins = 50, fill = "#7A7A7A", colour = "white", linewidth = 0.15) +
  scale_x_continuous(labels = label_number(accuracy = 0.1)) +
  labs(x = "Nominal P value", y = "Genes, n") +
  theme_publication(7.2)

s4 <- ((pS4A | pS4C) / (pS4B | pS4D)) +
  plot_layout(widths = c(1, 1), heights = c(0.72, 1.55), guides = "keep") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(family = "Arial", face = "bold", size = 12),
    plot.tag.position = c(0.00, 1.08),
    plot.margin = margin(7, 3, 3, 3, unit = "mm")
  )
save_publication(s4, file.path(DIR_SUPP, "Supplementary_Figure_4_bulk_signature_QC"),
         width_cm = 17.0, height_cm = 16.2, dpi = 600)

## Supplementary Figure 5: scoring robustness.
method_corr$m1 <- factor(method_corr$m1, levels = c("AddModuleScore", "AUCell", "UCell", "singscore", "ssGSEA"))
method_corr$m2 <- factor(method_corr$m2, levels = c("AddModuleScore", "AUCell", "UCell", "singscore", "ssGSEA"))
pS5A <- ggplot(method_corr, aes(m1, m2, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.1, fontface = "bold") +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B2182B",
                       midpoint = 0.50, limits = c(0, 1), name = "Spearman\nrho") +
  labs(x = NULL, y = NULL) +
  theme_publication(7.2) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

set.seed(20260531)
per_method_path <- file.path(NODE04, "source_data", "ED_Fig1_per_method_violin_source_data.csv")
per_method <- read.csv(per_method_path, check.names = FALSE)
per_method$celltype <- factor(per_method$celltype, levels = ct_levels)
per_method$method <- factor(per_method$method,
                            levels = c("AddModuleScore", "AUCell", "UCell", "singscore", "ssGSEA"))
per_method_long <- per_method %>%
  group_by(celltype, method) %>%
  group_modify(~ .x[sample(seq_len(nrow(.x)), min(nrow(.x), 3500)), , drop = FALSE]) %>%
  ungroup()
pS5B <- ggplot(per_method_long, aes(value, celltype, fill = celltype)) +
  geom_violin(scale = "width", linewidth = 0.18, colour = "black") +
  scale_fill_manual(values = ct_cols, guide = "none") +
  facet_wrap(~ method, ncol = 1) +
  labs(x = "Normalized score", y = NULL) +
  theme_publication(6.2) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

s5 <- ((panel_label(pS5A, "A") | panel_label(pD, "C")) /
       (panel_label(pS5B, "B") | panel_label(pE, "D"))) +
  plot_layout(heights = c(0.82, 1.18), widths = c(1.15, 0.85), guides = "keep") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_blank(),
    plot.margin = margin(3, 3, 3, 3, unit = "mm")
  )
save_publication(s5, file.path(DIR_SUPP, "Supplementary_Figure_5_AD_signature_scoring_robustness"),
         width_cm = 17.0, height_cm = 18.4, dpi = 600)

## Decision log.
decision_log <- c(
  "Figure 3 contract: demonstrate that the bulk AD-up transcriptional signature is not an isolated bulk result and can be localized to major cell compartments in the snRNA-seq atlas.",
  "Main Figure 3 panels: A bulk AD/CN volcano defining the AD-up gene set; B Fisher overlap of bulk DE genes with major-cell-type marker programs; C five-method score matrix across major cell classes; D UMAP projection of the combined AD-up signature.",
  "Excluded from main Figure 3: bulk PCA, top-DE heatmap, DE gene counts and P-value distribution, because these support bulk-data quality rather than the cell-type-mapping claim.",
  "Supplementary Figure 4: bulk sample structure and bulk DE diagnostics.",
  "Supplementary Figure 5: method-level score robustness, inter-method correlation, and the downshifted cell-level combined-score distribution and median-difference panels.",
  "Figure labels use harmonized AD and control terminology, with cerebrovascular cells represented as a major cell class.",
  "R backend used for all Figure 3 drawing and export, following the selected project plotting workflow."
)
writeLines(decision_log, file.path(DIR_LOG, "Figure3_bulk_signature_mapping_decision_log.txt"))

message("Figure 3 MANUSCRIPT outputs written to: ", OUT)
