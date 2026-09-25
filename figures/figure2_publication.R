################################################################################
# Rebuild MANUSCRIPT Figure 2 and related supplementary figures
# Claim: the integrated AD snRNA-seq atlas yields robust major-cell-type
# annotation and sample-level cell-type DE context without overextending into
# broad pathway analysis.
################################################################################

options(stringsAsFactors = FALSE)
set.seed(123456)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
})

ROOT <- "/path/to/project/results"
NODE01 <- file.path(ROOT, "01_atlas", "outputs")
NODE02 <- file.path(ROOT, "02_annotation")
OUT <- file.path(ROOT, "publication_outputs")
DIR_MAIN <- file.path(OUT, "main_figures")
DIR_SUPP <- file.path(OUT, "supplementary_figures")
DIR_TABLE <- file.path(OUT, "supplementary_tables")
DIR_LOG <- file.path(OUT, "decision_logs")
DIR_SRC <- file.path(OUT, "source_data")
DIR_SUPP_DRAFT <- file.path(OUT, "intermediate", "figure2_simplified_supp_drafts_not_for_submission")
for (d in c(DIR_MAIN, DIR_SUPP, DIR_TABLE, DIR_LOG, DIR_SRC, DIR_SUPP_DRAFT)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cm_to_in <- function(cm) cm / 2.54

theme_publication <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial", face = "bold", colour = "black"),
      axis.text = element_text(size = base_size, face = "bold", colour = "black"),
      axis.title = element_text(size = base_size, face = "bold", colour = "black"),
      axis.title.y = element_text(size = base_size, face = "bold", margin = margin(r = 1), vjust = 0.5),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      legend.text = element_text(size = base_size - 0.4, face = "bold", colour = "black"),
      legend.title = element_text(size = base_size - 0.2, face = "bold", colour = "black"),
      legend.key.size = unit(3.0, "mm"),
      strip.text = element_text(size = base_size, face = "bold", colour = "black"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.35),
      plot.title = element_text(size = base_size, face = "bold", hjust = 0.5, colour = "black"),
      panel.grid = element_blank(),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )
}

save_publication <- function(p, file, width_cm = 17, height_cm = 15, dpi = 600) {
  message("Saving: ", basename(file))
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
ct_short <- c(
  "Oligodendrocytes" = "Oligo",
  "Excitatory" = "Exc",
  "Inhibitory" = "Inh",
  "Astrocytes" = "Ast",
  "Microglia" = "Mic",
  "OPCs" = "OPC",
  "Cerebrovascular cells" = "Vasc"
)
group_cols <- c("CN" = "#4DBBD5", "AD" = "#E64B35")
dataset_cols <- c("D1" = "#00A087", "D2" = "#CC79A7")
de_cols <- c("AD-enriched" = "#D55E00", "CN-enriched" = "#0072B2", "Not significant" = "#C8C8C8")
de_labels <- c("AD-enriched" = "AD up", "CN-enriched" = "CN up", "Not significant" = "NS")

plot_points <- function(mapping, data, size = 0.03, alpha = 0.65) {
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::geom_point_rast(mapping = mapping, data = data, size = size, alpha = alpha, raster.dpi = 600)
  } else {
    geom_point(mapping = mapping, data = data, size = size, alpha = alpha)
  }
}

umap <- read_csv(file.path(NODE02, "source_data", "Fig1_UMAP_overview_source_data.csv"))
umap$celltype <- factor(umap$celltype, levels = ct_levels)
umap$group <- factor(umap$group, levels = c("CN", "AD"))
umap$dataset <- factor(umap$dataset, levels = c("D1", "D2"))

markers <- read_csv(file.path(NODE02, "source_data", "Fig2_Markers_dotplot_source_data.csv"))
markers$id <- paste0("g", markers$id)
markers$feature.groups <- factor(markers$feature.groups, levels = ct_levels)
markers$features.plot <- factor(markers$features.plot, levels = unique(markers$features.plot))

scores <- read_csv(file.path(NODE02, "source_data", "Fig3_Annotation_scores_source_data.csv"))
scores$cluster <- factor(scores$cluster, levels = unique(scores$cluster))
scores$celltype <- factor(scores$celltype, levels = ct_levels)
markers$id <- factor(markers$id, levels = levels(scores$cluster))

conf <- read_csv(file.path(NODE02, "source_data", "Fig3_Annotation_confidence_source_data.csv"))
conf$cluster <- factor(conf$cluster, levels = levels(scores$cluster))
conf$confidence <- factor(conf$confidence, levels = c("Low", "Medium", "High"))

comp_g <- read_csv(file.path(NODE02, "source_data", "Fig4_composition_group_source_data.csv"))
comp_d <- read_csv(file.path(NODE02, "source_data", "Fig4_composition_dataset_source_data.csv"))
comp_g$celltype <- factor(comp_g$celltype, levels = ct_levels)
comp_d$celltype <- factor(comp_d$celltype, levels = ct_levels)

de_path <- file.path(DIR_SRC, "Figure_2D_major_celltype_multivolcano_source_data.csv")
if (!file.exists(de_path)) {
  stop("Missing pseudo-bulk DE source data: ", de_path,
       "\nRun /path/to/code_workspace/run_figure2_major_celltype_pseudobulk_DE.R first.")
}
de <- read_csv(de_path)
de$celltype <- factor(de$celltype, levels = ct_levels)
LFC_CUTOFF <- 0.58
P_CUTOFF <- 0.05
de$direction2 <- ifelse(de$PValue < P_CUTOFF & de$logFC > LFC_CUTOFF, "AD-enriched",
                        ifelse(de$PValue < P_CUTOFF & de$logFC < -LFC_CUTOFF, "CN-enriched", "Not significant"))
de$direction2 <- factor(de$direction2, levels = names(de_cols))
de$neglog10P_capped <- pmin(-log10(pmax(de$PValue, .Machine$double.xmin)), 8)
de$ct_index <- as.numeric(de$celltype)
de$x_plot <- de$ct_index + pmax(pmin(de$logFC, 1.2), -1.2) * 0.23
de$is_fdr <- de$FDR < 0.05

is_labelable_gene <- function(x) {
  !grepl("^(AC[0-9]|AL[0-9]|AP[0-9]|LINC|MIR|SNOR|SNHG|RNU|RNA|LOC|CTD-|RP[0-9]|MT-)|(-AS[0-9]*$)|(^.*AS[0-9]+$)|\\.", x)
}
de$labelable_gene <- is_labelable_gene(de$gene)
label_candidates <- subset(de, PValue < P_CUTOFF & abs(logFC) > LFC_CUTOFF & labelable_gene)
label_list <- lapply(split(label_candidates, label_candidates$celltype), function(x) {
  x <- x[order(x$FDR, x$PValue, -abs(x$logFC)), ]
  head(x, 2)
})
label_de <- if (length(label_list)) do.call(rbind, label_list) else data.frame()
if (is.null(label_de) || !nrow(label_de)) {
  label_de <- do.call(rbind, lapply(split(de, de$celltype), function(x) head(x[order(x$PValue), ], 1)))
}

umap_labels <- aggregate(cbind(UMAP_1, UMAP_2) ~ celltype, umap, median, na.rm = TRUE)
umap_labels$label <- unname(ct_short[as.character(umap_labels$celltype)])
umap_labels$dx <- c(
  "Oligodendrocytes" = 0.00,
  "Excitatory" = -0.45,
  "Inhibitory" = 0.00,
  "Astrocytes" = 0.35,
  "Microglia" = 0.10,
  "OPCs" = 0.00,
  "Cerebrovascular cells" = 0.00
)[as.character(umap_labels$celltype)]
umap_labels$dy <- c(
  "Oligodendrocytes" = 0.00,
  "Excitatory" = 0.05,
  "Inhibitory" = 0.00,
  "Astrocytes" = 0.30,
  "Microglia" = 0.00,
  "OPCs" = 0.10,
  "Cerebrovascular cells" = 0.15
)[as.character(umap_labels$celltype)]
umap_labels$UMAP_1 <- umap_labels$UMAP_1 + umap_labels$dx
umap_labels$UMAP_2 <- umap_labels$UMAP_2 + umap_labels$dy

pA <- ggplot(umap, aes(UMAP_1, UMAP_2)) +
  plot_points(aes(colour = celltype), umap, size = 0.030, alpha = 0.70) +
  geom_label(
    data = umap_labels,
    aes(UMAP_1, UMAP_2, label = label),
    inherit.aes = FALSE,
    family = "Arial", fontface = "bold", size = 2.35,
    label.size = 0.28, label.padding = unit(1.15, "mm"),
    fill = scales::alpha("white", 0.82), colour = "black"
  ) +
  scale_colour_manual(values = ct_cols, drop = FALSE) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2", colour = NULL) +
  theme_publication(8) +
  theme(legend.position = "none",
        plot.margin = margin(1.5, 1.5, 1.5, 1.5, unit = "mm"))

comp_all <- rbind(
  data.frame(panel = "Group", label = comp_g$Group, celltype = comp_g$celltype, prop = as.numeric(comp_g$prop)),
  data.frame(panel = "Dataset", label = comp_d$Dataset, celltype = comp_d$celltype, prop = as.numeric(comp_d$prop))
)
comp_all$label <- factor(comp_all$label, levels = c("CN", "AD", "D1", "D2"))
comp_all$celltype <- factor(comp_all$celltype, levels = ct_levels)
pB <- ggplot(comp_all, aes(label, prop, fill = celltype)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = ct_cols, drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion", fill = NULL) +
  theme_publication(8) +
  theme(legend.position = "none", strip.text.x = element_text(size = 8, face = "bold"))

main_marker_genes <- c(
  "MBP", "PLP1", "SLC17A7", "CAMK2A", "GAD1", "GAD2", "AQP4", "GFAP",
  "C3", "P2RY12", "PDGFRA", "VCAN", "FLT1", "VWF"
)
markers_main <- markers[as.character(markers$features.plot) %in% main_marker_genes, ]
markers_main$features.plot <- factor(as.character(markers_main$features.plot), levels = main_marker_genes)
pC <- ggplot(markers_main, aes(features.plot, id)) +
  geom_point(aes(size = pct.exp, colour = avg.exp.scaled)) +
  facet_grid(. ~ feature.groups, scales = "free_x", space = "free_x",
             labeller = labeller(feature.groups = ct_short)) +
  scale_size(range = c(0.08, 2.15), breaks = c(0, 25, 50, 75, 100), name = "% expr.") +
  scale_colour_gradient2(low = "#F2F2F2", mid = "#FCAE91", high = "#CB181D",
                         midpoint = 0, name = "Avg. expr.") +
  labs(x = NULL, y = NULL) +
  theme_publication(8) +
  theme(strip.text.x = element_text(size = 6.7, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.2),
        axis.text.y = element_text(size = 6.7),
        legend.position = "right",
        panel.spacing.x = unit(0.8, "mm"))

pD <- ggplot(de, aes(x_plot, neglog10P_capped)) +
  geom_vline(xintercept = seq_along(ct_levels), colour = "grey82", linewidth = 0.25) +
  geom_hline(yintercept = -log10(P_CUTOFF), linetype = "dashed", colour = "grey45", linewidth = 0.3) +
  geom_point(aes(colour = direction2, alpha = direction2), size = 0.16, stroke = 0) +
  geom_point(data = subset(de, is_fdr), aes(x_plot, neglog10P_capped), shape = 21,
             fill = NA, colour = "black", size = 0.85, stroke = 0.25) +
  scale_colour_manual(values = de_cols, labels = de_labels, drop = FALSE) +
  scale_alpha_manual(values = c("AD-enriched" = 0.80, "CN-enriched" = 0.80, "Not significant" = 0.18), guide = "none") +
  scale_x_continuous(breaks = seq_along(ct_levels), labels = unname(ct_short[ct_levels]),
                     limits = c(0.62, length(ct_levels) + 0.38), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.08))) +
  labs(x = NULL, y = expression(-log[10](italic(P))), colour = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.8, alpha = 1))) +
  theme_publication(8) +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(3.2, "mm"),
        legend.text = element_text(size = 6.8, face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 7.0))
if (requireNamespace("ggrepel", quietly = TRUE) && nrow(label_de)) {
  pD <- pD +
    ggrepel::geom_text_repel(
      data = label_de,
      aes(x_plot, neglog10P_capped, label = gene),
      family = "Arial", fontface = "bold", size = 1.75,
      colour = "black", min.segment.length = 0,
      segment.size = 0.18, segment.colour = "grey45",
      max.overlaps = Inf, box.padding = 0.18, point.padding = 0.12,
      show.legend = FALSE
    )
}

fig2 <- (panel_label(pA, "A") | panel_label(pB, "B")) /
  (panel_label(pC, "C") | panel_label(pD, "D")) +
  plot_layout(widths = c(1.03, 0.97), heights = c(0.98, 1.02), guides = "keep")

save_publication(fig2, file.path(DIR_MAIN, "Figure_2_integrated_snRNA_atlas_annotation"),
         width_cm = 17, height_cm = 17.2)
save_publication(pA, file.path(DIR_MAIN, "Figure_2A_umap_annotation"), width_cm = 8.0, height_cm = 7.2)
save_publication(pB, file.path(DIR_MAIN, "Figure_2B_celltype_composition"), width_cm = 8.1, height_cm = 6.2)
save_publication(pC, file.path(DIR_MAIN, "Figure_2C_marker_support"), width_cm = 8.4, height_cm = 7.4)
save_publication(pD, file.path(DIR_MAIN, "Figure_2D_major_celltype_pseudobulk_multivolcano"), width_cm = 8.4, height_cm = 7.4)

# Supplementary Figure 1: atlas QC and integration support.
qc1 <- read_csv(file.path(NODE02, "source_data", "ED_Fig1_QC_per_cluster_source_data.csv"))
if (!("cluster" %in% colnames(qc1)) && "seurat_clusters" %in% colnames(qc1)) qc1$cluster <- qc1$seurat_clusters
qc_long <- reshape(qc1[, c("cluster", "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")],
                   varying = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
                   v.names = "value", timevar = "metric",
                   times = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
                   direction = "long")
qc_long$metric <- factor(qc_long$metric, levels = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"))
qc_long$cluster <- factor(qc_long$cluster, levels = sort(unique(qc_long$cluster)))
qc_summary <- aggregate(value ~ cluster + metric, qc_long, median, na.rm = TRUE)
qc_summary$metric <- factor(qc_summary$metric,
                            levels = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
                            labels = c("Genes", "UMIs", "Mito. %", "Ribo. %"))
qc_summary <- do.call(rbind, lapply(split(qc_summary, qc_summary$metric), function(x) {
  x$z <- as.numeric(scale(log1p(x$value)))
  x
}))
pS1A <- ggplot(qc_summary, aes(cluster, metric, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Median\n(z)") +
  labs(x = "Cluster", y = NULL) +
  theme_publication(8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.0),
        axis.text.y = element_text(size = 7.2),
        legend.position = "right")
pS1B <- ggplot(umap, aes(UMAP_1, UMAP_2)) +
  plot_points(aes(colour = group), umap, size = 0.018, alpha = 0.55) +
  scale_colour_manual(values = group_cols) +
  facet_wrap(~group, nrow = 1) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2", colour = NULL) +
  theme_publication(8) +
  theme(legend.position = "none")
pS1C <- ggplot(umap, aes(UMAP_1, UMAP_2)) +
  plot_points(aes(colour = dataset), umap, size = 0.018, alpha = 0.55) +
  scale_colour_manual(values = dataset_cols) +
  facet_wrap(~dataset, nrow = 1) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2", colour = NULL) +
  theme_publication(8) +
  theme(legend.position = "none")
supp1 <- pS1A | pS1B | pS1C +
  plot_layout(widths = c(1.05, 0.95, 0.95)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "black"))
save_publication(supp1, file.path(DIR_SUPP_DRAFT, "Supplementary_Figure_1_atlas_QC_integration"),
         width_cm = 17, height_cm = 6.8)

# Supplementary Figure 2: annotation score support, without duplicating main marker dot plot.
pS2A <- ggplot(scores, aes(celltype, cluster, fill = score)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#6A51A3", mid = "white", high = "#CB181D", midpoint = 0, name = "Score\n(z)") +
  scale_x_discrete(labels = ct_short) +
  labs(x = NULL, y = "Cluster") +
  theme_publication(8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.8),
        legend.position = "right")
conf$cluster <- factor(conf$cluster, levels = rev(levels(scores$cluster)))
pS2B <- ggplot(conf, aes(score_margin, cluster, colour = confidence)) +
  geom_segment(aes(x = 0, xend = score_margin, yend = cluster), linewidth = 0.35, colour = "grey70") +
  geom_point(size = 2.0) +
  scale_colour_manual(values = c("Low" = "#D73027", "Medium" = "#FDAE61", "High" = "#1A9850"), drop = FALSE) +
  labs(x = "Top score margin", y = NULL, colour = NULL) +
  theme_publication(8) +
  theme(axis.text.y = element_text(size = 6.8), legend.position = "top")
cons <- read_csv(file.path(NODE02, "source_data", "ED_Fig4_Annotation_consensus_source_data.csv"))
cons$celltype <- factor(cons$celltype, levels = ct_levels)
pS2C <- ggplot(cons, aes(mean_expr_z, module_z, colour = celltype)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.25) +
  geom_point(size = 1.25, alpha = 0.85) +
  scale_colour_manual(values = ct_cols, labels = ct_short, drop = FALSE) +
  labs(x = "Marker expression z", y = "Module score z", colour = NULL) +
  theme_publication(8) +
  theme(legend.position = "right")
supp2 <- (pS2A | (pS2B / pS2C)) +
  plot_layout(widths = c(1.05, 1.0)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "black"))
save_publication(supp2, file.path(DIR_SUPP_DRAFT, "Supplementary_Figure_2_annotation_scoring"),
         width_cm = 17, height_cm = 11.8)

# Supplementary Figure 3: pseudo-bulk DE details that complement, rather than repeat, main Figure 2D.
de_summary <- aggregate(gene ~ celltype, de, length)
colnames(de_summary)[2] <- "tested_genes"
de_summary$nominal_AD <- tapply(de$PValue < P_CUTOFF & de$logFC > LFC_CUTOFF, de$celltype, sum, na.rm = TRUE)[as.character(de_summary$celltype)]
de_summary$nominal_CN <- tapply(de$PValue < P_CUTOFF & de$logFC < -LFC_CUTOFF, de$celltype, sum, na.rm = TRUE)[as.character(de_summary$celltype)]
de_summary$FDR05_AD <- tapply(de$FDR < 0.05 & de$logFC > 0, de$celltype, sum, na.rm = TRUE)[as.character(de_summary$celltype)]
de_summary$FDR05_CN <- tapply(de$FDR < 0.05 & de$logFC < 0, de$celltype, sum, na.rm = TRUE)[as.character(de_summary$celltype)]
de_summary$celltype <- factor(de_summary$celltype, levels = ct_levels)
de_count <- reshape(
  de_summary[, c("celltype", "nominal_AD", "nominal_CN", "FDR05_AD", "FDR05_CN")],
  varying = c("nominal_AD", "nominal_CN", "FDR05_AD", "FDR05_CN"),
  v.names = "n_genes", timevar = "category",
  times = c("Nominal AD", "Nominal CN", "FDR AD", "FDR CN"),
  direction = "long"
)
de_count$category <- factor(de_count$category, levels = c("Nominal AD", "Nominal CN", "FDR AD", "FDR CN"))
pS3A <- ggplot(de_count, aes(celltype, n_genes, fill = category)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("Nominal AD" = "#F4A582", "Nominal CN" = "#92C5DE", "FDR AD" = "#D55E00", "FDR CN" = "#0072B2"),
                    labels = c("Nominal AD" = "AD nominal", "Nominal CN" = "CN nominal", "FDR AD" = "AD FDR", "FDR CN" = "CN FDR")) +
  scale_x_discrete(labels = ct_short) +
  labs(x = NULL, y = "Genes", fill = NULL) +
  theme_publication(8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")
top_fdr <- subset(de, FDR < 0.05 & labelable_gene)
top_fdr <- do.call(rbind, lapply(split(top_fdr, top_fdr$celltype), function(x) {
  x <- x[order(x$FDR, -abs(x$logFC)), ]
  head(x, 10)
}))
if (is.null(top_fdr) || !nrow(top_fdr)) {
  top_fdr <- subset(de, PValue < P_CUTOFF & abs(logFC) > LFC_CUTOFF & labelable_gene)
  top_fdr <- do.call(rbind, lapply(split(top_fdr, top_fdr$celltype), function(x) {
    x <- x[order(x$PValue, -abs(x$logFC)), ]
    head(x, 6)
  }))
  top_fdr$evidence <- "Nominal"
} else {
  top_fdr$evidence <- "FDR"
}
top_fdr$gene <- factor(top_fdr$gene, levels = rev(unique(top_fdr$gene[order(top_fdr$celltype, top_fdr$FDR)])))
top_fdr$neglog10FDR <- -log10(pmax(top_fdr$FDR, .Machine$double.xmin))
pS3B <- ggplot(top_fdr, aes(neglog10FDR, gene)) +
  geom_segment(aes(x = 0, xend = neglog10FDR, yend = gene), colour = "grey70", linewidth = 0.35) +
  geom_point(aes(fill = logFC), shape = 21, size = 2.2, colour = "black", stroke = 0.25) +
  facet_grid(celltype ~ ., scales = "free_y", space = "free_y", labeller = labeller(celltype = ct_short)) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, name = "log2FC") +
  labs(x = expression(-log[10](FDR)), y = NULL) +
  theme_publication(8) +
  theme(strip.text.y = element_text(angle = 0, size = 7.2),
        axis.text.y = element_text(size = 6.6),
        legend.position = "right")
supp3 <- pS3A / pS3B +
  plot_layout(heights = c(0.78, 1.22)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "black"))
save_publication(supp3, file.path(DIR_SUPP_DRAFT, "Supplementary_Figure_3_major_celltype_pseudobulk_DE"),
         width_cm = 17, height_cm = 14.5)
# Tables: descriptive cohort/QC and cell composition stay supplementary.
file.copy(file.path(NODE01, "main_tables", "Table_1_node01_cohort_QC_summary.csv"),
          file.path(DIR_TABLE, "Supplementary_Table_S1_dataset_QC_summary.csv"), overwrite = TRUE)
file.copy(file.path(NODE02, "outputs", "main_tables", "Table_2_node02_celltype_composition_summary.csv"),
          file.path(DIR_TABLE, "Supplementary_Table_S2_major_celltype_composition.csv"), overwrite = TRUE)
file.copy(file.path(OUT, "tables", "Table_Figure2_major_celltype_pseudobulk_DE_all.csv"),
          file.path(DIR_TABLE, "Supplementary_Table_S3_major_celltype_pseudobulk_DE_all.csv"), overwrite = TRUE)

writeLines(c(
  "Figure 2 MANUSCRIPT publication rebuild",
  "Main claim: integrated snRNA-seq atlas and major-cell-type annotation are technically defensible; cell-type DE is shown as sample-level pseudo-bulk context only.",
  "Panel A: integrated UMAP by seven major cell types.",
  "Panel B: major-cell-type composition by diagnostic group and source dataset.",
  "Panel C: canonical marker support, retained in the main figure as annotation evidence.",
  sprintf("Panel D: edgeR sample-level pseudo-bulk AD vs CN multi-volcano by seven major cell types; displayed coloured genes use P < %.2g and |log2FC| > %.2f, while black circles mark FDR < 0.05.", P_CUTOFF, LFC_CUTOFF),
  "Supplementary Figure 1: QC and annotation support, without duplicating the marker dot plot in main Figure 2C.",
  "Supplementary Figure 2: pseudo-bulk count summary and FDR-supported genes, without repeating the multi-volcano in main Figure 2D.",
  "Export-format details are retained here only as a production note and must not be included in the manuscript text."
), file.path(DIR_LOG, "Figure2_rebuild_notes.txt"))

