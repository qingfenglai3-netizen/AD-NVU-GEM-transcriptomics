#!/usr/bin/env Rscript
# node05 MiloR differential abundance analysis for AD-NVU
# publication-ready first-level cell-type analysis only.

.libPaths(unique(c(
  "/path/to/R_libraries",
  .libPaths()
)))

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(miloR)
  library(SingleCellExperiment)
  library(Matrix)
  library(scales)
  library(ggrepel)
  library(ggpubr)
})

set.seed(20260528)
options(stringsAsFactors = FALSE)

IN_QS <- "/path/to/project/results/02_annotation/rds/scRNA_annotated.qs"
OUT <- "/path/to/project/results/05_major_celltype_milor"

DIRS <- list(
  figs = file.path(OUT, "main_figures"),
  supp = file.path(OUT, "supplementary_figures"),
  tabs = file.path(OUT, "tables"),
  rds = file.path(OUT, "rds"),
  logs = file.path(OUT, "logs"),
  source = file.path(OUT, "source_data"),
  final = file.path(OUT, "outputs")
)
DIRS$final_main <- file.path(DIRS$final, "main_figures")
DIRS$final_supp <- file.path(DIRS$final, "supplementary_figures")
DIRS$final_edit <- file.path(DIRS$final, "editable_panels")
DIRS$final_tabs <- file.path(DIRS$final, "supplementary_tables")
DIRS$final_methods <- file.path(DIRS$final, "methods_results")
DIRS$final_source <- file.path(DIRS$final, "source_data")
DIRS$final_audit <- file.path(DIRS$final, "audit")
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

LOG_FILE <- file.path(DIRS$logs, "05_major_celltype_milor_final_regen.log")
log_msg <- function(...) {
  msg <- paste0(...)
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

save_plot_all <- function(plot, stem, width, height, dirs = c(DIRS$figs, DIRS$final_edit), dpi = 600) {
  for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(d, paste0(stem, ".pdf")), plot, width = width, height = height,
           device = cairo_pdf, limitsize = FALSE)
    ggsave(file.path(d, paste0(stem, ".png")), plot, width = width, height = height,
           dpi = dpi, bg = "white", limitsize = FALSE)
    ggsave(file.path(d, paste0(stem, ".tiff")), plot, width = width, height = height,
           dpi = dpi, compression = "lzw", bg = "white", limitsize = FALSE)
    try(suppressWarnings(ggsave(file.path(d, paste0(stem, ".svg")), plot,
                                width = width, height = height, limitsize = FALSE)), silent = TRUE)
  }
  log_msg("Saved plot: ", stem)
}

write_table <- function(x, stem, final = TRUE) {
  f1 <- file.path(DIRS$tabs, paste0(stem, ".csv"))
  write.csv(x, f1, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(x, file.path(DIRS$source, paste0(stem, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
  if (isTRUE(final)) {
    write.csv(x, file.path(DIRS$final_tabs, paste0(stem, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
    write.csv(x, file.path(DIRS$final_source, paste0(stem, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
  }
  log_msg("Saved table: ", stem)
}

theme_publication <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size - 1),
      panel.grid = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
}

standardize_group <- function(x) {
  x <- as.character(x)
  x[x %in% c("NC", "Normal", "Control", "CON")] <- "CN"
  x[x %in% c("UC", "AD", "Disease", "Case")] <- "AD"
  factor(x, levels = c("CN", "AD"))
}

standardize_celltype <- function(x) {
  x <- as.character(x)
  x[x %in% c("Endothelial", "Endothelial cells", "Vascular", "Vascular cells")] <- "Cerebrovascular cells"
  x
}

group_cols <- c(CN = "#4C78A8", AD = "#D55E00")
direction_cols <- c("AD-enriched" = "#B2182B", "CN-enriched" = "#2166AC", "Not significant" = "grey80")

log_msg("node05 MiloR final analysis started.")
log_msg("Input: ", IN_QS)

seu <- qread(IN_QS)
seu$Group <- standardize_group(seu$Group)
seu$celltype <- standardize_celltype(seu$celltype)
if (!"orig.ident" %in% colnames(seu@meta.data)) stop("orig.ident missing in metadata.")
if (!"umap" %in% names(seu@reductions)) stop("UMAP reduction missing.")
if (!"pca" %in% names(seu@reductions)) stop("PCA reduction missing.")

keep <- !is.na(seu$Group) & !is.na(seu$celltype)
seu <- subset(seu, cells = colnames(seu)[keep])
seu$celltype <- factor(seu$celltype)
seu$Group <- factor(as.character(seu$Group), levels = c("CN", "AD"))

log_msg("Loaded cells: ", ncol(seu), "; samples: ", length(unique(seu$orig.ident)),
        "; cell types: ", paste(levels(seu$celltype), collapse = ", "))

MAX_CELLS <- 50000
PROTECT_CT <- "Cerebrovascular cells"
FORCE_RECOMPUTE_DA <- TRUE
CKPT_SAMPLED <- file.path(DIRS$rds, "ckpt_sampled_ADCN_firstlevel_cerebro.qs")
CKPT_MILO <- file.path(DIRS$rds, "ckpt_milo_ADCN_firstlevel_cerebro.qs")
CKPT_DA <- file.path(DIRS$rds, "ckpt_da_ADCN_firstlevel_cerebro.qs")

if (file.exists(CKPT_SAMPLED)) {
  log_msg("Loading clean AD/CN sampled checkpoint.")
  seu_sub <- qread(CKPT_SAMPLED)
} else {
  log_msg("Sampling cells with protection for cerebrovascular cells.")
  meta <- seu@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(Group = as.character(Group), celltype = as.character(celltype))
  protected_cells <- meta$cell[meta$celltype %in% PROTECT_CT]
  remaining_budget <- max(1000, MAX_CELLS - length(protected_cells))
  non_protected <- meta %>% filter(!cell %in% protected_cells)
  strata <- non_protected %>%
    dplyr::count(orig.ident, Group, celltype, name = "n") %>%
    mutate(target = pmax(1, round(n / sum(n) * remaining_budget)))
  sampled <- protected_cells
  for (i in seq_len(nrow(strata))) {
    candidates <- non_protected %>%
      filter(orig.ident == strata$orig.ident[i], Group == strata$Group[i], celltype == strata$celltype[i]) %>%
      pull(cell)
    sampled <- c(sampled, sample(candidates, min(length(candidates), strata$target[i])))
  }
  sampled <- unique(sampled)
  if (length(sampled) > MAX_CELLS) {
    non_prot_sampled <- setdiff(sampled, protected_cells)
    sampled <- unique(c(protected_cells, sample(non_prot_sampled, MAX_CELLS - length(protected_cells))))
  }
  seu_sub <- subset(seu, cells = sampled)
  qsave(seu_sub, CKPT_SAMPLED, preset = "fast")
  sampling_info <- tibble(
    strategy = "stratified_by_sample_group_celltype_with_cerebrovascular_protection",
    total_cells = ncol(seu),
    sampled_cells = ncol(seu_sub),
    max_cells = MAX_CELLS,
    protected_celltype = PROTECT_CT,
    protected_cells = sum(seu$celltype %in% PROTECT_CT),
    sampled_protected_cells = sum(seu_sub$celltype %in% PROTECT_CT)
  )
  write_table(sampling_info, "Table_S18_MiloR_sampling_info_AD_CN")
}

seu_sub$Group <- factor(as.character(standardize_group(seu_sub$Group)), levels = c("CN", "AD"))
seu_sub$celltype <- factor(standardize_celltype(seu_sub$celltype))
sample_meta <- seu_sub@meta.data %>%
  distinct(orig.ident, Group) %>%
  mutate(Group = factor(as.character(Group), levels = c("CN", "AD")))
rownames(sample_meta) <- sample_meta$orig.ident

if (file.exists(CKPT_MILO)) {
  log_msg("Loading clean AD/CN Milo checkpoint.")
  milo <- qread(CKPT_MILO)
} else {
  log_msg("Building Milo object.")
  seu_sub <- JoinLayers(seu_sub)
  counts <- GetAssayData(seu_sub, assay = "RNA", layer = "counts")
  data <- GetAssayData(seu_sub, assay = "RNA", layer = "data")
  sce <- SingleCellExperiment(
    assays = list(counts = counts, logcounts = data),
    colData = seu_sub@meta.data
  )
  reducedDim(sce, "PCA") <- Embeddings(seu_sub, "pca")
  reducedDim(sce, "UMAP") <- Embeddings(seu_sub, "umap")
  milo <- Milo(sce)
  milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "PCA")
  milo <- makeNhoods(milo, prop = 0.10, k = 30, d = 30, refined = TRUE, reduced_dims = "PCA")
  milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample = "orig.ident")
  qsave(milo, CKPT_MILO, preset = "fast")
  log_msg("Milo neighborhoods: ", ncol(nhoods(milo)))
}

if (file.exists(CKPT_DA) && !FORCE_RECOMPUTE_DA) {
  log_msg("Loading clean AD/CN DA checkpoint.")
  ckpt <- qread(CKPT_DA)
  da <- ckpt$da
  nhood_meta <- ckpt$nhood_meta
  milo <- ckpt$milo
} else {
  log_msg("Testing differential abundance: AD versus CN.")
  design <- data.frame(
    Group = sample_meta[colnames(nhoodCounts(milo)), "Group"],
    row.names = colnames(nhoodCounts(milo))
  )
  design$Group <- factor(as.character(design$Group), levels = c("CN", "AD"))
  da <- testNhoods(milo, design = ~ Group, design.df = design, fdr.weighting = "none")
  da <- as.data.frame(da)
  da$Nhood <- seq_len(nrow(da))
  da$SpatialFDR_BH <- p.adjust(da$PValue, method = "BH")
  if (!"SpatialFDR" %in% colnames(da) || all(is.na(da$SpatialFDR))) {
    da$SpatialFDR <- da$SpatialFDR_BH
  }
  da$Direction <- case_when(
    da$SpatialFDR < 0.10 & da$logFC > 0 ~ "AD-enriched",
    da$SpatialFDR < 0.10 & da$logFC < 0 ~ "CN-enriched",
    TRUE ~ "Not significant"
  )
  da$FDR_bin <- cut(
    da$SpatialFDR,
    breaks = c(-Inf, 0.01, 0.05, 0.10, Inf),
    labels = c("FDR < 0.01", "FDR < 0.05", "FDR < 0.10", "Not significant")
  )

  nh <- nhoods(milo)
  cd <- as.data.frame(colData(milo))
  umap <- reducedDim(milo, "UMAP")
  dominant_ct <- character(ncol(nh))
  purity <- numeric(ncol(nh))
  size <- numeric(ncol(nh))
  nh_umap <- matrix(NA_real_, nrow = ncol(nh), ncol = 2)
  for (i in seq_len(ncol(nh))) {
    idx <- which(nh[, i] > 0)
    size[i] <- length(idx)
    cts <- as.character(cd$celltype[idx])
    tab <- sort(table(cts), decreasing = TRUE)
    dominant_ct[i] <- names(tab)[1]
    purity[i] <- as.numeric(tab[1]) / sum(tab)
    nh_umap[i, ] <- colMeans(umap[idx, , drop = FALSE])
  }
  da$dominant_celltype <- dominant_ct
  da$dominant_purity <- purity
  da$nhood_size <- size
  nhood_meta <- da %>%
    transmute(
      Nhood, UMAP1 = nh_umap[, 1], UMAP2 = nh_umap[, 2], nhood_size,
      dominant_celltype, dominant_purity, logFC, PValue, SpatialFDR,
      SpatialFDR_BH, Direction, FDR_bin
    )
  qsave(list(milo = milo, da = da, nhood_meta = nhood_meta), CKPT_DA, preset = "fast")
  write_table(da, "Table_S19_MiloR_neighborhood_DA_full_AD_vs_CN")
  write_table(nhood_meta, "Table_S20_MiloR_neighborhood_metadata_AD_vs_CN")
  log_msg("Significant neighborhoods FDR < 0.10: ", sum(da$SpatialFDR < 0.10, na.rm = TRUE),
          " / ", nrow(da))
}

sig_da <- da %>%
  filter(SpatialFDR < 0.10) %>%
  mutate(direction = if_else(logFC > 0, "AD-enriched", "CN-enriched"))

ct_summary <- da %>%
  mutate(direction = case_when(
    SpatialFDR < 0.10 & logFC > 0 ~ "AD-enriched",
    SpatialFDR < 0.10 & logFC < 0 ~ "CN-enriched",
    TRUE ~ "Not significant"
  )) %>%
  group_by(dominant_celltype) %>%
  summarise(
    total_neighborhoods = n(),
    significant_neighborhoods = sum(SpatialFDR < 0.10, na.rm = TRUE),
    AD_enriched = sum(direction == "AD-enriched", na.rm = TRUE),
    CN_enriched = sum(direction == "CN-enriched", na.rm = TRUE),
    mean_logFC_sig = if_else(significant_neighborhoods > 0, mean(logFC[SpatialFDR < 0.10], na.rm = TRUE), NA_real_),
    min_SpatialFDR = min(SpatialFDR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(significant_neighborhoods), min_SpatialFDR)
write_table(ct_summary, "Table_S21_MiloR_DA_by_major_celltype_AD_vs_CN")

cell_meta <- as.data.frame(colData(milo)) %>%
  rownames_to_column("cell") %>%
  mutate(Group = factor(as.character(Group), levels = c("CN", "AD")),
         celltype = factor(as.character(celltype)))

composition_sample <- cell_meta %>%
  dplyr::count(orig.ident, Group, celltype, name = "cells") %>%
  group_by(orig.ident) %>%
  mutate(prop = cells / sum(cells)) %>%
  ungroup()
write_table(composition_sample, "Table_S22_MiloR_sample_level_major_celltype_composition")

umap_cells <- as.data.frame(reducedDim(milo, "UMAP")) %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(celltype = as.character(colData(milo)$celltype), Group = as.character(colData(milo)$Group))

ct_labels <- umap_cells %>%
  group_by(celltype) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

top_labs <- da %>%
  filter(SpatialFDR < 0.10, Direction %in% c("AD-enriched", "CN-enriched")) %>%
  group_by(Direction, dominant_celltype) %>%
  arrange(SpatialFDR, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  group_by(Direction) %>%
  arrange(SpatialFDR, .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  ungroup() %>%
  mutate(label = paste0("N", Nhood, "\n", dominant_celltype))
top_labs <- top_labs %>%
  arrange(Direction, SpatialFDR) %>%
  group_by(Direction) %>%
  mutate(label_rank = row_number()) %>%
  ungroup() %>%
  mutate(
    label_x = case_when(
      Direction == "CN-enriched" ~ -3.65,
      Direction == "AD-enriched" ~ 2.15,
      TRUE ~ logFC
    ),
    label_y = case_when(
      Direction == "CN-enriched" ~ 7.95 - (label_rank - 1) * 0.90,
      Direction == "AD-enriched" ~ 6.65 - (label_rank - 1) * 0.90,
      TRUE ~ -log10(SpatialFDR)
    )
  )

panel_a <- ggplot() +
  geom_point(data = umap_cells, aes(UMAP1, UMAP2), colour = "grey88", size = 0.06, alpha = 0.55) +
  geom_point(data = nhood_meta %>% filter(Direction != "Not significant"),
             aes(UMAP1, UMAP2, colour = logFC, size = -log10(SpatialFDR)),
             alpha = 0.90) +
  geom_label_repel(data = ct_labels, aes(UMAP1, UMAP2, label = celltype),
                   seed = 20260528, size = 2.2, label.size = 0.15,
                   fill = "white", alpha = 0.92, max.overlaps = Inf, min.segment.length = 0) +
  scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "log2FC\nAD/CN") +
  scale_size_continuous(range = c(0.6, 2.8), name = "-log10\nMiloR FDR") +
  labs(title = "A  DA neighborhoods on major-cell UMAP", x = "UMAP 1", y = "UMAP 2") +
  theme_publication(8) +
  theme(legend.position = "right")

panel_b <- ggplot(da, aes(logFC, -log10(SpatialFDR))) +
  geom_point(aes(colour = Direction), size = 0.8, alpha = 0.72) +
  geom_hline(yintercept = -log10(0.10), linetype = "dashed", linewidth = 0.25, colour = "grey35") +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.25, colour = "grey50") +
  geom_segment(data = top_labs,
               aes(x = logFC, y = -log10(SpatialFDR), xend = label_x, yend = label_y),
               inherit.aes = FALSE, linewidth = 0.22, colour = "grey20") +
  geom_label(data = top_labs, aes(x = label_x, y = label_y, label = label),
             inherit.aes = FALSE, size = 1.55, label.size = 0.10,
             fill = "white", alpha = 0.95) +
  scale_colour_manual(values = direction_cols, name = NULL) +
  coord_cartesian(clip = "off") +
  labs(title = "B  Neighborhood-level DA statistics",
       x = "log2 fold-change (AD/CN)", y = "-log10 MiloR FDR") +
  theme_publication(8) +
  theme(plot.margin = margin(4, 14, 4, 4))

ct_long <- ct_summary %>%
  select(dominant_celltype, AD_enriched, CN_enriched) %>%
  pivot_longer(c(AD_enriched, CN_enriched), names_to = "direction", values_to = "n") %>%
  mutate(direction = recode(direction, AD_enriched = "AD-enriched", CN_enriched = "CN-enriched"),
         dominant_celltype = factor(dominant_celltype, levels = rev(ct_summary$dominant_celltype))) %>%
  filter(n > 0)

panel_c <- ggplot(ct_long, aes(dominant_celltype, n, fill = direction)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = "grey25", linewidth = 0.15) +
  coord_flip() +
  scale_fill_manual(values = direction_cols[c("AD-enriched", "CN-enriched")], name = NULL) +
  labs(title = "C  Direction of significant neighborhoods",
       x = NULL, y = "Significant neighborhoods (MiloR FDR < 0.10)") +
  theme_publication(8) +
  theme(legend.position = "bottom")

comp_stats <- composition_sample %>%
  group_by(celltype) %>%
  summarise(
    p = tryCatch(wilcox.test(prop ~ Group)$p.value, error = function(e) NA_real_),
    y = max(prop, na.rm = TRUE) * 1.12,
    .groups = "drop"
  ) %>%
  mutate(
    label = case_when(
      is.na(p) ~ "NA",
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    y_text = y * 1.045
  )

panel_d <- ggplot(composition_sample, aes(Group, prop, fill = Group)) +
  geom_boxplot(width = 0.48, outlier.shape = NA, linewidth = 0.25, alpha = 0.80) +
  geom_point(position = position_jitter(width = 0.10, height = 0), size = 0.65, alpha = 0.75) +
  geom_segment(data = comp_stats, aes(x = 1, xend = 2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.25) +
  geom_text(data = comp_stats, aes(x = 1.5, y = y_text, label = label),
            inherit.aes = FALSE, size = 2.3, family = "Arial") +
  facet_wrap(~celltype, ncol = 4, scales = "free_y") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = group_cols, guide = "none") +
  labs(title = "D  Sample-level composition in MiloR input",
       x = NULL, y = "Cells per sample (%)") +
  theme_publication(8)

save_plot_all(panel_a, "Fig6A_MiloR_DA_UMAP_major_celltypes", 6.8, 4.8)
save_plot_all(panel_b, "Fig6B_MiloR_DA_volcano_AD_vs_CN", 4.6, 3.7)
save_plot_all(panel_c, "Fig6C_MiloR_direction_by_major_celltype", 4.8, 3.8)
save_plot_all(panel_d, "Fig6D_MiloR_sample_composition_major_celltypes", 7.0, 4.4)

figure_6 <- (
  (panel_a + theme(legend.position = "right")) /
    ((panel_b + theme(legend.position = "none")) | panel_c) /
    (panel_d + theme(legend.position = "none"))
) +
  plot_layout(heights = c(1.20, 0.95, 1.02))
for (d in c(DIRS$final_main, DIRS$figs)) {
  ggsave(file.path(d, "Figure_6_node05_MiloR_major_celltype_DA.pdf"), figure_6,
         width = 8.3, height = 10.8, device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(d, "Figure_6_node05_MiloR_major_celltype_DA.png"), figure_6,
         width = 8.3, height = 10.8, dpi = 600, bg = "white", limitsize = FALSE)
  ggsave(file.path(d, "Figure_6_node05_MiloR_major_celltype_DA.tiff"), figure_6,
         width = 8.3, height = 10.8, dpi = 600, compression = "lzw", bg = "white", limitsize = FALSE)
  try(suppressWarnings(ggsave(file.path(d, "Figure_6_node05_MiloR_major_celltype_DA.svg"),
                              figure_6, width = 8.3, height = 10.8, limitsize = FALSE)), silent = TRUE)
}

panel_s1 <- ggplot(nhood_meta, aes(nhood_size)) +
  geom_histogram(bins = 45, fill = "#4C78A8", colour = "white", linewidth = 0.15) +
  scale_x_log10() +
  labs(title = "A  Neighborhood size distribution", x = "Cells per neighborhood (log10)", y = "Neighborhoods") +
  theme_publication(8)

panel_s2 <- ggplot(da, aes(dominant_purity, -log10(SpatialFDR), colour = Direction)) +
  geom_point(size = 0.75, alpha = 0.65) +
  geom_hline(yintercept = -log10(0.10), linetype = "dashed", linewidth = 0.25) +
  scale_colour_manual(values = direction_cols, name = NULL) +
  labs(title = "B  Dominant-cell purity of neighborhoods", x = "Dominant-cell fraction", y = "-log10 MiloR FDR") +
  theme_publication(8)

panel_s3 <- ggplot() +
  geom_point(data = umap_cells, aes(UMAP1, UMAP2, colour = celltype), size = 0.05, alpha = 0.62) +
  labs(title = "C  Major-cell landscape used for MiloR", x = "UMAP 1", y = "UMAP 2") +
  theme_publication(8) +
  guides(colour = guide_legend(title = NULL, override.aes = list(size = 2, alpha = 1))) +
  theme(legend.position = "bottom")

figure_s15 <- (panel_s1 | panel_s2) / panel_s3 + plot_layout(heights = c(0.82, 1.18))
for (d in c(DIRS$final_supp, DIRS$supp)) {
  ggsave(file.path(d, "Figure_S15_node05_MiloR_quality_controls.pdf"), figure_s15,
         width = 8.3, height = 7.8, device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(d, "Figure_S15_node05_MiloR_quality_controls.png"), figure_s15,
         width = 8.3, height = 7.8, dpi = 600, bg = "white", limitsize = FALSE)
  ggsave(file.path(d, "Figure_S15_node05_MiloR_quality_controls.tiff"), figure_s15,
         width = 8.3, height = 7.8, dpi = 600, compression = "lzw", bg = "white", limitsize = FALSE)
  try(suppressWarnings(ggsave(file.path(d, "Figure_S15_node05_MiloR_quality_controls.svg"),
                              figure_s15, width = 8.3, height = 7.8, limitsize = FALSE)), silent = TRUE)
}

audit <- tibble(
  item = c(
    "Analysis scope",
    "Group labels",
    "Vascular label",
    "Checkpoint policy",
    "Main figure",
    "Supplementary figure",
    "Interpretation"
  ),
  status = c(
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS_PENDING_VISUAL_REVIEW",
    "PASS_PENDING_VISUAL_REVIEW",
    "PASS_PENDING_VISUAL_REVIEW"
  ),
  note = c(
    "node05 is restricted to first-level major cell types; subtype-level MiloR is reserved for 05_nvu_subtype_milor.R.",
    "NC/UC aliases were standardized to CN/AD and all DA labels use AD/CN.",
    "Endothelial aliases were standardized to Cerebrovascular cells.",
    "New AD/CN checkpoint names prevent reuse of stale NC/UC outputs.",
    "Figure 6 combines DA UMAP, volcano, cell-type directional summary and sample-level composition.",
    "Figure S15 provides neighborhood size, purity and input UMAP quality controls.",
    "The panel set tests neighborhood abundance shifts without forcing sharp subtype boundaries."
  )
)
write_table(audit, "Table_S23_node05_MiloR_quality_control", final = TRUE)

legend_text <- c(
  "Figure 6. Major cell-type differential abundance in AD by MiloR.",
  "A, UMAP of sampled single nuclei with significant MiloR neighborhoods overlaid and colored by log2 fold-change for AD versus CN. Labels indicate the first-level cell-type annotations used for the node05 analysis. B, Volcano plot of neighborhood-level differential abundance statistics; the dashed line indicates MiloR graph FDR = 0.10. C, Number of significant neighborhoods per dominant major cell type, separated by AD-enriched and CN-enriched directions. D, Sample-level composition of major cell types in the MiloR input, shown to separate neighborhood-based inference from raw abundance summaries. CN, cognitively normal control; AD, Alzheimer disease; DA, differential abundance; FDR, false discovery rate.",
  "",
  "Figure S15. MiloR quality-control summaries.",
  "A, Distribution of cells per neighborhood. B, Relationship between dominant-cell purity and neighborhood-level MiloR graph FDR. C, Major-cell UMAP landscape used as the input for neighborhood construction. The analysis is intentionally restricted to first-level annotations; subtype-level MiloR is addressed in the downstream NVU-focused MiloR node."
)
writeLines(legend_text, file.path(DIRS$final_methods, "node05_MiloR_figure_legends.txt"), useBytes = TRUE)

methods_text <- c(
  "Methods - MiloR major cell-type differential abundance",
  "",
  "Neighborhood-based differential abundance testing was performed with MiloR on the node02 annotated single-nucleus object. Group labels were harmonized to CN and AD, and the first-level vascular annotation was harmonized to Cerebrovascular cells. To preserve rare vascular populations while maintaining tractable runtime, cells were sampled in a stratified manner by sample, group and first-level cell type, with all Cerebrovascular cells retained when possible. A Milo object was constructed from RNA counts and log-normalized expression, using PCA coordinates for k-nearest-neighbor graph construction (k = 30, d = 30) and refined neighborhood sampling (prop = 0.10). Neighborhood cell counts were modeled by sample with group as the design term, testing AD versus CN. Neighborhoods were annotated by their dominant first-level cell type, and significance was defined at MiloR graph FDR < 0.10. This node was restricted to major cell classes; subtype-level MiloR was reserved for the downstream NVU-focused analysis.",
  "",
  "Results - major cell-type abundance remodeling",
  "",
  paste0("MiloR identified ", sum(da$SpatialFDR < 0.10, na.rm = TRUE), " significant neighborhoods among ", nrow(da), " tested neighborhoods at MiloR graph FDR < 0.10. Significant neighborhoods were summarized by dominant first-level cell type to identify the major cell classes contributing to AD-associated abundance shifts. This analysis complements the node04 subtype prioritization by providing a topology-aware, first-level differential-abundance view, without forcing subtype boundaries into this upstream node."),
  "",
  "Supplementary note",
  "",
  "The input UMAP, neighborhood size distribution and dominant-cell purity summaries are provided as quality controls. These outputs are intended to document that the reported DA signal reflects local graph neighborhoods and major cell-type context, not merely a relabeling of raw sample composition."
)
writeLines(methods_text, file.path(DIRS$final_methods, "node05_MiloR_methods_results_draft.txt"), useBytes = TRUE)

writeLines(capture.output(sessionInfo()), file.path(DIRS$logs, "sessionInfo_05_major_celltype_milor_final.txt"))
log_msg("node05 MiloR final analysis completed.")
