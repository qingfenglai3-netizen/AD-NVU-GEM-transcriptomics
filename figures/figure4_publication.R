################################################################################
# Rebuild MANUSCRIPT Figure 4 and related supplementary figures
# Claim: astrocytes, microglia and cerebrovascular cells can be decomposed into
# marker-supported second-level subtypes for downstream mechanism analyses.
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

theme_publication <- function(base_size = 6.0) {
  base_size <- 6.0
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial", face = "bold", colour = "black", size = base_size),
      axis.text = element_text(size = base_size, face = "bold", colour = "black"),
      axis.title = element_text(size = base_size, face = "bold", colour = "black"),
      axis.title.y = element_text(size = base_size, face = "bold", margin = margin(r = 3), vjust = 0.5),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      legend.text = element_text(size = base_size, face = "bold", colour = "black"),
      legend.title = element_text(size = base_size, face = "bold", colour = "black"),
      legend.key.size = unit(2.4, "mm"),
      strip.text = element_text(size = base_size, face = "bold", colour = "black"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.35),
      plot.title = element_text(size = base_size, face = "bold", hjust = 0.5, colour = "black"),
      panel.grid = element_blank(),
      legend.spacing.y = unit(0.6, "mm"),
      legend.box.spacing = unit(0.6, "mm"),
      plot.margin = margin(1.2, 1.2, 1.2, 1.2, unit = "mm")
    )
}

save_publication <- function(p, file, width_cm = 17.4, height_cm = 18.0, dpi = 600) {
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

copy_source <- function(from, name) {
  to <- file.path(DIR_SRC, name)
  file.copy(from, to, overwrite = TRUE)
  to
}

umap_file <- function(compartment) {
  file.path(NODE04, "source_data", paste0("Fig_", compartment, "_UMAP_source_data.csv"))
}
dot_file <- function(compartment) {
  file.path(NODE04, "source_data", paste0("Fig_", compartment, "_DotPlot_source_data.csv"))
}
comp_file <- function(compartment) {
  file.path(NODE04, "source_data", paste0("Fig_", compartment, "_Composition_source_data.csv"))
}
heat_file <- function(compartment) {
  file.path(NODE04, "source_data", paste0("Fig_", compartment, "_Heatmap_source_data.csv"))
}

cv_levels <- c("Arterial", "Capillary", "Pericyte", "SMC", "Venous")
astro_levels <- c("Homeostatic", "Intermediate", "Reactive")
micro_levels <- c("Homeostatic", "DAM")

## Manuscript-wide subtype palette.
## Low-saturation, colorblind-aware hues are grouped by biological family:
## vascular subtypes span ochre/teal/blue/violet; astrocyte states use a rose-lilac
## progression; microglial states use pale-to-deep orange. This avoids the previous
## red/green contrast and keeps subtype colors stable across main and supplementary figures.
cv_cols <- c("Arterial" = "#CC6677", "Capillary" = "#B79F00", "Pericyte" = "#44AA99",
             "SMC" = "#4477AA", "Venous" = "#AA4499")
astro_cols <- c("Homeostatic" = "#D8A7B1", "Intermediate" = "#B07AA1", "Reactive" = "#CC79A7")
micro_cols <- c("Homeostatic" = "#E6AB6F", "DAM" = "#D55E00")
group_cols <- c("CN" = "#4DBBD5", "AD" = "#E64B35")

make_umap <- function(df, levels, cols, title, max_n = 65000) {
  df$subtype <- factor(df$subtype, levels = levels)
  if (nrow(df) > max_n) {
    df <- df %>%
      group_by(subtype) %>%
      group_modify(~ .x[sample(seq_len(nrow(.x)), min(nrow(.x), ceiling(max_n / length(levels)))), , drop = FALSE]) %>%
      ungroup()
  }
  lab <- aggregate(cbind(UMAP_1, UMAP_2) ~ subtype, df, median, na.rm = TRUE)
  lab$label <- as.character(lab$subtype)
  if (title == "Cerebrovascular cells") {
    lab$dx <- c("Arterial" = 0.15, "Capillary" = -0.25, "Pericyte" = 0.00,
                "SMC" = 0.00, "Venous" = -0.80)[as.character(lab$subtype)]
    lab$dy <- c("Arterial" = 0.35, "Capillary" = 0.25, "Pericyte" = 0.10,
                "SMC" = 0.30, "Venous" = 0.15)[as.character(lab$subtype)]
  } else if (title == "Astrocytes") {
    lab$dx <- c("Homeostatic" = 0.50, "Intermediate" = 0.25, "Reactive" = -0.95)[as.character(lab$subtype)]
    lab$dy <- c("Homeostatic" = -0.75, "Intermediate" = 0.70, "Reactive" = 0.55)[as.character(lab$subtype)]
  } else {
    lab$dx <- c("Homeostatic" = 1.00, "DAM" = 0.25)[as.character(lab$subtype)]
    lab$dy <- c("Homeostatic" = -0.10, "DAM" = 0.10)[as.character(lab$subtype)]
  }
  lab$UMAP_1 <- lab$UMAP_1 + lab$dx
  lab$UMAP_2 <- lab$UMAP_2 + lab$dy
  ggplot(df, aes(UMAP_1, UMAP_2, colour = subtype)) +
    geom_point(size = 0.11, alpha = 0.72) +
    geom_label(
      data = lab,
      aes(UMAP_1, UMAP_2, label = label),
      inherit.aes = FALSE,
      family = "Arial", fontface = "bold", size = 2.0,
      label.size = 0.24, label.padding = unit(0.9, "mm"),
      fill = scales::alpha("white", 0.84), colour = "black"
    ) +
    scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
    scale_y_continuous(labels = function(x) sprintf("%3.0f", x)) +
    coord_cartesian(clip = "off") +
    labs(x = "UMAP 1", y = "UMAP 2", title = title) +
    theme_publication() +
    theme(legend.position = "none",
          aspect.ratio = 1,
          plot.margin = margin(1.0, 1.0, 1.0, 1.0, unit = "mm"))
}

make_dot <- function(df, levels, title, show_legend = FALSE) {
  feature_order <- df %>%
    mutate(feature.groups = factor(feature.groups, levels = levels)) %>%
    arrange(feature.groups, features.plot) %>%
    pull(features.plot) %>%
    unique()
  df$id <- factor(df$id, levels = levels)
  df$feature.groups <- factor(df$feature.groups, levels = levels)
  df$features.plot <- factor(df$features.plot, levels = feature_order)
  ggplot(df, aes(features.plot, id)) +
    geom_point(aes(size = pct.exp, fill = avg.exp.scaled),
               shape = 21, colour = "black", stroke = 0.15) +
    scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                         midpoint = 0, limits = c(-1.2, 1.6), oob = squish,
                         breaks = c(-1, 0, 1), name = "Scaled\nexpression",
                         guide = if (show_legend) guide_colorbar(barheight = unit(11, "mm"),
                                                                 barwidth = unit(3, "mm")) else "none") +
    scale_size_area(max_size = 3.4, name = "% cells", breaks = c(20, 50, 80),
                    guide = if (show_legend) guide_legend(override.aes = list(fill = "white")) else "none") +
    labs(x = NULL, y = NULL, title = title) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = if (show_legend) "right" else "none")
}

cv_umap <- read_csv(umap_file("Cerebrovascular cells"))
astro_umap <- read_csv(umap_file("Astrocytes"))
micro_umap <- read_csv(umap_file("Microglia"))

cv_dot <- read_csv(dot_file("Cerebrovascular cells"))
astro_dot <- read_csv(dot_file("Astrocytes"))
micro_dot <- read_csv(dot_file("Microglia"))

copy_source(umap_file("Cerebrovascular cells"), "Figure_4A_cerebrovascular_subtype_UMAP_source_data.csv")
copy_source(umap_file("Astrocytes"), "Figure_4B_astrocyte_subtype_UMAP_source_data.csv")
copy_source(umap_file("Microglia"), "Figure_4C_microglia_subtype_UMAP_source_data.csv")
copy_source(dot_file("Cerebrovascular cells"), "Figure_4D_cerebrovascular_marker_dotplot_source_data.csv")
copy_source(dot_file("Astrocytes"), "Figure_4E_astrocyte_marker_dotplot_source_data.csv")
copy_source(dot_file("Microglia"), "Figure_4F_microglia_marker_dotplot_source_data.csv")

make_comp <- function(df, levels, cols, title) {
  df$Group <- factor(df$Group, levels = c("CN", "AD"))
  df$subtype <- factor(df$subtype, levels = levels)
  ggplot(df, aes(Group, prop, fill = subtype)) +
    geom_col(width = 0.62, colour = "black", linewidth = 0.18) +
    scale_fill_manual(values = cols, name = NULL, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.03))) +
    labs(x = NULL, y = "Proportion", title = title) +
    theme_publication() +
    theme(legend.position = "none")
}

cv_comp <- read_csv(comp_file("Cerebrovascular cells"))
astro_comp <- read_csv(comp_file("Astrocytes"))
micro_comp <- read_csv(comp_file("Microglia"))
copy_source(comp_file("Cerebrovascular cells"), "Figure_4C_cerebrovascular_composition_source_data.csv")
copy_source(comp_file("Astrocytes"), "Figure_4F_astrocyte_composition_source_data.csv")
copy_source(comp_file("Microglia"), "Figure_4I_microglia_composition_source_data.csv")

pA <- make_umap(cv_umap, cv_levels, cv_cols, "Cerebrovascular cells")
pB <- make_dot(cv_dot, cv_levels, "Cerebrovascular markers")
pC <- make_comp(cv_comp, cv_levels, cv_cols, "Cerebrovascular composition")
pD <- make_umap(astro_umap, astro_levels, astro_cols, "Astrocytes")
pE <- make_dot(astro_dot, astro_levels, "Astrocyte markers")
pF <- make_comp(astro_comp, astro_levels, astro_cols, "Astrocyte composition")
pG <- make_umap(micro_umap, micro_levels, micro_cols, "Microglia")
pH <- make_dot(micro_dot, micro_levels, "Microglial markers", show_legend = TRUE)
pI <- make_comp(micro_comp, micro_levels, micro_cols, "Microglial composition")

main_fig <- (panel_label(pA, "A") | panel_label(pB, "B") | panel_label(pC, "C")) /
  (panel_label(pD, "D") | panel_label(pE, "E") | panel_label(pF, "F")) /
  (panel_label(pG, "G") | panel_label(pH, "H") | panel_label(pI, "I")) +
  plot_layout(widths = c(1.14, 1.22, 0.74), heights = c(1, 1, 0.92), guides = "keep")

save_publication(main_fig, file.path(DIR_MAIN, "Figure_4_AD_NVU_compartment_subtype_atlas"),
         width_cm = 17.0, height_cm = 17.4, dpi = 600)

## Supplementary Figure 6: cluster-level assignment and marker heat-map evidence.
make_assignment <- function(path, levels, title) {
  df <- read_csv(path)
  long <- df %>%
    select(cluster, assigned, confidence, all_of(levels)) %>%
    pivot_longer(cols = all_of(levels), names_to = "module", values_to = "score")
  long$cluster <- factor(long$cluster, levels = sort(unique(long$cluster)))
  long$module <- factor(long$module, levels = levels)
  ggplot(long, aes(cluster, module, fill = score)) +
    geom_tile(colour = "white", linewidth = 0.22) +
    geom_point(data = subset(long, as.character(module) == assigned),
               shape = 21, size = 1.6, stroke = 0.25, colour = "black", fill = "white") +
    scale_fill_gradient(low = "#F3F3F3", high = "#B2182B", name = "Module\nscore") +
    labs(x = "Cluster", y = NULL, title = title) +
    theme_publication(6.5) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}

copy_source(file.path(NODE04, "tables", "subtype_assignment_Cerebrovascular cells.csv"),
            "Supplementary_Figure_6A_cerebrovascular_assignment_source_data.csv")
copy_source(file.path(NODE04, "tables", "subtype_assignment_Astrocytes.csv"),
            "Supplementary_Figure_6B_astrocyte_assignment_source_data.csv")
copy_source(file.path(NODE04, "tables", "subtype_assignment_Microglia.csv"),
            "Supplementary_Figure_6C_microglia_assignment_source_data.csv")

s6A <- make_assignment(file.path(NODE04, "tables", "subtype_assignment_Cerebrovascular cells.csv"),
                       cv_levels, "Cerebrovascular assignment")
s6B <- make_assignment(file.path(NODE04, "tables", "subtype_assignment_Astrocytes.csv"),
                       astro_levels, "Astrocyte assignment")
s6C <- make_assignment(file.path(NODE04, "tables", "subtype_assignment_Microglia.csv"),
                       micro_levels, "Microglial assignment")

## Average marker heat maps by subtype.
make_heat <- function(path, levels, title) {
  df <- read_csv(path) %>%
    pivot_longer(cols = all_of(levels), names_to = "subtype", values_to = "z")
  df$subtype <- factor(df$subtype, levels = levels)
  df$gene <- factor(df$gene, levels = rev(unique(df$gene)))
  ggplot(df, aes(subtype, gene, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.20) +
    scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                         midpoint = 0, name = "z-score") +
    labs(x = NULL, y = NULL, title = title) +
    theme_publication(6.2) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

copy_source(heat_file("Cerebrovascular cells"), "Supplementary_Figure_7A_cerebrovascular_marker_heatmap_source_data.csv")
copy_source(heat_file("Astrocytes"), "Supplementary_Figure_7B_astrocyte_marker_heatmap_source_data.csv")
copy_source(heat_file("Microglia"), "Supplementary_Figure_7C_microglia_marker_heatmap_source_data.csv")

s7A <- make_heat(heat_file("Cerebrovascular cells"), cv_levels, "Cerebrovascular cells")
s7B <- make_heat(heat_file("Astrocytes"), astro_levels, "Astrocytes")
s7C <- make_heat(heat_file("Microglia"), c("DAM", "Homeostatic"), "Microglia")
s6 <- (panel_label(s6A, "A") | panel_label(s6B, "B") | panel_label(s6C, "C")) /
  (panel_label(s7A, "D") | panel_label(s7B, "E") | panel_label(s7C, "F")) +
  plot_layout(heights = c(0.78, 1.05), guides = "keep")
save_publication(s6, file.path(DIR_SUPP, "Supplementary_Figure_6_subtype_assignment_marker_heatmaps"),
         width_cm = 17.0, height_cm = 17.6, dpi = 600)

decision_log <- c(
  "Figure 4 contract: decompose astrocytes, microglia and cerebrovascular cells into marker-supported second-level subtypes that define objects for later mechanism analyses.",
  "Node boundary: AD signature mapping is Figure 3; differential abundance, MiloR, and bulk subtype validation are reserved for Figure 5. Figure 4 therefore avoids prioritization claims beyond subtype construction and marker support.",
  "Main Figure 4 panels: each row presents one compartment with subtype UMAP, canonical marker dotplot and descriptive AD/CN subtype composition.",
  "Supplementary Figure 6: cluster-level assignment module scores and average marker z-score heat maps by subtype, retained as supporting evidence for subtype annotation.",
  "Excluded from main Figure 4: bulk subtype validation and AD-vs-CN subtype boxplots, because these answer the next node's question about which subtypes/programs are AD-associated.",
  "Figure labels use harmonized AD and control terminology; formal titles are provided in the legends.",
  "R backend used for all Figure 4 drawing and export."
)
writeLines(decision_log, file.path(DIR_LOG, "Figure4_subtype_atlas_decision_log.txt"))

message("Figure 4 MANUSCRIPT outputs written to: ", OUT)
