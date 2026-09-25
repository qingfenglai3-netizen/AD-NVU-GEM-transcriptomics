suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
})

root <- "/path/to/project/results"
manuscript_dir <- file.path(root, "publication_outputs")
out_main <- file.path(manuscript_dir, "main_figures")
out_source <- file.path(manuscript_dir, "source_data")
dir.create(out_main, recursive = TRUE, showWarnings = FALSE)
dir.create(out_source, recursive = TRUE, showWarnings = FALSE)

milo_path <- file.path(root, "05_nvu_subtype_milor", "outputs", "source_data", "Table_S27_NVU_MiloR_all_subtypes_including_non_significant.csv")
major_milo_path <- file.path(root, "05_major_celltype_milor", "outputs", "source_data", "Table_S21_MiloR_DA_by_major_celltype_AD_vs_CN.csv")
bulk_stats_path <- file.path(root, "04_subtypes", "tables", "bulk_subtype_AD_vs_CN_stats.csv")
bulk_long_path <- file.path(out_source, "Figure_5_bulk_ssGSEA_subtype_scores_long_source_data.csv")

theme_publication <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = min(base_size + 1, 8.0), hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = max(base_size - 0.5, 6.0), colour = "black"),
      legend.title = element_text(size = max(base_size - 0.5, 6.0)),
      legend.text = element_text(size = max(base_size - 1, 6.0)),
      axis.title.y = element_text(size = base_size, margin = margin(r = 1), vjust = 0.5),
      strip.background = element_rect(fill = "#F1F1F1", colour = "#D0D0D0", linewidth = 0.25),
      strip.text = element_text(size = max(base_size - 0.5, 6.0), face = "bold"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

panel_label <- function(p, lab) {
  label <- ggplot() +
    theme_void(base_family = "Arial") +
    annotate("text", x = 0, y = 0.5, label = lab, hjust = 0, vjust = 0.5,
             family = "Arial", fontface = "bold", size = 12 / .pt) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme(plot.margin = margin(0, 0, 0, 0))
  wrap_plots(label, wrap_elements(full = p), ncol = 1, heights = c(0.105, 1))
}

label_only <- function(lab) {
  ggplot() +
    theme_void() +
    labs(title = lab) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0, margin = margin(t = 1)))
}

ad_col <- "#B2182B"
cn_col <- "#2166AC"
ns_col <- "#BDBDBD"
pos_col <- "#D6604D"
neg_col <- "#4393C3"

milo <- read_csv(milo_path, show_col_types = FALSE) %>%
  mutate(
    subtype_program = paste0(Compartment, " | ", Subtype),
    display = case_when(
      Compartment == "Cerebrovascular cells" ~ paste0("Cerebrovascular | ", Subtype),
      TRUE ~ paste0(Compartment, " | ", Subtype)
    ),
    ad_n = AD_enriched,
    cn_n = CN_enriched,
    milo_direction = case_when(
      ad_n > cn_n & Significant_neighborhoods_FDR10 > 0 ~ "AD-enriched",
      cn_n > ad_n & Significant_neighborhoods_FDR10 > 0 ~ "CN-enriched",
      Significant_neighborhoods_FDR10 > 0 ~ "Mixed",
      TRUE ~ "Not significant"
    )
  )

major_milo <- read_csv(major_milo_path, show_col_types = FALSE) %>%
  mutate(
    dominant_celltype = if_else(dominant_celltype == "Cerebrovascular cells", "Cerebrovascular", dominant_celltype),
    dominant_celltype = factor(dominant_celltype,
                               levels = c("Oligodendrocytes", "Excitatory", "Astrocytes",
                                          "Cerebrovascular", "Microglia", "Inhibitory", "OPCs"))
  )

bulk_stats <- read_csv(bulk_stats_path, show_col_types = FALSE) %>%
  mutate(
    subtype_program = str_replace(geneset, "_", " | "),
    subtype_program = str_replace(subtype_program, "Cerebrovascular cells", "Cerebrovascular"),
    display = subtype_program,
    bulk_direction = case_when(
      padj < 0.05 & delta > 0 ~ "AD-enriched",
      padj < 0.05 & delta < 0 ~ "CN-enriched",
      padj < 0.10 & delta > 0 ~ "AD trend",
      padj < 0.10 & delta < 0 ~ "CN trend",
      TRUE ~ "Not significant"
    ),
    fdr_class = case_when(
      padj < 0.05 ~ "FDR < 0.05",
      padj < 0.10 ~ "FDR < 0.10",
      TRUE ~ "NS"
    )
  )

priority <- milo %>%
  select(Compartment, Subtype, display, Significant_neighborhoods_FDR10, AD_enriched, CN_enriched,
         Mean_logFC_significant, Min_SpatialFDR, milo_direction) %>%
  left_join(
    bulk_stats %>% mutate(join_display = display) %>%
      transmute(join_display, bulk_delta = delta, bulk_padj = padj, bulk_direction),
    by = c("display" = "join_display")
  ) %>%
  mutate(
    bulk_direction = replace_na(bulk_direction, "Not significant"),
    evidence_class = case_when(
      milo_direction == "AD-enriched" & bulk_direction == "AD-enriched" ~ "Convergent AD",
      milo_direction == "CN-enriched" & bulk_direction == "CN-enriched" ~ "Convergent CN",
      milo_direction == "AD-enriched" | bulk_direction == "AD-enriched" ~ "AD-supported",
      milo_direction == "CN-enriched" | bulk_direction == "CN-enriched" ~ "CN-supported",
      TRUE ~ "Not prioritized"
    ),
    priority_rank = case_when(
      evidence_class == "Convergent AD" ~ 1,
      evidence_class == "AD-supported" ~ 2,
      evidence_class == "Convergent CN" ~ 3,
      evidence_class == "CN-supported" ~ 4,
      TRUE ~ 5
    )
  ) %>%
  arrange(priority_rank, desc(abs(coalesce(Mean_logFC_significant, 0))), desc(abs(coalesce(bulk_delta, 0))))

write_csv(milo, file.path(out_source, "Figure_5A-C_MiloR_subtype_summary_source_data.csv"))
write_csv(major_milo, file.path(out_source, "Figure_5A_major_celltype_MiloR_summary_source_data.csv"))
write_csv(bulk_stats, file.path(out_source, "Figure_5D_bulk_subtype_ssGSEA_stats_source_data.csv"))
write_csv(priority, file.path(out_source, "Figure_5F_convergent_priority_source_data.csv"))

major_long <- major_milo %>%
  transmute(dominant_celltype, `AD-enriched` = AD_enriched, `CN-enriched` = CN_enriched) %>%
  pivot_longer(c(`AD-enriched`, `CN-enriched`), names_to = "Direction", values_to = "Neighborhoods")

pA <- ggplot(major_long, aes(dominant_celltype, Neighborhoods, fill = Direction)) +
  geom_col(width = 0.72, colour = "grey25", linewidth = 0.15) +
  scale_fill_manual(values = c("AD-enriched" = ad_col, "CN-enriched" = cn_col)) +
  labs(x = NULL, y = "Major-cell DA\nneighborhoods") +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 1), vjust = 0.5),
    legend.position = "none",
    plot.margin = margin(2, 1, 2, 0)
  )

sub_long <- milo %>%
  filter(Significant_neighborhoods_FDR10 > 0) %>%
  select(display, `AD-enriched` = AD_enriched, `CN-enriched` = CN_enriched) %>%
  pivot_longer(c(`AD-enriched`, `CN-enriched`), names_to = "Direction", values_to = "Neighborhoods") %>%
  mutate(display = factor(display, levels = rev(unique(priority$display[priority$Significant_neighborhoods_FDR10 > 0]))))

pB <- ggplot(sub_long, aes(display, Neighborhoods, fill = Direction)) +
  geom_col(width = 0.62, colour = "grey25", linewidth = 0.12, position = "stack") +
  coord_flip() +
  scale_fill_manual(values = c("AD-enriched" = ad_col, "CN-enriched" = cn_col)) +
  guides(fill = guide_legend(title = NULL, nrow = 1, byrow = TRUE)) +
  labs(x = NULL, y = NULL) +
  theme_publication() +
  theme(
    legend.position = "none"
  )

umap_meta <- bind_rows(
  read_csv(file.path(root, "05_nvu_subtype_milor", "tables", "Table_Cerebrovascular_nhood_metadata_full.csv"), show_col_types = FALSE) %>%
    mutate(Compartment = "Cerebrovascular cells"),
  read_csv(file.path(root, "05_nvu_subtype_milor", "tables", "Table_Astrocytes_nhood_metadata_full.csv"), show_col_types = FALSE) %>%
    mutate(Compartment = "Astrocytes"),
  read_csv(file.path(root, "05_nvu_subtype_milor", "tables", "Table_Microglia_nhood_metadata_full.csv"), show_col_types = FALSE) %>%
    mutate(Compartment = "Microglia")
) %>%
  mutate(
    Compartment = factor(Compartment, levels = c("Cerebrovascular cells", "Astrocytes", "Microglia")),
    DA_direction = case_when(
      SpatialFDR < 0.10 & logFC > 0 ~ "AD-enriched",
      SpatialFDR < 0.10 & logFC < 0 ~ "CN-enriched",
      TRUE ~ "Not significant"
    ),
    DA_direction = factor(DA_direction, levels = c("AD-enriched", "CN-enriched", "Not significant")),
    point_size = if_else(SpatialFDR < 0.10, pmin(-log10(pmax(SpatialFDR, 1e-12)), 4), 0.45)
  )

write_csv(umap_meta, file.path(out_source, "Figure_5C_subtype_MiloR_UMAP_source_data.csv"))

plot_umap_compartment <- function(compartment_name, show_y = TRUE, show_legend = FALSE) {
  df <- filter(umap_meta, Compartment == compartment_name)
  ggplot(df, aes(UMAP1, UMAP2)) +
    geom_point(data = filter(df, DA_direction == "Not significant"),
               colour = "#D8D8D8", size = 0.18, alpha = 0.55) +
    geom_point(data = filter(df, DA_direction != "Not significant"),
               aes(fill = DA_direction, size = point_size),
               shape = 21, colour = "white", stroke = 0.08, alpha = 0.92) +
    scale_fill_manual(values = c("AD-enriched" = ad_col, "CN-enriched" = cn_col, "Not significant" = "#D8D8D8"),
                      name = NULL, drop = FALSE) +
    scale_size_continuous(range = c(1.0, 3.4), name = "-log10(FDR)") +
    labs(x = "UMAP 1", y = if (show_y) "UMAP 2" else NULL, title = as.character(compartment_name)) +
    coord_fixed() +
    theme_publication(base_size = 6.5) +
    theme(
      plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5),
      legend.position = if (show_legend) "right" else "none",
      legend.key.size = unit(3.0, "mm"),
      axis.text = element_text(size = 6.0),
      axis.title = element_text(size = 6.0),
      axis.text.y = if (show_y) element_text(size = 6.0) else element_blank(),
      axis.ticks.y = if (show_y) element_line(linewidth = 0.25) else element_blank()
    )
}

pB_labeled <- label_only("C") |
  (plot_umap_compartment("Cerebrovascular cells", TRUE, FALSE) |
     plot_umap_compartment("Astrocytes", FALSE, FALSE) |
     plot_umap_compartment("Microglia", FALSE, TRUE))
pB_labeled <- pB_labeled + plot_layout(widths = c(0.028, 1, 1, 1.12))

bubble_df <- milo %>%
  filter(Significant_neighborhoods_FDR10 > 0) %>%
  mutate(display = factor(display, levels = rev(unique(priority$display[priority$Significant_neighborhoods_FDR10 > 0]))))

pC <- ggplot(bubble_df, aes(Mean_logFC_significant, display)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey45") +
  geom_point(aes(size = Significant_neighborhoods_FDR10, fill = Mean_logFC_significant),
             shape = 21, colour = "white", stroke = 0.2, alpha = 0.95) +
  scale_fill_gradient2(low = cn_col, mid = "white", high = ad_col, midpoint = 0, name = "Mean log2FC\nAD/CN") +
  scale_size_area(max_size = 6.2, name = "DA neighborhoods") +
  labs(x = "Mean log2FC (AD/CN)", y = NULL) +
  theme_publication() +
  theme(legend.position = "none")

bulk_order <- bulk_stats %>% arrange(delta) %>% pull(display)
pD <- ggplot(bulk_stats %>% mutate(display = factor(display, levels = bulk_order)),
             aes(delta, display)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey45") +
  geom_segment(aes(x = 0, xend = delta, yend = display), colour = "grey70", linewidth = 0.35) +
  geom_point(aes(fill = delta, size = -log10(padj), shape = fdr_class),
             colour = "grey25", stroke = 0.2) +
  scale_fill_gradient2(low = cn_col, mid = "white", high = ad_col, midpoint = 0, name = "Bulk delta\nAD-CN") +
  scale_size_continuous(range = c(1.7, 4.8), name = "-log10(FDR)") +
  scale_shape_manual(values = c("FDR < 0.05" = 21, "FDR < 0.10" = 24, "NS" = 22), name = "Bulk ssGSEA") +
  labs(x = "Bulk ssGSEA score difference", y = NULL) +
  theme_publication() +
  theme(legend.position = "none")

bulk_long <- read_csv(bulk_long_path, show_col_types = FALSE) %>%
  mutate(display = str_replace(geneset, "_", " | "),
         display = str_replace(display, "Cerebrovascular cells", "Cerebrovascular"))

selected_bulk <- c("Astrocytes | Reactive", "Astrocytes | Intermediate", "Microglia | DAM",
                   "Cerebrovascular | Venous", "Cerebrovascular | Capillary", "Cerebrovascular | Pericyte")
short_bulk_labs <- c(
  "Astrocytes | Reactive" = "Reactive\nastrocytes",
  "Astrocytes | Intermediate" = "Intermediate\nastrocytes",
  "Microglia | DAM" = "DAM\nmicroglia",
  "Cerebrovascular | Venous" = "Venous\nvascular",
  "Cerebrovascular | Capillary" = "Capillary\nvascular",
  "Cerebrovascular | Pericyte" = "Pericyte\nvascular"
)

box_df <- bulk_long %>% filter(display %in% selected_bulk) %>%
  mutate(display = factor(short_bulk_labs[display], levels = unname(short_bulk_labs[selected_bulk])),
         group = factor(group, levels = c("CN", "AD")))

ann <- bulk_stats %>% filter(display %in% selected_bulk) %>%
  mutate(display = factor(short_bulk_labs[display], levels = unname(short_bulk_labs[selected_bulk])),
         label = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ "ns")) %>%
  group_by(display) %>% summarise(label = first(label), y = max(box_df$score[box_df$display == first(display)], na.rm = TRUE) * 1.04, .groups = "drop")

pE <- ggplot(box_df, aes(group, score, fill = group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.28) +
  geom_jitter(width = 0.11, size = 0.35, alpha = 0.65, colour = "grey10") +
  geom_text(data = ann, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 2.2) +
  facet_wrap(~display, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("CN" = "#4DBBD5", "AD" = "#E64B35")) +
  labs(x = NULL, y = "Bulk ssGSEA score") +
  theme_publication() +
  theme(legend.position = "none", strip.text = element_text(size = 6.0, lineheight = 0.88),
        axis.text.x = element_text(size = 6),
        axis.title.y = element_text(margin = margin(r = 1), vjust = 0.5),
        plot.margin = margin(2, 2, 2, 0))

matrix_df <- priority %>%
  transmute(display, evidence_class,
            MiloR = milo_direction,
            `Bulk ssGSEA` = bulk_direction) %>%
  pivot_longer(c(MiloR, `Bulk ssGSEA`), names_to = "Evidence", values_to = "Direction") %>%
  mutate(
    display = factor(display, levels = rev(priority$display)),
    Evidence_direction = case_when(
      Direction %in% c("AD-enriched", "AD trend") ~ "AD-associated",
      Direction %in% c("CN-enriched", "CN trend") ~ "CN-associated",
      Direction == "Mixed" ~ "Mixed/discordant",
      TRUE ~ "Not significant"
    ),
    Evidence_direction = factor(
      Evidence_direction,
      levels = c("AD-associated", "CN-associated", "Mixed/discordant", "Not significant")
    )
  )

pF <- ggplot(matrix_df, aes(Evidence, display, fill = Evidence_direction)) +
  geom_tile(colour = "white", linewidth = 0.4, width = 0.92, height = 0.82) +
  scale_fill_manual(values = c("AD-associated" = ad_col,
                               "CN-associated" = cn_col,
                               "Mixed/discordant" = "#B2ABD2",
                               "Not significant" = "#E5E5E5"),
                    drop = TRUE) +
  guides(fill = guide_legend(title = NULL, nrow = 1, byrow = TRUE,
                             override.aes = list(colour = NA))) +
  labs(x = NULL, y = NULL) +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "bottom",
    legend.key.size = unit(3.0, "mm"),
    legend.text = element_text(size = 6.3),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

final <- (
  (panel_label(pA, "A") | panel_label(pB, "B")) /
    (panel_label(pC, "C") | panel_label(pF, "D")) /
    panel_label(pE, "E")
) +
  plot_layout(heights = c(0.92, 0.96, 0.92), widths = c(0.84, 1.16))

base <- file.path(out_main, "Figure_5_AD_associated_subtype_program_prioritization")
ggsave(paste0(base, ".pdf"), final, width = 17.0, height = 17.2, units = "cm", device = cairo_pdf)
ggsave(paste0(base, ".svg"), final, width = 17.0, height = 17.2, units = "cm", device = svglite)
ggsave(paste0(base, ".png"), final, width = 17.0, height = 17.2, units = "cm", dpi = 600, device = ragg::agg_png)
ggsave(paste0(base, ".tiff"), final, width = 17.0, height = 17.2, units = "cm", dpi = 600, compression = "lzw", device = ragg::agg_tiff)

supp_pD <- panel_label(pD, "J")
supp_base <- file.path(manuscript_dir, "supplementary_figures", "Supplementary_Figure_8J_bulk_ssGSEA_delta_for_merge")
ggsave(paste0(supp_base, ".pdf"), supp_pD, width = 8.5, height = 8.0, units = "cm", device = cairo_pdf)
ggsave(paste0(supp_base, ".svg"), supp_pD, width = 8.5, height = 8.0, units = "cm", device = svglite)
ggsave(paste0(supp_base, ".png"), supp_pD, width = 8.5, height = 8.0, units = "cm", dpi = 600, device = ragg::agg_png)
ggsave(paste0(supp_base, ".tiff"), supp_pD, width = 8.5, height = 8.0, units = "cm", dpi = 600, compression = "lzw", device = ragg::agg_tiff)

message("Saved Figure 5 to: ", base)

