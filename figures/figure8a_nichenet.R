#!/usr/bin/env Rscript
# Figure 8A: glia-to-vascular NicheNet ligand activity and vascular target priority.

set.seed(20260608)
options(stringsAsFactors = FALSE)

paths_to_add <- c(
  "<R_SITE_LIBRARY>",
  "<R_USER_LIBRARY>",
  "<R_BASE_LIBRARY>"
)
for (p in paths_to_add) {
  if (dir.exists(p) && !(p %in% .libPaths())) .libPaths(c(p, .libPaths()))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(grid)
})

ROOT <- "<ANALYSIS_ROOT>"
JTM <- file.path(ROOT, "fuxian_test", "JTM_manuscript")
SRC_RESEL <- file.path(JTM, "sandbox", "figure8_target_reselection")
OUT <- "<WORKDIR>/figures/Fig8A_source"
OUT_FIG <- file.path(OUT, "figures")
OUT_SRC <- file.path(OUT, "source_data")
OUT_LOG <- file.path(OUT, "logs")
for (d in c(OUT, OUT_FIG, OUT_SRC, OUT_LOG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

safe_gene <- function(x) toupper(trimws(as.character(x)))

nature_theme <- function(base = 7) {
  theme_void(base_size = base, base_family = "Arial") +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 4, 5),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = base),
      legend.title = element_text(size = base, face = "bold"),
      legend.key.size = unit(0.32, "cm")
    )
}

save_multi <- function(plot, stem, width = 3.15, height = 2.65, dpi = 600) {
  ggsave(file.path(OUT_FIG, paste0(stem, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf,
         bg = "white")
  ggsave(file.path(OUT_FIG, paste0(stem, ".svg")), plot,
         width = width, height = height, units = "in", bg = "white")
  ggsave(file.path(OUT_FIG, paste0(stem, ".png")), plot,
         width = width, height = height, units = "in", dpi = dpi,
         bg = "white")
  ggsave(file.path(OUT_FIG, paste0(stem, ".tiff")), plot,
         width = width, height = height, units = "in", dpi = dpi,
         compression = "lzw", device = "tiff", bg = "white")
}

parse_support <- function(x) {
  if (is.na(x) || !nzchar(x)) return(tibble(ligand_clean = character(), regulatory_potential = numeric()))
  parts <- unlist(strsplit(x, ";", fixed = TRUE))
  tibble(raw = parts) |>
    mutate(
      ligand_clean = safe_gene(sub(":.*$", "", raw)),
      regulatory_potential = as.numeric(sub("^.*:", "", raw))
    ) |>
    filter(!is.na(regulatory_potential)) |>
    arrange(desc(regulatory_potential)) |>
    slice_head(n = 5)
}

ligand_activities <- read_csv(file.path(SRC_RESEL, "glia_to_vascular_ligand_activities.csv"),
                              show_col_types = FALSE) |>
  mutate(ligand_clean = safe_gene(test_ligand))

target_priority <- read_csv(file.path(SRC_RESEL, "glia_to_vascular_target_priority.csv"),
                            show_col_types = FALSE) |>
  mutate(target_clean = safe_gene(target_clean))

top_targets <- target_priority |>
  arrange(rank_in_vascular_receiver) |>
  slice_head(n = 8) |>
  mutate(
    target_label = target_clean,
    target_role = if_else(target_clean == "GEM", "Selected target", "Vascular AD-up target")
  ) |>
  arrange(desc(logFC), PValue) |>
  mutate(display_rank_by_logFC = row_number()) |>
  arrange(rank_in_vascular_receiver)

edge_tbl <- bind_rows(lapply(seq_len(nrow(top_targets)), function(i) {
  parse_support(top_targets$supporting_ligands_top5[i]) |>
    mutate(target_clean = top_targets$target_clean[i])
})) |>
  filter(target_clean %in% top_targets$target_clean)

ligands_from_edges <- edge_tbl |>
  group_by(ligand_clean) |>
  summarise(max_regulatory_potential = max(regulatory_potential, na.rm = TRUE),
            n_targets = n_distinct(target_clean), .groups = "drop") |>
  arrange(desc(n_targets), desc(max_regulatory_potential)) |>
  slice_head(n = 14)

top_ligands <- ligand_activities |>
  distinct(ligand_clean, .keep_all = TRUE) |>
  left_join(ligands_from_edges, by = "ligand_clean") |>
  mutate(
    source_class = case_when(
      source_families == "Astrocyte" ~ "Astrocyte ligand",
      source_families == "Microglia" ~ "Microglial ligand",
      TRUE ~ "Glial ligand"
    )
  ) |>
  filter(source_class %in% c("Astrocyte ligand", "Microglial ligand")) |>
  group_by(source_class) |>
  arrange(desc(aupr), desc(max_regulatory_potential), .by_group = TRUE) |>
  slice_head(n = 5) |>
  ungroup()

edge_tbl <- edge_tbl |>
  semi_join(top_ligands, by = "ligand_clean") |>
  semi_join(top_targets, by = "target_clean")

target_order <- top_targets |>
  arrange(display_rank_by_logFC) |>
  pull(target_clean)

astro_nodes <- top_ligands |>
  filter(source_class == "Astrocyte ligand") |>
  mutate(
    node = ligand_clean,
    label = ligand_clean,
    x = 0.18,
    y = seq(0.80, 0.58, length.out = n()),
    node_type = source_class,
    node_size = rescale(aupr, to = c(1.7, 3.0))
  )

micro_nodes <- top_ligands |>
  filter(source_class == "Microglial ligand") |>
  mutate(
    node = ligand_clean,
    label = ligand_clean,
    x = 0.18,
    y = seq(0.47, 0.22, length.out = n()),
    node_type = source_class,
    node_size = rescale(aupr, to = c(1.7, 3.0))
  )

ligand_nodes <- bind_rows(astro_nodes, micro_nodes)

target_nodes <- top_targets |>
  arrange(display_rank_by_logFC) |>
  mutate(
    node = target_clean,
    label = target_clean,
    x = 0.86,
    y = seq(0.80, 0.28, length.out = n()),
    node_type = target_role,
    node_size = rescale(abs(logFC), to = c(1.9, 3.3))
  )

edges_plot <- edge_tbl |>
  left_join(ligand_nodes |> select(ligand_clean = node, x1 = x, y1 = y, source_class = node_type), by = "ligand_clean") |>
  left_join(target_nodes |> select(target_clean = node, x2 = x, y2 = y), by = "target_clean") |>
  mutate(
    edge_width = rescale(regulatory_potential, to = c(0.08, 0.55)),
    edge_alpha = rescale(regulatory_potential, to = c(0.12, 0.45))
  )

write_csv(top_ligands, file.path(OUT_SRC, "Figure8A_glia_to_vascular_ligands_source_data.csv"))
write_csv(top_targets, file.path(OUT_SRC, "Figure8A_glia_to_vascular_targets_source_data.csv"))
write_csv(edges_plot, file.path(OUT_SRC, "Figure8A_glia_to_vascular_edges_source_data.csv"))

ligand_pal <- c(
  "Astrocyte ligand" = "#00A087",
  "Microglial ligand" = "#925E9F",
  "Glial ligand" = "#8A8A8A"
)
target_pal <- c(
  "Vascular AD-up target" = "#F39B7F",
  "Selected target" = "#E64B35"
)
target_shape <- c(
  "Vascular AD-up target" = 21,
  "Selected target" = 23
)

panel_body <- ggplot() +
  geom_curve(
    data = edges_plot,
    aes(x = x1 + 0.035, y = y1, xend = x2 - 0.035, yend = y2,
        linewidth = edge_width, alpha = edge_alpha, color = source_class),
    curvature = 0.12,
    lineend = "round"
  ) +
  scale_linewidth_identity() +
  scale_alpha_identity() +
  geom_point(
    data = ligand_nodes,
    aes(x = x, y = y, size = node_size, fill = node_type),
    shape = 21, color = "white", stroke = 0.30
  ) +
  geom_point(
    data = target_nodes,
    aes(x = x, y = y, size = node_size, fill = node_type, shape = node_type),
    color = "#4A4A4A", stroke = 0.36
  ) +
  geom_text(
    data = ligand_nodes,
    aes(x = x - 0.050, y = y, label = label),
    hjust = 1, vjust = 0.5, family = "Arial", fontface = "bold", size = 2.12, color = "#222222"
  ) +
  geom_text(
    data = target_nodes,
    aes(x = x + 0.045, y = y, label = label,
        fontface = if_else(node == "GEM", "bold", "plain")),
    hjust = 0, vjust = 0.5, family = "Arial", size = 2.12, color = "#222222"
  ) +
  geom_text(aes(x = 0.23, y = 0.94, label = "Prioritized glial ligands"),
            family = "Arial", fontface = "bold", size = 2.12, hjust = 0.5) +
  geom_text(aes(x = 0.86, y = 0.945, label = "Vascular targets"),
            family = "Arial", fontface = "bold", size = 2.12, hjust = 0.5) +
  geom_segment(aes(x = -0.125, xend = -0.125, y = 0.55, yend = 0.83),
               linewidth = 0.42, color = "#00A087", lineend = "round") +
  geom_segment(aes(x = -0.125, xend = -0.125, y = 0.19, yend = 0.50),
               linewidth = 0.42, color = "#925E9F", lineend = "round") +
  geom_text(aes(x = -0.160, y = 0.69, label = "Astrocyte"),
            angle = 90, family = "Arial", fontface = "bold", size = 2.12,
            color = "#00A087") +
  geom_text(aes(x = -0.160, y = 0.345, label = "Microglia"),
            angle = 90, family = "Arial", fontface = "bold", size = 2.12,
            color = "#925E9F") +
  geom_segment(aes(x = 0.38, xend = 0.46, y = 0.120, yend = 0.120),
               linewidth = 0.25, color = "#777777", alpha = 0.30) +
  geom_segment(aes(x = 0.48, xend = 0.56, y = 0.120, yend = 0.120),
               linewidth = 0.65, color = "#777777", alpha = 0.70) +
  geom_text(aes(x = 0.585, y = 0.120, label = "Regulatory potential"),
            family = "Arial", fontface = "bold", size = 2.12, hjust = 0, vjust = 0.5) +
  scale_fill_manual(values = c(ligand_pal, target_pal), guide = "none") +
  scale_color_manual(values = c(
    "Astrocyte ligand" = "#00A087",
    "Microglial ligand" = "#925E9F"
  ), guide = "none") +
  scale_shape_manual(values = target_shape, guide = "none") +
  scale_size_identity() +
  coord_cartesian(xlim = c(-0.135, 1.0), ylim = c(0.0, 1.0), clip = "off") +
  nature_theme(7)

panel <- panel_body +
  geom_text(aes(x = 0.005, y = 0.975, label = "A"),
            family = "Arial", fontface = "bold", size = 4.0, hjust = 0, vjust = 0.5)
panel_no_tag <- panel_body

save_multi(panel, "Figure8A_glia_to_vascular_NicheNet")
save_multi(panel_no_tag, "Figure8A_glia_to_vascular_NicheNet_no_panel_tag")


message(file.path(OUT_FIG, "Figure8A_glia_to_vascular_NicheNet.png"))
