#!/usr/bin/env Rscript

set.seed(20260529)
options(stringsAsFactors = FALSE)

paths_to_add <- "/path/to/R_libraries"
for (p in paths_to_add) {
  if (dir.exists(p) && !(p %in% .libPaths())) .libPaths(c(p, .libPaths()))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(readr)
  library(forcats)
})

OUT <- "/path/to/project/results/06_communication/outputs"
TAB <- file.path(OUT, "tables")
MAIN <- file.path(OUT, "main_figures")
SUPP <- file.path(OUT, "supplementary_figures")
SRC <- file.path(OUT, "source_data")
DOC <- file.path(OUT, "methods_results")
LEGACY <- file.path(OUT, "legacy_uncombined_panels_removed_from_submission")
for (d in c(MAIN, SUPP, SRC, DOC, LEGACY)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

archive_nonfinal_figures <- function(dir_path, expected_basenames, legacy_subdir) {
  if (!dir.exists(dir_path)) return(invisible(NULL))
  files <- list.files(dir_path, full.names = TRUE)
  if (length(files) == 0) return(invisible(NULL))
  keep <- basename(files) %in% expected_basenames
  to_archive <- files[!keep]
  if (length(to_archive) == 0) return(invisible(NULL))
  archive_dir <- file.path(LEGACY, legacy_subdir)
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  for (f in to_archive) {
    dest <- file.path(archive_dir, basename(f))
    if (file.exists(dest)) {
      ext <- tools::file_ext(dest)
      stem <- sub(paste0("\\.", ext, "$"), "", basename(dest))
      dest <- file.path(archive_dir, paste0(stem, "_archived_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext))
    }
    ok <- file.rename(f, dest)
    if (!ok) warning("Could not archive non-final figure: ", f)
  }
  invisible(NULL)
}

archive_nonfinal_figures(
  MAIN,
  paste0("Figure_8_node06_NVU_cell_cell_communication.", c("png", "pdf", "svg", "tiff")),
  "main_figures"
)
archive_nonfinal_figures(
  SUPP,
  c(
    paste0("Figure_S20_node06_cerebrovascular_functional_enrichment.", c("png", "pdf", "svg", "tiff")),
    paste0("Figure_S21_node06_glial_functional_enrichment.", c("png", "pdf", "svg", "tiff")),
    paste0("Figure_S22_node06_pericyte_SMC_module_statistics.", c("png", "pdf", "svg", "tiff"))
  ),
  "supplementary_figures"
)

cn_col <- "#4DBBD5"
ad_col <- "#E64B35"
grey_col <- "#8A8A8A"
okabe <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#F0E442")

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text = element_text(face = "bold", size = base_size),
      legend.key.size = unit(0.35, "cm"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.5),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      plot.margin = margin(4, 4, 4, 4)
    )
}

clean_gene <- function(x) {
  x <- gsub("_", "/", as.character(x))
  x
}

clean_label <- function(x) {
  x <- as.character(x)
  x <- gsub("^Cerebrovascular_", "", x)
  x <- gsub("^CapEC_", "Capillary ", x)
  x <- gsub("^Astro_", "Astro ", x)
  x <- gsub("^Micro_", "Micro ", x)
  x <- gsub("_", " ", x)
  x <- gsub("ActHigh", "ActHigh", x)
  x <- gsub("ActMid", "ActMid", x)
  x <- gsub("ActLow", "ActLow", x)
  trimws(x)
}

p_lab <- function(p, q = NA_real_) {
  p_txt <- ifelse(is.na(p), "p=NA", ifelse(p < 0.001, "p<0.001", paste0("p=", signif(p, 2))))
  q_txt <- ifelse(is.na(q), NA_character_, ifelse(q < 0.001, "q<0.001", paste0("q=", signif(q, 2))))
  ifelse(is.na(q_txt), p_txt, paste0(p_txt, "\n", q_txt))
}

save_all <- function(plot, prefix, width, height, dpi = 600) {
  ggsave(paste0(prefix, ".png"), plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(prefix, ".pdf"), plot, width = width, height = height, bg = "white")
  ggsave(paste0(prefix, ".svg"), plot, width = width, height = height, bg = "white")
  ggsave(paste0(prefix, ".tiff"), plot, width = width, height = height, dpi = dpi, bg = "white", compression = "lzw")
}

diff_ccc <- read_csv(file.path(TAB, "node06_Differential_CCC_Cerebrovascular_to_NVU_AD_vs_CN.csv"), show_col_types = FALSE)
sample_scores <- read_csv(file.path(TAB, "node06_CapEC_ModuleScore_SampleLevel.csv"), show_col_types = FALSE)
sample_stats <- read_csv(file.path(TAB, "node06_CapEC_ModuleScore_SampleLevel_Stats.csv"), show_col_types = FALSE)
state_comp <- read_csv(file.path(TAB, "node06_CapEC_ActivationState_Composition.csv"), show_col_types = FALSE)
state_lr <- read_csv(file.path(TAB, "node06_CapEC_State_LIANA_TargetFocus.csv"), show_col_types = FALSE)
act_lr <- read_csv(file.path(TAB, "node06_CapEC_ActHigh_Curated_LR_Candidates_AD.csv"), show_col_types = FALSE)
go_df <- read_csv(file.path(TAB, "node06_GO_BP_enrichment_all.csv"), show_col_types = FALSE)
react_df <- read_csv(file.path(TAB, "node06_Reactome_enrichment_all.csv"), show_col_types = FALSE)

write_csv(diff_ccc, file.path(TAB, "Table_S28_node06_cerebrovascular_to_NVU_differential_LIANA_AD_vs_CN.csv"))
write_csv(sample_stats, file.path(TAB, "Table_S29_node06_capillary_module_and_state_statistics.csv"))
write_csv(act_lr, file.path(TAB, "Table_S30_node06_activated_capillary_state_LR_candidates_AD.csv"))
write_csv(bind_rows(go_df, react_df), file.path(TAB, "Table_S31_node06_GO_Reactome_enrichment_all.csv"))

# Figure 8A: balanced differential cerebrovascular-to-NVU ligand-receptor pairs.
diff_plot <- bind_rows(
  diff_ccc %>% filter(delta_AD_vs_CN > 0) %>% arrange(desc(delta_AD_vs_CN)) %>% slice_head(n = 8),
  diff_ccc %>% filter(delta_AD_vs_CN < 0) %>% arrange(delta_AD_vs_CN) %>% slice_head(n = 8)
) %>%
  mutate(
    direction = ifelse(delta_AD_vs_CN > 0, "AD-enriched", "CN-enriched"),
    lr_pair = paste(clean_gene(ligand), clean_gene(receptor), sep = " - "),
    interaction = paste0(clean_label(source), " -> ", clean_label(target), "\n", lr_pair),
    interaction = fct_reorder(interaction, delta_AD_vs_CN)
  )
write_csv(diff_plot, file.path(SRC, "Figure_8A_source_data.csv"))

pA <- ggplot(diff_plot, aes(x = delta_AD_vs_CN, y = interaction, fill = direction)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey35") +
  scale_fill_manual(values = c("AD-enriched" = ad_col, "CN-enriched" = cn_col)) +
  theme_pub(6.5) +
  theme(legend.position = "top", axis.text.y = element_text(size = 5.2), legend.justification = "left") +
  labs(x = "Delta communication strength (AD - CN)", y = NULL, fill = NULL, title = "Differential cerebrovascular-to-NVU signaling")

# Figure 8B: sample-level capillary endothelial modules, selected for NVU relevance.
module_labels <- c(
  CapEC_BBB_maintenance_score_z = "BBB maintenance",
  CapEC_Transport_score_z = "Transport",
  CapEC_ECM_remodeling_score_z = "ECM remodeling",
  CapEC_Inflammatory_score_z = "Inflammatory",
  CapEC_Permeability_angiogenic_score_z = "Permeability",
  CapEC_dysfunction_composite_score_z = "Dysfunction"
)
score_long <- sample_scores %>%
  select(Group, sample_id, all_of(names(module_labels))) %>%
  pivot_longer(all_of(names(module_labels)), names_to = "score", values_to = "value") %>%
  mutate(
    Group = factor(Group, levels = c("CN", "AD")),
    module = factor(module_labels[score], levels = module_labels)
  )
score_ann <- sample_stats %>%
  filter(score %in% names(module_labels)) %>%
  mutate(module = factor(module_labels[score], levels = module_labels)) %>%
  select(module, p_value, p_adj_BH)
score_y <- score_long %>% group_by(module) %>% summarise(y = max(value, na.rm = TRUE) + 0.15 * diff(range(value, na.rm = TRUE)), .groups = "drop")
score_ann <- left_join(score_ann, score_y, by = "module") %>% mutate(label = p_lab(p_value, p_adj_BH))
write_csv(score_long, file.path(SRC, "Figure_8B_source_data.csv"))

pB <- ggplot(score_long, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.22, alpha = 0.75) +
  geom_point(position = position_jitter(width = 0.11, height = 0), size = 0.75, alpha = 0.75, shape = 21, stroke = 0.12) +
  geom_text(data = score_ann, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 1.75, lineheight = 0.85) +
  facet_wrap(~ module, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(CN = cn_col, AD = ad_col)) +
  theme_pub(6.5) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 5.8),
    panel.spacing = unit(0.8, "lines")
  ) +
  labs(x = NULL, y = "Mean z-score per sample", title = "Capillary endothelial programs")

# Figure 8C: capillary activation-state composition.
state_stats <- sample_stats %>%
  filter(score %in% c("frac_ActLow", "frac_ActMid", "frac_ActHigh")) %>%
  mutate(
    CapEC_activation_state = recode(score, frac_ActLow = "Capillary ActLow", frac_ActMid = "Capillary ActMid", frac_ActHigh = "Capillary ActHigh")
  ) %>%
  select(CapEC_activation_state, p_value, p_adj_BH)
state_long <- state_comp %>%
  mutate(
    Group = factor(Group, levels = c("CN", "AD")),
    CapEC_activation_state = recode(CapEC_activation_state,
      CapEC_ActLow = "Capillary ActLow",
      CapEC_ActMid = "Capillary ActMid",
      CapEC_ActHigh = "Capillary ActHigh"
    ),
    CapEC_activation_state = factor(CapEC_activation_state, levels = c("Capillary ActLow", "Capillary ActMid", "Capillary ActHigh"))
  )
state_y <- state_long %>% group_by(CapEC_activation_state) %>% summarise(y = max(frac, na.rm = TRUE) + 0.07, .groups = "drop")
state_stats <- left_join(state_stats, state_y, by = "CapEC_activation_state") %>% mutate(label = p_lab(p_value, p_adj_BH))
write_csv(state_long, file.path(SRC, "Figure_8C_source_data.csv"))

pC <- ggplot(state_long, aes(x = Group, y = frac, fill = Group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.22, alpha = 0.75) +
  geom_point(position = position_jitter(width = 0.11, height = 0), size = 0.85, alpha = 0.8, shape = 21, stroke = 0.12) +
  geom_text(data = state_stats, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 1.8, lineheight = 0.85) +
  facet_wrap(~ CapEC_activation_state, nrow = 1) +
  scale_fill_manual(values = c(CN = cn_col, AD = ad_col)) +
  theme_pub(6.5) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, max(state_stats$y, na.rm = TRUE) + 0.05)) +
  labs(x = NULL, y = "Fraction per sample", title = "Capillary activation-state composition")

# Figure 8D: AD capillary state-to-NVU subtype communication.
heat_df <- state_lr %>%
  filter(group == "AD", source %in% c("CapEC_ActLow", "CapEC_ActMid", "CapEC_ActHigh")) %>%
  mutate(
    source = recode(source, CapEC_ActLow = "ActLow", CapEC_ActMid = "ActMid", CapEC_ActHigh = "ActHigh"),
    target_clean = clean_label(target)
  ) %>%
  group_by(source, target_clean) %>%
  summarise(mean_strength = mean(communication_strength, na.rm = TRUE), .groups = "drop")
top_targets <- heat_df %>%
  group_by(target_clean) %>%
  summarise(max_strength = max(mean_strength, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_strength)) %>%
  slice_head(n = 10) %>%
  pull(target_clean)
heat_df <- heat_df %>%
  filter(target_clean %in% top_targets) %>%
  mutate(
    source = factor(source, levels = c("ActLow", "ActMid", "ActHigh")),
    target_clean = factor(target_clean, levels = top_targets)
  )
write_csv(heat_df, file.path(SRC, "Figure_8D_source_data.csv"))

pD <- ggplot(heat_df, aes(x = target_clean, y = source, fill = mean_strength)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "#F7F7F7", high = ad_col) +
  theme_pub(6.5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.3),
    legend.position = "right",
    legend.key.height = unit(0.32, "cm")
  ) +
  labs(x = NULL, y = "Capillary state", fill = "Mean\nstrength", title = "AD capillary-to-NVU communication")

# Figure 8E: ActHigh-enriched ligand-receptor candidates, no stacking.
act_plot <- act_lr %>%
  mutate(
    target_clean = clean_label(target),
    lr_pair = paste(clean_gene(ligand), clean_gene(receptor), sep = " - "),
    lr_target = paste0(lr_pair, " -> ", target_clean)
  ) %>%
  arrange(desc(delta_high_vs_low)) %>%
  slice_head(n = 16) %>%
  mutate(lr_target = fct_reorder(lr_target, delta_high_vs_low))
write_csv(act_plot, file.path(SRC, "Figure_8E_source_data.csv"))

target_fam <- function(x) {
  case_when(
    str_detect(x, "^Astro") ~ "Astrocyte",
    str_detect(x, "^Micro") ~ "Microglia",
    str_detect(x, "^Capillary|^Pericyte|^Venous|^SMC|^Arterial") ~ "Vascular",
    str_detect(x, "OPC") ~ "OPC",
    str_detect(x, "Oligodendrocyte") ~ "Oligodendrocyte",
    str_detect(x, "Excitatory|Inhibitory") ~ "Neuron",
    TRUE ~ "Other"
  )
}
act_plot <- act_plot %>% mutate(target_family = target_fam(target_clean))

pE <- ggplot(act_plot, aes(x = delta_high_vs_low, y = lr_target, fill = target_family)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c(Astrocyte = "#E69F00", Microglia = "#56B4E9", Vascular = "#009E73", OPC = "#CC79A7", Oligodendrocyte = "#D55E00", Neuron = "#0072B2", Other = grey_col)) +
  theme_pub(6.5) +
  theme(axis.text.y = element_text(size = 5.6), legend.position = "top", legend.justification = "left") +
  labs(x = "ActHigh - ActLow communication strength", y = NULL, fill = "Target family", title = "Activated capillary ligand-receptor candidates in AD")

fig8 <- ((pA | pB) / (pC | pD) / pE) +
  plot_layout(widths = c(1.15, 1), heights = c(1.15, 0.95, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11), plot.tag.position = c(0, 1))

save_all(fig8, file.path(MAIN, "Figure_8_node06_NVU_cell_cell_communication"), width = 10.8, height = 12.2)

# Supplementary Figure S20: GO and Reactome programs across three NVU families.
make_enrich_plot <- function(df, fam, db_label, top_n = 6) {
  d <- df %>%
    filter(family == fam, !is.na(p.adjust), p.adjust < 0.05) %>%
    group_by(subtype) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(
      score = -log10(p.adjust + 1e-12),
      subtype = clean_label(subtype),
      Description = str_wrap(Description, width = 42),
      Description = fct_reorder(Description, score)
    )
  write_csv(d, file.path(SRC, paste0("Figure_S20_", fam, "_", db_label, "_source_data.csv")))
  ggplot(d, aes(x = subtype, y = Description, size = Count, colour = score)) +
    geom_point(alpha = 0.95) +
    scale_colour_gradient(low = "#3C5488", high = ad_col) +
    scale_size(range = c(1.4, 4.2)) +
    theme_pub(6) +
    theme(axis.text.y = element_text(size = 4.6, lineheight = 0.9), axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right") +
    labs(x = NULL, y = NULL, colour = "-log10\nFDR", size = "Genes", title = paste0(clean_label(fam), " ", db_label))
}

s20 <- (
  make_enrich_plot(go_df, "Cerebrovascular", "GO BP") |
    make_enrich_plot(react_df, "Cerebrovascular", "Reactome")
) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))
save_all(s20, file.path(SUPP, "Figure_S20_node06_cerebrovascular_functional_enrichment"), width = 12.2, height = 6.2)

s21 <- (
  make_enrich_plot(go_df, "Astro", "GO BP") |
    make_enrich_plot(react_df, "Astro", "Reactome")
) / (
  make_enrich_plot(go_df, "Micro", "GO BP") |
    make_enrich_plot(react_df, "Micro", "Reactome")
) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))
save_all(s21, file.path(SUPP, "Figure_S21_node06_glial_functional_enrichment"), width = 12.2, height = 9.4)

# Supplementary Figure S21: vascular support module checks.
peri <- read_csv(file.path(TAB, "node06_Pericyte_ModuleScore_Stats.csv"), show_col_types = FALSE) %>% mutate(family = "Pericyte")
vsmc <- read_csv(file.path(TAB, "node06_VSMC_ModuleScore_Stats.csv"), show_col_types = FALSE) %>% mutate(family = "Arterial/SMC")
vascular_stats <- bind_rows(peri, vsmc) %>%
  mutate(
    Module = gsub("Pericyte_|VSMC_", "", Module),
    Module = gsub("_", " ", Module),
    delta_AD_vs_CN = AD_mean - CN_mean
  )
write_csv(vascular_stats, file.path(TAB, "Table_S32_node06_pericyte_vsmc_module_statistics.csv"))
write_csv(vascular_stats, file.path(SRC, "Figure_S21_source_data.csv"))

pS21 <- ggplot(vascular_stats, aes(x = delta_AD_vs_CN, y = fct_reorder(Module, delta_AD_vs_CN), fill = family)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey35") +
  facet_wrap(~ family, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("Pericyte" = "#009E73", "Arterial/SMC" = "#0072B2")) +
  theme_pub(7) +
  theme(legend.position = "none") +
  labs(x = "Delta mean module score (AD - CN)", y = NULL, title = "Pericyte and arterial/SMC module statistics")
save_all(pS21, file.path(SUPP, "Figure_S22_node06_pericyte_SMC_module_statistics"), width = 5.2, height = 5.8)

# Audit table for generated final figure package.
audit <- tibble::tibble(
  item = c("Figure_8", "Figure_S20", "Figure_S21", "Figure_S22", "Table_S28", "Table_S29", "Table_S30", "Table_S31", "Table_S32"),
  role = c(
    "Main figure: NVU communication and capillary state",
    "Supplementary figure: cerebrovascular subtype marker GO/Reactome enrichment",
    "Supplementary figure: astrocyte and microglial subtype marker GO/Reactome enrichment",
    "Supplementary figure: pericyte and arterial/SMC module checks",
    "Supplementary table: differential cerebrovascular-to-NVU LIANA",
    "Supplementary table: capillary module/state statistics",
    "Supplementary table: activated capillary LR candidates",
    "Supplementary table: GO/Reactome enrichment",
    "Supplementary table: pericyte and arterial/SMC modules"
  ),
  selection_status = c(
    "Pass after replotting; state LIANA rerun with astro/micro subtypes",
    "Pass as supplement; split from glial enrichment for readability",
    "Pass as supplement; split from vascular enrichment for readability",
    "Pass as supplement; supports vascular specificity but not central",
    "Pass",
    "Pass with caution: only nominal ActMid shift, FDR not significant",
    "Pass; hypothesis-generating mechanistic candidates",
    "Pass",
    "Pass; no significant vascular support shift"
  )
)
write_csv(audit, file.path(TAB, "node06_final_publication_audit.csv"))

message("node06 final publication package generated.")
