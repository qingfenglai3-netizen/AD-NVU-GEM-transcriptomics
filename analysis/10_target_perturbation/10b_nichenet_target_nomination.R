library(readr)
library(dplyr)
library(nichenetr)

root <- "<ANALYSIS_ROOT>/fuxian_test/JTM_manuscript"
src <- file.path(root, "source_data")
out_dir <- file.path(root, "sandbox", "figure8_target_reselection")
prior <- file.path(root, "external_data", "nichenet_v2", "ligand_target_matrix_nsga2r_final.rds")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_gene <- function(x) toupper(trimws(as.character(x)))

ltm <- readRDS(prior)
rownames(ltm) <- safe_gene(rownames(ltm))
colnames(ltm) <- safe_gene(colnames(ltm))

final_lr <- read_csv(file.path(src, "Figure_6_v3_full_filtered_LIANA_LR_source_data.csv"),
                     show_col_types = FALSE) |>
  mutate(
    ligand_clean = safe_gene(ligand),
    receptor_clean = safe_gene(receptor),
    source_family = as.character(source_family),
    target_family = as.character(target_family)
  )

major_de <- read_csv(file.path(src, "Figure_2D_major_celltype_multivolcano_source_data.csv"),
                     show_col_types = FALSE) |>
  mutate(gene_clean = safe_gene(gene))

candidate_ligands <- final_lr |>
  filter(
    delta_AD_vs_CN > 0,
    source_family %in% c("Astrocyte", "Microglia"),
    target_family == "Vascular",
    ligand_clean %in% colnames(ltm)
  ) |>
  group_by(ligand_clean, ligand) |>
  summarise(
    n_lr_edges = n(),
    max_delta_AD_vs_CN = max(delta_AD_vs_CN, na.rm = TRUE),
    mean_delta_AD_vs_CN = mean(delta_AD_vs_CN, na.rm = TRUE),
    source_families = paste(sort(unique(source_family)), collapse = ";"),
    source_cellstates = paste(sort(unique(source)), collapse = ";"),
    target_cellstates = paste(sort(unique(target)), collapse = ";"),
    lr_pair_examples = paste(sort(unique(lr_pair))[seq_len(min(5, n_distinct(lr_pair)))], collapse = ";"),
    .groups = "drop"
  ) |>
  arrange(desc(max_delta_AD_vs_CN), desc(n_lr_edges))

receiver_de <- major_de |>
  filter(celltype == "Cerebrovascular cells")

background <- receiver_de |>
  filter(!is.na(PValue), gene_clean %in% rownames(ltm)) |>
  pull(gene_clean) |>
  unique()

receiver_targets <- receiver_de |>
  filter(logFC > 0, PValue < 0.05, gene_clean %in% background) |>
  arrange(PValue, FDR, desc(logFC)) |>
  pull(gene_clean) |>
  unique()

top_ligands <- candidate_ligands |>
  pull(ligand_clean) |>
  unique() |>
  intersect(colnames(ltm))

ligand_activities <- predict_ligand_activities(
  geneset = receiver_targets,
  background_expressed_genes = background,
  ligand_target_matrix = ltm,
  potential_ligands = top_ligands
) |>
  as_tibble() |>
  mutate(ligand_clean = safe_gene(test_ligand)) |>
  left_join(candidate_ligands, by = "ligand_clean") |>
  arrange(desc(aupr), desc(pearson), desc(auroc))

top_ligands_by_aupr <- ligand_activities |>
  slice_head(n = 20) |>
  pull(ligand_clean) |>
  intersect(colnames(ltm))

target_rows <- list()
for (g in receiver_targets) {
  vals <- ltm[g, top_ligands_by_aupr, drop = TRUE]
  vals <- vals[order(vals, decreasing = TRUE)]
  vals_pos <- vals[vals > 0]
  if (length(vals_pos) == 0) next
  target_rows[[length(target_rows) + 1]] <- tibble(
    receiver_context = "Cerebrovascular cells",
    target_clean = g,
    max_regulatory_potential = max(vals_pos),
    sum_regulatory_potential = sum(vals_pos),
    n_supporting_glial_ligands = length(vals_pos),
    best_ligand = names(vals_pos)[1],
    supporting_ligands_top5 = paste(names(vals_pos)[seq_len(min(5, length(vals_pos)))],
                                    round(vals_pos[seq_len(min(5, length(vals_pos)))], 4),
                                    sep = ":", collapse = ";")
  )
}

target_priority <- bind_rows(target_rows) |>
  left_join(
    receiver_de |>
      select(celltype, gene_clean, logFC, PValue, FDR, direction),
    by = c("target_clean" = "gene_clean")
  ) |>
  arrange(desc(max_regulatory_potential), desc(n_supporting_glial_ligands), PValue) |>
  mutate(rank_in_vascular_receiver = row_number())

write_csv(candidate_ligands, file.path(out_dir, "glia_to_vascular_candidate_ligands.csv"))
write_csv(ligand_activities, file.path(out_dir, "glia_to_vascular_ligand_activities.csv"))
write_csv(target_priority, file.path(out_dir, "glia_to_vascular_target_priority.csv"))

print("Candidate ligand summary")
print(candidate_ligands |> summarise(n_ligands = n_distinct(ligand_clean),
                                     n_edges = sum(n_lr_edges)))
print("Top 30 target-side vascular genes from glial ligands")
print(target_priority |>
        select(rank_in_vascular_receiver, target_clean, max_regulatory_potential,
               n_supporting_glial_ligands, best_ligand, logFC, PValue, FDR) |>
        slice_head(n = 30))
print("GEM/NAMPT/PLAT")
print(target_priority |>
        filter(target_clean %in% c("GEM", "NAMPT", "PLAT")) |>
        select(rank_in_vascular_receiver, target_clean, max_regulatory_potential,
               n_supporting_glial_ligands, best_ligand, supporting_ligands_top5,
               logFC, PValue, FDR))
