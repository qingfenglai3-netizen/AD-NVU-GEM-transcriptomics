#!/usr/bin/env Rscript

# Generate the GSE5281 external bulk-expression support figure and table.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
  library(svglite)
  library(ragg)
})

PROJECT_ROOT <- Sys.getenv("AD_NVU_PROJECT_ROOT", "/path/to/project")
in_dir <- file.path(PROJECT_ROOT, "results", "12_bulk_external_support")
out_dir <- file.path(PROJECT_ROOT, "results", "12_gse5281_support")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

all_res <- read.csv(
  file.path(in_dir, "AD_bulk_GEO_GEM_module_screen_results_v2_logscale_checked.csv"),
  stringsAsFactors = FALSE
)
reg <- read.csv(
  file.path(in_dir, "GSE5281_region_stratified_GEM_module_screen_v2.csv"),
  stringsAsFactors = FALSE
)

gse5281 <- all_res[all_res$GSE == "GSE5281" & all_res$status == "ok", ]
gse5281$gene <- factor(gse5281$gene, levels = gse5281$gene[order(gse5281$logFC)])
gse5281$sig <- ifelse(
  gse5281$adj.P.Val < 0.05 & gse5281$logFC > 0,
  "AD-up FDR < 0.05",
  ifelse(gse5281$adj.P.Val < 0.05, "FDR < 0.05 opposite", "NS")
)

reg_gem <- reg[reg$gene == "GEM" & reg$status == "ok", ]
reg_gem$subset <- factor(reg_gem$subset, levels = reg_gem$subset[order(reg_gem$logFC)])
reg_gem$sig <- ifelse(
  reg_gem$adj.P.Val < 0.05,
  "FDR < 0.05",
  ifelse(reg_gem$P.Value < 0.05, "nominal P < 0.05", "NS")
)

wb <- createWorkbook()
addWorksheet(wb, "GSE5281 module genes")
writeData(
  wb,
  "GSE5281 module genes",
  gse5281[, c("GSE", "platform", "n_AD", "n_Control", "gene", "probe", "logFC", "P.Value", "adj.P.Val", "transformed")]
)
addWorksheet(wb, "GSE5281 GEM by region")
writeData(
  wb,
  "GSE5281 GEM by region",
  reg_gem[, c("subset", "n_AD", "n_Control", "gene", "probe", "logFC", "P.Value", "adj.P.Val", "transformed")]
)
saveWorkbook(
  wb,
  file.path(out_dir, "Supplementary_Table_external_GEM_validation_GSE5281.xlsx"),
  overwrite = TRUE
)

base_theme <- theme_classic(base_family = "Arial", base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

pB <- ggplot(gse5281, aes(x = gene, y = logFC, fill = sig)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c(
    "AD-up FDR < 0.05" = "#D55E00",
    "FDR < 0.05 opposite" = "#0072B2",
    "NS" = "#BDBDBD"
  )) +
  labs(
    title = "GSE5281 bulk AD brain: GEM module genes",
    x = NULL,
    y = "AD vs control log2 fold change",
    fill = "Status"
  ) +
  base_theme

pC <- ggplot(reg_gem, aes(x = subset, y = logFC, fill = sig)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c(
    "FDR < 0.05" = "#D55E00",
    "nominal P < 0.05" = "#E69F00",
    "NS" = "#BDBDBD"
  )) +
  labs(
    title = "GSE5281 GEM by brain region",
    x = NULL,
    y = "AD vs control log2 fold change",
    fill = "Status"
  ) +
  base_theme

fig <- pB | pC
ggsave(file.path(out_dir, "Supplementary_Figure_external_GEM_validation_GSE5281.pdf"), fig, width = 180, height = 170, units = "mm")
ggsave(file.path(out_dir, "Supplementary_Figure_external_GEM_validation_GSE5281.svg"), fig, width = 180, height = 170, units = "mm", device = svglite::svglite)
ggsave(file.path(out_dir, "Supplementary_Figure_external_GEM_validation_GSE5281.png"), fig, width = 180, height = 170, units = "mm", dpi = 600, device = ragg::agg_png)
ggsave(file.path(out_dir, "Supplementary_Figure_external_GEM_validation_GSE5281.tiff"), fig, width = 180, height = 170, units = "mm", dpi = 600, compression = "lzw", device = ragg::agg_tiff)

cat("Completed. Output directory:", out_dir, "
")
