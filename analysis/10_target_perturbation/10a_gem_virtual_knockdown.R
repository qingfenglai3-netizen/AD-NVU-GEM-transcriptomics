#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scTenifoldKnk)
  library(qs)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

SEED <- 20260607
set.seed(SEED)

PROJECT_ROOT <- "/path/to/project"
TARGET_GENE <- Sys.getenv("SCTK_TARGET", "GEM")
OUT <- file.path(PROJECT_ROOT, "results", sprintf("10_gem_virtual_knockdown_%s_LiteratureAlignedProteinCodingHVG", TARGET_GENE))
OUT_FIG <- file.path(OUT, "figures")
OUT_TAB <- file.path(OUT, "tables")
OUT_RDS <- file.path(OUT, "rds")
OUT_LOG <- file.path(OUT, "logs")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_RDS, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_LOG, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUT_LOG, sprintf("10_gem_virtual_knockdown_%s_literature_aligned_protein_coding_HVG.log", TARGET_GENE))
if (file.exists(LOG_FILE)) file.remove(LOG_FILE)

GROUP_VALUE <- "AD"
lineages_env <- Sys.getenv("SCTK_LINEAGES", "")
LINEAGES <- if (nzchar(lineages_env)) strsplit(lineages_env, "\\|")[[1]] else c("Astrocytes", "Cerebrovascular cells")
HVG_N <- as.integer(Sys.getenv("SCTK_HVG_N", "2000"))
MAX_CELLS <- as.integer(Sys.getenv("SCTK_MAX_CELLS", "900"))
NC_NCELLS <- as.integer(Sys.getenv("SCTK_NC_NCELLS", "350"))
MA_NDIM <- as.integer(Sys.getenv("SCTK_MA_NDIM", "2"))

msg <- function(x) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), x)
  message(line)
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

safe_name <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

get_counts <- function(obj) {
  assay <- DefaultAssay(obj)
  out <- tryCatch(GetAssayData(obj, assay = assay, layer = "counts"), error = function(e) NULL)
  if (is.null(out)) {
    out <- GetAssayData(obj, assay = assay, slot = "counts")
  }
  out
}

is_technical_gene <- function(genes) {
  g <- toupper(genes)
  grepl("^MT-", g) |
    grepl("^RPS[0-9A-Z]", g) |
    grepl("^RPL[0-9A-Z]", g) |
    grepl("^HIST[0-9H]", g) |
    grepl("^H[1234][A-Z0-9]", g)
}

balanced_cells <- function(meta, lineage, group_value = "AD", max_cells = 900) {
  hit <- rownames(meta)[meta$celltype == lineage & meta$Group == group_value]
  if (length(hit) <= max_cells) return(hit)
  sample(hit, max_cells)
}

annotate_gene_type <- function(genes) {
  genes_up <- toupper(genes)
  ann <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(genes_up),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "GENETYPE", "ENTREZID", "GENENAME")
  ) |>
    as_tibble() |>
    group_by(SYMBOL) |>
    summarise(
      GENETYPE = paste(sort(unique(na.omit(GENETYPE))), collapse = ";"),
      ENTREZID = paste(sort(unique(na.omit(ENTREZID))), collapse = ";"),
      GENENAME = paste(sort(unique(na.omit(GENENAME))), collapse = ";"),
      .groups = "drop"
    ) |>
    mutate(
      GENETYPE = ifelse(GENETYPE == "", NA_character_, GENETYPE),
      ENTREZID = ifelse(ENTREZID == "", NA_character_, ENTREZID),
      GENENAME = ifelse(GENENAME == "", NA_character_, GENENAME)
    )
  tibble(gene = genes, SYMBOL = genes_up) |>
    left_join(ann, by = "SYMBOL")
}

select_literature_aligned_protein_coding_hvg <- function(counts_all_cells, target_gene, hvg_n) {
  genes0 <- rownames(counts_all_cells)
  if (!target_gene %in% genes0) stop(sprintf("Target gene %s is absent from the input count matrix", target_gene))

  gene_annot <- annotate_gene_type(genes0)
  tech <- is_technical_gene(genes0)
  protein_coding <- !is.na(gene_annot$GENETYPE) & gene_annot$GENETYPE == "protein-coding"
  keep <- (!tech & protein_coding) | genes0 == target_gene
  keep_genes <- genes0[keep]
  filtered_counts <- counts_all_cells[keep_genes, , drop = FALSE]

  seu <- CreateSeuratObject(counts = filtered_counts, project = paste0(target_gene, "_scTenifoldKnk_ProteinCodingHVG"))
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = hvg_n, verbose = FALSE)
  hvg <- VariableFeatures(seu)
  selected <- unique(c(hvg, target_gene))
  selected <- selected[selected %in% rownames(filtered_counts)]

  metrics <- gene_annot |>
    mutate(
      technical_gene = tech,
      protein_coding_gene = protein_coding,
      retained_for_hvg_selection = keep,
      selected_hvg = gene %in% hvg,
      force_included_target = gene == target_gene,
      selected_for_sctenifold = gene %in% selected
    )

  list(
    selected = selected,
    metrics = metrics,
    target_force_included = target_gene %in% selected && !target_gene %in% hvg,
    n_technical_removed = sum(tech),
    n_non_protein_coding_removed = sum(!protein_coding & genes0 != target_gene, na.rm = TRUE),
    n_protein_coding_available = sum(protein_coding, na.rm = TRUE)
  )
}

run_one <- function(counts, lineage, group_value, gko) {
  tag <- paste(safe_name(lineage), group_value, gko, paste0("ProteinCodingHVG", HVG_N), paste0("nCells", NC_NCELLS), paste0("d", MA_NDIM), sep = "__")
  rds_file <- file.path(OUT_RDS, paste0("scTenifoldKnk__", tag, ".qs"))
  tab_file <- file.path(OUT_TAB, paste0("diffRegulation__", tag, ".csv"))
  if (file.exists(rds_file) && file.exists(tab_file)) {
    msg(sprintf("Reuse existing protein-coding HVG result: %s", tag))
    dr <- read_csv(tab_file, show_col_types = FALSE)
    return(list(tag = tag, diff = dr, status = "reused"))
  }

  msg(sprintf("Running scTenifoldKnk: lineage=%s group=%s gKO=%s genes=%d cells=%d",
              lineage, group_value, gko, nrow(counts), ncol(counts)))
  res <- scTenifoldKnk(
    countMatrix = as.matrix(counts),
    qc = FALSE,
    gKO = gko,
    nc_nNet = 10,
    nc_nCells = min(NC_NCELLS, ncol(counts)),
    nc_nComp = 3,
    nc_q = 0.9,
    td_K = 3,
    td_maxIter = 1000,
    ma_nDim = MA_NDIM,
    nCores = 2
  )
  dr <- as_tibble(res$diffRegulation) |>
    mutate(lineage = lineage, group = group_value, gKO = gko,
           hvg_n = HVG_N, gene_universe_size = nrow(counts), gene_universe = "technical-filtered protein-coding HVG", .before = 1) |>
    arrange(p.adj, desc(distance))
  write_csv(dr, tab_file)
  qs::qsave(res, rds_file)
  list(tag = tag, diff = dr, status = "new")
}

msg(sprintf("=== Literature-aligned protein-coding HVG %s scTenifoldKnk started ===", TARGET_GENE))
msg(sprintf("Seed=%d; target=%s; group=%s; lineages=%s; hvg_n=%d; max_cells=%d; nc_nCells=%d; ma_nDim=%d",
            SEED, TARGET_GENE, GROUP_VALUE, paste(LINEAGES, collapse = "|"), HVG_N, MAX_CELLS, NC_NCELLS, MA_NDIM))

sc_file <- file.path(PROJECT_ROOT, "results", "02_annotation", "rds", "scRNA_annotated.qs")
obj <- qs::qread(sc_file)
obj$Group <- dplyr::recode(as.character(obj$Group), "NC" = "CN", "UC" = "AD", .default = as.character(obj$Group))
all_counts <- get_counts(obj)
meta <- obj@meta.data

all_summary <- list()
all_diff <- list()
all_gene_metrics <- list()

for (lineage in LINEAGES) {
  all_lineage_cells <- rownames(meta)[meta$celltype == lineage & meta$Group == GROUP_VALUE]
  cells <- balanced_cells(meta, lineage, GROUP_VALUE, max_cells = MAX_CELLS)
  if (length(cells) < 120) {
    msg(sprintf("Skip %s: too few %s cells (%d)", lineage, GROUP_VALUE, length(cells)))
    next
  }

  lineage_counts_all_cells <- all_counts[, all_lineage_cells, drop = FALSE]
  selection <- select_literature_aligned_protein_coding_hvg(lineage_counts_all_cells, TARGET_GENE, HVG_N)
  genes <- selection$selected
  lineage_counts <- all_counts[genes, cells, drop = FALSE]
  gene_metrics <- selection$metrics |>
    mutate(lineage = lineage, group = GROUP_VALUE, hvg_n = HVG_N, .before = 1)
  all_gene_metrics[[lineage]] <- gene_metrics

  msg(sprintf("%s %s protein-coding aligned matrix: %d genes x %d sampled cells; %d available cells; protein-coding available=%d; technical removed=%d; non-protein-coding/unknown removed=%d",
              lineage, GROUP_VALUE, nrow(lineage_counts), ncol(lineage_counts), length(all_lineage_cells),
              selection$n_protein_coding_available, selection$n_technical_removed, selection$n_non_protein_coding_removed))

  out <- tryCatch(
    run_one(lineage_counts, lineage, GROUP_VALUE, TARGET_GENE),
    error = function(e) {
      msg(sprintf("FAILED lineage=%s gKO=%s error=%s", lineage, TARGET_GENE, e$message))
      NULL
    }
  )
  if (is.null(out)) next

  dr <- out$diff
  sig_fdr <- dr |> filter(gene != TARGET_GENE, p.adj < 0.05, abs(Z) > 2)
  sig_fdr_only <- dr |> filter(gene != TARGET_GENE, p.adj < 0.05)
  sig_nominal <- dr |> filter(gene != TARGET_GENE, p.value < 0.05, abs(Z) > 2)
  top <- dr |> filter(gene != TARGET_GENE) |> arrange(p.adj, desc(abs(Z))) |> slice_head(n = 20)
  all_summary[[lineage]] <- tibble(
    lineage = lineage,
    group = GROUP_VALUE,
    gKO = TARGET_GENE,
    status = out$status,
    total_available_cells = length(all_lineage_cells),
    sampled_cells = ncol(lineage_counts),
    technical_genes_removed = selection$n_technical_removed,
    non_protein_coding_or_unknown_removed = selection$n_non_protein_coding_removed,
    protein_coding_genes_available = selection$n_protein_coding_available,
    hvg_n_requested = HVG_N,
    final_input_genes = nrow(lineage_counts),
    gene_universe = sprintf("technical-filtered protein-coding HVG plus forced %s", TARGET_GENE),
    target_force_included = selection$target_force_included,
    n_fdr_lt_0_05 = nrow(sig_fdr_only),
    n_nominal_p_lt_0_05_absZ_gt_2 = nrow(sig_nominal),
    n_fdr_lt_0_05_absZ_gt_2 = nrow(sig_fdr),
    median_distance = median(dr$distance, na.rm = TRUE),
    max_distance = max(dr$distance, na.rm = TRUE),
    min_pvalue = suppressWarnings(min(dr$p.value[dr$gene != TARGET_GENE], na.rm = TRUE)),
    min_fdr = suppressWarnings(min(dr$p.adj[dr$gene != TARGET_GENE], na.rm = TRUE)),
    top_affected_genes = paste(top$gene, collapse = ";"),
    analysis_role = sprintf("Post-nomination %s virtual knockout; protein-coding HVG universe; GSEA should use the complete diffRegulation ranked list from this input matrix", TARGET_GENE)
  )
  all_diff[[lineage]] <- dr
}

summary_tbl <- bind_rows(all_summary)
diff_tbl <- bind_rows(all_diff)
gene_metric_tbl <- bind_rows(all_gene_metrics)

write_csv(summary_tbl, file.path(OUT_TAB, sprintf("TableS1_%s_literature_aligned_protein_coding_HVG_context_summary.csv", TARGET_GENE)))
write_csv(diff_tbl, file.path(OUT_TAB, sprintf("TableS2_%s_literature_aligned_protein_coding_HVG_all_diffRegulation.csv", TARGET_GENE)))
write_csv(gene_metric_tbl, file.path(OUT_TAB, sprintf("TableS3_%s_literature_aligned_protein_coding_HVG_selection_metrics.csv", TARGET_GENE)))

if (nrow(summary_tbl) > 0) {
  p <- summary_tbl |>
    mutate(lineage = factor(lineage, levels = LINEAGES)) |>
    ggplot(aes(lineage, final_input_genes, fill = lineage)) +
    geom_col(width = 0.65, color = "black", linewidth = 0.25) +
    geom_text(aes(label = final_input_genes), vjust = -0.4, size = 2.5, fontface = "bold") +
    scale_fill_manual(values = c("Astrocytes" = "#E69F00", "Cerebrovascular cells" = "#4C78A8")) +
    labs(x = NULL, y = "Genes entering scTenifoldKnk", title = sprintf("%s virtual KO protein-coding HVG input", TARGET_GENE)) +
    theme_classic(base_size = 7, base_family = "Arial") +
    theme(text = element_text(face = "bold"), legend.position = "none")
  ggsave(file.path(OUT_FIG, sprintf("FigQC_%s_literature_aligned_protein_coding_HVG_input_counts.pdf", TARGET_GENE)), p,
         width = 4.4, height = 3.2, units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(OUT_FIG, sprintf("FigQC_%s_literature_aligned_protein_coding_HVG_input_counts.png", TARGET_GENE)), p,
         width = 4.4, height = 3.2, units = "in", dpi = 600, bg = "white")
}

manifest <- tibble(
  output = list.files(OUT, recursive = TRUE, full.names = FALSE)
) |>
  mutate(size_bytes = file.info(file.path(OUT, output))$size)
write_csv(manifest, file.path(OUT_TAB, sprintf("output_manifest_%s_literature_aligned_protein_coding_HVG.csv", TARGET_GENE)))
writeLines(capture.output(sessionInfo()), file.path(OUT_LOG, "sessionInfo.txt"))

msg(sprintf("=== Literature-aligned protein-coding HVG %s scTenifoldKnk complete ===", TARGET_GENE))
