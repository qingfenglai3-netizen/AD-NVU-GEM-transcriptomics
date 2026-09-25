#!/usr/bin/env Rscript
# Figure 8B-E and S12 Fig: GEM scTenifoldKnk virtual-knockdown network, affected genes, GSEA and GO.

set.seed(20260608)
options(stringsAsFactors = FALSE)
if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
}

paths_to_add <- c(
  "<R_SITE_LIBRARY>",
  "<R_USER_LIBRARY>",
  "<R_BASE_LIBRARY>"
)
for (p in paths_to_add) {
  if (dir.exists(p) && !(p %in% .libPaths())) .libPaths(c(p, .libPaths()))
}

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(clusterProfiler)
  library(patchwork)
  library(scales)
  library(grid)
  library(scatterpie)
  library(png)
})

ROOT <- "<ANALYSIS_ROOT>"
JTM <- file.path(ROOT, "fuxian_test", "JTM_manuscript")
SCTK <- file.path(ROOT, "fuxian_test", "day5_3_scTenifoldKnk_GEM_QCMinCells25_ProteinCodingHVG")
GSEA_DIR <- file.path(JTM, "sandbox", "figure8_GEM_GSEA_literature_alignment_audit")
OUT <- "<WORKDIR>/figures/Fig8_source"
OUT_FIG <- file.path(OUT, "figures")
OUT_SRC <- file.path(OUT, "source_data")
OUT_LOG <- file.path(OUT, "logs")
for (d in c(OUT, OUT_FIG, OUT_SRC, OUT_LOG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

TARGET <- "GEM"
TAG <- "Cerebrovascular_cells__AD__GEM__ProteinCodingHVG2000__nCells350__d2"
RDS_FILE <- file.path(SCTK, "rds", paste0("scTenifoldKnk__", TAG, ".qs"))
DIFF_FILE <- file.path(SCTK, "tables", paste0("diffRegulation__", TAG, ".csv"))
GMT_GO <- file.path(JTM, "literature", "figure8_virtual_perturbation_benchmark", "gene_sets", "GO_Biological_Process_2023.gmt")
GMT_REACTOME <- file.path(JTM, "literature", "figure8_virtual_perturbation_benchmark", "gene_sets", "Reactome_2022.gmt")
GSEA_FILE <- file.path(GSEA_DIR, "GEM_scTenifoldKnk_GSEA_all_rank_library_results.csv")
A_SCRIPT <- "<WORKDIR>/figures/figure8a_nichenet.R"
stopifnot(
  file.exists(RDS_FILE),
  file.exists(DIFF_FILE),
  file.exists(GMT_GO),
  file.exists(GMT_REACTOME),
  file.exists(GSEA_FILE),
  file.exists(A_SCRIPT)
)

safe_gene <- function(x) toupper(trimws(as.character(x)))

theme_lit <- function(base = 6.8) {
  theme_classic(base_size = base, base_family = "Arial") +
    theme(
      text = element_text(color = "black", face = "bold", size = base),
      axis.text = element_text(color = "black", face = "bold", size = base),
      axis.title = element_text(color = "black", face = "bold", size = base),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.30, color = "black"),
      plot.title = element_text(face = "bold", size = base, hjust = 0),
      legend.title = element_text(size = base, face = "bold"),
      legend.text = element_text(size = base, face = "bold"),
      legend.key.size = unit(0.28, "cm"),
      plot.margin = margin(3, 4, 3, 4)
    )
}

clean_term <- function(x) {
  x <- sub("^.*::", "", x)
  x <- sub(" \\(GO:.*$", "", x)
  x <- sub(" R-HSA-.*$", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

short_term <- function(x) {
  x <- clean_term(x)
  x <- gsub("^Regulation Of ", "", x, ignore.case = TRUE)
  x <- gsub("^Positive Regulation Of ", "Pos. ", x, ignore.case = TRUE)
  x <- gsub("G Protein Coupled Receptor", "GPCR", x, ignore.case = TRUE)
  x <- gsub("Monoatomic", "Ion", x, ignore.case = TRUE)
  x <- gsub("Blood-Brain Barrier", "BBB", x, ignore.case = TRUE)
  x <- gsub("^Ion Cation Transmembrane Transport$", "Cation Transport", x, ignore.case = TRUE)
  x <- gsub("^Inorganic Cation Transmembrane Transport$", "Inorganic Cation Transport", x, ignore.case = TRUE)
  x <- gsub("^Cgmp-Mediated Signaling$", "cGMP Signaling", x, ignore.case = TRUE)
  x <- gsub("^Cardiac Muscle Contraction By Calcium Ion Signaling$", "Calcium-Linked Contraction", x, ignore.case = TRUE)
  x <- gsub("^Cardiac Muscle Contraction By Regulation Of The Release Of Sequestered Calcium Ion$", "Calcium Release Signaling", x, ignore.case = TRUE)
  x <- gsub("^Release Of Sequestered Calcium Ion Into Cytosol By Sarcoplasmic Reticulum$", "Calcium Release", x, ignore.case = TRUE)
  x <- gsub("^Ion Cation Import Across Plasma Membrane$", "Cation Import", x, ignore.case = TRUE)
  x
}

wrap_text <- function(x, width = 28) {
  vapply(strwrap(x, width = width, simplify = FALSE), paste, collapse = "\n", FUN.VALUE = character(1))
}

read_gmt <- function(files) {
  rows <- list()
  for (f in files) {
    lib <- tools::file_path_sans_ext(basename(f))
    lines <- readLines(f, warn = FALSE)
    rows[[basename(f)]] <- bind_rows(lapply(lines, function(line) {
      parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(parts) < 3) return(NULL)
      tibble(term = paste(lib, parts[1], sep = "::"), gene = unique(toupper(parts[-c(1, 2)])))
    }))
  }
  bind_rows(rows) |> distinct(term, gene)
}

obj <- qs::qread(RDS_FILE)
reg <- read_csv(DIFF_FILE, show_col_types = FALSE) |>
  mutate(
    gene = safe_gene(gene),
    neglog10_fdr = -log10(pmax(p.adj, 1e-300)),
    sig_dr = p.adj < 0.05 & gene != TARGET,
    sig_strict = p.adj < 0.05 & abs(Z) > 2 & gene != TARGET
  )

wt <- as.matrix(obj$tensorNetworks$WT)
if (is.null(rownames(wt))) rownames(wt) <- obj$diffRegulation$gene
if (is.null(colnames(wt))) colnames(wt) <- obj$diffRegulation$gene
rownames(wt) <- safe_gene(rownames(wt))
colnames(wt) <- safe_gene(colnames(wt))

term2gene_all <- read_gmt(c(GMT_GO, GMT_REACTOME))
term2gene_go <- read_gmt(GMT_GO)

sig_genes <- reg |>
  filter(sig_dr, distance > 1e-10) |>
  filter(!grepl("^(MT-|RPL|RPS)", gene, ignore.case = TRUE)) |>
  arrange(p.adj, desc(distance)) |>
  pull(gene)
glist <- unique(c(TARGET, sig_genes))
glist <- intersect(glist, rownames(wt))

plot_ko_network <- function(q = 0.99) {
  display_genes <- reg |>
    filter(gene %in% sig_genes, gene != TARGET) |>
    arrange(desc(FC), p.adj, desc(distance)) |>
    slice_head(n = 38) |>
    pull(gene)
  display_glist <- intersect(unique(c(TARGET, display_genes)), rownames(wt))

  scluster <- wt[display_glist, display_glist, drop = FALSE]
  ko_info <- scluster[TARGET, ]
  scluster[abs(scluster) <= stats::quantile(abs(scluster), q, na.rm = TRUE)] <- 0
  scluster[TARGET, ] <- ko_info
  diag(scluster) <- 0

  edge_tbl <- as.data.frame(as.table(scluster), stringsAsFactors = FALSE)
  colnames(edge_tbl) <- c("from", "to", "W")
  edge_tbl <- edge_tbl |>
    filter(W != 0, from != to) |>
    mutate(
      edge_sign = if_else(from == TARGET, "GEM-centered perturbation link", "peripheral WT subnetwork link"),
      absW = abs(W)
    )
  edge_tbl <- bind_rows(
    edge_tbl |> filter(from == TARGET),
    edge_tbl |> filter(from != TARGET) |> arrange(desc(absW)) |> slice_head(n = 12)
  )

  ego <- enricher(
    gene = sig_genes,
    universe = reg$gene,
    TERM2GENE = term2gene_all,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 500
  )
  pathway_tbl <- as.data.frame(ego@result) |>
    mutate(
      Description_clean = clean_term(Description),
      term_short = short_term(Description_clean)
    ) |>
    filter(!grepl("disease|cancer|infection|viral|covid|organismal|nonsense mediated decay|translation initiation|ribosome",
                  Description_clean, ignore.case = TRUE)) |>
    arrange(p.adjust, pvalue, desc(Count)) |>
    slice_head(n = 5) |>
    mutate(term_key = paste0("T", row_number()))
  if (nrow(pathway_tbl) == 0) {
    pathway_tbl <- tibble(
      ID = character(), Description = character(), Description_clean = character(),
      term_short = character(), p.adjust = numeric(), Count = integer(),
      geneID = character(), term_key = character()
    )
  }

  term_genes <- setNames(vector("list", nrow(pathway_tbl)), pathway_tbl$term_key)
  for (i in seq_len(nrow(pathway_tbl))) {
    term_genes[[pathway_tbl$term_key[[i]]]] <- unique(toupper(strsplit(pathway_tbl$geneID[[i]], "/")[[1]]))
  }

  net_plot <- igraph::graph_from_data_frame(edge_tbl, directed = TRUE)
  degree_vec <- igraph::centr_degree(net_plot)$res
  edge_weights <- rep(1, nrow(edge_tbl))
  secondary_genes <- (names(igraph::V(net_plot))[degree_vec > 1])[-1]
  edge_weights[edge_tbl$from %in% secondary_genes] <- 0.2
  edge_weights[edge_tbl$to %in% secondary_genes] <- 0.2
  edge_weights[edge_tbl$from %in% TARGET] <- 1
  edge_weights[edge_tbl$from %in% TARGET & edge_tbl$to %in% secondary_genes] <- 0.8
  set.seed(1)
  layout_tbl <- as.data.frame(igraph::layout_with_fr(net_plot, weights = edge_weights))
  colnames(layout_tbl) <- c("x", "y")
  layout_tbl$name <- names(igraph::V(net_plot))
  layout_tbl <- layout_tbl |>
    mutate(
      x = x - x[name == TARGET][1],
      y = y - y[name == TARGET][1],
      scale_denom = max(sqrt(x^2 + y^2), na.rm = TRUE),
      x = x / scale_denom * 1.48 - 1.30,
      y = y / scale_denom * 1.48
    ) |>
    select(name, x, y)

  node_tbl <- tibble(name = names(igraph::V(net_plot))) |>
    left_join(reg |> select(name = gene, distance, Z, FC, p.adj, neglog10_fdr, sig_dr), by = "name") |>
    left_join(layout_tbl, by = "name") |>
    mutate(
      degree_scaled = degree_vec[match(name, names(degree_vec))] / max(degree_vec, na.rm = TRUE),
      degree_scaled = if_else(is.finite(degree_scaled), degree_scaled, 0),
      node_size = if_else(name == TARGET, 22, 8 + 10 * degree_scaled)
    )

  for (tk in pathway_tbl$term_key) {
    node_tbl[[tk]] <- ifelse(node_tbl$name %in% term_genes[[tk]], 1, 0)
  }
  if (nrow(pathway_tbl) > 0) {
    node_tbl <- node_tbl |>
      mutate(n_terms = rowSums(across(all_of(pathway_tbl$term_key)))) |>
      mutate(across(all_of(pathway_tbl$term_key), ~ ifelse(n_terms > 0, .x / n_terms, 0)))
  } else {
    node_tbl$n_terms <- 0
  }

  edge_plot <- edge_tbl |>
    left_join(node_tbl |> select(from = name, x_from = x, y_from = y), by = "from") |>
    left_join(node_tbl |> select(to = name, x_to = x, y_to = y), by = "to") |>
    filter(!is.na(x_from), !is.na(x_to))

  write_csv(pathway_tbl, file.path(OUT_SRC, "GEM_DR_network_pathway_legend.csv"))
  write_csv(edge_plot, file.path(OUT_SRC, "GEM_DR_network_plotKO_edges.csv"))
  write_csv(node_tbl, file.path(OUT_SRC, "GEM_DR_network_nodes.csv"))

  pal <- c("#77B7C8", "#9BCDB3", "#F2D65C", "#E7A348", "#F05A5A")
  names(pal) <- pathway_tbl$term_key

  base_png <- file.path(OUT_SRC, "GEM_plotKO_tutorial_style_base_render.png")
  png(base_png, width = 4800, height = 3200, res = 600, type = "cairo", family = "Arial")
  oldpar <- par(no.readonly = TRUE)
  par(fig = c(0.00, 0.72, 0.00, 1.00), mar = c(0.02, 0.02, 0.02, 0.02), xpd = TRUE, family = "Arial")

  vertex_names <- names(igraph::V(net_plot))
  p_mat <- matrix(0, nrow = length(vertex_names), ncol = max(nrow(pathway_tbl), 1))
  rownames(p_mat) <- toupper(vertex_names)
  if (nrow(pathway_tbl) > 0) {
    for (i in seq_len(nrow(pathway_tbl))) {
      genes_i <- unique(toupper(strsplit(pathway_tbl$geneID[[i]], "/")[[1]]))
      p_mat[intersect(rownames(p_mat), genes_i), i] <- 1
    }
  }
  p_vec <- lapply(seq_len(nrow(p_mat)), function(i) as.vector(p_mat[i, ]))
  names(p_vec) <- vertex_names
  enriched_genes <- if (nrow(pathway_tbl) > 0) unique(unlist(lapply(pathway_tbl$geneID, function(x) strsplit(x, "/")[[1]]))) else character()
  is_enriched <- toupper(vertex_names) %in% toupper(enriched_genes)
  vertex_color <- rgb(195 / 255, 199 / 255, 198 / 255, 0.30)
  vertex_size <- (7 + (degree_vec / max(degree_vec, na.rm = TRUE)) * 12) * 1.35
  vertex_size[is_enriched] <- pmax(vertex_size[is_enriched] * 1.16, 18)
  vertex_size[vertex_names == TARGET] <- 34
  pie_colors <- if (nrow(pathway_tbl) > 0) pal[pathway_tbl$term_key] else "#77B7C8"
  suppressWarnings(plot(
    net_plot,
    layout = as.matrix(layout_tbl[match(vertex_names, layout_tbl$name), c("x", "y")]),
    rescale = FALSE,
    xlim = c(-3.02, 0.26),
    ylim = c(-1.62, 1.62),
    edge.arrow.size = 0.18,
    edge.width = ifelse(edge_tbl$from == TARGET, 0.75, 0.45),
    vertex.label.color = "black",
    vertex.label.family = "Arial",
    vertex.shape = ifelse(is_enriched, "pie", "circle"),
    vertex.pie = p_vec,
    vertex.pie.color = list(pie_colors),
    vertex.size = vertex_size,
    vertex.label.cex = ifelse(vertex_names == TARGET, 1.34, 0.98),
    vertex.label.font = ifelse(is_enriched, 2, 1),
    edge.color = ifelse(edge_tbl$from == TARGET, "#152BFF", "#E41A1C"),
    edge.curved = ifelse(edge_tbl$from == TARGET, 0.04, 0.12),
    vertex.color = vertex_color,
    vertex.frame.color = NA
  ))

  if (nrow(pathway_tbl) > 0) {
    par(fig = c(0.70, 1.00, 0.04, 0.98), mar = c(0.05, 0.05, 0.05, 0.05), new = TRUE, xpd = TRUE, family = "Arial")
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    sig_level <- formatC(pathway_tbl$p.adjust, digits = 2, format = "g", width = 0, drop0trailing = TRUE)
    legend_terms <- gsub("Gpcr Downstream Signaling", "GPCR downstream", pathway_tbl$term_short, ignore.case = TRUE)
    legend_terms <- gsub("Signaling by Gpcr", "GPCR signaling", legend_terms, ignore.case = TRUE)
    legend_terms <- gsub("Inorganic Cation Transport", "Inorganic cation", legend_terms, ignore.case = TRUE)
    legend_terms <- gsub("Cation Transport", "Cation transport", legend_terms, ignore.case = TRUE)
    legend_terms <- gsub("cGMP Signaling", "cGMP signaling", legend_terms, ignore.case = TRUE)
    legend_terms <- vapply(legend_terms, function(x) paste(strwrap(x, width=18), collapse="\n"), character(1))
    labels <- sprintf("(%d) %s\nFDR=%s", pathway_tbl$Count, legend_terms, sig_level)
    legend_y <- seq(0.83, 0.18, length.out = nrow(pathway_tbl))
    points(rep(0.055, length(legend_y)), legend_y,
           pch = 16, col = pie_colors, cex = 1.72)
    text(rep(0.125, length(legend_y)), legend_y, labels = labels,
         adj = c(0, 0.5), cex = 1.02, font = 2, family = "Arial")
  }

  par(oldpar)
  dev.off()

  img_full <- png::readPNG(base_png)
  content_mask <- apply(img_full[, , 1:3, drop = FALSE] < 0.992, c(1, 2), any)
  content_rows <- which(rowSums(content_mask) > 0)
  content_cols <- which(colSums(content_mask) > 0)
  pad_y <- 36
  pad_x <- 42
  r1 <- max(1, min(content_rows) - pad_y)
  r2 <- min(dim(img_full)[1], max(content_rows) + pad_y)
  c1 <- max(1, min(content_cols) - pad_x)
  c2 <- min(dim(img_full)[2], max(content_cols) + pad_x)
  img <- img_full[r1:r2, c1:c2, , drop = FALSE]
  png::writePNG(img, file.path(OUT_SRC, "GEM_plotKO_tutorial_style_base_render_tight.png"))
  ggplot() +
    annotation_custom(
      grid::rasterGrob(img, interpolate = TRUE),
      xmin = -0.08, xmax = 0.98, ymin = -0.06, ymax = 1.18
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

plot_shift_scatter <- function() {
  d <- reg |> filter(gene != TARGET)
  lab <- d |> filter(sig_dr) |> arrange(p.adj, desc(distance)) |> slice_head(n = 8)
  ggplot(d, aes(Z, neglog10_fdr)) +
    geom_point(aes(color = sig_dr), size = 0.75, alpha = 0.76) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.30, color = "grey50") +
    geom_vline(xintercept = c(0, 2), linetype = c("solid", "dashed"),
               linewidth = c(0.25, 0.30), color = c("grey75", "grey50")) +
    geom_text_repel(data = lab, aes(label = gene),
                    size = 2.12, family = "Arial", fontface = "bold",
                    min.segment.length = 0, segment.size = 0.14,
                    box.padding = 0.18, point.padding = 0.12,
                    force = 4.0, force_pull = 0.15,
                    max.overlaps = Inf, show.legend = FALSE) +
    scale_color_manual(values = c("FALSE" = "#BDBDBD", "TRUE" = "#E64B35"),
                       labels = c("Not significant", "DR gene"), name = NULL) +
    labs(title = "Differentially regulated genes", x = "DR Z-score", y = expression(-log[10](FDR))) +
    theme_lit(6.0) +
    theme(legend.position = "bottom")
}

plot_top_bar <- function() {
  d <- reg |>
    filter(gene != TARGET) |>
    arrange(desc(FC), p.adj) |>
    slice_head(n = 20) |>
    mutate(gene = factor(gene, levels = rev(gene)))
  ggplot(d, aes(FC, gene)) +
    geom_col(width = 0.72, fill = "#E64B35") +
    labs(title = "Top 20 affected genes", x = "scTenifoldKnk FC", y = NULL) +
    theme_lit(6.0) +
    theme(axis.text.y = element_text(face = "bold.italic", size = 6.0),
          legend.position = "none")
}

plot_go <- function() {
  top_genes <- sig_genes
  ego <- enricher(gene = top_genes, universe = reg$gene, TERM2GENE = term2gene_go,
                  pAdjustMethod = "BH", pvalueCutoff = 1, minGSSize = 5, maxGSSize = 500)
  ora <- as.data.frame(ego@result) |>
    mutate(Description_clean = clean_term(Description)) |>
    filter(!grepl("disease|cancer|infection|viral|covid|organismal", Description_clean, ignore.case = TRUE)) |>
    arrange(p.adjust, pvalue, desc(Count))
  write_csv(ora, file.path(OUT_SRC, "GEM_GO_ORA_significant_DR_genes.csv"))
  show <- ora |>
    distinct(Description_clean, .keep_all = TRUE) |>
    slice_head(n = 6) |>
    mutate(term_label = short_term(Description_clean),
           term = factor(wrap_text(term_label, 18),
                         levels = rev(wrap_text(term_label, 18))),
           neglog10_fdr = -log10(pmax(p.adjust, 1e-300)))
  ggplot(show, aes(Count, term)) +
    geom_col(aes(fill = neglog10_fdr), width = 0.68) +
    scale_fill_gradient(low = "#A6CEE3", high = "#E64B35", name = expression(-log[10](FDR))) +
    labs(title = "GO enrichment of DR genes", x = "Gene count", y = NULL) +
    theme_lit(6.0) +
    theme(axis.text.y = element_text(size = 6.0, face = "bold", lineheight = 0.84),
          legend.position = "right")
}

plot_gsea_distance <- function() {
  gsea_all <- read_csv(GSEA_FILE, show_col_types = FALSE) |>
    filter(rank_metric == "distance_desc", library == "local_GO_Reactome") |>
    arrange(p.adjust, pvalue, desc(abs(NES)))
  selected <- gsea_all |> slice_head(n = 1)
  term_id <- selected$ID[[1]]
  write_csv(selected, file.path(OUT_SRC, "GEM_GSEA_selected_distance_ranked_term.csv"))

  gene_list <- reg$distance
  names(gene_list) <- reg$gene
  gene_list <- sort(gene_list[is.finite(gene_list)], decreasing = TRUE)
  hit_genes <- term2gene_all |> filter(term == term_id) |> pull(gene) |> unique()
  genes <- names(gene_list)
  scores <- as.numeric(gene_list)
  hits <- genes %in% hit_genes
  weights <- abs(scores)
  running <- cumsum(ifelse(hits, weights / sum(weights[hits]), 0) -
                      ifelse(!hits, 1 / sum(!hits), 0))
  eps <- min(scores[scores > 0], na.rm = TRUE) * 0.5
  metric_display <- as.numeric(scale(log10(scores + eps)))
  df <- tibble(rank = seq_along(scores), score = scores, metric_display = metric_display,
               hit = hits, running_es = running)
  write_csv(df, file.path(OUT_SRC, "GEM_GSEA_distance_ranked_running_score.csv"))

  stat_label <- sprintf("NES %.2f, FDR %s", selected$NES[[1]], formatC(selected$p.adjust[[1]], format = "e", digits = 1))
  title <- wrap_text(clean_term(selected$Description_clean[[1]]), 32)
  p_es <- ggplot(df, aes(rank, running_es)) +
    geom_vline(xintercept = pretty(df$rank, n = 5), color = "#E6E6E6", linewidth = 0.30) +
    geom_line(color = "#D62728", linewidth = 0.82) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.28) +
    annotate("label", x = Inf, y = Inf, label = stat_label, hjust = 1.02, vjust = 1.08,
             size = 2.12, family = "Arial", fontface = "bold", linewidth = 0.12, fill = "#EFEFEF") +
    labs(title = title, x = NULL, y = "Running ES") +
    theme_lit(6.0) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(size = 6.0, face = "bold", lineheight = 0.90))
  p_hits <- ggplot(df, aes(rank, 0)) +
    geom_segment(data = df |> filter(hit), aes(x = rank, xend = rank, y = -0.40, yend = 0.40),
                 linewidth = 0.31, color = "black") +
    geom_tile(aes(y = -0.68, fill = score), height = 0.28, width = 1) +
    scale_fill_gradient(low = "#F7F7F7", high = "#E64B35", guide = "none") +
    coord_cartesian(ylim = c(-0.85, 0.45), expand = FALSE) +
    theme_void(base_family = "Arial", base_size = 6.0) +
    theme(text = element_text(face = "bold", size = 6.0))
  p_metric <- ggplot(df, aes(rank, metric_display)) +
    geom_vline(xintercept = pretty(df$rank, n = 5), color = "#E6E6E6", linewidth = 0.30) +
    geom_area(fill = "#BDBDBD", color = NA, alpha = 0.92) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.25) +
    labs(x = "Rank in ordered dataset", y = "Ranked\nmetric") +
    theme_lit(6.0)
  p_es / p_hits / p_metric + plot_layout(heights = c(1.45, 0.32, 0.70))
}

tight_raster_plot <- function(plot, stem, width_mm = 94, height_mm = 62,
                              pad_x = 18, pad_y = 20,
                              xmin = -0.02, xmax = 1.02,
                              ymin = -0.02, ymax = 1.02) {
  full_png <- file.path(OUT_SRC, paste0(stem, "_full_render.png"))
  tight_png <- file.path(OUT_SRC, paste0(stem, "_tight_render.png"))
  ggsave(full_png, plot, width = width_mm, height = height_mm,
         units = "mm", dpi = 600, bg = "white")
  img_full <- png::readPNG(full_png)
  content_mask <- apply(img_full[, , 1:3, drop = FALSE] < 0.992, c(1, 2), any)
  content_rows <- which(rowSums(content_mask) > 0)
  content_cols <- which(colSums(content_mask) > 0)
  r1 <- max(1, min(content_rows) - pad_y)
  r2 <- min(dim(img_full)[1], max(content_rows) + pad_y)
  c1 <- max(1, min(content_cols) - pad_x)
  c2 <- min(dim(img_full)[2], max(content_cols) + pad_x)
  img <- img_full[r1:r2, c1:c2, , drop = FALSE]
  png::writePNG(img, tight_png)
  ggplot() +
    annotation_custom(
      grid::rasterGrob(img, interpolate = TRUE),
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

a_env <- new.env(parent = globalenv())
source(A_SCRIPT, local = a_env)
p_a <- a_env$panel_no_tag
p_a_main <- p_a + theme(plot.margin = margin(0, -18, 0, -128))
p_net <- plot_ko_network(q = 0.99)
p_net_main <- p_net + theme(plot.margin = margin(-18, -42, -18, 10))
p_scatter <- plot_shift_scatter()
p_bar <- plot_top_bar()
p_go_supp <- plot_go()
p_gsea <- plot_gsea_distance()

combined_bcde <- ((p_net_main | p_scatter) / (p_bar | p_gsea)) +
  plot_layout(widths = c(1.05, 0.95), heights = c(1.05, 1.00)) +
  plot_annotation(tag_levels = list(c("B", "C", "D", "E"))) &
  theme(plot.tag = element_text(face = "bold", size = 12, family = "Arial"))

figure8_main <- wrap_plots(
  A = p_a_main,
  B = p_net_main,
  C = p_bar,
  D = p_gsea,
  design = "
  AB
  CD
  "
) +
  plot_layout(widths = c(0.96, 1.04), heights = c(1.05, 0.95)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) &
  theme(plot.tag = element_text(face = "bold", size = 12, family = "Arial"))

supp_dr_go <- (p_scatter | p_go_supp) +
  plot_layout(widths = c(1.02, 0.98)) +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(plot.tag = element_text(face = "bold", size = 12, family = "Arial"))

save_pub <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  ggsave(file.path(OUT_FIG, paste0(stem, ".png")), plot, width = width_mm, height = height_mm,
         units = "mm", dpi = dpi, bg = "white")
  ggsave(file.path(OUT_FIG, paste0(stem, ".tiff")), plot, width = width_mm, height = height_mm,
         units = "mm", dpi = dpi, compression = "lzw", device = "tiff", bg = "white")
  ggsave(file.path(OUT_FIG, paste0(stem, ".pdf")), plot, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(OUT_FIG, paste0(stem, ".svg")), plot, width = width_mm, height = height_mm,
         units = "mm", bg = "white")
}

save_pub(p_net, "Figure8B_GEM_official_plotKO_style_egocentric_network", 86, 66)
save_pub(p_scatter, "Figure8C_GEM_scTenifoldKnk_DR_scatter", 80, 58)
save_pub(p_bar, "Figure8D_GEM_top20_affected_genes", 80, 64)
save_pub(p_gsea, "Figure8E_GEM_distance_ranked_GSEA", 80, 64)
save_pub(p_go_supp, "Supplementary_Figure8_GEM_GO_ORA_significant_DR_genes", 80, 64)
save_pub(supp_dr_go, "Supplementary_Figure8_GEM_scTenifoldKnk_DR_scatter_GO", 170, 74)
save_pub(combined_bcde, "Figure8_GEM_scTenifoldKnk_BCDE", 170, 126)
save_pub(figure8_main, "Figure8_GEM_main_ABCD_no_volcano_layout", 170, 132)
save_pub(figure8_main, "Figure8_GEM_main_ABCDE_halfwidth_GSEA_layout", 170, 132)


message(file.path(OUT_FIG, "Figure8_GEM_main_ABCDE_halfwidth_GSEA_layout.png"))
