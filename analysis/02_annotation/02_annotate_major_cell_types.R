################################################################################
#
# Hybrid: Integrated analysis structure with a curated seven-class marker panel
#
# Implementation notes
#   - AverageExpression forces layer = "data"; JoinLayers first
#   - Module scoring uses only genes present (no silent drop)
#   - z-score normalization across cell types per cluster
#   - AddModuleScore column rename robust
# [Annotation]
#   - Dual evidence: mean_expr_z + module_z consensus
#   - Confidence: High / Medium / Low by score margin
#   - Cerebrovascular cells defined as the coarse brain vascular compartment (refinement removed per user decision)
# [publication-required outputs]
#   - Main Fig 1: UMAP by celltype + group + dataset
#   - Main Fig 2: Canonical marker dot plot + stacked violin
#   - Main Fig 3: Annotation score heatmap + confidence bar
#   - Main Fig 4: Composition (group / dataset / sample)
#   - ED Fig 1: QC violins per cluster
#   - ED Fig 2: QC metrics on UMAP
#   - ED Fig 3: Full marker dot plot (cluster level)
#   - ED Fig 4: Annotation consensus scatter
#   - ED Fig 5: Optional doublet diagnostics
#   - ED Fig 6: Inhibitory vs Excitatory UMAP detail
#   - Source data CSV for every figure (publication output)
################################################################################

set.seed(123456)
options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 200 * 1024^3)
paths_to_add <- "/path/to/R_libraries"
for (p in paths_to_add) {
  if (dir.exists(p) && !(p %in% .libPaths())) {
    .libPaths(c(p, .libPaths()))
  }
}
if (!is.null(Sys.getenv("R_LIBS_USER", unset = "")) &&
    dir.exists(Sys.getenv("R_LIBS_USER")) &&
    !(Sys.getenv("R_LIBS_USER") %in% .libPaths())) {
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
  library(grid)
})
has_qs <- requireNamespace("qs", quietly = TRUE)
if (has_qs) suppressPackageStartupMessages(library(qs))

read_annot_object <- function(path) {
  if (has_qs && endsWith(path, ".qs")) {
    qread(path)
  } else {
    readRDS(path)
  }
}

save_annot_object <- function(object, path) {
  if (has_qs && endsWith(path, ".qs")) {
    qs::qsave(object, path, preset = "fast")
  } else {
    saveRDS(object, path)
  }
}

safe_copy_file <- function(src, dst) {
  tryCatch(
    {
      if (file.exists(src)) {
        dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
        file.copy(src, dst, overwrite = TRUE)
      }
    },
    error = function(e) {
      log_msg(sprintf("[WARN] Copy failed: %s -> %s (%s)", src, dst, e$message))
      return(FALSE)
    }
  )
}

# Optional packages (graceful degradation)
has_ggsci   <- requireNamespace("ggsci",   quietly = TRUE)
has_pals    <- requireNamespace("pals",    quietly = TRUE)
has_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)

# ==============================================================================
# USER PARAMETERS
# ==============================================================================
TEST_MODE <- FALSE
if (TEST_MODE) {
  IN_RDS  <- "/path/to/project/results/01_atlas/rds"
  IN_FILE <- file.path(IN_RDS, "scRNA_node1_final.qs")  # from node01
  OUT     <- "/path/to/project/results/02_annotation"
} else {
  IN_FILE <- "scRNA_node1_final.qs"
  OUT     <- "02_annotation"
}

OUT_FIG <- file.path(OUT, "figures_main")
OUT_ED  <- file.path(OUT, "figures_extended_data")
OUT_TAB <- file.path(OUT, "tables")
OUT_RDS <- file.path(OUT, "rds")
OUT_LOG <- file.path(OUT, "logs")
OUT_SRC <- file.path(OUT, "source_data")
for (d in c(OUT_FIG, OUT_ED, OUT_TAB, OUT_RDS, OUT_LOG, OUT_SRC))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

LOG_FILE <- file.path(OUT_LOG, "node02_run.log")
log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", ts, msg)
  cat(line, "\n")
  tryCatch(
    cat(line, "\n", file = LOG_FILE, append = TRUE),
    error = function(e) {
      cat(sprintf("[WARN] log write failed, continue without file log: %s\n", e$message))
    }
  )
}

# ==============================================================================
# publication THEME & PALETTES
# ==============================================================================
publication_theme <- function(base = 7) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      axis.text       = element_text(color = "black", size = base),
      axis.title      = element_text(color = "black", size = base + 0.5),
      axis.line       = element_line(color = "black", linewidth = 0.35),
      axis.ticks      = element_line(color = "black", linewidth = 0.35),
      legend.text     = element_text(size = base - 0.5),
      legend.title    = element_text(size = base),
      legend.key.size = unit(0.32, "cm"),
      strip.background= element_blank(),
      strip.text      = element_text(size = base, face = "bold"),
      plot.title      = element_text(size = base + 1, face = "bold", hjust = 0.5),
      panel.grid      = element_blank(),
      plot.margin     = margin(2, 2, 2, 2)
    )
}

# publication palette: NPG for <=10 levels
publication_palette <- function(n) {
  if (n <= 10 && has_ggsci) {
    ggsci::pal_npg("nrc")(n)
  } else if (has_pals) {
    as.character(pals::glasbey(n))
  } else {
    colorRampPalette(brewer.pal(8, "Set1"))(n)
  }
}

group_palette <- c(
  "CN" = "#4DBBD5", "AD" = "#E64B35",
  "NC" = "#4DBBD5", "UC" = "#E64B35",
  "Ctrl" = "#4DBBD5", "Control" = "#4DBBD5", "Disease" = "#E64B35"
)

# 7 cell type colors (colorblind-friendly, Wong-style)
ct_order <- c("Oligodendrocytes", "OPCs", "Astrocytes", "Microglia",
              "Cerebrovascular cells", "Inhibitory", "Excitatory")
ct_colors_default <- c(
  "Oligodendrocytes" = "#0072B2",
  "OPCs"             = "#009E73",
  "Astrocytes"       = "#CC79A7",
  "Microglia"        = "#D55E00",
  "Cerebrovascular cells"      = "#B79F00",
  "Inhibitory"       = "#56B4E9",
  "Excitatory"       = "#999999"
)

# ==============================================================================
# HELPERS
# ==============================================================================
save_fig <- function(p, name, dir, w = 6, h = 4) {
  pdf_f  <- file.path(dir, paste0(name, ".pdf"))
  png_f  <- file.path(dir, paste0(name, ".png"))
  svg_f  <- file.path(dir, paste0(name, ".svg"))
  tiff_f <- file.path(dir, paste0(name, ".tiff"))
  tryCatch(
    ggsave(pdf_f,  p, width = w, height = h, units = "in", device = cairo_pdf),
    error = function(e) log_msg(sprintf("  PDF skipped (%s): %s", name, e$message))
  )
  tryCatch(
    ggsave(png_f,  p, width = w, height = h, units = "in", dpi = 300),
    error = function(e) log_msg(sprintf("  PNG skipped (%s): %s", name, e$message))
  )
  tryCatch(
    ggsave(svg_f,  p, width = w, height = h, units = "in"),
    error = function(e) log_msg(sprintf("  SVG skipped (%s): %s", name, e$message))
  )
  tryCatch(
    ggsave(tiff_f, p, width = w, height = h, units = "in", dpi = 600,
           compression = "lzw", device = "tiff"),
    error = function(e) log_msg(sprintf("  TIFF skipped (%s): %s", name, e$message))
  )
  log_msg(sprintf("  Saved: %s.{pdf,png,svg,tiff}", name))
}

safe_write_csv <- function(x, f) {
  tryCatch(write.csv(x, f, row.names = FALSE),
           error = function(e) log_msg(sprintf("  CSV save failed: %s", f)))
}

save_source <- function(df, fig_name) {
  safe_write_csv(df, file.path(OUT_SRC, paste0(fig_name, "_source_data.csv")))
}

standardize_group <- function(x) {
  v <- tolower(trimws(as.character(x)))
  v[v %in% c("cn", "nc", "ctrl", "control", "normal", "normal_control", "healthy")] <- "CN"
  v[v %in% c("ad", "uc", "alz", "alzheimers", "alzheimersdisease", "disease", "case")] <- "AD"
  v
}

find_optional_doublet_cols <- function(meta) {
  cand <- c("DF.classifications", "doublet_finder", "scDblFinder.class",
            "scDblFinder_class", "doublet", "is_doublet")
  hits <- intersect(cand, colnames(meta))
  hits <- c(hits, grep("(?i)doublet|dblfinder", colnames(meta), value = TRUE))
  unique(hits)
}

# ==============================================================================
# LOAD DATA
# ==============================================================================
log_msg("===== node02 annotation (7 types x 3 markers, publication-ready) =====")
log_msg(sprintf("Input: %s", IN_FILE))
stopifnot(file.exists(IN_FILE))
sc <- read_annot_object(IN_FILE)
log_msg(sprintf("Loaded: %d cells x %d features", ncol(sc), nrow(sc)))

DefaultAssay(sc) <- "RNA"

# Seurat v5: join layers and ensure data layer exists
if (inherits(sc[["RNA"]], "Assay5")) {
  log_msg("Detected Seurat v5 Assay5 - joining layers ...")
  sc[["RNA"]] <- JoinLayers(sc[["RNA"]])
}

# Ensure normalized data
data_layers <- Layers(sc[["RNA"]])
if (!"data" %in% data_layers || all(GetAssayData(sc, assay = "RNA", layer = "data")@x == 0)) {
  log_msg("No normalized data found - running NormalizeData ...")
  sc <- NormalizeData(sc, verbose = FALSE)
}

# Ensure clusters present
if (!"seurat_clusters" %in% colnames(sc[[]])) {
  stop("seurat_clusters not found in metadata. Did node01 finish clustering?")
}
Idents(sc) <- "seurat_clusters"
n_clusters <- length(unique(as.character(sc$seurat_clusters)))
log_msg(sprintf("Clusters: %d", n_clusters))

# Standardize group labels for manuscript consistency (NC/UC -> CN/AD)
if ("Group" %in% colnames(sc[[]])) {
  sc$Group <- standardize_group(sc$Group)
} else if ("group" %in% colnames(sc[[]])) {
  sc$group <- standardize_group(sc$group)
  sc$Group <- sc$group
} else if ("Condition" %in% colnames(sc[[]])) {
  sc$Group <- standardize_group(sc$Condition)
}

# ==============================================================================
# MARKER PANEL: 7 cell types x 3 markers
# (Cerebrovascular cells defined as the coarse brain vascular compartment per user decision)
# ==============================================================================
markers_list <- list(
  Oligodendrocytes = c("MBP", "PLP1", "MOBP"),
  OPCs             = c("PDGFRA", "VCAN", "CSPG4"),
  Astrocytes       = c("AQP4", "GFAP", "ADGRV1"),
  Microglia        = c("C3", "P2RY12", "TREM2"),
  "Cerebrovascular cells" = c("CLDN5", "FLT1", "VWF"),
  Inhibitory       = c("GAD1", "GAD2", "PCDH15"),
  Excitatory       = c("CAMK2A", "SLC17A7", "SATB2")
)

# Filter to genes present
present <- rownames(sc)
markers_present <- lapply(markers_list, function(x) intersect(x, present))
marker_avail_tbl <- data.frame(
  cell_type = names(markers_list),
  total_markers = sapply(markers_list, length),
  available     = sapply(markers_present, length),
  available_genes = sapply(markers_present, paste, collapse = ","),
  missing_genes   = sapply(seq_along(markers_list), function(i)
    paste(setdiff(markers_list[[i]], markers_present[[i]]), collapse = ","))
)
safe_write_csv(marker_avail_tbl, file.path(OUT_TAB, "marker_availability.csv"))
log_msg(sprintf("Marker availability: %d/%d cell types have >=2 markers",
                sum(marker_avail_tbl$available >= 2), nrow(marker_avail_tbl)))

# Drop cell types with <2 markers
markers_use <- markers_present[sapply(markers_present, length) >= 2]
if (length(markers_use) < 3) stop("Too few cell types have markers.")

# ==============================================================================
# [1] ANNOTATION SCORING - dual evidence
# ==============================================================================
log_msg("\n[1] Computing annotation scores (mean expression + module score) ...")

## (a) Cluster-level mean normalized expression per cell type
avg_expr <- AverageExpression(
  sc, assays = "RNA", layer = "data",
  features = unique(unlist(markers_use)),
  group.by = "seurat_clusters", verbose = FALSE
)$RNA

# z-score per gene across clusters, then average per cell type
z_expr <- t(scale(t(as.matrix(avg_expr))))
z_expr[is.na(z_expr)] <- 0

ct_score_mean <- sapply(markers_use, function(g) {
  g <- intersect(g, rownames(z_expr))
  if (length(g) == 0) return(rep(0, ncol(z_expr)))
  colMeans(z_expr[g, , drop = FALSE])
})
ct_score_mean <- as.data.frame(ct_score_mean)
rownames(ct_score_mean) <- colnames(z_expr)

## (b) Module score per cell, then aggregate by cluster
log_msg("  Running AddModuleScore for each cell type ...")
mod_names <- paste0("Mod_", names(markers_use))
sc <- AddModuleScore(sc, features = markers_use, name = "Mod_", seed = 123,
                     assay = "RNA", search = FALSE)
mod_cols <- paste0("Mod_", seq_along(markers_use))
colnames(sc[[]])[match(mod_cols, colnames(sc[[]]))] <- mod_names

ct_score_mod <- sc[[]] %>%
  as.data.frame() %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(mod_names), mean), .groups = "drop") %>%
  as.data.frame()
rownames(ct_score_mod) <- as.character(ct_score_mod$seurat_clusters)
ct_score_mod$seurat_clusters <- NULL
colnames(ct_score_mod) <- names(markers_use)
# z-score across cell types per cluster (column-wise)
ct_score_mod_z <- as.data.frame(t(scale(t(as.matrix(ct_score_mod)))))
ct_score_mod_z[is.na(ct_score_mod_z)] <- 0

# Align rows: AverageExpression appends "g" prefix
rownames(ct_score_mod_z) <- paste0("g", rownames(ct_score_mod_z))
common_rows <- intersect(rownames(ct_score_mean), rownames(ct_score_mod_z))
ct_score_mean   <- ct_score_mean[common_rows, , drop = FALSE]
ct_score_mod_z  <- ct_score_mod_z[common_rows, , drop = FALSE]

## (c) Consensus = average of two scores
ct_score_consensus <- (as.matrix(ct_score_mean) + as.matrix(ct_score_mod_z)) / 2

safe_write_csv(cbind(cluster = rownames(ct_score_consensus), as.data.frame(ct_score_consensus)),
               file.path(OUT_TAB, "cluster_celltype_score_matrix.csv"))

# ==============================================================================
# [2] INITIAL ANNOTATION + CONFIDENCE
# ==============================================================================
log_msg("\n[2] Assigning initial annotation per cluster ...")

assign_top <- function(scores) {
  apply(scores, 1, function(x) {
    o <- order(x, decreasing = TRUE)
    list(top1 = colnames(scores)[o[1]], top2 = colnames(scores)[o[2]],
         s1 = x[o[1]], s2 = x[o[2]], margin = x[o[1]] - x[o[2]])
  })
}
ann <- assign_top(ct_score_consensus)
ann_df <- data.frame(
  cluster = rownames(ct_score_consensus),
  initial_celltype = sapply(ann, `[[`, "top1"),
  second_choice    = sapply(ann, `[[`, "top2"),
  top_score        = sapply(ann, `[[`, "s1"),
  second_score     = sapply(ann, `[[`, "s2"),
  score_margin     = as.numeric(sapply(ann, `[[`, "margin")),
  stringsAsFactors = FALSE
)
ann_df$score_margin[is.na(ann_df$score_margin)] <- 0
cat("\nScore margins:", paste(round(ann_df$score_margin, 4), collapse = ", "), "\n")
ann_df$confidence <- cut(as.numeric(ann_df$score_margin),
                         breaks = c(-Inf, 0.15, 0.4, Inf),
                         labels = c("Low", "Medium", "High"))
safe_write_csv(ann_df, file.path(OUT_TAB, "cluster_annotation_initial.csv"))
log_msg(sprintf("  Confidence: High=%d  Medium=%d  Low=%d",
                sum(ann_df$confidence == "High"),
                sum(ann_df$confidence == "Medium"),
                sum(ann_df$confidence == "Low")))

# Map to cells
cl2ct <- setNames(ann_df$initial_celltype, ann_df$cluster)
sc@meta.data[["celltype"]] <- cl2ct[paste0("g", as.character(sc$seurat_clusters))]

# Log annotation
cat("\n=== INITIAL ANNOTATION ===\n")
for(i in seq_len(nrow(ann_df))) {
  cat(sprintf("  %-4s -> %-16s margin=%+.3f [%s]\n",
      ann_df$cluster[i], ann_df$initial_celltype[i],
      as.numeric(ann_df$score_margin[i]), as.character(ann_df$confidence[i])))
}
cat(sprintf("  NAs in celltype: %d\n", sum(is.na(sc@meta.data[["celltype"]]))))

# No Pericyte refinement (Cerebrovascular cells defined as the coarse brain vascular compartment)

# ==============================================================================
# ANNOTATION SUMMARY
# ==============================================================================
final_counts <- table(sc@meta.data[["celltype"]])
log_msg("Final cell type counts:")
for (n in names(final_counts)) log_msg(sprintf("  %-16s %d", n, final_counts[[n]]))

ct_levels <- names(sort(final_counts, decreasing = TRUE))
sc@meta.data[["celltype"]] <- factor(sc@meta.data[["celltype"]], levels = ct_levels)
ct_cols <- ct_colors_default[intersect(names(ct_colors_default), ct_levels)]

# ===== Reorder clusters by cell type =====
cluster_ct_map <- unique(data.frame(
  cluster = as.character(sc$seurat_clusters),
  celltype = as.character(sc$celltype),
  stringsAsFactors = FALSE
))
cluster_ct_map <- cluster_ct_map[order(
  match(cluster_ct_map$celltype, ct_levels),
  as.numeric(cluster_ct_map$cluster)
), ]
sc$cluster_by_ct <- factor(sc$seurat_clusters,
  levels = cluster_ct_map$cluster)
log_msg(sprintf("Cluster order by cell type: %s",
  paste(paste0(cluster_ct_map$cluster, "(", cluster_ct_map$celltype, ")"),
        collapse = " -> ")))

# Also reorder the consensus matrix rows for Fig 3 heatmap
ct_score_consensus <- ct_score_consensus[paste0("g", cluster_ct_map$cluster), , drop = FALSE]
# Reorder ann_df too
ann_df <- ann_df[match(paste0("g", cluster_ct_map$cluster), ann_df$cluster), ]

# ==============================================================================
# [3] MAIN FIGURES
# ==============================================================================
log_msg("\n[3] Building main figures ...")

## Fig 1: UMAP triptych (celltype / group / dataset)
umap_df <- as.data.frame(Embeddings(sc, "umap"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$celltype <- sc@meta.data[["celltype"]]
  if ("Group" %in% colnames(sc[[]])) {
    umap_df$group <- standardize_group(sc$Group)
  } else if ("group" %in% colnames(sc[[]])) {
    umap_df$group <- standardize_group(sc$group)
  }
if ("Dataset" %in% colnames(sc[[]])) umap_df$dataset <- sc$Dataset

geom_pt <- function(...) {
  if (has_ggrastr) ggrastr::geom_point_rast(..., raster.dpi = 300)
  else geom_point(...)
}

p1a <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = celltype)) +
  geom_pt(size = 0.15, alpha = 0.75) +
  scale_color_manual(values = ct_cols) +
  guides(color = guide_legend(override.aes = list(size = 2.2), ncol = 1)) +
  labs(title = "Cell types", x = "UMAP 1", y = "UMAP 2", color = NULL) +
  publication_theme(7) + theme(aspect.ratio = 1)

panels <- list(p1a)
if ("group" %in% colnames(umap_df)) {
  p1b <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = group)) +
    geom_pt(size = 0.08, alpha = 0.55) +
    scale_color_manual(values = group_palette) +
    facet_wrap(~ group, nrow = 1) +
    guides(color = "none") +
    labs(title = "Condition", x = "UMAP 1", y = "UMAP 2", color = NULL) +
    publication_theme(7) + theme(aspect.ratio = 1, strip.text = element_text(size = 7, face = "bold"))
  panels[[length(panels) + 1]] <- p1b
}
if ("dataset" %in% colnames(umap_df)) {
  ds_n <- length(unique(umap_df$dataset))
  p1c <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = dataset)) +
    geom_pt(size = 0.08, alpha = 0.55) +
    scale_color_manual(values = publication_palette(ds_n)) +
    facet_wrap(~ dataset, nrow = 1) +
    guides(color = "none") +
    labs(title = "Dataset", x = "UMAP 1", y = "UMAP 2", color = NULL) +
    publication_theme(7) + theme(aspect.ratio = 1, strip.text = element_text(size = 7, face = "bold"))
  panels[[length(panels) + 1]] <- p1c
}

fig1 <- wrap_plots(panels, ncol = length(panels)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9))
save_fig(fig1, "Fig1_UMAP_overview", OUT_FIG,
         w = 4.8 * length(panels) + 1.2, h = 4.6)
save_source(umap_df, "Fig1_UMAP_overview")

## Fig 2: Canonical marker dot plot + stacked violin
markers_by_ct <- markers_use[intersect(names(markers_use), ct_levels)]

# Fig2a: DotPlot (x=markers grouped by ct, y=clusters by ct, ct names as headers)
# Extract DotPlot data, reorder axes manually
sc$clust_ordered <- factor(sc$seurat_clusters, levels = cluster_ct_map$cluster)
p_tmp <- DotPlot(sc, features = markers_by_ct, group.by = "clust_ordered",
                 cols = c("grey92", "#E64B35"))
df <- p_tmp$data
correct_order <- unlist(markers_by_ct[intersect(names(markers_by_ct), ct_levels)])
df$features.plot <- factor(df$features.plot, levels = correct_order)
df$feature.groups <- factor(df$feature.groups, levels = ct_levels)
df$id <- factor(df$id, levels = rev(levels(df$id)))

p2a <- ggplot(df, aes(features.plot, id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 2.5), name = "% expr") +
  scale_color_gradient(low = "grey92", high = "#E64B35", name = "Avg. expr") +
  facet_grid(~ feature.groups, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL, title = "Canonical markers") +
  publication_theme(5.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4.5),
        axis.text.y = element_text(size = 4.5),
        strip.text = element_text(size = 5, face = "bold"),
        strip.background = element_rect(fill = "grey95"),
        panel.spacing = unit(0.1, "lines"),
        plot.margin = margin(1, 1, 1, 1))

# Top-1 marker per cell type for stacked violin, in ct_levels order
top1_markers <- sapply(ct_levels, function(ct) {
  if (ct %in% names(markers_use)) intersect(markers_use[[ct]], rownames(sc))[1] else NA
})
top1_markers <- top1_markers[!is.na(top1_markers)]
top1_markers <- unname(top1_markers)

p2b <- VlnPlot(sc, features = unname(top1_markers),
               stack = TRUE, flip = TRUE, fill.by = "feature") +
  NoLegend() + labs(title = "Top markers", x = NULL, y = NULL) +
  publication_theme(6.5)

fig2 <- (p2a / p2b) + plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9))
save_fig(fig2, "Fig2_Markers", OUT_FIG, w = 8.5, h = 7)
save_source(p2a$data, "Fig2_Markers_dotplot")

## Fig 3: Annotation score heatmap
# Y-axis: clusters ordered by cell type (Oligo at bottom, Endo at top)
# X-axis: cell types in ct_levels order
cluster_order_y <- rev(paste0("g", cluster_ct_map$cluster))  # bottom = first in ct_levels

score_long <- ct_score_consensus %>%
  as.data.frame() %>%
  rownames_to_column("cluster") %>%
  pivot_longer(-cluster, names_to = "celltype", values_to = "score")
score_long$cluster  <- factor(score_long$cluster, levels = cluster_order_y)
score_long$celltype <- factor(score_long$celltype, levels = ct_levels)

p3 <- ggplot(score_long, aes(celltype, cluster, fill = score)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                       midpoint = 0, name = "Consensus\nscore (z)") +
  labs(x = NULL, y = "Cluster", title = "Cluster x cell type scoring") +
  publication_theme(6.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ann_df$cluster <- factor(ann_df$cluster, levels = levels(score_long$cluster))
p3b <- ggplot(ann_df, aes(x = 1, y = cluster, fill = confidence)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(High = "#1B7837", Medium = "#FDB863", Low = "#B2182B")) +
  labs(x = "Conf.", y = NULL, fill = NULL) +
  publication_theme(6.5) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank())

fig3 <- p3 + p3b + plot_layout(widths = c(8, 1))
save_fig(fig3, "Fig3_Annotation_scores", OUT_FIG, w = 7.5, h = 5.5)
save_source(score_long, "Fig3_Annotation_scores")
save_source(ann_df, "Fig3_Annotation_confidence")

## Fig 4: Composition
comp_panels <- list()
meta_df <- sc[[]] %>% as.data.frame()

if ("Group" %in% colnames(meta_df)) {
  comp_g <- meta_df %>%
    dplyr::count(Group, celltype) %>%
    group_by(Group) %>% mutate(prop = n / sum(n)) %>% ungroup()
  comp_panels$group <- ggplot(comp_g, aes(Group, prop, fill = celltype)) +
    geom_col(color = "white", linewidth = 0.2) +
    scale_fill_manual(values = ct_cols) +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
    labs(x = NULL, y = "Proportion", title = "By group", fill = NULL) +
    publication_theme(6.5)
  save_source(comp_g, "Fig4_composition_group")
}
if ("Dataset" %in% colnames(meta_df)) {
  comp_d <- meta_df %>%
    dplyr::count(Dataset, celltype) %>%
    group_by(Dataset) %>% mutate(prop = n / sum(n)) %>% ungroup()
  comp_panels$dataset <- ggplot(comp_d, aes(Dataset, prop, fill = celltype)) +
    geom_col(color = "white", linewidth = 0.2) +
    scale_fill_manual(values = ct_cols) +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
    labs(x = NULL, y = "Proportion", title = "By dataset", fill = NULL) +
    publication_theme(6.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_source(comp_d, "Fig4_composition_dataset")
}
if ("SampleID" %in% colnames(meta_df)) {
  comp_s <- meta_df %>%
    dplyr::count(SampleID, celltype) %>%
    group_by(SampleID) %>% mutate(prop = n / sum(n)) %>% ungroup()
  p_comp_sample <- ggplot(comp_s, aes(SampleID, prop, fill = celltype)) +
    geom_col(color = "white", linewidth = 0.15) +
    scale_fill_manual(values = ct_cols) +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
    labs(x = NULL, y = "Proportion", title = "Cell type composition by sample", fill = NULL) +
    publication_theme(6) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))
  save_fig(p_comp_sample, "ED_Fig7_Composition_by_sample", OUT_ED, w = 12, h = 4.8)
  save_source(comp_s, "Fig4_composition_sample")
}

if (length(comp_panels) > 0) {
  fig4_main <- wrap_plots(comp_panels, ncol = length(comp_panels), guides = "collect") &
    theme(legend.position = "right")
  fig4 <- fig4_main +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 9), legend.position = "right")
  save_fig(fig4, "Fig4_Composition", OUT_FIG,
           w = 3.2 * length(comp_panels) + 1.5, h = 4)
}

# ==============================================================================
# [4] EXTENDED DATA FIGURES
# ==============================================================================
log_msg("\n[4] Building Extended Data figures ...")

## ED Fig 1: QC violins per cluster
qc_metrics <- intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
                        colnames(sc[[]]))
if (length(qc_metrics) > 0) {
  Idents(sc) <- "seurat_clusters"
  p_ed1 <- VlnPlot(sc, features = qc_metrics, pt.size = 0,
                   ncol = 1, log = FALSE) &
    publication_theme(6) &
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))
  save_fig(p_ed1, "ED_Fig1_QC_per_cluster", OUT_ED,
           w = 8, h = 1.8 * length(qc_metrics))
  qc_long <- sc[[]] %>% as.data.frame() %>%
    select(seurat_clusters, all_of(qc_metrics))
  save_source(qc_long, "ED_Fig1_QC_per_cluster")
}

## ED Fig 2: UMAP per QC metric
if (length(qc_metrics) > 0) {
p_ed2 <- FeaturePlot(sc, features = qc_metrics, reduction = "umap",
                       raster = TRUE, order = TRUE, ncol = min(3, length(qc_metrics))) &
    scale_color_gradientn(colors = c("#0d1a50", "#4b6db8", "#fdae61", "#d7191c")) &
    publication_theme(6) &
    theme(aspect.ratio = 1)
  save_fig(p_ed2, "ED_Fig2_QC_UMAP", OUT_ED,
           w = 3.5 * min(3, length(qc_metrics)), h = 3.5 * ceiling(length(qc_metrics)/3))
}

## ED Fig 3: Full marker DotPlot (x=markers by ct, y=clusters by ct)
ed3_markers <- markers_use[intersect(names(markers_use), ct_levels)]
p_tmp3 <- DotPlot(sc, features = ed3_markers, group.by = "clust_ordered",
                   cols = c("grey92", "#E64B35"), dot.scale = 3)
df3 <- p_tmp3$data
df3$features.plot <- factor(df3$features.plot, levels = unlist(ed3_markers[intersect(names(ed3_markers), ct_levels)]))
df3$feature.groups <- factor(df3$feature.groups, levels = ct_levels)
df3$id <- factor(df3$id, levels = rev(levels(df3$id)))

p_ed3 <- ggplot(df3, aes(features.plot, id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 2), name = "% expr") +
  scale_color_gradient(low = "grey92", high = "#E64B35", name = "Avg. expr") +
  facet_grid(~ feature.groups, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL, title = NULL) +
  publication_theme(5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
        axis.text.y = element_text(size = 4),
        strip.text = element_text(size = 4.5, face = "bold"),
        strip.background = element_rect(fill = "grey95"),
        panel.spacing = unit(0.05, "lines"),
        legend.key.size = unit(0.2, "cm"),
        plot.margin = margin(1, 1, 1, 1))
save_fig(p_ed3, "ED_Fig3_Marker_dotplot_full", OUT_ED,
         w = 0.3 * length(unique(df3$features.plot)) + 1.2,
         h = 0.18 * length(unique(df3$id)) + 0.8)

## ED Fig 4: Annotation consensus scatter
consensus_df <- data.frame(
  cluster = rep(rownames(ct_score_mean), ncol(ct_score_mean)),
  celltype = rep(colnames(ct_score_mean), each = nrow(ct_score_mean)),
  mean_expr_z = as.vector(as.matrix(ct_score_mean)),
  module_z    = as.vector(as.matrix(ct_score_mod_z[rownames(ct_score_mean), colnames(ct_score_mean)]))
)
p_ed4 <- ggplot(consensus_df, aes(mean_expr_z, module_z, color = celltype)) +
  geom_point(size = 0.9, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = publication_palette(length(unique(consensus_df$celltype)))) +
  labs(x = "Mean-expression z-score", y = "Module z-score",
       title = "Annotation method consensus", color = NULL) +
  publication_theme(7) +
  theme(legend.position = "right", plot.title = element_text(size = 8, face = "bold"))
save_fig(p_ed4, "ED_Fig4_Annotation_consensus", OUT_ED, w = 6, h = 4.5)
save_source(consensus_df, "ED_Fig4_Annotation_consensus")

## ED Fig 5: Inhibitory detail UMAP
inh_cells <- which(sc$celltype %in% c("Inhibitory", "Excitatory"))
p_ed5 <- DimPlot(sc, cells = inh_cells, group.by = "celltype", reduction = "umap",
                 cols = c(Inhibitory = "#D55E00", Excitatory = "#0072B2"),
                 pt.size = 0.2, label = TRUE, label.size = 3, repel = TRUE, raster = TRUE) +
  labs(title = "Excitatory and inhibitory neurons", x = "UMAP 1", y = "UMAP 2") +
  publication_theme(7) + theme(aspect.ratio = 1, legend.position = "right")
save_fig(p_ed5, "ED_Fig5_Inh_vs_Ex", OUT_ED, w = 7, h = 5.5)
save_source(umap_df[inh_cells, ], "ED_Fig5_Inh_vs_Ex")

## ED Fig 6: UMAP colored by cell type
p_ed6 <- DimPlot(sc, group.by = "celltype", reduction = "umap",
                 cols = ct_cols, raster = TRUE) &
  publication_theme(6) & theme(aspect.ratio = 1) &
  labs(title = "Cell type")
save_fig(p_ed6, "ED_Fig6_CellType_UMAP", OUT_ED, w = 5.2, h = 4.8)
log_msg("  ED Fig 6: UMAP colored by cell type")

# ==============================================================================
# [5] CELL-LEVEL METADATA EXPORT
# ==============================================================================
log_msg("\n[5] Exporting cell-level metadata ...")
cell_meta <- sc[[]] %>% as.data.frame() %>% rownames_to_column("cell_id")
keep_cols <- c("cell_id", "SampleID", "Group", "Dataset",
               "seurat_clusters", "celltype",
               "nFeature_RNA", "nCount_RNA", "percent.mt",
               grep("Mod_", colnames(cell_meta), value = TRUE))
keep_cols <- intersect(keep_cols, colnames(cell_meta))
safe_write_csv(cell_meta[, keep_cols], file.path(OUT_TAB, "cell_metadata_annotated.csv"))

if (exists("comp_g")) safe_write_csv(comp_g, file.path(OUT_TAB, "composition_by_group.csv"))
if (exists("comp_d")) safe_write_csv(comp_d, file.path(OUT_TAB, "composition_by_dataset.csv"))

# ==============================================================================
# [6] OUTPUT MANIFEST
# ==============================================================================
manifest <- data.frame(
  file = c(
    "figures_main/Fig1_UMAP_overview.{pdf,png,tiff}",
    "figures_main/Fig2_Markers.{pdf,png,tiff}",
    "figures_main/Fig3_Annotation_scores.{pdf,png,tiff}",
    "figures_main/Fig4_Composition.{pdf,png,tiff}",
    "figures_extended_data/ED_Fig1_QC_per_cluster.*",
    "figures_extended_data/ED_Fig2_QC_UMAP.*",
    "figures_extended_data/ED_Fig3_Marker_dotplot_full.*",
    "figures_extended_data/ED_Fig4_Annotation_consensus.*",
    "figures_extended_data/ED_Fig5_Inh_vs_Ex.*",
    "figures_extended_data/ED_Fig6_CellType_UMAP.{pdf,png,svg,tiff}",
    "tables/marker_availability.csv",
    "tables/cluster_celltype_score_matrix.csv",
    "tables/cluster_annotation_initial.csv",
    "tables/composition_by_group.csv",
    "tables/composition_by_dataset.csv",
    "tables/cell_metadata_annotated.csv",
    "source_data/*_source_data.csv",
    "rds/scRNA_annotated.qs"
  ),
  description = c(
    "Main Fig 1: UMAPs by celltype/group/dataset",
    "Main Fig 2: Canonical marker dot plot + stacked violin",
    "Main Fig 3: Cluster x celltype consensus heatmap + confidence",
    "Main Fig 4: Cell type composition by group/dataset/sample",
    "ED Fig 1: QC violins per cluster",
    "ED Fig 2: QC metrics on UMAP",
    "ED Fig 3: Full marker dot plot (cluster level)",
    "ED Fig 4: Annotation consensus (mean expr vs module score)",
    "ED Fig 5: Inhibitory vs Excitatory UMAP detail",
    "ED Fig 6: Doublet diagnostics (if available)",
    "Marker availability per cell type",
    "Cluster x celltype consensus score matrix",
    "Initial cluster annotation with confidence",
    "Composition by group",
    "Composition per dataset",
    "Per-cell annotated metadata",
    "Per-figure source data (publication requirement)",
    "Annotated Seurat object"
  ),
  stringsAsFactors = FALSE
)
safe_write_csv(manifest, file.path(OUT_TAB, "output_manifest.csv"))

# ==============================================================================
# [7] SAVE OBJECT + SESSION INFO
# ==============================================================================
log_msg("\n[7] Saving annotated object ...")
Idents(sc) <- "celltype"
save_annot_object(sc, file.path(OUT_RDS, "scRNA_annotated.qs"))
capture.output(sessionInfo(), file = file.path(OUT_LOG, "sessionInfo_node02.txt"))

summary_txt <- c(
  "===== node02 annotation summary =====",
  sprintf("Input:    %s", IN_FILE),
  sprintf("Output:   %s", OUT),
  sprintf("Cells:    %d", ncol(sc)),
  sprintf("Features: %d", nrow(sc)),
  sprintf("Clusters: %d", n_clusters),
  sprintf("Cell types (7): %s", paste(names(final_counts), collapse = ", ")),
  "",
  sprintf("High-confidence clusters:   %d", sum(ann_df$confidence == "High")),
  sprintf("Medium-confidence clusters: %d", sum(ann_df$confidence == "Medium")),
  sprintf("Low-confidence clusters:    %d", sum(ann_df$confidence == "Low")),
  "",
  "Marker panel: MBP/PLP1/MOBP (Oligo), PDGFRA/VCAN/CSPG4 (OPC),",
  "             AQP4/GFAP/ADGRV1 (Astro), C3/P2RY12/TREM2 (Micro),",
  "             CLDN5/FLT1/VWF (cerebrovascular), GAD1/GAD2/PCDH15 (Inhib),",
  "             CAMK2A/SLC17A7/SATB2 (Excit)",
  "References: Mathys 2019 publication, Lau et al. 2020 PNAS",
  "",
  "Key outputs:",
  "  Main Fig 1-4 + ED Fig 1-6",
  "  Per-figure source_data CSVs (publication output)",
  "",
  "Notes:",
  "  - Annotation uses dual evidence (mean-expr z + module score z) consensus.",
  "  - Cerebrovascular cells defined as the coarse brain vascular compartment (S/N ratio < 2x, no independent cluster).",
  "  - Low-confidence clusters should be re-examined in subcluster analysis."
)
writeLines(summary_txt, con = file.path(OUT_LOG, "run_summary_node02.txt"))

log_msg("\n===== ALL DONE =====")
log_msg("Main: Fig1-4 | Extended Data: ED1-6 | Source data + manifest written.")
# ==============================================================================
# [6.5] BUILD FINAL PUBLICATION-BUNDLE (publication-style)
# ============================================================================== 
OUT_FINAL <- file.path(OUT, "outputs")
OUT_F_MAIN <- file.path(OUT_FINAL, "main_figures")
OUT_F_SUPP <- file.path(OUT_FINAL, "supplementary_figures")
OUT_F_EXTD <- file.path(OUT_FINAL, "extended_data_figures")
OUT_T_MAIN <- file.path(OUT_FINAL, "main_tables")
OUT_T_SUPP <- file.path(OUT_FINAL, "supplementary_tables")
OUT_METH <- file.path(OUT_FINAL, "methods_results")
OUT_FIG_EDIT <- file.path(OUT_FINAL, "editable_panels")
for (d in c(OUT_F_MAIN, OUT_F_SUPP, OUT_F_EXTD, OUT_T_MAIN, OUT_T_SUPP, OUT_METH, OUT_FIG_EDIT)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

copy_multi_format <- function(src_base, dst_base, src_dir, dst_dir) {
  for (ext in c("pdf", "png", "svg", "tiff")) {
    src <- file.path(src_dir, paste0(src_base, ".", ext))
    dst <- file.path(dst_dir, paste0(dst_base, ".", ext))
    safe_copy_file(src, dst)
  }
}

copy_editable <- function(src_base, dst_base, src_dir, dst_dir) {
  for (ext in c("pdf", "svg")) {
    src <- file.path(src_dir, paste0(src_base, ".", ext))
    dst <- file.path(dst_dir, paste0(dst_base, ".", ext))
    safe_copy_file(src, dst)
  }
}

copy_multi_format("Fig1_UMAP_overview", "Figure_2_node02_panelA_umap_annotation", OUT_FIG, OUT_F_MAIN)
copy_multi_format("Fig2_Markers", "Figure_2_node02_panelB_marker_validation", OUT_FIG, OUT_F_MAIN)
copy_multi_format("Fig3_Annotation_scores", "Figure_2_node02_panelC_annotation_scoring", OUT_FIG, OUT_F_MAIN)
if (file.exists(file.path(OUT_FIG, "Fig4_Composition.pdf"))) {
  copy_multi_format("Fig4_Composition", "Figure_2_node02_panelD_celltype_composition", OUT_FIG, OUT_F_MAIN)
}

if (exists("p1a") && exists("p2a") && exists("fig3") && exists("fig4_main")) {
  final_main <- (wrap_elements(full = p1a) | wrap_elements(full = p2a)) /
    (wrap_elements(full = fig3) | wrap_elements(full = fig4_main)) +
    plot_layout(heights = c(1, 1.05)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 11),
          plot.title = element_text(size = 8, face = "bold"))
  save_fig(final_main, "Figure_2_node02_celltype_annotation_FINAL", OUT_F_MAIN, w = 14, h = 11)
}

copy_multi_format("ED_Fig1_QC_per_cluster", "Extended_Data_Figure_1_node02_qc_per_cluster", OUT_ED, OUT_F_EXTD)
copy_multi_format("ED_Fig2_QC_UMAP", "Extended_Data_Figure_2_node02_qc_embedding", OUT_ED, OUT_F_EXTD)
copy_multi_format("ED_Fig3_Marker_dotplot_full", "Extended_Data_Figure_3_node02_marker_dotplot_full", OUT_ED, OUT_F_EXTD)
copy_multi_format("ED_Fig4_Annotation_consensus", "Extended_Data_Figure_4_02_annotation_consensus", OUT_ED, OUT_F_EXTD)
copy_multi_format("ED_Fig6_CellType_UMAP", "Extended_Data_Figure_6_node02_celltype_umap", OUT_ED, OUT_F_EXTD)
copy_multi_format("ED_Fig7_Composition_by_sample", "Extended_Data_Figure_7_node02_composition_by_sample", OUT_ED, OUT_F_EXTD)

copy_multi_format("ED_Fig1_QC_per_cluster", "Figure_S2_node02_qc_per_cluster", OUT_ED, OUT_F_SUPP)
copy_multi_format("ED_Fig2_QC_UMAP", "Figure_S3_node02_qc_embedding", OUT_ED, OUT_F_SUPP)
copy_multi_format("ED_Fig3_Marker_dotplot_full", "Figure_S4_node02_marker_dotplot_full", OUT_ED, OUT_F_SUPP)
copy_multi_format("ED_Fig4_Annotation_consensus", "Figure_S5_02_annotation_consensus", OUT_ED, OUT_F_SUPP)
copy_multi_format("ED_Fig5_Inh_vs_Ex", "Figure_S6_node02_excitatory_inhibitory_umap", OUT_ED, OUT_F_SUPP)
copy_multi_format("ED_Fig6_CellType_UMAP", "Figure_S7_node02_celltype_umap", OUT_ED, OUT_F_SUPP)
copy_multi_format("ED_Fig7_Composition_by_sample", "Figure_S8_node02_composition_by_sample", OUT_ED, OUT_F_SUPP)

copy_editable("Fig1_UMAP_overview", "node02_Fig1_UMAP_overview_edit", OUT_FIG, OUT_FIG_EDIT)
copy_editable("Fig2_Markers", "node02_Fig2_Markers_edit", OUT_FIG, OUT_FIG_EDIT)
copy_editable("Fig3_Annotation_scores", "node02_Fig3_Annotation_scores_edit", OUT_FIG, OUT_FIG_EDIT)
copy_editable("Fig4_Composition", "node02_Fig4_Composition_edit", OUT_FIG, OUT_FIG_EDIT)
copy_editable("ED_Fig1_QC_per_cluster", "node02_ED_Fig1_QC_per_cluster_edit", OUT_ED, OUT_FIG_EDIT)
copy_editable("ED_Fig2_QC_UMAP", "node02_ED_Fig2_QC_UMAP_edit", OUT_ED, OUT_FIG_EDIT)
copy_editable("ED_Fig3_Marker_dotplot_full", "node02_ED_Fig3_Marker_dotplot_full_edit", OUT_ED, OUT_FIG_EDIT)
copy_editable("ED_Fig4_Annotation_consensus", "node02_ED_Fig4_Annotation_consensus_edit", OUT_ED, OUT_FIG_EDIT)
copy_editable("ED_Fig5_Inh_vs_Ex", "node02_ED_Fig5_Inh_vs_Ex_edit", OUT_ED, OUT_FIG_EDIT)
copy_editable("ED_Fig6_CellType_UMAP", "node02_ED_Fig6_CellType_UMAP_edit", OUT_ED, OUT_FIG_EDIT)
copy_editable("ED_Fig7_Composition_by_sample", "node02_ED_Fig7_Composition_by_sample_edit", OUT_ED, OUT_FIG_EDIT)

if (exists("comp_g")) {
  comp_summary <- comp_g %>%
    dplyr::mutate(group = as.character(Group)) %>%
    dplyr::filter(group %in% c("CN", "AD")) %>%
    dplyr::select(group, celltype, n, prop)
  if (nrow(comp_summary) > 0) {
    safe_write_csv(comp_summary, file.path(OUT_T_MAIN, "Table_2_node02_celltype_composition_summary.csv"))
  }
}

if (file.exists(file.path(OUT_TAB, "marker_availability.csv"))) {
  safe_copy_file(
    file.path(OUT_TAB, "marker_availability.csv"),
    file.path(OUT_T_SUPP, "Supplementary_Table_1_node02_marker_availability.csv")
  )
}
if (file.exists(file.path(OUT_TAB, "cluster_annotation_initial.csv"))) {
  safe_copy_file(
    file.path(OUT_TAB, "cluster_annotation_initial.csv"),
    file.path(OUT_T_SUPP, "Supplementary_Table_2_node02_cluster_annotation_initial.csv")
  )
}
if (file.exists(file.path(OUT_TAB, "cluster_celltype_score_matrix.csv"))) {
  safe_copy_file(
    file.path(OUT_TAB, "cluster_celltype_score_matrix.csv"),
    file.path(OUT_T_SUPP, "Supplementary_Table_3_node02_cluster_celltype_score_matrix.csv")
  )
}
if (file.exists(file.path(OUT_TAB, "composition_by_dataset.csv"))) {
  safe_copy_file(
    file.path(OUT_TAB, "composition_by_dataset.csv"),
    file.path(OUT_T_SUPP, "Supplementary_Table_4_node02_composition_by_dataset.csv")
  )
}
if (file.exists(file.path(OUT_TAB, "cell_metadata_annotated.csv"))) {
  safe_copy_file(
    file.path(OUT_TAB, "cell_metadata_annotated.csv"),
    file.path(OUT_T_SUPP, "Supplementary_Table_5_node02_cell_metadata_annotated.csv")
  )
}

safe_copy_file(file.path(OUT_LOG, "run_summary_node02.txt"), file.path(OUT_METH, "node02_methods_and_results.md"))
supp_note <- c(
  "## node02 Cell type annotation (publication-ready summary)",
  "",
  "### Figure legends",
  "- Figure 2A. UMAP overview annotated by cerebrovascular-aware cell types and sample-level projections (cell type / CN-AD group / dataset).",
  "- Figure 2B. Canonical marker validation with dotplot + stacked top marker violin.",
  "- Figure 2C. Dual-evidence annotation confidence map (cluster-level consensus score + confidence bar).",
  "- Figure 2D. Cell type composition by group (CN vs AD) and dataset.",
  "- Extended Data Figure 1. QC distribution by cluster (nFeature/nCount/percent.mt).",
  "- Extended Data Figure 2. QC metrics on UMAP embeddings.",
  "- Extended Data Figure 3. Full cluster-level marker dotplot.",
  "- Extended Data Figure 4. Annotation consensus between mean-expression and module-score evidence.",
  "- Extended Data Figure 6. Cell-type UMAP.",
  "",
  "### Table legends",
  "- Table 2. Grouped cell-type composition summary by CN/AD (counts and proportions).",
  "- Supplementary Table 1. Marker availability by cell type.",
  "- Supplementary Table 2. Initial cluster-cell type annotations with confidence.",
  "- Supplementary Table 3. Cluster-cell type consensus score matrix.",
  "- Supplementary Table 4. Dataset-level cell type composition.",
  "- Supplementary Table 5. Per-cell annotated metadata."
)
writeLines(supp_note, con = file.path(OUT_METH, "node02_final_legends_methods_results.md"))
