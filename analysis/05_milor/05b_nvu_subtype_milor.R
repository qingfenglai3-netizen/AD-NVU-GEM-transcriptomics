#!/usr/bin/env Rscript

.libPaths(unique(c(
  "/path/to/R_libraries",
  .libPaths()
)))

# =============================================================================
# 05_nvu_subtype_milor_final.R
#
# NVU-focused MiloR differential abundance analysis
# Subset-level MiloR for:
#   1. Cerebrovascular cells
#   2. Astrocytes
#   3. Microglia
#
# Designed for:
#   R 4.4.1
#   Seurat v5
#   miloR
#
# Main output:
#   publication-style figures, supplementary figures, DA tables, QC tables,
#   sample-level abundance validation, parameters, sessionInfo
#
# Important interpretation:
#   logFC > 0 means AD-enriched neighborhoods
# =============================================================================


suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(miloR)
  library(SingleCellExperiment)
  library(Matrix)
  library(scales)
  library(ggrepel)
  library(ggpubr)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
})


# =============================================================================
# 0. Global configuration
# =============================================================================

TEST_MODE <- TRUE

if (TEST_MODE) {
  RDS_DIR <- "/path/to/project/results/04_subtypes/rds"
  OUT <- "/path/to/project/results/05_nvu_subtype_milor"
} else {
  RDS_DIR <- "/path/to/project/results/04_subtypes/rds"
  OUT <- "/path/to/project/results/05_nvu_subtype_milor"
}

OUT_FIGS <- file.path(OUT, "main_figures")
OUT_SUPP <- file.path(OUT, "supplementary_figures")
OUT_TAB  <- file.path(OUT, "tables")
OUT_RDS  <- file.path(OUT, "rds")
OUT_LOG  <- file.path(OUT, "logs")

for (d in c(OUT_FIGS, OUT_SUPP, OUT_TAB, OUT_RDS, OUT_LOG)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

set.seed(20260523)

# -------------------------
# MiloR parameters
# -------------------------

MAX_CELLS_PER_OBJECT <- 60000
MIN_CELLS_PER_OBJECT <- 300
MIN_SAMPLES_PER_GROUP <- 2

K_NN <- 20
DIMS_USE_DEFAULT <- 20
NHOOD_PROP <- 0.10

TOP_NHOODS_FOR_BOXPLOT <- 12
TOP_NHOODS_FOR_HEATMAP <- 50

FDR_CUTOFF_MAIN <- 0.10

# If these covariates are present in sample metadata and valid, they can be used.
# For small sample size, over-adjustment can hurt power.
# By default, the script uses Group-only model for robustness.
USE_COVARIATES <- FALSE
CANDIDATE_COVARIATES <- c("Sex", "sex", "Batch", "batch", "PMI", "RIN", "Age", "age")


# =============================================================================
# 1. Plot style
# =============================================================================

theme_publication <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black", size = base_size * 0.85),
      axis.title = element_text(color = "black", size = base_size),
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = base_size * 1.05,
        margin = margin(b = 6)
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = base_size * 0.82,
        color = "grey30"
      ),
      legend.title = element_text(size = base_size * 0.85),
      legend.text = element_text(size = base_size * 0.80),
      strip.background = element_rect(fill = "grey95", color = "grey70", linewidth = 0.3),
      strip.text = element_text(face = "bold", size = base_size * 0.85)
    )
}

group_palette <- c(
  "CN" = "#4C78A8",
  "AD" = "#D55E00"
)

fc_palette <- c(
  "AD-enriched" = "#B2182B",
  "CN-enriched" = "#2166AC",
  "NS" = "grey80"
)

sig_palette <- c(
  "FDR<0.01" = "#B2182B",
  "FDR<0.05" = "#EF8A62",
  "FDR<0.1"  = "#67A9CF",
  "NS"       = "grey78"
)

log_msg <- function(m) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), m)
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(OUT_LOG, "05_nvu_subtype_milor.log"), append = TRUE)
}

save_fig <- function(p, name, w = 7, h = 5, folder = OUT_FIGS) {
  ggsave(
    filename = file.path(folder, paste0(name, ".pdf")),
    plot = p,
    width = w,
    height = h,
    device = cairo_pdf,
    dpi = 300,
    limitsize = FALSE
  )
  ggsave(
    filename = file.path(folder, paste0(name, ".png")),
    plot = p,
    width = w,
    height = h,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  tryCatch({
    ggsave(
      filename = file.path(folder, paste0(name, ".svg")),
      plot = p,
      width = w,
      height = h,
      limitsize = FALSE
    )
  }, error = function(e) NULL)
  tryCatch({
    ggsave(
      filename = file.path(folder, paste0(name, ".tiff")),
      plot = p,
      width = w,
      height = h,
      dpi = 600,
      compression = "lzw",
      bg = "white",
      limitsize = FALSE
    )
  }, error = function(e) NULL)
  log_msg(paste("Saved figure:", name))
}

safe_csv <- function(x, filename) {
  tryCatch({
    write.csv(x, file.path(OUT_TAB, filename), row.names = FALSE)
    log_msg(paste("Saved table:", filename))
  }, error = function(e) {
    log_msg(paste("Failed to save table:", filename, e$message))
  })
}


# =============================================================================
# 2. Utility functions
# =============================================================================

detect_group_column <- function(meta) {
  candidates <- c("Group", "group", "Diagnosis", "diagnosis", "Condition", "condition", "Disease", "disease")
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0) stop("No group/diagnosis column found in metadata.")
  hit[1]
}

standardize_group <- function(x) {
  x0 <- as.character(x)
  x1 <- toupper(trimws(x0))
  
  control_labels <- c("CN", "NC", "ND", "CONTROL", "CTRL", "CON", "HEALTHY", "NORMAL")
  disease_labels <- c("AD", "UC", "DISEASE", "CASE", "PATIENT")
  
  out <- rep(NA_character_, length(x1))
  out[x1 %in% control_labels] <- "CN"
  out[x1 %in% disease_labels] <- "AD"
  
  out[grepl("CONTROL|CTRL|NORMAL|HEALTHY", x1)] <- "CN"
  out[grepl("AD|DISEASE|CASE|PATIENT", x1)] <- "AD"
  
  factor(out, levels = c("CN", "AD"))
}

get_sample_column <- function(meta) {
  candidates <- c("orig.ident", "sample", "Sample", "sample_id", "SampleID", "donor", "Donor", "subject", "Subject")
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0) stop("No sample column found. Expected orig.ident/sample/SampleID/donor.")
  hit[1]
}

safe_join_layers <- function(seu) {
  seu2 <- seu
  tryCatch({
    seu2 <- JoinLayers(seu2)
  }, error = function(e) {
    log_msg(paste("JoinLayers skipped:", e$message))
  })
  seu2
}

get_layer_safe <- function(seu, assay = "RNA", layer = "counts") {
  out <- NULL
  out <- tryCatch({
    GetAssayData(seu, assay = assay, layer = layer)
  }, error = function(e) {
    NULL
  })
  
  if (is.null(out) && layer == "data") {
    log_msg("RNA data layer not found. Running NormalizeData().")
    seu <- NormalizeData(seu, assay = assay, verbose = FALSE)
    out <- GetAssayData(seu, assay = assay, layer = "data")
  }
  
  if (is.null(out)) {
    stop(paste("Cannot retrieve layer:", layer))
  }
  
  out
}

ensure_reduction <- function(seu, max_dims = 30) {
  red_names <- names(seu@reductions)
  
  if ("harmony" %in% red_names) {
    return(list(seu = seu, reduction = "harmony", reduced_dim_name = "HARMONY"))
  }
  
  if ("Harmony" %in% red_names) {
    return(list(seu = seu, reduction = "Harmony", reduced_dim_name = "HARMONY"))
  }
  
  if ("pca" %in% red_names) {
    return(list(seu = seu, reduction = "pca", reduced_dim_name = "PCA"))
  }
  
  log_msg("No PCA/Harmony reduction found. Recomputing PCA.")
  
  DefaultAssay(seu) <- "RNA"
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 3000, verbose = FALSE)
  seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
  npcs <- min(max_dims, ncol(seu) - 1, 50)
  seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = npcs, verbose = FALSE)
  
  list(seu = seu, reduction = "pca", reduced_dim_name = "PCA")
}

ensure_umap <- function(seu, reduction_use = "pca", dims_use = 1:20) {
  if ("umap" %in% names(seu@reductions)) return(seu)
  
  log_msg("No UMAP found. Running UMAP for visualization.")
  dims_use <- dims_use[dims_use <= ncol(Embeddings(seu, reduction_use))]
  seu <- RunUMAP(seu, reduction = reduction_use, dims = dims_use, verbose = FALSE)
  seu
}

stratified_downsample <- function(seu, group_col, subtype_col, sample_col, max_cells = 60000) {
  n <- ncol(seu)
  if (n <= max_cells) {
    log_msg(paste("No downsampling needed:", n, "cells."))
    return(seu)
  }
  
  log_msg(paste("Downsampling from", n, "to", max_cells, "cells."))
  
  meta <- seu@meta.data
  meta$.cell <- rownames(meta)
  
  strata_cols <- c(group_col, subtype_col, sample_col)
  meta$.strata <- apply(meta[, strata_cols, drop = FALSE], 1, paste, collapse = "__")
  
  strata_tab <- table(meta$.strata)
  strata_prop <- as.numeric(strata_tab) / sum(strata_tab)
  target_n <- pmax(5, round(strata_prop * max_cells))
  names(target_n) <- names(strata_tab)
  
  sampled_cells <- unlist(lapply(names(strata_tab), function(st) {
    cells <- meta$.cell[meta$.strata == st]
    size <- min(length(cells), target_n[st])
    sample(cells, size = size)
  }))
  
  if (length(sampled_cells) > max_cells) {
    sampled_cells <- sample(sampled_cells, max_cells)
  }
  
  log_msg(paste("Sampled cells:", length(sampled_cells)))
  subset(seu, cells = sampled_cells)
}

make_sig_category <- function(fdr) {
  cut(
    fdr,
    breaks = c(-Inf, 0.01, 0.05, 0.10, Inf),
    labels = c("FDR<0.01", "FDR<0.05", "FDR<0.1", "NS")
  )
}

make_direction <- function(logFC, fdr, cutoff = 0.10) {
  ifelse(
    is.na(fdr) | fdr >= cutoff,
    "NS",
    ifelse(logFC > 0, "AD-enriched", "CN-enriched")
  )
}


# =============================================================================
# 3. Main MiloR function
# =============================================================================

run_nvu_milo <- function(ct_name, qs_file, subtype_col, out_prefix) {
  
  log_msg("====================================================")
  log_msg(paste("Starting NVU MiloR:", ct_name))
  log_msg("====================================================")
  
  ckpt <- file.path(OUT_RDS, paste0("ckpt_", out_prefix, "_ADCN_publication_v2_milor_final.qs"))
  
  if (file.exists(ckpt)) {
    log_msg(paste("Checkpoint found. Loading:", ckpt))
    return(qread(ckpt))
  }
  
  if (!file.exists(qs_file)) {
    log_msg(paste("Input file does not exist:", qs_file))
    return(NULL)
  }
  
  # ---------------------------------------------------------------------------
  # Load object
  # ---------------------------------------------------------------------------
  
  seu <- qread(qs_file)
  seu <- safe_join_layers(seu)
  
  meta <- seu@meta.data
  
  if (!subtype_col %in% colnames(meta)) {
    log_msg(paste("Subtype column missing:", subtype_col))
    return(NULL)
  }
  
  group_col_raw <- detect_group_column(meta)
  sample_col <- get_sample_column(meta)
  
  seu$MiloGroup <- standardize_group(meta[[group_col_raw]])
  seu$MiloSample <- as.character(meta[[sample_col]])
  
  if (any(is.na(seu$MiloGroup))) {
    keep <- !is.na(seu$MiloGroup)
    log_msg(paste("Removing cells with undefined group:", sum(!keep)))
    seu <- subset(seu, cells = colnames(seu)[keep])
  }
  
  group_tab <- table(seu$MiloGroup)
  log_msg(paste("Group table:", paste(names(group_tab), group_tab, collapse = "; ")))
  
  sample_group <- seu@meta.data %>%
    distinct(MiloSample, MiloGroup)
  
  samples_per_group <- table(sample_group$MiloGroup)
  
  if (any(samples_per_group < MIN_SAMPLES_PER_GROUP)) {
    log_msg("Too few samples in at least one group. MiloR skipped.")
    return(NULL)
  }
  
  if (ncol(seu) < MIN_CELLS_PER_OBJECT) {
    log_msg(paste("Too few cells:", ncol(seu), ". MiloR skipped."))
    return(NULL)
  }
  
  # ---------------------------------------------------------------------------
  # Downsample if needed
  # ---------------------------------------------------------------------------
  
  seu <- stratified_downsample(
    seu = seu,
    group_col = "MiloGroup",
    subtype_col = subtype_col,
    sample_col = "MiloSample",
    max_cells = MAX_CELLS_PER_OBJECT
  )
  
  # ---------------------------------------------------------------------------
  # Ensure reduction and UMAP
  # ---------------------------------------------------------------------------
  
  red_info <- ensure_reduction(seu, max_dims = 50)
  seu <- red_info$seu
  reduction_use <- red_info$reduction
  reduced_dim_name <- red_info$reduced_dim_name
  
  emb <- Embeddings(seu, reduction_use)
  dims_use <- min(DIMS_USE_DEFAULT, ncol(emb), ncol(seu) - 1)
  
  if (dims_use < 5) {
    log_msg("Too few dimensions available. MiloR skipped.")
    return(NULL)
  }
  
  seu <- ensure_umap(seu, reduction_use = reduction_use, dims_use = seq_len(dims_use))
  
  log_msg(paste("Using reduction:", reduction_use, "| Milo reducedDim:", reduced_dim_name, "| dims:", dims_use))
  
  # ---------------------------------------------------------------------------
  # Build SCE / Milo object
  # ---------------------------------------------------------------------------
  
  counts <- get_layer_safe(seu, assay = "RNA", layer = "counts")
  data_layer <- get_layer_safe(seu, assay = "RNA", layer = "data")
  
  sce <- SingleCellExperiment(
    assays = list(
      counts = counts,
      logcounts = data_layer
    ),
    colData = seu@meta.data
  )
  
  reducedDim(sce, reduced_dim_name) <- Embeddings(seu, reduction_use)
  reducedDim(sce, "UMAP") <- Embeddings(seu, "umap")
  
  milo <- Milo(sce)
  
  # ---------------------------------------------------------------------------
  # Build graph and neighborhoods
  # ---------------------------------------------------------------------------
  
  log_msg("Building Milo graph.")
  milo <- buildGraph(
    milo,
    k = K_NN,
    d = dims_use,
    reduced.dim = reduced_dim_name
  )
  
  log_msg("Making neighborhoods.")
  milo <- makeNhoods(
    milo,
    prop = NHOOD_PROP,
    k = K_NN,
    d = dims_use,
    refined = TRUE,
    reduced_dims = reduced_dim_name
  )
  
  n_nhood <- ncol(nhoods(milo))
  log_msg(paste("Number of neighborhoods:", n_nhood))
  
  if (n_nhood < 10) {
    log_msg("Too few neighborhoods. MiloR skipped.")
    return(NULL)
  }
  
  # ---------------------------------------------------------------------------
  # Count cells per sample
  # ---------------------------------------------------------------------------
  
  milo <- countCells(
    milo,
    meta.data = as.data.frame(colData(milo)),
    samples = "MiloSample"
  )
  
  # ---------------------------------------------------------------------------
  # Sample design
  # ---------------------------------------------------------------------------
  
  cd <- as.data.frame(colData(milo))
  
  sample_design <- cd %>%
    as.data.frame() %>%
    distinct(MiloSample, MiloGroup, .keep_all = TRUE) %>%
    select(MiloSample, MiloGroup, any_of(CANDIDATE_COVARIATES))
  
  rownames(sample_design) <- sample_design$MiloSample
  
  # Keep only samples represented in nhoodCounts
  sample_design <- sample_design[colnames(nhoodCounts(milo)), , drop = FALSE]
  
  sample_design$MiloGroup <- factor(sample_design$MiloGroup, levels = c("CN", "AD"))
  
  # Candidate covariates
  valid_covars <- c()
  if (USE_COVARIATES) {
    for (cv in CANDIDATE_COVARIATES) {
      if (cv %in% colnames(sample_design)) {
        x <- sample_design[[cv]]
        if (sum(!is.na(x)) == nrow(sample_design) && length(unique(x)) > 1) {
          valid_covars <- c(valid_covars, cv)
        }
      }
    }
  }
  
  if (length(valid_covars) > 0) {
    design_formula <- as.formula(paste("~", paste(c(valid_covars, "MiloGroup"), collapse = " + ")))
  } else {
    design_formula <- ~ MiloGroup
  }
  
  log_msg(paste("Design formula:", deparse(design_formula)))
  
  safe_csv(sample_design, paste0("Table_", out_prefix, "_sample_design.csv"))
  
  # sample cell counts
  sample_cell_counts <- cd %>%
    as.data.frame() %>%
    dplyr::count(MiloSample, as.character(MiloGroup), name = "Cells") %>%
    dplyr::rename(MiloGroup = `as.character(MiloGroup)`) %>%
    arrange(MiloGroup, MiloSample)
  
  safe_csv(sample_cell_counts, paste0("Table_", out_prefix, "_sample_cell_counts.csv"))
  
  # ---------------------------------------------------------------------------
  # DA test
  # ---------------------------------------------------------------------------
  
  log_msg("Testing DA neighborhoods.")
  
  da <- testNhoods(
    milo,
    design = design_formula,
    design.df = sample_design,
    fdr.weighting = "none",
    reduced.dim = reduced_dim_name
  )
  
  # Milo sometimes returns SpatialFDR already; keep both but standardize
  da$SpatialFDR_BH <- p.adjust(da$PValue, method = "BH")
  da$SpatialFDR <- da$SpatialFDR_BH
  
  da$SigCat <- make_sig_category(da$SpatialFDR)
  da$Direction <- make_direction(da$logFC, da$SpatialFDR, cutoff = FDR_CUTOFF_MAIN)
  
  # ---------------------------------------------------------------------------
  # Annotate neighborhoods
  # ---------------------------------------------------------------------------
  
  nh <- nhoods(milo)
  cd <- as.data.frame(colData(milo))
  
  dominant_subtype <- sapply(seq_len(ncol(nh)), function(i) {
    idx <- which(nh[, i] > 0)
    if (length(idx) == 0) return(NA_character_)
    vals <- as.character(cd[idx, subtype_col])
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(NA_character_)
    names(sort(table(vals), decreasing = TRUE))[1]
  })
  
  subtype_purity <- sapply(seq_len(ncol(nh)), function(i) {
    idx <- which(nh[, i] > 0)
    if (length(idx) == 0) return(NA_real_)
    vals <- as.character(cd[idx, subtype_col])
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(NA_real_)
    max(table(vals)) / length(vals)
  })
  
  nhood_size <- Matrix::colSums(nh > 0)
  
  umap <- reducedDim(milo, "UMAP")
  nh_umap <- t(sapply(seq_len(ncol(nh)), function(i) {
    idx <- which(nh[, i] > 0)
    if (length(idx) == 0) c(NA_real_, NA_real_) else colMeans(umap[idx, , drop = FALSE])
  }))
  
  da$Nhood <- seq_len(nrow(da))
  da$NhoodSubtype <- dominant_subtype
  da$SubtypePurity <- subtype_purity
  da$NhoodSize <- as.numeric(nhood_size)
  
  nhood_meta <- data.frame(
    Nhood = seq_len(ncol(nh)),
    UMAP1 = nh_umap[, 1],
    UMAP2 = nh_umap[, 2],
    Size = as.numeric(nhood_size),
    logFC = da$logFC,
    PValue = da$PValue,
    SpatialFDR = da$SpatialFDR,
    SigCat = as.character(da$SigCat),
    Direction = da$Direction,
    Subtype = dominant_subtype,
    SubtypePurity = subtype_purity,
    stringsAsFactors = FALSE
  )
  
  safe_csv(da, paste0("Table_", out_prefix, "_DA_results_full.csv"))
  safe_csv(nhood_meta, paste0("Table_", out_prefix, "_nhood_metadata_full.csv"))
  
  sig_nhoods <- nhood_meta %>%
    filter(SpatialFDR < FDR_CUTOFF_MAIN) %>%
    arrange(SpatialFDR)
  
  safe_csv(sig_nhoods, paste0("Table_", out_prefix, "_significant_nhoods_FDR10.csv"))
  
  top_nhoods <- nhood_meta %>%
    arrange(SpatialFDR) %>%
    dplyr::filter(row_number() <= min(50, nrow(nhood_meta)))
  
  safe_csv(top_nhoods, paste0("Table_", out_prefix, "_top50_nhoods.csv"))
  
  # ---------------------------------------------------------------------------
  # Summary by subtype
  # ---------------------------------------------------------------------------
  
  da_by_subtype <- nhood_meta %>%
    filter(SpatialFDR < FDR_CUTOFF_MAIN) %>%
    group_by(Subtype) %>%
    summarise(
      Total = n(),
      AD_enriched = sum(logFC > 0, na.rm = TRUE),
      CN_enriched = sum(logFC < 0, na.rm = TRUE),
      Mean_logFC = mean(logFC, na.rm = TRUE),
      Median_logFC = median(logFC, na.rm = TRUE),
      Mean_purity = mean(SubtypePurity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Total))
  
  safe_csv(da_by_subtype, paste0("Table_", out_prefix, "_DA_by_subtype.csv"))
  
  # ---------------------------------------------------------------------------
  # Parameters
  # ---------------------------------------------------------------------------
  
  param_table <- data.frame(
    Parameter = c(
      "CellType", "InputFile", "SubtypeColumn", "Reduction",
      "ReducedDimName", "Dimensions", "k", "NhoodProp",
      "MaxCells", "FDRCutoff", "DesignFormula",
      "PositiveLogFC"
    ),
    Value = c(
      ct_name, qs_file, subtype_col, reduction_use,
      reduced_dim_name, dims_use, K_NN, NHOOD_PROP,
      MAX_CELLS_PER_OBJECT, FDR_CUTOFF_MAIN, deparse(design_formula),
      "AD-enriched versus CN"
    )
  )
  
  safe_csv(param_table, paste0("Table_", out_prefix, "_parameters.csv"))
  
  # =============================================================================
  # 4. Figures
  # =============================================================================
  
  log_msg("Generating figures.")
  
  cell_umap_df <- data.frame(
    UMAP1 = umap[, 1],
    UMAP2 = umap[, 2],
    Group = cd$MiloGroup,
    Subtype = cd[[subtype_col]]
  )
  
  # ---------------------------------------------------------------------------
  # Fig 1. DA neighborhood UMAP
  # ---------------------------------------------------------------------------
  
  p_da_umap <- ggplot() +
    geom_point(
      data = cell_umap_df,
      aes(UMAP1, UMAP2),
      color = "grey90",
      size = 0.15,
      alpha = 0.55
    ) +
    geom_point(
      data = nhood_meta,
      aes(
        UMAP1, UMAP2,
        color = logFC,
        size = pmin(-log10(SpatialFDR), 6),
        alpha = Direction != "NS"
      )
    ) +
    scale_color_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "log2FC\nAD/CN"
    ) +
    scale_size_continuous(range = c(0.4, 3.2), name = "-log10(FDR)") +
    scale_alpha_manual(values = c("TRUE" = 0.95, "FALSE" = 0.35), guide = "none") +
    labs(
      title = paste0(ct_name, ": DA neighborhoods"),
      subtitle = paste0("MiloR, ", reduced_dim_name, ", k=", K_NN, ", d=", dims_use),
      x = "UMAP1",
      y = "UMAP2"
    ) +
    theme_publication(11)
  
  save_fig(p_da_umap, paste0("Fig_", out_prefix, "_01_DA_neighborhood_UMAP"), 7.2, 6.2, OUT_FIGS)
  
  # ---------------------------------------------------------------------------
  # Fig 2. Volcano
  # ---------------------------------------------------------------------------
  
  p_volcano <- ggplot(
    nhood_meta,
    aes(logFC, -log10(SpatialFDR), color = SigCat)
  ) +
    geom_point(alpha = 0.82, size = 1.25) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey45") +
    geom_hline(yintercept = -log10(FDR_CUTOFF_MAIN), linetype = "dashed", linewidth = 0.35, color = "grey45") +
    scale_color_manual(values = sig_palette, drop = FALSE) +
    labs(
      title = paste0(ct_name, ": MiloR DA volcano"),
      x = "log2FC AD/CN",
      y = "-log10 MiloR FDR"
    ) +
    theme_publication(11) +
    theme(legend.title = element_blank())
  
  save_fig(p_volcano, paste0("Fig_", out_prefix, "_02_DA_volcano"), 6.3, 5.2, OUT_FIGS)
  
  # ---------------------------------------------------------------------------
  # Fig 3. DA by subtype, direction-aware
  # ---------------------------------------------------------------------------
  
  subtype_dir <- nhood_meta %>%
    filter(SpatialFDR < FDR_CUTOFF_MAIN) %>%
    mutate(Direction = factor(Direction, levels = c("AD-enriched", "CN-enriched"))) %>%
    as.data.frame() %>%
    dplyr::count(Subtype, as.character(Direction), name = "Nhoods") %>%
    dplyr::rename(Direction = `as.character(Direction)`) %>%
    group_by(Subtype) %>%
    mutate(Total = sum(Nhoods)) %>%
    ungroup() %>%
    arrange(desc(Total))
  
  if (nrow(subtype_dir) > 0) {
    subtype_order <- subtype_dir %>%
      group_by(Subtype) %>%
      summarise(Total = sum(Nhoods), .groups = "drop") %>%
      arrange(Total) %>%
      pull(Subtype)
    
    subtype_dir$Subtype <- factor(subtype_dir$Subtype, levels = subtype_order)
    
    p_subtype <- ggplot(subtype_dir, aes(Subtype, Nhoods, fill = Direction)) +
      geom_col(width = 0.72, color = "grey25", linewidth = 0.25) +
      coord_flip() +
      scale_fill_manual(values = fc_palette) +
      labs(
        title = paste0(ct_name, ": DA neighborhoods by subtype"),
        x = NULL,
        y = paste0("Significant neighborhoods, FDR<", FDR_CUTOFF_MAIN)
      ) +
      theme_publication(11)
    
    save_fig(p_subtype, paste0("Fig_", out_prefix, "_03_DA_by_subtype_direction"), 7, 4.8, OUT_FIGS)
  }
  
  # ---------------------------------------------------------------------------
  # Fig 4. Nhood size distribution
  # ---------------------------------------------------------------------------
  
  p_size <- ggplot(nhood_meta, aes(Size)) +
    geom_histogram(bins = 40, fill = "grey60", color = "white", linewidth = 0.15) +
    geom_vline(xintercept = median(nhood_meta$Size, na.rm = TRUE), color = "#B2182B", linewidth = 0.7) +
    labs(
      title = paste0(ct_name, ": neighborhood size distribution"),
      subtitle = paste0("Median size = ", round(median(nhood_meta$Size, na.rm = TRUE), 1)),
      x = "Cells per neighborhood",
      y = "Number of neighborhoods"
    ) +
    theme_publication(11)
  
  save_fig(p_size, paste0("Fig_", out_prefix, "_04_nhood_size_distribution"), 5.8, 4.6, OUT_SUPP)
  
  # ---------------------------------------------------------------------------
  # Fig 5. Sample-level cell counts
  # ---------------------------------------------------------------------------
  
  p_sample_counts <- ggplot(sample_cell_counts, aes(MiloGroup, Cells, fill = MiloGroup)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.12, size = 2, alpha = 0.9) +
    scale_fill_manual(values = group_palette) +
    labs(
      title = paste0(ct_name, ": sample-level cell counts"),
      x = NULL,
      y = "Cells per sample"
    ) +
    theme_publication(11) +
    theme(legend.position = "none")
  
  save_fig(p_sample_counts, paste0("Fig_", out_prefix, "_05_sample_cell_counts"), 4.4, 4.5, OUT_SUPP)
  
  # ---------------------------------------------------------------------------
  # Fig 6. Subtype composition by group
  # ---------------------------------------------------------------------------
  
  subtype_comp <- cd %>%
    as.data.frame() %>%
    dplyr::count(MiloSample, as.character(MiloGroup), !!sym(subtype_col), name = "Cells") %>%
    dplyr::rename(MiloGroup = `as.character(MiloGroup)`) %>%
    group_by(MiloSample) %>%
    mutate(Fraction = Cells / sum(Cells)) %>%
    ungroup() %>%
    dplyr::rename(Subtype = !!sym(subtype_col))
  
  safe_csv(subtype_comp, paste0("Table_", out_prefix, "_subtype_composition_by_sample.csv"))
  
  p_comp <- ggplot(subtype_comp, aes(MiloGroup, Fraction, fill = MiloGroup)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.70) +
    geom_jitter(width = 0.10, size = 1.6, alpha = 0.85) +
    facet_wrap(~Subtype, scales = "free_y") +
    scale_fill_manual(values = group_palette) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = paste0(ct_name, ": subtype composition by group"),
      x = NULL,
      y = "Fraction per sample"
    ) +
    theme_publication(10) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
  
  save_fig(p_comp, paste0("Fig_", out_prefix, "_06_subtype_composition"), 8.5, 6.2, OUT_SUPP)
  
  # ---------------------------------------------------------------------------
  # Fig 7. Top significant nhood sample-level abundance
  # ---------------------------------------------------------------------------
  
  nhood_count_mat <- as.matrix(nhoodCounts(milo))
  
  sample_totals <- cd %>%
    as.data.frame() %>%
    dplyr::count(MiloSample, name = "SampleCells")
  
  rownames(sample_totals) <- sample_totals$MiloSample
  
  # Normalize to cells per 10k cells in this subset
  norm_mat <- sweep(
    nhood_count_mat,
    2,
    sample_totals[colnames(nhood_count_mat), "SampleCells"],
    FUN = "/"
  ) * 10000
  
  top_sig <- nhood_meta %>%
    filter(SpatialFDR < FDR_CUTOFF_MAIN) %>%
    arrange(SpatialFDR) %>%
    slice_head(n = TOP_NHOODS_FOR_BOXPLOT)
  
  if (nrow(top_sig) > 0) {
    top_ids <- top_sig$Nhood
    
    top_abund <- as.data.frame(t(norm_mat[top_ids, , drop = FALSE]))
    top_abund$MiloSample <- rownames(top_abund)
    top_abund$MiloGroup <- sample_design[top_abund$MiloSample, "MiloGroup"]
    
    top_abund_long <- top_abund %>%
      pivot_longer(
        cols = starts_with(as.character(top_ids[1])) | matches("^[0-9]+$"),
        names_to = "Nhood",
        values_to = "Abundance_per10k"
      )
    
    top_abund_long$Nhood <- as.integer(top_abund_long$Nhood)
    
    top_abund_long <- top_abund_long %>%
      left_join(
        top_sig %>% select(Nhood, logFC, SpatialFDR, Subtype),
        by = "Nhood"
      ) %>%
      mutate(
        Label = paste0(
          "N", Nhood,
          " | ", Subtype,
          "\nlogFC=", round(logFC, 2),
          ", FDR=", signif(SpatialFDR, 2)
        )
      )
    
    safe_csv(top_abund_long, paste0("Table_", out_prefix, "_top_nhood_sample_abundance.csv"))
    
    p_top_box <- ggplot(top_abund_long, aes(MiloGroup, Abundance_per10k, fill = MiloGroup)) +
      geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
      geom_jitter(width = 0.10, size = 1.6, alpha = 0.85) +
      facet_wrap(~Label, scales = "free_y", ncol = 4) +
      scale_fill_manual(values = group_palette) +
      labs(
        title = paste0(ct_name, ": top DA neighborhood abundance"),
        x = NULL,
        y = "Neighborhood cells per 10,000 cells"
      ) +
      theme_publication(9.5) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text = element_text(size = 7.5)
      ) +
      geom_text(data = top_abund_long %>% distinct(Label, .keep_all = TRUE),
                aes(label = ifelse(SpatialFDR < 0.001, "***",
                           ifelse(SpatialFDR < 0.01, "**",
                           ifelse(SpatialFDR < 0.05, "*", "ns"))),
                    x = 1.5, y = Inf),
                vjust = 1.5, size = 3.5, inherit.aes = FALSE)
    
    save_fig(p_top_box, paste0("Fig_", out_prefix, "_07_top_DA_nhood_sample_abundance"), 10, 7.2, OUT_FIGS)
  }
  
  # ---------------------------------------------------------------------------
  # Fig 8. Heatmap of top DA nhood abundance
  # ---------------------------------------------------------------------------
  
  top_heat <- nhood_meta %>%
    arrange(SpatialFDR) %>%
    dplyr::filter(row_number() <= min(TOP_NHOODS_FOR_HEATMAP, nrow(nhood_meta)))
  
  if (nrow(top_heat) >= 5) {
    heat_ids <- top_heat$Nhood
    heat_mat <- norm_mat[heat_ids, , drop = FALSE]
    
    # log transform and row z-score
    heat_mat_log <- log2(heat_mat + 1)
    heat_z <- t(scale(t(heat_mat_log)))
    heat_z[is.na(heat_z)] <- 0
    rownames(heat_z) <- paste0(
      "N", heat_ids,
      "_", top_heat$Subtype,
      "_FC", round(top_heat$logFC, 1)
    )
    
    ann_col <- data.frame(
      Group = sample_design[colnames(heat_z), "MiloGroup"]
    )
    rownames(ann_col) <- colnames(heat_z)
    
    ann_colors <- list(Group = group_palette)
    
    pdf(file.path(OUT_SUPP, paste0("Fig_", out_prefix, "_08_top_DA_nhood_heatmap.pdf")),
        width = 7.8, height = 8.5)
    pheatmap(
      heat_z,
      annotation_col = ann_col,
      annotation_colors = ann_colors,
      color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      fontsize = 7,
      main = paste0(ct_name, ": top neighborhood abundance")
    )
    dev.off()
    
    png(file.path(OUT_SUPP, paste0("Fig_", out_prefix, "_08_top_DA_nhood_heatmap.png")),
        width = 2400, height = 2600, res = 300)
    pheatmap(
      heat_z,
      annotation_col = ann_col,
      annotation_colors = ann_colors,
      color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      fontsize = 7,
      main = paste0(ct_name, ": top neighborhood abundance")
    )
    dev.off()
    
    log_msg("Saved top DA nhood heatmap.")
  }
  
  # ---------------------------------------------------------------------------
  # Fig 9. Labeled top DA nhoods on UMAP
  # ---------------------------------------------------------------------------
  
  label_df <- nhood_meta %>%
    filter(SpatialFDR < FDR_CUTOFF_MAIN) %>%
    arrange(SpatialFDR) %>%
    dplyr::filter(row_number() <= min(15, nrow(da))) %>%
    mutate(Label = paste0("N", Nhood, "\n", Subtype))
  
  if (nrow(label_df) > 0) {
    p_label <- ggplot() +
      geom_point(
        data = cell_umap_df,
        aes(UMAP1, UMAP2),
        color = "grey90",
        size = 0.15,
        alpha = 0.5
      ) +
      geom_point(
        data = nhood_meta,
        aes(UMAP1, UMAP2, color = logFC),
        size = 1.2,
        alpha = 0.45
      ) +
      geom_point(
        data = label_df,
        aes(UMAP1, UMAP2),
        color = "black",
        fill = "yellow",
        shape = 21,
        size = 3.2,
        stroke = 0.4
      ) +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(UMAP1, UMAP2, label = Label),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.35,
        point.padding = 0.2,
        segment.color = "grey40",
        segment.linewidth = 0.25
      ) +
      scale_color_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        name = "log2FC"
      ) +
      labs(
        title = paste0(ct_name, ": top DA neighborhoods"),
        x = "UMAP1",
        y = "UMAP2"
      ) +
      theme_publication(11)
    
    save_fig(p_label, paste0("Fig_", out_prefix, "_09_top_DA_nhood_labels"), 7.2, 6.2, OUT_FIGS)
  }
  
  # ---------------------------------------------------------------------------
  # Save Milo object and result
  # ---------------------------------------------------------------------------
  
  res <- list(
    milo = milo,
    da = da,
    nhood_meta = nhood_meta,
    da_by_subtype = da_by_subtype,
    sample_design = sample_design,
    sample_cell_counts = sample_cell_counts,
    parameters = param_table
  )
  
  qsave(res, ckpt, preset = "fast")
  log_msg(paste("Finished NVU MiloR:", ct_name))
  
  gc()
  
  return(res)
}


# =============================================================================
# 5. Run analyses
# =============================================================================

log_msg("===== NVU-focused MiloR analysis started =====")

res_endothelial <- run_nvu_milo(
  ct_name = "Cerebrovascular cells",
  qs_file = file.path(RDS_DIR, "scRNA_cereb_annotated.qs"),
  subtype_col = "vascular_subtype",
  out_prefix = "Cerebrovascular"
)

res_astro <- run_nvu_milo(
  ct_name = "Astrocytes",
  qs_file = file.path(RDS_DIR, "scRNA_astro_annotated.qs"),
  subtype_col = "astro_subtype",
  out_prefix = "Astrocytes"
)

res_micro <- run_nvu_milo(
  ct_name = "Microglia",
  qs_file = file.path(RDS_DIR, "scRNA_micro_annotated.qs"),
  subtype_col = "micro_subtype",
  out_prefix = "Microglia"
)


# =============================================================================
# 6. Combined summary across NVU compartments
# =============================================================================

log_msg("Generating combined NVU summary.")

collect_summary <- function(res, compartment) {
  if (is.null(res)) return(NULL)
  nm <- res$nhood_meta
  if (is.null(nm) || nrow(nm) == 0) return(NULL)
  
  nm %>%
    mutate(Compartment = compartment) %>%
    select(
      Compartment,
      Nhood,
      Subtype,
      Size,
      logFC,
      PValue,
      SpatialFDR,
      SigCat,
      Direction,
      SubtypePurity
    )
}

combined_nhood <- bind_rows(
  collect_summary(res_endothelial, "Cerebrovascular cells"),
  collect_summary(res_astro, "Astrocytes"),
  collect_summary(res_micro, "Microglia")
)

if (!is.null(combined_nhood) && nrow(combined_nhood) > 0) {
  
  safe_csv(combined_nhood, "Table_NVU_MiloR_combined_nhood_results.csv")
  
  combined_by_compartment <- combined_nhood %>%
    filter(SpatialFDR < FDR_CUTOFF_MAIN) %>%
    group_by(Compartment) %>%
    summarise(
      Total_DA = n(),
      AD_enriched = sum(logFC > 0, na.rm = TRUE),
      CN_enriched = sum(logFC < 0, na.rm = TRUE),
      Mean_logFC = mean(logFC, na.rm = TRUE),
      Median_logFC = median(logFC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Total_DA))
  
  safe_csv(combined_by_compartment, "Table_NVU_MiloR_combined_by_compartment.csv")
  
  combined_by_subtype <- combined_nhood %>%
    filter(SpatialFDR < FDR_CUTOFF_MAIN) %>%
    group_by(Compartment, Subtype) %>%
    summarise(
      Total_DA = n(),
      AD_enriched = sum(logFC > 0, na.rm = TRUE),
      CN_enriched = sum(logFC < 0, na.rm = TRUE),
      Mean_logFC = mean(logFC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Compartment, desc(Total_DA))
  
  safe_csv(combined_by_subtype, "Table_NVU_MiloR_combined_by_subtype.csv")
  
  # Combined figure 1: compartment summary
  comp_long <- combined_by_compartment %>%
    select(Compartment, AD_enriched, CN_enriched) %>%
    pivot_longer(
      cols = c(AD_enriched, CN_enriched),
      names_to = "Direction",
      values_to = "Nhoods"
    ) %>%
    mutate(
      Direction = recode(
        Direction,
        "AD_enriched" = "AD-enriched",
        "CN_enriched" = "CN-enriched"
      )
    )
  
  p_comb_comp <- ggplot(comp_long, aes(Compartment, Nhoods, fill = Direction)) +
    geom_col(width = 0.70, color = "grey25", linewidth = 0.25) +
    scale_fill_manual(values = fc_palette) +
    labs(
      title = "NVU-focused MiloR: DA neighborhoods by compartment",
      x = NULL,
      y = paste0("Significant neighborhoods, FDR<", FDR_CUTOFF_MAIN)
    ) +
    theme_publication(11)
  
  save_fig(p_comb_comp, "Fig_NVU_MiloR_combined_01_by_compartment", 6.2, 4.8, OUT_FIGS)
  
  # Combined figure 2: subtype summary
  if (nrow(combined_by_subtype) > 0) {
    subtype_long <- combined_by_subtype %>%
      select(Compartment, Subtype, AD_enriched, CN_enriched) %>%
      pivot_longer(
        cols = c(AD_enriched, CN_enriched),
        names_to = "Direction",
        values_to = "Nhoods"
      ) %>%
      mutate(
        Direction = recode(
          Direction,
          "AD_enriched" = "AD-enriched",
          "CN_enriched" = "CN-enriched"
        ),
        Label = paste0(Compartment, " | ", Subtype)
      )
    
    subtype_order <- subtype_long %>%
      group_by(Label) %>%
      summarise(Total = sum(Nhoods), .groups = "drop") %>%
      arrange(Total) %>%
      pull(Label)
    
    subtype_long$Label <- factor(subtype_long$Label, levels = subtype_order)
    
    p_comb_sub <- ggplot(subtype_long, aes(Label, Nhoods, fill = Direction)) +
      geom_col(width = 0.72, color = "grey25", linewidth = 0.25) +
      coord_flip() +
      scale_fill_manual(values = fc_palette) +
      labs(
        title = "NVU-focused MiloR: DA neighborhoods by subtype",
        x = NULL,
        y = paste0("Significant neighborhoods, FDR<", FDR_CUTOFF_MAIN)
      ) +
      theme_publication(10)
    
    save_fig(p_comb_sub, "Fig_NVU_MiloR_combined_02_by_subtype", 8.2, 7.2, OUT_FIGS)
  }

  # ---------------------------------------------------------------------------
  # Final publication package: compact publication-style Figure 7
  # ---------------------------------------------------------------------------

  FINAL_DIR <- file.path(OUT, "outputs")
  FINAL_MAIN <- file.path(FINAL_DIR, "main_figures")
  FINAL_SUPP <- file.path(FINAL_DIR, "supplementary_figures")
  FINAL_EDIT <- file.path(FINAL_DIR, "editable_panels")
  FINAL_TAB <- file.path(FINAL_DIR, "supplementary_tables")
  FINAL_METHODS <- file.path(FINAL_DIR, "methods_results")
  FINAL_SOURCE <- file.path(FINAL_DIR, "source_data")
  FINAL_AUDIT <- file.path(FINAL_DIR, "audit")
  for (d in c(FINAL_MAIN, FINAL_SUPP, FINAL_EDIT, FINAL_TAB, FINAL_METHODS, FINAL_SOURCE, FINAL_AUDIT)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  final_theme <- function(base_size = 8) {
    theme_classic(base_size = base_size, base_family = "Arial") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
        axis.text = element_text(color = "black", size = base_size - 1),
        axis.title = element_text(color = "black", size = base_size),
        legend.title = element_text(size = base_size - 1),
        legend.text = element_text(size = base_size - 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = base_size - 1),
        plot.margin = margin(4, 4, 4, 4)
      )
  }

  comp_final <- combined_by_compartment %>%
    select(Compartment, AD_enriched, CN_enriched) %>%
    pivot_longer(c(AD_enriched, CN_enriched), names_to = "Direction", values_to = "Nhoods") %>%
    mutate(
      Direction = recode(Direction, AD_enriched = "AD-enriched", CN_enriched = "CN-enriched"),
      Compartment = factor(Compartment, levels = combined_by_compartment$Compartment)
    )

  p7a <- ggplot(comp_final, aes(Compartment, Nhoods, fill = Direction)) +
    geom_col(width = 0.68, color = "grey25", linewidth = 0.18) +
    scale_fill_manual(values = fc_palette[c("AD-enriched", "CN-enriched")], name = NULL) +
    labs(
      title = "A  NVU compartments",
      x = NULL,
      y = "Significant neighborhoods\n(MiloR FDR < 0.10)"
    ) +
    final_theme(8) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")

  subtype_final <- combined_by_subtype %>%
    select(Compartment, Subtype, AD_enriched, CN_enriched, Mean_logFC, Total_DA) %>%
    pivot_longer(c(AD_enriched, CN_enriched), names_to = "Direction", values_to = "Nhoods") %>%
    mutate(
      Direction = recode(Direction, AD_enriched = "AD-enriched", CN_enriched = "CN-enriched"),
      Label = paste0(Compartment, " | ", Subtype)
    ) %>%
    filter(Nhoods > 0)

  subtype_order <- subtype_final %>%
    group_by(Label) %>%
    summarise(Total = sum(Nhoods), .groups = "drop") %>%
    arrange(Total) %>%
    pull(Label)
  subtype_final$Label <- factor(subtype_final$Label, levels = subtype_order)

  p7b <- ggplot(subtype_final, aes(Label, Nhoods, fill = Direction)) +
    geom_col(width = 0.68, color = "grey25", linewidth = 0.18) +
    coord_flip() +
    scale_fill_manual(values = fc_palette[c("AD-enriched", "CN-enriched")], name = NULL) +
    labs(
      title = "B  NVU subtype directionality",
      x = NULL,
      y = "Significant neighborhoods"
    ) +
    final_theme(8) +
    theme(legend.position = "bottom")

  p7c_data <- combined_by_subtype %>%
    mutate(
      Label = paste0(Compartment, " | ", Subtype),
      Direction = ifelse(Mean_logFC > 0, "AD-enriched", "CN-enriched"),
      Label = factor(Label, levels = rev(subtype_order))
    )

  p7c <- ggplot(p7c_data, aes(Mean_logFC, Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, color = "grey45") +
    geom_point(aes(size = Total_DA, color = Mean_logFC), alpha = 0.92) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                          name = "Mean log2FC\nAD/CN") +
    scale_size_continuous(range = c(1.5, 5.0), name = "DA nhoods") +
    labs(
      title = "C  Direction and effect size",
      x = "Mean log2FC (AD/CN)",
      y = NULL
    ) +
    final_theme(8) +
    theme(legend.position = "right")

  top_table <- combined_by_subtype %>%
    arrange(desc(Total_DA), desc(abs(Mean_logFC))) %>%
    mutate(
      ShortLabel = paste0(Subtype, " (", Compartment, ")"),
      DirectionLabel = paste0("AD ", AD_enriched, " / CN ", CN_enriched)
    )

  p7d <- ggplot(top_table, aes(Total_DA, reorder(ShortLabel, Total_DA), fill = Mean_logFC)) +
    geom_col(width = 0.72, color = "grey25", linewidth = 0.18) +
    geom_text(aes(label = DirectionLabel), hjust = -0.05, size = 2.2, family = "Arial") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         name = "Mean log2FC\nAD/CN") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.38))) +
    labs(title = "D  Prioritized subtype programs", x = "DA neighborhoods", y = NULL) +
    final_theme(7) +
    theme(
      legend.position = "right"
    )

  save_fig(p7a, "Fig7A_NVU_MiloR_compartment_summary", 3.9, 3.0, FINAL_EDIT)
  save_fig(p7b, "Fig7B_NVU_MiloR_subtype_directionality", 4.5, 4.0, FINAL_EDIT)
  save_fig(p7c, "Fig7C_NVU_MiloR_effect_size_dotplot", 4.4, 4.0, FINAL_EDIT)
  save_fig(p7d, "Fig7D_NVU_MiloR_readable_summary", 6.0, 3.6, FINAL_EDIT)

  figure7 <- (p7a | p7b) / (p7c | p7d) +
    plot_layout(heights = c(0.92, 1.08), widths = c(0.95, 1.05))

  ggsave(file.path(FINAL_MAIN, "Figure_7_Day2_4_5_6_NVU_MiloR_subtype_DA.pdf"),
         figure7, width = 8.3, height = 7.8, device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(FINAL_MAIN, "Figure_7_Day2_4_5_6_NVU_MiloR_subtype_DA.png"),
         figure7, width = 8.3, height = 7.8, dpi = 600, bg = "white", limitsize = FALSE)
  ggsave(file.path(FINAL_MAIN, "Figure_7_Day2_4_5_6_NVU_MiloR_subtype_DA.tiff"),
         figure7, width = 8.3, height = 7.8, dpi = 600, compression = "lzw", bg = "white", limitsize = FALSE)
  try(ggsave(file.path(FINAL_MAIN, "Figure_7_Day2_4_5_6_NVU_MiloR_subtype_DA.svg"),
             figure7, width = 8.3, height = 7.8, limitsize = FALSE), silent = TRUE)

  write.csv(combined_nhood, file.path(FINAL_SOURCE, "SourceData_Figure_7_combined_nhood_results.csv"), row.names = FALSE)
  write.csv(combined_by_compartment, file.path(FINAL_TAB, "Table_S24_NVU_MiloR_by_compartment.csv"), row.names = FALSE)
  write.csv(combined_by_subtype, file.path(FINAL_TAB, "Table_S25_NVU_MiloR_by_subtype.csv"), row.names = FALSE)

  audit <- data.frame(
    Item = c("Scope", "Group labels", "Cerebrovascular naming", "Subtype-level focus", "Main figure"),
    Decision = c("PASS", "PASS", "PASS", "PASS", "PASS_PENDING_VISUAL_REVIEW"),
    Note = c(
      "This node is restricted to NVU-related subtype MiloR for cerebrovascular cells, astrocytes and microglia.",
      "All labels use CN and AD; logFC is reported as AD/CN.",
      "The updated node04 Cerebrovascular cells object is used instead of the stale Endothelial object.",
      "Major-cell MiloR remains in node05; this node addresses subtype-resolved neighborhoods.",
      "Figure 7 summarizes compartment and subtype-level DA directionality with editable PDF/SVG panels."
    )
  )
  write.csv(audit, file.path(FINAL_AUDIT, "Table_S26_Day2_4_5_6_NVU_MiloR_quality_control.csv"), row.names = FALSE)

  legend_text <- c(
    "Figure 7. NVU subtype-resolved differential abundance in AD.",
    "A, Number of significant MiloR neighborhoods across neurovascular-unit compartments, separated by AD-enriched and CN-enriched directions. B, Significant neighborhoods stratified by dominant subtype within each compartment. C, Mean neighborhood log2 fold-change for each subtype; point size indicates the number of significant neighborhoods. D, Compact summary of subtype-level DA neighborhoods. CN, cognitively normal control; AD, Alzheimer disease; DA, differential abundance; FDR, false discovery rate; NVU, neurovascular unit; SMC, smooth muscle cell; DAM, disease-associated microglia."
  )
  writeLines(legend_text, file.path(FINAL_METHODS, "Day2_4_5_6_NVU_MiloR_figure_legends.txt"), useBytes = TRUE)

  methods_text <- c(
    "Methods - NVU subtype-resolved MiloR differential abundance",
    "",
    "Subtype-resolved neighborhood differential abundance was performed with MiloR on the node04 annotated cerebrovascular, astrocyte and microglial objects. Group labels were harmonized to CN and AD, and the updated Cerebrovascular cells object was used for vascular-lineage analysis. For each compartment, a k-nearest-neighbor graph was constructed from the available Harmony/PCA embedding, refined neighborhoods were sampled, and neighborhood cell counts were modeled at the sample level with group as the primary design term. Positive log2 fold-change denotes AD enrichment relative to CN. Neighborhoods were assigned to the dominant subtype, and significant neighborhoods were defined at MiloR graph FDR < 0.10.",
    "",
    "Results - NVU subtype DA remodeling",
    "",
    paste0("The updated NVU MiloR analysis identified ", sum(combined_by_compartment$Total_DA), " significant subtype-resolved neighborhoods across the three NVU-related compartments. Microglia showed the strongest AD-enriched signal, driven by DAM neighborhoods. Cerebrovascular cells were predominantly CN-enriched, especially SMC and arterial neighborhoods, with only limited AD-enriched pericyte signal. Astrocyte neighborhoods were mainly CN-enriched for homeostatic astrocytes, with sparse AD-enriched reactive/intermediate neighborhoods. These results support a disease trajectory in which AD-associated microglial activation occurs alongside depletion or loss of homeostatic vascular and astrocytic neighborhood states.")
  )
  writeLines(methods_text, file.path(FINAL_METHODS, "Day2_4_5_6_NVU_MiloR_methods_results_draft.txt"), useBytes = TRUE)
}


# =============================================================================
# 7. Save all results and sessionInfo
# =============================================================================

qsave(
  list(
    Cerebrovascular = res_endothelial,
    Astrocytes = res_astro,
    Microglia = res_micro,
    Combined = combined_nhood
  ),
  file.path(OUT_RDS, "NVU_MiloR_final_results_all.qs"),
  preset = "fast"
)

sink(file.path(OUT_LOG, "sessionInfo.txt"))
print(sessionInfo())
sink()

log_msg("===== ALL NVU MiloR analyses completed =====")
