################################################################################
# node01_complete_final_v5_publication.R
#
# PURPOSE
#   End-to-end node1 preprocessing pipeline for large multi-sample snRNA-seq/
#   scRNA-seq datasets under R 4.4.1 + Seurat 5.5.0 / SeuratObject 5.4.0.
#
# KEY DESIGN CHOICES
#   1) Full pipeline from raw 10x input OR resume from node1_merged_after_sample_qc.rds
#   2) Sample-wise scDblFinder to avoid Windows native crash on huge merged objects
#   3) Sample-wise decontX retained (required by user)
#   4) Explicit Seurat v5 layer handling (LayerData / JoinLayers / Assay5 creation)
#   5) publication-style figure export defaults and quality-control QC outputs
#
# NOTE
#   This script is written to be robust and general. If you already have
#   node1_merged_after_sample_qc.rds, set PARAMS$RESUME_MODE <- "merged_qc".
#
# EXPECTED INPUT MODES
#   A) RAW MODE
#      - Provide BASE_INPUT_DIR containing 10x folders or .h5 files, OR
#      - Provide SAMPLE_SHEET CSV with at least: SampleID,input_path
#      - Optional columns in SAMPLE_SHEET: Dataset, Group, SubjectID, Sex, Batch
#
#   B) RESUME MODE
#      - Provide MERGED_QC_RDS = path/to/node1_merged_after_sample_qc.rds
#
# AUTHORING TARGET
#   publication / quality-control-ready preprocessing outputs
################################################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 120 * 1024^3)
set.seed(123456)

# Limit threads for Windows stability during Bioconductor native code.
Sys.setenv(
  TMPDIR = "/path/to/tmp",
  TEMP = "/path/to/tmp",
  TMP = "/path/to/tmp",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(4)

# ==============================================================================
# 0. USER PARAMETERS
# ==============================================================================

PARAMS <- list(
  # ------------------------- RUN MODE -----------------------------------------
  RESUME_MODE = "raw",              # "raw" or "merged_qc"
  TEST_MODE = FALSE,

  # ------------------------- INPUTS -------------------------------------------
  SAMPLE_SHEET = "/path/to/project/results/sample_sheet.csv",
  BASE_INPUT_DIR = NA_character_,
  MERGED_QC_RDS = NA_character_,

  # ------------------------- OUTPUT ROOT --------------------------------------
  OUT = "/path/to/project/results/01_atlas",

  # ------------------------- GENERAL ------------------------------------------
  SPECIES = "human",                # "human" or "mouse"
  PROJECT = "AD_brain_snRNA",
  MIN_CELLS_PER_GENE = 3,
  MIN_GENES_PER_CELL = 200,
  SAVE_SAMPLE_CHECKPOINTS = TRUE,
  SAVE_INTERMEDIATE_RDS = TRUE,

  # ------------------------- QC FILTERS ---------------------------------------
  USE_ADAPTIVE_QC = TRUE,
  MAD_NCOUNT_LOW = 3,
  MAD_NCOUNT_HIGH = 4,
  MAD_NFEATURE_LOW = 3,
  MAD_NFEATURE_HIGH = 4,
  DEFAULT_MAX_MT = 5,
  DEFAULT_MAX_RIBO = 20,
  DEFAULT_MAX_HB = 5,
  ABS_MIN_COUNTS = 300,
  ABS_MIN_FEATURES = 200,
  ABS_MAX_FEATURES = 10000,
  ABS_MAX_COUNTS = 100000,

  # ------------------------- DOUBLETS -----------------------------------------
  DOUBLET_RATE = 0.04,
  DOUBLET_DBR_SD = 1.0,
  SCDBL_CLUSTERS = FALSE,
  MIN_ARTIFICIAL_DOUBLETS = 3000,
  MAX_ARTIFICIAL_DOUBLETS = 15000,
  SCDBL_NFEATURES = 2000,

  # ------------------------- AMBIENT RNA --------------------------------------
  RUN_DECONTX = TRUE,
  CONTAM_CUTOFF = 0.25,

  # ------------------------- ANALYSIS -----------------------------------------
  ANALYSIS_ASSAY = "decontX",       # if unavailable automatically falls back to RNA
  N_HVG = 3000,
  N_PCS = 50,
  N_DIMS_USE = 15,
  CLUST_RES = 0.3,
  UMAP_N_NEIGHBORS = 30,
  UMAP_MIN_DIST = 0.30,
  REGRESS_VARS = character(0),       # e.g. c("percent.mt", "decontX_contamination")

  # ------------------------- MARKERS ------------------------------------------
  RUN_MARKERS = TRUE,
  MARKER_ONLY_POS = TRUE,
  MARKER_TEST = "wilcox",
  MARKER_MIN_PCT = 0.25,
  MARKER_LOGFC = 0.25,
  MARKER_MAX_CELLS_PER_IDENT = NULL, # set e.g. 3000 if runtime is too long
  TOP_N_MARKERS = c(5, 10, 20),

  # ------------------------- EXPORT -------------------------------------------
  OUT_DPI_RASTER = 450,
  MAIN_PDF_WIDTH = 5.0,
  MAIN_PDF_HEIGHT = 4.0,
  FONT_FAMILY = "Arial"
)

# Optional lightweight test overrides.
if (isTRUE(PARAMS$TEST_MODE)) {
  PARAMS$OUT <- file.path(PARAMS$OUT, "test_run")
}

# Normalize some parameters.
PARAMS$N_DIMS_USE <- seq_len(PARAMS$N_DIMS_USE)

# Reproducibility
set.seed(42)

# ==============================================================================
# 1. OUTPUT DIRECTORIES
# ==============================================================================

OUT <- PARAMS$OUT
OUT_FIGS    <- file.path(OUT, "main_figures")
OUT_EDFIGS  <- file.path(OUT, "extended_data_figures")
OUT_SUPPFIGS <- file.path(OUT, "supplementary_figures")
OUT_TABLES  <- file.path(OUT, "tables")
OUT_RDS     <- file.path(OUT, "rds")
OUT_LOGS    <- file.path(OUT, "logs")
OUT_METHODS <- file.path(OUT, "methods")
OUT_SAMPLES <- file.path(OUT, "sample_checkpoints")

for (d in c(OUT, OUT_FIGS, OUT_EDFIGS, OUT_SUPPFIGS, OUT_TABLES, OUT_RDS, OUT_LOGS, OUT_METHODS, OUT_SAMPLES)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

GROUP_PAL <- c("CN" = "#4DBBD5", "AD" = "#E64B35",
               "NC" = "#4DBBD5", "UC" = "#E64B35",
               "D1" = "#00A087", "D2" = "#E64B35", "Case" = "#999999")

LOG_FILE <- file.path(OUT_LOGS, "node01_complete.log")
log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

write_progress <- function(step, msg) {
  line <- paste0("[STEP ", step, "] ", msg)
  log_msg("\n", strrep("=", nchar(line)), "\n", line, "\n", strrep("=", nchar(line)))
}

# ==============================================================================
# 2. PACKAGES
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
  library(harmony)
  library(celda)
  library(qs)
})

# ==============================================================================
# 3. FIGURE STYLE HELPERS
# ==============================================================================

publication_palette <- c(
  "blue" = "#0072B2",
  "orange" = "#E69F00",
  "green" = "#009E73",
  "vermillion" = "#D55E00",
  "purple" = "#CC79A7",
  "sky" = "#56B4E9",
  "yellow" = "#F0E442",
  "black" = "#000000"
)

theme_publication <- function(base_size = 7, family = PARAMS$FONT_FAMILY) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = family, colour = "black"),
      axis.title = ggplot2::element_text(size = base_size + 0.5),
      axis.text = ggplot2::element_text(size = base_size - 0.5),
      axis.line = ggplot2::element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, colour = "black"),
      legend.title = ggplot2::element_text(size = base_size - 0.2),
      legend.text = ggplot2::element_text(size = base_size - 0.8),
      plot.title = ggplot2::element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.margin = ggplot2::margin(5, 10, 5, 5),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size, face = "bold"),
      panel.grid = ggplot2::element_blank()
    )
}

save_fig <- function(p, filename, outdir, w = 5, h = 4, dpi = PARAMS$OUT_DPI_RASTER) {
  pdf_file <- file.path(outdir, paste0(filename, ".pdf"))
  png_file <- file.path(outdir, paste0(filename, ".png"))
  svg_file <- file.path(outdir, paste0(filename, ".svg"))
  ggplot2::ggsave(pdf_file, plot = p, width = w, height = h, device = cairo_pdf, bg = "white")
  ggplot2::ggsave(png_file, plot = p, width = w, height = h, dpi = dpi, bg = "white")
  ggplot2::ggsave(svg_file, plot = p, width = w, height = h, device = "svg", bg = "white")
  log_msg("  [FIG] ", filename)
}

# ==============================================================================
# 4. SEURAT V5 SAFE ACCESSORS / HELPERS
# ==============================================================================

safe_create_assay <- function(counts_mat) {
  if (!inherits(counts_mat, "dgCMatrix")) counts_mat <- as(counts_mat, "dgCMatrix")
  obj <- tryCatch(
    SeuratObject::CreateAssay5Object(counts = counts_mat),
    error = function(e) SeuratObject::CreateAssayObject(counts = counts_mat)
  )
  obj
}

safe_join_layers <- function(seu, assay = NULL) {
  if (is.null(assay)) assay <- DefaultAssay(seu)
  DefaultAssay(seu) <- assay

  lyr <- tryCatch(SeuratObject::Layers(seu[[assay]]), error = function(e) NULL)
  if (is.null(lyr) || length(lyr) <= 1) return(seu)

  need_join <- FALSE
  if (sum(grepl("^counts", lyr)) > 1) need_join <- TRUE
  if (sum(grepl("^data", lyr)) > 1) need_join <- TRUE
  if (sum(grepl("^scale", lyr)) > 1) need_join <- TRUE

  if (need_join) {
    seu <- tryCatch(
      SeuratObject::JoinLayers(seu, assay = assay),
      error = function(e) {
        log_msg("  [WARN] JoinLayers failed on assay ", assay, ": ", conditionMessage(e))
        seu
      }
    )
  }
  seu
}

safe_get_layer <- function(seu, assay = NULL, layer = "counts") {
  if (is.null(assay)) assay <- DefaultAssay(seu)
  mat <- tryCatch(
    SeuratObject::LayerData(seu, assay = assay, layer = layer),
    error = function(e1) {
      tryCatch(
        Seurat::GetAssayData(seu, assay = assay, layer = layer),
        error = function(e2) {
          slot_name <- switch(layer,
                              counts = "counts",
                              data = "data",
                              `scale.data` = "scale.data",
                              layer)
          tryCatch(
            Seurat::GetAssayData(seu, assay = assay, slot = slot_name),
            error = function(e3) NULL
          )
        }
      )
    }
  )
  mat
}

safe_get_counts <- function(seu, assay = NULL) {
  mat <- safe_get_layer(seu, assay = assay, layer = "counts")
  if (is.null(mat)) stop("Cannot retrieve counts for assay: ", ifelse(is.null(assay), DefaultAssay(seu), assay))
  if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
  mat
}

safe_get_data <- function(seu, assay = NULL) {
  mat <- safe_get_layer(seu, assay = assay, layer = "data")
  if (!is.null(mat) && !inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
  }
  mat
}

# ==============================================================================
# 5. IO HELPERS
# ==============================================================================

strip_suffixes <- function(x) {
  x <- gsub("\\.h5$", "", x, ignore.case = TRUE)
  x <- gsub("_filtered_feature_bc_matrix$", "", x, ignore.case = TRUE)
  x <- gsub("_feature_bc_matrix$", "", x, ignore.case = TRUE)
  x <- gsub("filtered_feature_bc_matrix$", "", x, ignore.case = TRUE)
  x <- gsub("[\\/]+$", "", x)
  basename(x)
}

is_display_gene <- function(g) {
  if (grepl("^(LINC|AC[0-9]|AL[0-9]|AP[0-9]|MIR[0-9]|SNORD|SNORA|SCARNA|RNU[0-9]|RN7S|RPPH1|RMRP|VTRNA|MT-|MTRNR)", g, ignore.case = TRUE)) return(FALSE)
  if (grepl("-AS[0-9]?$|-DT$|^RP[0-9]+-", g, ignore.case = TRUE)) return(FALSE)
  if (grepl("^(RPS|RPL|MRPS|MRPL)[0-9]", g, ignore.case = TRUE)) return(FALSE)
  TRUE
}

infer_group_from_name <- function(x) {
  xl <- tolower(x)
  if (grepl("(^|[_-])(nc|cn|ctrl|control|normal|wt)([_-]|$)", xl)) return("CN")
  if (grepl("(^|[_-])(ad|alz|alzheimer|uc)([_-]|$)", xl)) return("AD")
  return("Case")
}

standardize_group_names <- function(meta) {
  if ("Group" %in% colnames(meta)) {
    g <- tolower(as.character(meta$Group))
    meta$Group <- ifelse(g %in% c("control", "ctrl", "normal", "wt", "nc", "cn"), "CN",
                  ifelse(g %in% c("ad", "alzheimer", "alz", "uc"), "AD", as.character(meta$Group)))
  }
  if ("Dataset" %in% colnames(meta)) {
    meta$Dataset <- ifelse(meta$Dataset == "GSE237718", "D1", meta$Dataset)
  }
  if ("SubjectID" %in% colnames(meta) && "Dataset" %in% colnames(meta)) {
    meta$SubjectID <- ifelse(grepl("^GSE237718_", meta$SubjectID),
                             gsub("^GSE237718_", "D1_", meta$SubjectID),
                             meta$SubjectID)
  }
  meta
}

infer_dataset_from_path <- function(x) {
  parts <- strsplit(normalizePath(x, winslash = "/", mustWork = FALSE), "/", fixed = TRUE)[[1]]
  hits <- parts[grepl("^(D\\d+|GSE\\d+)$", parts, ignore.case = TRUE)]
  if (length(hits) > 0) return(hits[[1]])
  return("Dataset1")
}

infer_subject_from_sample <- function(x) {
  x2 <- gsub("(_NC|_AD|_UC|_CASE)$", "", x, ignore.case = TRUE)
  x2 <- gsub("(_rep\\d+)$", "", x2, ignore.case = TRUE)
  x2
}

read_10x_any <- function(path) {
  if (!file.exists(path)) stop("Input path does not exist: ", path)

  obj <- NULL
  if (dir.exists(path)) {
    obj <- Seurat::Read10X(path)
  } else if (grepl("\\.h5$", path, ignore.case = TRUE)) {
    obj <- Seurat::Read10X_h5(path)
  } else {
    stop("Unsupported input path: ", path)
  }

  if (is.list(obj)) {
    if ("Gene Expression" %in% names(obj)) {
      obj <- obj[["Gene Expression"]]
    } else {
      obj <- obj[[1]]
    }
  }
  if (!inherits(obj, "dgCMatrix")) obj <- as(obj, "dgCMatrix")
  obj
}

discover_inputs <- function(base_dir) {
  if (is.na(base_dir) || !nzchar(base_dir)) stop("BASE_INPUT_DIR is empty.")
  if (!dir.exists(base_dir)) stop("BASE_INPUT_DIR does not exist: ", base_dir)

  dirs <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
  tenx_dirs <- dirs[grepl("filtered_feature_bc_matrix$", dirs, ignore.case = TRUE)]
  h5_files <- list.files(base_dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

  all_paths <- unique(c(tenx_dirs, h5_files))
  if (length(all_paths) == 0) stop("No 10x directories or .h5 files found under: ", base_dir)

  df <- data.frame(
    input_path = all_paths,
    stringsAsFactors = FALSE
  )
  df$SampleID <- vapply(df$input_path, strip_suffixes, character(1))
  df$Dataset <- vapply(df$input_path, infer_dataset_from_path, character(1))
  df$Group <- vapply(df$SampleID, infer_group_from_name, character(1))
  df$SubjectID <- vapply(df$SampleID, infer_subject_from_sample, character(1))
  df <- standardize_group_names(df)
  df
}

load_sample_sheet <- function(path_csv) {
  df <- data.table::fread(path_csv) |> as.data.frame()
  req <- c("SampleID", "input_path")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop("Sample sheet is missing columns: ", paste(miss, collapse = ", "))

  if (!"Dataset" %in% colnames(df)) df$Dataset <- vapply(df$input_path, infer_dataset_from_path, character(1))
  if (!"Group" %in% colnames(df)) df$Group <- vapply(df$SampleID, infer_group_from_name, character(1))
  if (!"SubjectID" %in% colnames(df)) df$SubjectID <- vapply(df$SampleID, infer_subject_from_sample, character(1))
  df <- standardize_group_names(df)
  df
}

# ==============================================================================
# 6. QC HELPERS
# ==============================================================================

get_species_patterns <- function(species = "human") {
  species <- tolower(species)
  if (species == "mouse") {
    list(
      mt = "^mt-",
      ribo = "^Rpl|^Rps",
      hb = "^Hb[ab]"
    )
  } else {
    list(
      mt = "^MT-",
      ribo = "^RPL|^RPS",
      hb = "^HB[ABDEGQZ]"
    )
  }
}

calc_qc_metrics <- function(seu, species = PARAMS$SPECIES) {
  pats <- get_species_patterns(species)
  seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = pats$mt)
  seu[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seu, pattern = pats$ribo)
  seu[["percent.hb"]] <- Seurat::PercentageFeatureSet(seu, pattern = pats$hb)
  seu
}

mad_cutoffs <- function(x, nmads_low = 3, nmads_high = 3, lower_bound = -Inf, upper_bound = Inf, log1p_transform = TRUE) {
  y <- x
  if (log1p_transform) y <- log10(y + 1)
  med <- stats::median(y, na.rm = TRUE)
  madv <- stats::mad(y, center = med, constant = 1.4826, na.rm = TRUE)
  lo <- med - nmads_low * madv
  hi <- med + nmads_high * madv
  if (log1p_transform) {
    lo <- 10^lo - 1
    hi <- 10^hi - 1
  }
  lo <- max(lo, lower_bound)
  hi <- min(hi, upper_bound)
  c(lower = lo, upper = hi)
}

get_sample_qc_thresholds <- function(meta) {
  if (!PARAMS$USE_ADAPTIVE_QC) {
    return(list(
      nCount_low = PARAMS$ABS_MIN_COUNTS,
      nCount_high = PARAMS$ABS_MAX_COUNTS,
      nFeature_low = PARAMS$ABS_MIN_FEATURES,
      nFeature_high = PARAMS$ABS_MAX_FEATURES,
      mt_high = PARAMS$DEFAULT_MAX_MT,
      ribo_high = PARAMS$DEFAULT_MAX_RIBO,
      hb_high = PARAMS$DEFAULT_MAX_HB
    ))
  }

  cnt <- mad_cutoffs(meta$nCount_RNA,
                     nmads_low = PARAMS$MAD_NCOUNT_LOW,
                     nmads_high = PARAMS$MAD_NCOUNT_HIGH,
                     lower_bound = PARAMS$ABS_MIN_COUNTS,
                     upper_bound = PARAMS$ABS_MAX_COUNTS,
                     log1p_transform = TRUE)
  feat <- mad_cutoffs(meta$nFeature_RNA,
                      nmads_low = PARAMS$MAD_NFEATURE_LOW,
                      nmads_high = PARAMS$MAD_NFEATURE_HIGH,
                      lower_bound = PARAMS$ABS_MIN_FEATURES,
                      upper_bound = PARAMS$ABS_MAX_FEATURES,
                      log1p_transform = TRUE)

  list(
    nCount_low = cnt["lower"],
    nCount_high = cnt["upper"],
    nFeature_low = feat["lower"],
    nFeature_high = feat["upper"],
    mt_high = PARAMS$DEFAULT_MAX_MT,
    ribo_high = PARAMS$DEFAULT_MAX_RIBO,
    hb_high = PARAMS$DEFAULT_MAX_HB
  )
}

apply_sample_qc_filter <- function(seu) {
  meta <- seu@meta.data
  th <- get_sample_qc_thresholds(meta)

  keep <- (
    meta$nCount_RNA >= th$nCount_low &
      meta$nCount_RNA <= th$nCount_high &
      meta$nFeature_RNA >= th$nFeature_low &
      meta$nFeature_RNA <= th$nFeature_high &
      meta$percent.mt <= th$mt_high &
      meta$percent.ribo <= th$ribo_high &
      meta$percent.hb <= th$hb_high
  )
  keep[is.na(keep)] <- FALSE

  list(
    object = subset(seu, cells = rownames(meta)[keep]),
    thresholds = th,
    keep = keep
  )
}

make_qc_plots_before_after <- function(meta_before, meta_after, sid) {
  dummy_before <- Matrix::Matrix(0, nrow = 1, ncol = nrow(meta_before), sparse = TRUE)
  rownames(dummy_before) <- "DUMMY"
  colnames(dummy_before) <- rownames(meta_before)
  obj_before <- Seurat::CreateSeuratObject(dummy_before, meta.data = meta_before)

  dummy_after <- Matrix::Matrix(0, nrow = 1, ncol = nrow(meta_after), sparse = TRUE)
  rownames(dummy_after) <- "DUMMY"
  colnames(dummy_after) <- rownames(meta_after)
  obj_after <- Seurat::CreateSeuratObject(dummy_after, meta.data = meta_after)

  p_vln_before <- Seurat::VlnPlot(
    object = obj_before,
    features = intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), colnames(meta_before)),
    pt.size = 0,
    ncol = 5
  ) + theme_publication(7)

  p_vln_after <- Seurat::VlnPlot(
    object = obj_after,
    features = intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), colnames(meta_after)),
    pt.size = 0,
    ncol = 5
  ) + theme_publication(7)

  p_scatter_before <- ggplot2::ggplot(meta_before, ggplot2::aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
    ggplot2::geom_point(size = 0.15, alpha = 0.5) +
    ggplot2::scale_color_viridis_c(option = "D") +
    theme_publication(7) +
    ggplot2::labs(title = paste0("QC before filtering: ", sid), x = "nCount_RNA", y = "nFeature_RNA", color = "%MT")

  p_scatter_after <- ggplot2::ggplot(meta_after, ggplot2::aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
    ggplot2::geom_point(size = 0.15, alpha = 0.5) +
    ggplot2::scale_color_viridis_c(option = "D") +
    theme_publication(7) +
    ggplot2::labs(title = paste0("QC after filtering: ", sid), x = "nCount_RNA", y = "nFeature_RNA", color = "%MT")

  list(
    vln_before = p_vln_before,
    vln_after = p_vln_after,
    scatter_before = p_scatter_before,
    scatter_after = p_scatter_after
  )
}

# ==============================================================================
# 7. SINGLE-CELL CONVERSION HELPERS
# ==============================================================================

make_sce_from_seurat <- function(seu, assay = NULL) {
  if (is.null(assay)) assay <- DefaultAssay(seu)
  seu <- safe_join_layers(seu, assay = assay)
  counts_mat <- safe_get_counts(seu, assay = assay)

  keep_nonzero <- Matrix::colSums(counts_mat) > 0
  if (!all(keep_nonzero)) {
    seu <- subset(seu, cells = colnames(seu)[keep_nonzero])
    counts_mat <- counts_mat[, keep_nonzero, drop = FALSE]
  }

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts_mat),
    colData = S4Vectors::DataFrame(seu@meta.data)
  )
  sce
}

run_scDblFinder_samplewise <- function(seu, sid) {
  DefaultAssay(seu) <- "RNA"
  seu <- safe_join_layers(seu, assay = "RNA")
  counts_mat <- safe_get_counts(seu, assay = "RNA")

  keep_nonzero <- Matrix::colSums(counts_mat) > 0
  if (!all(keep_nonzero)) {
    seu <- subset(seu, cells = colnames(seu)[keep_nonzero])
    counts_mat <- counts_mat[, keep_nonzero, drop = FALSE]
  }

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts_mat),
    colData = S4Vectors::DataFrame(seu@meta.data)
  )

  artificial_n <- max(
    PARAMS$MIN_ARTIFICIAL_DOUBLETS,
    min(PARAMS$MAX_ARTIFICIAL_DOUBLETS, as.integer(round(ncol(sce) * 0.75)))
  )

  bp <- BiocParallel::SerialParam(RNGseed = 123456, progressbar = FALSE)

  sce <- scDblFinder::scDblFinder(
    sce,
    samples = NULL,
    dbr = PARAMS$DOUBLET_RATE,
    dbr.sd = PARAMS$DOUBLET_DBR_SD,
    clusters = PARAMS$SCDBL_CLUSTERS,
    nfeatures = PARAMS$SCDBL_NFEATURES,
    artificialDoublets = artificial_n,
    BPPARAM = bp,
    verbose = TRUE
  )

  seu$scDblFinder.score <- sce$scDblFinder.score
  seu$scDblFinder.class <- sce$scDblFinder.class
  if ("scDblFinder.weighted" %in% colnames(SummarizedExperiment::colData(sce))) {
    seu$scDblFinder.weighted <- sce$scDblFinder.weighted
  }

  meta <- seu@meta.data
  p_hist <- ggplot2::ggplot(meta, ggplot2::aes(x = scDblFinder.score, fill = scDblFinder.class)) +
    ggplot2::geom_histogram(position = "identity", bins = 60, alpha = 0.85) +
    ggplot2::scale_fill_manual(values = c("singlet" = publication_palette["sky"], "doublet" = publication_palette["vermillion"])) +
    theme_publication(7) +
    ggplot2::labs(title = paste0("scDblFinder score: ", sid), x = "scDblFinder score", y = "Cells", fill = NULL)

  dbl_rate <- mean(seu$scDblFinder.class == "doublet", na.rm = TRUE)
  class_tab <- as.data.frame(table(seu$scDblFinder.class), stringsAsFactors = FALSE)
  colnames(class_tab) <- c("Class", "Cells")

  seu_singlet <- subset(seu, subset = scDblFinder.class == "singlet")
  list(object = seu_singlet, doublet_rate = dbl_rate, class_table = class_tab)
}

run_decontx_samplewise <- function(seu, sid) {
  DefaultAssay(seu) <- "RNA"
  seu <- safe_join_layers(seu, assay = "RNA")
  sce <- make_sce_from_seurat(seu, assay = "RNA")

  set.seed(123456)
  sce <- celda::decontX(sce, verbose = FALSE)

  contam_col <- NULL
  cd <- as.data.frame(SingleCellExperiment::colData(sce))
  if ("decontX_contamination" %in% colnames(cd)) {
    contam_col <- "decontX_contamination"
  } else if ("contamination" %in% colnames(cd)) {
    contam_col <- "contamination"
  }

  if (!is.null(contam_col)) {
    seu$decontX_contamination <- cd[[contam_col]]
  } else {
    seu$decontX_contamination <- NA_real_
  }

  corrected_counts <- celda::decontXcounts(sce)
  corrected_counts <- round(corrected_counts)
  if (!inherits(corrected_counts, "dgCMatrix")) corrected_counts <- as(corrected_counts, "dgCMatrix")

  seu[["decontX"]] <- safe_create_assay(corrected_counts)

  p_contam <- ggplot2::ggplot(seu@meta.data, ggplot2::aes(x = decontX_contamination)) +
    ggplot2::geom_histogram(bins = 60, fill = publication_palette["blue"], color = "white", linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = PARAMS$CONTAM_CUTOFF, color = publication_palette["vermillion"], linetype = 2, linewidth = 0.4) +
    theme_publication(7) +
    ggplot2::labs(title = paste0("decontX contamination: ", sid), x = "decontX contamination", y = "Cells")

  keep <- is.na(seu$decontX_contamination) | seu$decontX_contamination < PARAMS$CONTAM_CUTOFF
  seu <- subset(seu, cells = colnames(seu)[keep])
  seu
}

# ==============================================================================
# 8. LOAD OR BUILD SAMPLE OBJECTS AFTER SAMPLE QC
# ==============================================================================

t0 <- Sys.time()
log_msg("============================================================")
log_msg("node01_complete_final_v5_publication_gptV5.R")
log_msg("R version             : ", R.version.string)
log_msg("Seurat version        : ", as.character(packageVersion("Seurat")))
log_msg("SeuratObject version  : ", as.character(packageVersion("SeuratObject")))
log_msg("Harmony version       : ", as.character(packageVersion("harmony")))
log_msg("scDblFinder version   : ", as.character(packageVersion("scDblFinder")))
log_msg("celda version         : ", as.character(packageVersion("celda")))
log_msg("Run mode              : ", PARAMS$RESUME_MODE)
log_msg("Output directory      : ", OUT)
log_msg("============================================================")

seu_list_filt <- NULL
qc_summary_list <- list()

if (identical(PARAMS$RESUME_MODE, "raw")) {
  write_progress(1, "Loading raw sample inputs and running sample-level QC")

  sample_table <- if (!is.na(PARAMS$SAMPLE_SHEET) && nzchar(PARAMS$SAMPLE_SHEET) && file.exists(PARAMS$SAMPLE_SHEET)) {
    load_sample_sheet(PARAMS$SAMPLE_SHEET)
  } else {
    discover_inputs(PARAMS$BASE_INPUT_DIR)
  }

  sample_table <- sample_table |>
    dplyr::distinct(SampleID, .keep_all = TRUE) |>
    dplyr::arrange(Dataset, SampleID)

  data.table::fwrite(sample_table, file.path(OUT_TABLES, "sample_input_manifest_used.csv"))

  seu_list_filt <- list()
  sample_qc_rows <- list()

  for (i in seq_len(nrow(sample_table))) {
    info <- sample_table[i, , drop = FALSE]
    sid <- info$SampleID[[1]]
    in_path <- info$input_path[[1]]
    ds <- info$Dataset[[1]]
    grp <- info$Group[[1]]
    subj <- info$SubjectID[[1]]

    log_msg("\n  [LOAD] ", sid, " <- ", in_path)

    counts <- read_10x_any(in_path)
    seu <- Seurat::CreateSeuratObject(
      counts = counts,
      min.cells = PARAMS$MIN_CELLS_PER_GENE,
      min.features = PARAMS$MIN_GENES_PER_CELL,
      project = PARAMS$PROJECT
    )

    seu$SampleID <- sid
    seu$orig.ident <- sid
    seu$Dataset <- ds
    seu$Group <- grp
    seu$SubjectID <- subj

    extra_cols <- setdiff(colnames(info), c("SampleID", "input_path", "Dataset", "Group", "SubjectID"))
    if (length(extra_cols) > 0) {
      for (cc in extra_cols) seu[[cc]] <- info[[cc]][[1]]
    }

    seu <- calc_qc_metrics(seu, species = PARAMS$SPECIES)

    meta_before <- seu@meta.data
    n_before <- ncol(seu)

    qc_res <- apply_sample_qc_filter(seu)
    seu2 <- qc_res$object
    th <- qc_res$thresholds
    meta_after <- seu2@meta.data
    n_after <- ncol(seu2)

    plots <- make_qc_plots_before_after(meta_before, meta_after, sid)

    sample_qc_rows[[sid]] <- data.frame(
      SampleID = sid,
      Dataset = ds,
      Group = grp,
      SubjectID = subj,
      Cells_before_QC = n_before,
      Cells_after_QC = n_after,
      nCount_low = th$nCount_low,
      nCount_high = th$nCount_high,
      nFeature_low = th$nFeature_low,
      nFeature_high = th$nFeature_high,
      mt_high = th$mt_high,
      ribo_high = th$ribo_high,
      hb_high = th$hb_high,
      median_nCount_before = median(meta_before$nCount_RNA, na.rm = TRUE),
      median_nCount_after = median(meta_after$nCount_RNA, na.rm = TRUE),
      median_nFeature_before = median(meta_before$nFeature_RNA, na.rm = TRUE),
      median_nFeature_after = median(meta_after$nFeature_RNA, na.rm = TRUE),
      median_percent_mt_before = median(meta_before$percent.mt, na.rm = TRUE),
      median_percent_mt_after = median(meta_after$percent.mt, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    if (ncol(seu2) == 0) {
      log_msg("  [WARN] ", sid, " has 0 cells after QC and will be skipped.")
    } else {
      seu_list_filt[[sid]] <- seu2
      if (PARAMS$SAVE_SAMPLE_CHECKPOINTS) {
        qsave(seu2, file.path(OUT_SAMPLES, paste0("sample_after_sample_QC_", sid, ".qs")))
      }
    }

    rm(counts, seu, seu2, meta_before, meta_after, plots, qc_res)
    gc()
  }

  qc_tbl <- dplyr::bind_rows(sample_qc_rows)
  data.table::fwrite(qc_tbl, file.path(OUT_TABLES, "sample_qc_filtering_summary.csv"))

  if (length(seu_list_filt) == 0) stop("No samples survived sample-level QC.")

  # --- Exclude low-quality samples (after QC documentation, before merge) ---
  # quality-control: all QC plots/tables above document WHY these are excluded.
  SAMPLES_TO_EXCLUDE <- c("D2_UC1", "D2_UC6")  # MT=7.81%, nCount=591, nFeature=419 (dead/dying cells)
  for (bad in SAMPLES_TO_EXCLUDE) {
    if (bad %in% names(seu_list_filt)) {
      log_msg("  [EXCLUDE] ", bad, ": removed after QC review (see QC figures above for justification)")
      seu_list_filt[[bad]] <- NULL
    }
  }

  # Save a merge checkpoint compatible with resume mode.
  write_progress(2, "Saving post-QC merged checkpoint")
  if (length(seu_list_filt) == 1) {
    scRNA_qc <- seu_list_filt[[1]]
  } else {
    scRNA_qc <- merge(
      x = seu_list_filt[[1]],
      y = seu_list_filt[-1],
      add.cell.ids = names(seu_list_filt),
      project = PARAMS$PROJECT
    )
  }
  qsave(scRNA_qc, file.path(OUT_RDS, "node1_merged_after_sample_qc.qs"))
  if (PARAMS$SAVE_INTERMEDIATE_RDS) {
    qsave(seu_list_filt, file.path(OUT_RDS, "node1_samples_after_sample_qc_list.qs"))
  }
  log_msg("  Merged QC-passed cells: ", ncol(scRNA_qc))

  rm(scRNA_qc)
  gc()

} else if (identical(PARAMS$RESUME_MODE, "merged_qc")) {
  write_progress(1, "Resuming from node1_merged_after_sample_qc.qs")

  if (!file.exists(PARAMS$MERGED_QC_RDS)) stop("MERGED_QC_RDS not found: ", PARAMS$MERGED_QC_RDS)
  scRNA_qc <- qread(PARAMS$MERGED_QC_RDS)

  if (!"SampleID" %in% colnames(scRNA_qc@meta.data)) {
    stop("MERGED_QC_RDS is missing SampleID in metadata.")
  }
  if (!"orig.ident" %in% colnames(scRNA_qc@meta.data)) {
    scRNA_qc$orig.ident <- scRNA_qc$SampleID
  }
  if (!"Dataset" %in% colnames(scRNA_qc@meta.data)) {
    scRNA_qc$Dataset <- vapply(scRNA_qc$SampleID, infer_dataset_from_path, character(1))
  }
  if (!"Group" %in% colnames(scRNA_qc@meta.data)) {
    scRNA_qc$Group <- vapply(scRNA_qc$SampleID, infer_group_from_name, character(1))
  }
  if (!"SubjectID" %in% colnames(scRNA_qc@meta.data)) {
    scRNA_qc$SubjectID <- vapply(scRNA_qc$SampleID, infer_subject_from_sample, character(1))
  }

  for (col in c("Group", "Dataset")) {
    if (col %in% colnames(scRNA_qc@meta.data)) {
      scRNA_qc[[col]] <- standardize_group_names(scRNA_qc@meta.data)[[col]]
    }
  }

  # --- Exclude low-quality samples (after metadata validation, before processing) ---
  SAMPLES_TO_EXCLUDE <- c("D2_UC1", "D2_UC6")  # MT=7.81%, nCount=591, nFeature=419 (dead/dying cells)
  for (bad in SAMPLES_TO_EXCLUDE) {
    if (bad %in% unique(scRNA_qc$SampleID)) {
      n_before <- ncol(scRNA_qc)
      scRNA_qc <- subset(scRNA_qc, subset = SampleID != bad)
      log_msg("  [EXCLUDE] ", bad, ": removed after QC review (", n_before, " -> ", ncol(scRNA_qc), " cells)")
    }
  }

  scRNA_qc <- safe_join_layers(scRNA_qc, assay = "RNA")
  seu_list_filt <- SplitObject(scRNA_qc, split.by = "SampleID")
  if (PARAMS$SAVE_INTERMEDIATE_RDS) {
    qsave(seu_list_filt, file.path(OUT_RDS, "node1_samples_after_sample_qc_list_from_resume.qs"))
  }
  log_msg("  Samples recovered from merged checkpoint: ", length(seu_list_filt))

  rm(scRNA_qc)
  gc()

} else {
  stop("PARAMS$RESUME_MODE must be 'raw' or 'merged_qc'.")
}

# ==============================================================================
# 9. SAMPLE-WISE DOUBLETS + DECONTX
# ==============================================================================

write_progress(3, "Running sample-wise scDblFinder and decontX")

sample_post_list <- list()
sample_flow_rows <- list()
doublet_class_rows <- list()
decontx_rows <- list()

sample_ids <- names(seu_list_filt)
for (sid in sample_ids) {
  log_msg("\n  [SAMPLE] ", sid)

  seu <- seu_list_filt[[sid]]
  DefaultAssay(seu) <- "RNA"

  n_before_doublet <- ncol(seu)
  dbl_res <- run_scDblFinder_samplewise(seu, sid)
  seu_singlet <- dbl_res$object
  dbl_rate <- dbl_res$doublet_rate
  dbl_class_tab <- dbl_res$class_table
  n_after_doublet <- ncol(seu_singlet)

  doublet_tab <- dbl_class_tab
  doublet_tab$SampleID <- sid
  doublet_class_rows[[sid]] <- doublet_tab[, c("SampleID", "Class", "Cells")]

  seu_final <- seu_singlet
  n_after_decontx <- n_after_doublet
  median_decontx <- NA_real_

  if (isTRUE(PARAMS$RUN_DECONTX)) {
    seu_final <- run_decontx_samplewise(seu_singlet, sid)
    n_after_decontx <- ncol(seu_final)
    if ("decontX_contamination" %in% colnames(seu_final@meta.data)) {
      median_decontx <- median(seu_final$decontX_contamination, na.rm = TRUE)
    }
  }

  sample_post_list[[sid]] <- seu_final
  sample_flow_rows[[sid]] <- data.frame(
    SampleID = sid,
    Dataset = unique(seu_final$Dataset)[1],
    Group = unique(seu_final$Group)[1],
    Cells_after_sample_QC = n_before_doublet,
    Cells_after_doublet = n_after_doublet,
    Cells_after_decontX = n_after_decontx,
    Doublet_rate = dbl_rate,
    Median_decontX = median_decontx,
    stringsAsFactors = FALSE
  )

  if ("decontX_contamination" %in% colnames(seu_final@meta.data)) {
    decontx_rows[[sid]] <- data.frame(
      SampleID = sid,
      Mean_contamination = mean(seu_final$decontX_contamination, na.rm = TRUE),
      Median_contamination = median(seu_final$decontX_contamination, na.rm = TRUE),
      Max_contamination = max(seu_final$decontX_contamination, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  if (PARAMS$SAVE_SAMPLE_CHECKPOINTS) {
    qsave(seu_final, file.path(OUT_SAMPLES, paste0("sample_after_doublet_decontX_", sid, ".qs")))
  }

  rm(seu, seu_singlet, seu_final, dbl_res)
  gc()
}

flow_tbl <- dplyr::bind_rows(sample_flow_rows)
doublet_tbl <- dplyr::bind_rows(doublet_class_rows)
decontx_tbl <- if (length(decontx_rows) > 0) dplyr::bind_rows(decontx_rows) else data.frame()

data.table::fwrite(flow_tbl, file.path(OUT_TABLES, "sample_flow_doublet_decontx_summary.csv"))
data.table::fwrite(doublet_tbl, file.path(OUT_TABLES, "doublet_class_by_sample.csv"))
data.table::fwrite(decontx_tbl, file.path(OUT_TABLES, "decontX_summary_by_sample.csv"))

# ==============================================================================
# 10. MERGE FILTERED SAMPLES
# ==============================================================================

write_progress(4, "Merging samples after sample-wise scDblFinder / decontX")

if (length(sample_post_list) == 0) stop("No samples available after doublet/decontX filtering.")

if (length(sample_post_list) == 1) {
  scRNA <- sample_post_list[[1]]
} else {
  scRNA <- merge(
    x = sample_post_list[[1]],
    y = sample_post_list[-1],
    add.cell.ids = names(sample_post_list),
    project = PARAMS$PROJECT
  )
}

qsave(scRNA, file.path(OUT_RDS, "node1_after_samplewise_doublet_decontx_merge.qs"))
log_msg("  Merged filtered object cells: ", ncol(scRNA))

# ==============================================================================
# 11. SAMPLE FLOW / QC SUMMARY FIGURES
# ==============================================================================

write_progress(5, "Generating sample flow and preprocessing summary figures")

p_flow_long <- flow_tbl |>
  tidyr::pivot_longer(
    cols = c("Cells_after_sample_QC", "Cells_after_doublet", "Cells_after_decontX"),
    names_to = "Step",
    values_to = "Cells"
  ) |>
  dplyr::mutate(
    Step = factor(Step,
                  levels = c("Cells_after_sample_QC", "Cells_after_doublet", "Cells_after_decontX"),
                  labels = c("After sample QC", "After doublet", "After decontX"))
  )

p_flow <- p_flow_long |>
  ggplot2::ggplot(ggplot2::aes(x = Step, y = Cells, group = SampleID, color = SampleID)) +
  ggplot2::geom_line(linewidth = 0.45, alpha = 0.8) +
  ggplot2::geom_point(size = 1.4) +
  theme_publication(7) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 5), legend.position = "none") +
  ggplot2::labs(title = "Cell number flow across preprocessing steps", x = NULL, y = "Cells")
save_fig(p_flow, "supp_Node1_sample_cell_flow_all_samples", OUT_SUPPFIGS, w = 6.2, h = 4.4)

p_flow_summary <- p_flow_long |>
  dplyr::group_by(Step) |>
  dplyr::summarise(Cells = sum(Cells), .groups = "drop") |>
  ggplot2::ggplot(ggplot2::aes(x = Step, y = Cells / 1000, group = 1)) +
  ggplot2::geom_area(fill = "grey85", alpha = 0.8) +
  ggplot2::geom_line(color = "grey20", linewidth = 0.6) +
  ggplot2::geom_point(color = "grey20", size = 1.8) +
  ggplot2::geom_text(ggplot2::aes(label = format(Cells, big.mark = ",")), vjust = -0.8, size = 2.4) +
  theme_publication(7) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)) +
  ggplot2::labs(title = "Cell retention through preprocessing", x = NULL, y = "Nuclei/cells (x10^3)")
save_fig(p_flow_summary, "Fig1A_Node1_cell_retention_summary", OUT_FIGS, w = 4.0, h = 3.0)

p_doublet_rate <- ggplot2::ggplot(flow_tbl, ggplot2::aes(x = reorder(SampleID, Doublet_rate), y = 100 * Doublet_rate)) +
  ggplot2::geom_col(fill = publication_palette["vermillion"], width = 0.72) +
  ggplot2::coord_flip() +
  theme_publication(7) +
  ggplot2::labs(title = "Doublet rate by sample", x = NULL, y = "Doublet rate (%)")
save_fig(p_doublet_rate, "ED_Node1_doublet_rate_by_sample", OUT_EDFIGS, w = 4.8, h = max(4.5, 0.18 * nrow(flow_tbl) + 1.5))

if ("decontX_contamination" %in% colnames(scRNA@meta.data)) {
  p_decontx_all <- ggplot2::ggplot(scRNA@meta.data, ggplot2::aes(x = decontX_contamination)) +
    ggplot2::geom_histogram(bins = 60, fill = publication_palette["blue"], color = "white", linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = PARAMS$CONTAM_CUTOFF, color = publication_palette["vermillion"], linetype = 2, linewidth = 0.4) +
    theme_publication(7) +
    ggplot2::labs(title = "decontX contamination distribution", x = "decontX contamination", y = "Cells")
  save_fig(p_decontx_all, "ED_Node1_decontX_contamination_hist", OUT_EDFIGS, w = 4.2, h = 3.0)
}

if (all(c("nCount_RNA", "nFeature_RNA", "percent.mt") %in% colnames(scRNA@meta.data))) {
  suppressMessages({
    qc_long <- scRNA@meta.data |>
      dplyr::select(SampleID, nCount_RNA, nFeature_RNA, percent.mt) |>
      tidyr::pivot_longer(cols = -SampleID, names_to = "Metric", values_to = "Value")
  })
  for (m in c("nCount_RNA", "nFeature_RNA", "percent.mt")) {
    sub <- qc_long |> dplyr::filter(Metric == m)
    p_qc <- ggplot2::ggplot(sub, ggplot2::aes(x = SampleID, y = Value, fill = SampleID)) +
      ggplot2::geom_violin(scale = "width", trim = TRUE) +
      theme_publication(7) + ggplot2::theme(legend.position = "none",
                                        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 5)) +
      ggplot2::labs(title = m, x = NULL, y = m)
    save_fig(p_qc, paste0("supp_QC_violin_", m), OUT_SUPPFIGS, w = 11, h = 4.8)
  }
  p_qc_scatter <- ggplot2::ggplot(scRNA@meta.data, ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
    ggplot2::geom_point(size = 0.3, alpha = 0.5) +
    ggplot2::scale_color_gradient(low = "grey80", high = "red") +
    theme_publication(7) + ggplot2::labs(title = "QC: nCount vs nFeature (all samples)", color = "MT%")
  save_fig(p_qc_scatter, "supp_QC_scatter_nCount_vs_nFeature", OUT_SUPPFIGS, w = 12, h = 8)
  log_msg("  Combined QC violin/scatter plots added to supplementary_figures")
}

# ==============================================================================
# 12. NORMALIZATION / HVG / PCA / HARMONY / UMAP / CLUSTERING
# ==============================================================================

write_progress(6, "Normalization, HVG, PCA, Harmony, clustering, UMAP")

analysis_assay <- PARAMS$ANALYSIS_ASSAY
if (!(analysis_assay %in% Seurat::Assays(scRNA))) analysis_assay <- "RNA"
DefaultAssay(scRNA) <- analysis_assay
scRNA <- safe_join_layers(scRNA, assay = analysis_assay)

# Always (re)normalize to ensure 'data' layer exists in decontX assay
scRNA <- Seurat::NormalizeData(
    scRNA,
    assay = analysis_assay,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

scRNA <- Seurat::FindVariableFeatures(
  scRNA,
  assay = analysis_assay,
  selection.method = "vst",
  nfeatures = PARAMS$N_HVG,
  verbose = FALSE
)

vars_to_regress <- intersect(PARAMS$REGRESS_VARS, colnames(scRNA@meta.data))
if (length(vars_to_regress) > 0) {
  scRNA <- Seurat::ScaleData(
    scRNA,
    assay = analysis_assay,
    features = VariableFeatures(scRNA),
    vars.to.regress = vars_to_regress,
    verbose = FALSE
  )
} else {
  scRNA <- Seurat::ScaleData(
    scRNA,
    assay = analysis_assay,
    features = VariableFeatures(scRNA),
    verbose = FALSE
  )
}

scRNA <- Seurat::RunPCA(
  scRNA,
  assay = analysis_assay,
  features = VariableFeatures(scRNA),
  npcs = PARAMS$N_PCS,
  verbose = FALSE
)

p_elbow <- Seurat::ElbowPlot(scRNA, ndims = min(PARAMS$N_PCS, 50)) +
  theme_publication(7) +
  ggplot2::labs(title = "PCA elbow plot")
save_fig(p_elbow, "ED_Node1_PCA_elbow", OUT_EDFIGS, w = 3.8, h = 3.0)

p_pca_s <- Seurat::DimPlot(scRNA, reduction = "pca", group.by = "SampleID", raster = TRUE) +
  theme_publication(7) + ggplot2::labs(title = "PCA before Harmony (by sample)")
save_fig(p_pca_s, "ED_Node1_PCA_by_sample_before_Harmony", OUT_EDFIGS, w = 8, h = 6)

if ("Group" %in% colnames(scRNA@meta.data)) {
  p_pca_g <- Seurat::DimPlot(scRNA, reduction = "pca", group.by = "Group", cols = GROUP_PAL, raster = TRUE) +
    theme_publication(7) + ggplot2::labs(title = "PCA before Harmony (by group)")
  save_fig(p_pca_g, "supp_PCA_by_group_before_Harmony", OUT_SUPPFIGS, w = 8, h = 6)
}

scRNA <- harmony::RunHarmony(
  object = scRNA,
  group.by.vars = "SampleID",
  reduction = "pca",
  dims.use = PARAMS$N_DIMS_USE,
  reduction.save = "harmony",
  project.dim = TRUE,
  verbose = FALSE
)

p_har_s <- Seurat::DimPlot(scRNA, reduction = "harmony", group.by = "SampleID", raster = TRUE) +
  theme_publication(7) + ggplot2::labs(title = "Harmony embedding (by sample)")
save_fig(p_har_s, "ED_Node1_Harmony_by_sample", OUT_EDFIGS, w = 8, h = 6)

if ("Group" %in% colnames(scRNA@meta.data)) {
  p_har_g <- Seurat::DimPlot(scRNA, reduction = "harmony", group.by = "Group", cols = GROUP_PAL, raster = TRUE) +
    theme_publication(7) + ggplot2::labs(title = "Harmony embedding (by group)")
  save_fig(p_har_g, "supp_Harmony_by_group", OUT_SUPPFIGS, w = 8, h = 6)
}

scRNA <- Seurat::FindNeighbors(
  scRNA,
  reduction = "harmony",
  dims = PARAMS$N_DIMS_USE,
  verbose = FALSE
)

scRNA <- Seurat::FindClusters(
  scRNA,
  resolution = PARAMS$CLUST_RES,
  verbose = FALSE
)

scRNA <- Seurat::RunUMAP(
  scRNA,
  reduction = "harmony",
  dims = PARAMS$N_DIMS_USE,
  n.neighbors = PARAMS$UMAP_N_NEIGHBORS,
  min.dist = PARAMS$UMAP_MIN_DIST,
  verbose = FALSE
)

# ==============================================================================
# 13. MAIN FIGURES / EXTENDED DATA FIGURES
# ==============================================================================

write_progress(7, "Generating integrated embedding, QC, and composition figures")

p_umap_cluster <- Seurat::DimPlot(
  scRNA,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 0.15
) + theme_publication(7) + ggplot2::theme(aspect.ratio = 1) + 
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(title = "Clusters")
save_fig(p_umap_cluster, "Fig2_Node1_UMAP_clusters", OUT_FIGS, w = 5.2, h = 4.8)

p_umap_sample <- Seurat::DimPlot(
  scRNA,
  reduction = "umap",
  group.by = "SampleID",
  raster = TRUE,
  pt.size = 0.15
) + theme_publication(7) + ggplot2::theme(aspect.ratio = 1) + ggplot2::labs(title = "Sample")
save_fig(p_umap_sample, "Fig3_Node1_UMAP_sample", OUT_FIGS, w = 5.6, h = 4.8)

if ("Group" %in% colnames(scRNA@meta.data)) {
  scRNA$Group <- factor(as.character(scRNA$Group), levels = intersect(c("CN", "AD"), unique(as.character(scRNA$Group))))
  group_levels <- levels(droplevels(scRNA$Group))
  group_cols <- GROUP_PAL[group_levels]

  p_umap_group <- Seurat::DimPlot(
    scRNA,
    reduction = "umap",
    group.by = "Group",
    cols = group_cols,
    raster = TRUE,
    pt.size = 0.25
  ) + theme_publication(7) + ggplot2::theme(aspect.ratio = 1) + ggplot2::labs(title = "Group")
  save_fig(p_umap_group, "Fig4_Node1_UMAP_group_CN_AD", OUT_FIGS, w = 4.8, h = 4.3)
}

if ("Dataset" %in% colnames(scRNA@meta.data)) {
  p_umap_dataset <- Seurat::DimPlot(
    scRNA,
    reduction = "umap",
    group.by = "Dataset",
    raster = TRUE,
    pt.size = 0.15
  ) + theme_publication(7) + ggplot2::theme(aspect.ratio = 1) + ggplot2::labs(title = "Dataset")
  save_fig(p_umap_dataset, "ED_Node1_UMAP_dataset", OUT_EDFIGS, w = 4.8, h = 4.2)
}

if ("Group" %in% colnames(scRNA@meta.data)) {
  p_umap_split <- Seurat::DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters",
                                   split.by = "Group", ncol = 2, raster = TRUE) +
    theme_publication(7) + ggplot2::labs(title = "UMAP split by group")
  save_fig(p_umap_split, "supp_UMAP_split_by_group", OUT_SUPPFIGS, w = 12, h = 6)
}

qc_features <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hb", "scDblFinder.score", "decontX_contamination"), colnames(scRNA@meta.data))
if (length(qc_features) > 0) {
  fp_list <- lapply(qc_features, function(feat) {
    Seurat::FeaturePlot(
      scRNA,
      features = feat,
      reduction = "umap",
      raster = TRUE,
      order = TRUE,
      pt.size = 0.15
    ) + theme_publication(7) + ggplot2::theme(aspect.ratio = 1) + ggplot2::labs(title = feat)
  })
  p_qc_overlay <- patchwork::wrap_plots(fp_list, ncol = 2)
  save_fig(p_qc_overlay, "ED_Node1_UMAP_QC_feature_overlays", OUT_EDFIGS, w = 7.2, h = 3.4 * ceiling(length(fp_list) / 2))
}

qc_violin_features <- intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt", "decontX_contamination", "scDblFinder.score"), colnames(scRNA@meta.data))
if (length(qc_violin_features) > 0) {
  p_qc_vln <- Seurat::VlnPlot(scRNA, features = qc_violin_features, group.by = "seurat_clusters", pt.size = 0, ncol = min(3, length(qc_violin_features))) +
    theme_publication(7) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_fig(p_qc_vln, "ED_Node1_QC_violin_by_cluster", OUT_EDFIGS, w = 8.0, h = 2.8 * ceiling(length(qc_violin_features) / min(3, length(qc_violin_features))))
}

if (all(c("nCount_RNA", "nFeature_RNA", "percent.mt", "SampleID") %in% colnames(scRNA@meta.data))) {
  qc_after_df <- scRNA@meta.data
  p_scatter_after <- ggplot2::ggplot(qc_after_df, ggplot2::aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
    ggplot2::geom_point(size = 0.15, alpha = 0.5) +
    ggplot2::scale_color_viridis_c(option = "D") +
    ggplot2::facet_wrap(~SampleID, scales = "free") +
    theme_publication(7) +
    ggplot2::labs(title = "QC scatter after filtering", x = "nCount_RNA", y = "nFeature_RNA", color = "%MT")
  save_fig(p_scatter_after, "ED_Node1_QC_scatter_after_filtering_all_samples", OUT_EDFIGS, w = 11.0, h = 7.8)
}

# Cell composition summaries.
comp_sample_cluster <- scRNA@meta.data |>
  dplyr::count(SampleID, seurat_clusters, name = "Cells") |>
  dplyr::group_by(SampleID) |>
  dplyr::mutate(Prop = Cells / sum(Cells)) |>
  dplyr::ungroup()
data.table::fwrite(comp_sample_cluster, file.path(OUT_TABLES, "composition_sample_by_cluster.csv"))

comp_group_cluster <- scRNA@meta.data |>
  dplyr::count(Group, seurat_clusters, name = "Cells") |>
  dplyr::group_by(Group) |>
  dplyr::mutate(Prop = Cells / sum(Cells)) |>
  dplyr::ungroup()
data.table::fwrite(comp_group_cluster, file.path(OUT_TABLES, "composition_group_by_cluster.csv"))

p_comp_group <- ggplot2::ggplot(comp_group_cluster, ggplot2::aes(x = Group, y = Prop, fill = seurat_clusters)) +
  ggplot2::geom_col(width = 0.65) +
  theme_publication(7) +
  ggplot2::theme(legend.position = "right") +
  ggplot2::labs(title = "Unannotated cluster composition by group", x = NULL, y = "Cluster fraction", fill = "Cluster")
save_fig(p_comp_group, "Fig1B_Node1_cluster_composition_by_group_CN_AD", OUT_FIGS, w = 4.6, h = 3.2)

p_comp <- ggplot2::ggplot(comp_sample_cluster, ggplot2::aes(x = SampleID, y = Prop, fill = seurat_clusters)) +
  ggplot2::geom_col(width = 0.9) +
  theme_publication(7) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(title = "Cluster composition by sample", x = NULL, y = "Proportion")
save_fig(p_comp, "ED_Node1_cluster_composition_by_sample", OUT_EDFIGS, w = 8.0, h = 4.2)

# ==============================================================================
# 14. MARKERS
# ==============================================================================

write_progress(8, "Running marker detection and exporting top marker figures")

markers_all <- data.frame()
if (isTRUE(PARAMS$RUN_MARKERS)) {
  Seurat::Idents(scRNA) <- "seurat_clusters"
  rna <- scRNA[["RNA"]]
  rna_lyr <- SeuratObject::Layers(rna)
  if (length(rna_lyr) > 1 && sum(grepl("^counts", rna_lyr)) > 1) {
    cat("  [MARKERS] Joining RNA layers...\n")
    scRNA <- SeuratObject::JoinLayers(scRNA, assay = "RNA")
  }
  cat("  [MARKERS] Normalizing RNA assay...\n")
  scRNA <- Seurat::NormalizeData(scRNA, assay = "RNA",
                                  normalization.method = "LogNormalize",
                                  scale.factor = 10000, verbose = TRUE)
  DefaultAssay(scRNA) <- "RNA"
  markers_all <- Seurat::FindAllMarkers(
    scRNA,
    assay = "RNA",
    only.pos = PARAMS$MARKER_ONLY_POS,
    test.use = PARAMS$MARKER_TEST,
    min.pct = PARAMS$MARKER_MIN_PCT,
    logfc.threshold = PARAMS$MARKER_LOGFC,
    max.cells.per.ident = PARAMS$MARKER_MAX_CELLS_PER_IDENT,
    verbose = TRUE
  )
}

data.table::fwrite(markers_all, file.path(OUT_TABLES, "Markers_all_clusters.csv"))
if (isTRUE(PARAMS$RUN_MARKERS) && nrow(markers_all) == 0) {
  stop("Marker detection returned zero rows. This is not acceptable for the final reproducible node01 output; inspect RNA layers/identities or set a deliberate marker fallback before publication export.")
}

if (nrow(markers_all) > 0) {
  markers_display <- markers_all |>
    dplyr::filter(is_display_gene(gene))
  
  for (nn in PARAMS$TOP_N_MARKERS) {
    topn <- markers_display |>
      dplyr::group_by(cluster) |>
      dplyr::slice_max(order_by = avg_log2FC, n = nn, with_ties = FALSE) |>
      dplyr::ungroup()
    data.table::fwrite(topn, file.path(OUT_TABLES, paste0("Markers_top", nn, "_each_cluster.csv")))
  }

  top5 <- markers_display |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) |>
    dplyr::ungroup()

  genes_dot <- unique(top5$gene)
  if (length(genes_dot) > 0) {
    p_dot <- Seurat::DotPlot(
      scRNA,
      features = rev(genes_dot),
      assay = analysis_assay,
      group.by = "seurat_clusters",
      cols = c("lightgrey", publication_palette["vermillion"]),
      dot.scale = 4
    ) +
      theme_publication(7) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 5),
        axis.text.y = ggplot2::element_text(size = 6)
      ) +
      ggplot2::labs(title = "Top 5 marker genes per cluster", x = "Gene", y = "Cluster")
    save_fig(p_dot, "Fig5_Node1_marker_dotplot_top5", OUT_FIGS, w = 8.0, h = 4.6)
  }
  
  # ---- Canonical brain marker dotplot ----
  canonical_genes <- c("SLC17A7", "SATB2", "GAD1", "GAD2", "PVALB", "SST", "VIP",
                       "GFAP", "AQP4", "MOG", "PLP1", "MBP", "PDGFRA", "CSPG4",
                       "AIF1", "CD68", "P2RY12", "CLDN5", "VWF", "PDGFRB", "RGS5",
                       "SNAP25", "SYT1", "NRGN")
  canonical_avail <- intersect(canonical_genes, rownames(scRNA))
  if (length(canonical_avail) >= 5) {
    scRNA_tmp <- suppressWarnings(tryCatch({
      rna_cts <- Seurat::GetAssayData(scRNA, assay = "RNA", layer = "counts")
      obj <- Seurat::CreateSeuratObject(rna_cts, meta.data = scRNA@meta.data)
      obj$seurat_clusters <- scRNA$seurat_clusters
      Seurat::NormalizeData(obj, verbose = FALSE)
    }, error = function(e) NULL))
    if (!is.null(scRNA_tmp)) {
      Seurat::Idents(scRNA_tmp) <- "seurat_clusters"
      p_canon <- Seurat::DotPlot(scRNA_tmp, features = rev(canonical_avail), assay = "RNA",
                                  group.by = "seurat_clusters",
                                  cols = c("lightgrey", publication_palette["blue"])) +
        theme_publication(7) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7)) +
        ggplot2::labs(title = "Canonical brain cell marker dotplot", x = NULL)
      save_fig(p_canon, "ED_Node1_canonical_brain_marker_dotplot", OUT_EDFIGS, w = 13, h = 5.5)
      rm(scRNA_tmp); gc()
    }
  }
  
  # ---- Marker heatmap (top5) ----
  genes_hm <- unique(top5$gene)
  if (length(genes_hm) > 5 && length(genes_hm) <= 300) {
    cluster_levels <- sort(unique(as.character(scRNA$seurat_clusters)))
    cluster_cols <- setNames(rep(publication_palette, length.out = length(cluster_levels)), cluster_levels)
    p_heat <- Seurat::DoHeatmap(scRNA, features = genes_hm, assay = analysis_assay,
                                 group.by = "seurat_clusters", group.colors = cluster_cols,
                                 raster = TRUE, size = 3) +
      theme_publication(7) + Seurat::NoLegend() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5)) +
      ggplot2::labs(title = "Top 5 marker heatmap")
    save_fig(p_heat, "supp_Marker_heatmap_top5", OUT_SUPPFIGS, w = 10, h = 12)
  }
  
} else {
  for (nn in PARAMS$TOP_N_MARKERS) {
    data.table::fwrite(data.frame(), file.path(OUT_TABLES, paste0("Markers_top", nn, "_each_cluster.csv")))
  }
}

# ==============================================================================
# 15. FINAL METADATA / TABLES / METHODS SUMMARY
# ==============================================================================

write_progress(9, "Exporting final metadata, summaries, and methods text")

final_meta <- scRNA@meta.data
final_meta$CellID <- rownames(final_meta)
data.table::fwrite(final_meta, file.path(OUT_TABLES, "metadata_final.csv"))

has_percent_mt <- "percent.mt" %in% colnames(final_meta)
has_scdbl <- "scDblFinder.score" %in% colnames(final_meta)
has_decontx <- "decontX_contamination" %in% colnames(final_meta)

final_sample_summary <- final_meta |>
  dplyr::group_by(SampleID, Dataset, Group, SubjectID) |>
  dplyr::summarise(
    Cells = dplyr::n(),
    Median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    Median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    Median_percent_mt = if (has_percent_mt) median(percent.mt, na.rm = TRUE) else NA_real_,
    Median_scDblFinder_score = if (has_scdbl) median(scDblFinder.score, na.rm = TRUE) else NA_real_,
    Median_decontX = if (has_decontx) median(decontX_contamination, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )
data.table::fwrite(final_sample_summary, file.path(OUT_TABLES, "sample_summary_final.csv"))

cluster_summary <- final_meta |>
  dplyr::count(seurat_clusters, SampleID, Group, Dataset, name = "Cells")
data.table::fwrite(cluster_summary, file.path(OUT_TABLES, "cluster_sample_group_summary.csv"))

cluster_overview <- final_meta |>
  dplyr::count(seurat_clusters, name = "Cells") |>
  dplyr::mutate(Proportion = Cells / sum(Cells))
data.table::fwrite(cluster_overview, file.path(OUT_TABLES, "cluster_overview_summary.csv"))

method_lines <- c(
  "Node1 preprocessing summary (Seurat v5 final complete version)",
  paste0("- Project: ", PARAMS$PROJECT),
  paste0("- Run mode: ", PARAMS$RESUME_MODE),
  if (identical(PARAMS$RESUME_MODE, "merged_qc")) paste0("- Resume checkpoint: ", PARAMS$MERGED_QC_RDS) else paste0("- Raw input root: ", PARAMS$BASE_INPUT_DIR),
  "",
  "Sample loading and QC:",
  paste0("- Species mode: ", PARAMS$SPECIES),
  paste0("- CreateSeuratObject min.cells = ", PARAMS$MIN_CELLS_PER_GENE, "; min.features = ", PARAMS$MIN_GENES_PER_CELL),
  paste0("- Adaptive QC enabled: ", PARAMS$USE_ADAPTIVE_QC),
  paste0("- Mitochondrial percentage cutoff: ", PARAMS$DEFAULT_MAX_MT, "%"),
  paste0("- Ribosomal percentage cutoff: ", PARAMS$DEFAULT_MAX_RIBO, "%"),
  paste0("- Hemoglobin percentage cutoff: ", PARAMS$DEFAULT_MAX_HB, "%"),
  "",
  "Doublet and ambient RNA filtering:",
  "- scDblFinder was run sample-wise rather than on the full merged multi-sample object.",
  paste0("- scDblFinder expected doublet rate dbr = ", PARAMS$DOUBLET_RATE, "; dbr.sd = ", PARAMS$DOUBLET_DBR_SD),
  paste0("- scDblFinder nfeatures = ", PARAMS$SCDBL_NFEATURES, "; artificial doublets capped between ", PARAMS$MIN_ARTIFICIAL_DOUBLETS, " and ", PARAMS$MAX_ARTIFICIAL_DOUBLETS),
  paste0("- decontX enabled: ", PARAMS$RUN_DECONTX),
  paste0("- decontX contamination cutoff: ", PARAMS$CONTAM_CUTOFF),
  "",
  "Dimensionality reduction and clustering:",
  paste0("- Analysis assay: ", analysis_assay),
  paste0("- Highly variable genes: ", PARAMS$N_HVG),
  paste0("- PCA components computed: ", PARAMS$N_PCS),
  paste0("- Harmony integration variable: SampleID"),
  paste0("- Harmony dimensions used: 1:", max(PARAMS$N_DIMS_USE)),
  paste0("- Clustering resolution: ", PARAMS$CLUST_RES),
  paste0("- UMAP n.neighbors: ", PARAMS$UMAP_N_NEIGHBORS),
  paste0("- UMAP min.dist: ", PARAMS$UMAP_MIN_DIST),
  "",
  "Marker analysis:",
  paste0("- FindAllMarkers test = ", PARAMS$MARKER_TEST),
  paste0("- only.pos = ", PARAMS$MARKER_ONLY_POS),
  paste0("- min.pct = ", PARAMS$MARKER_MIN_PCT),
  paste0("- logfc.threshold = ", PARAMS$MARKER_LOGFC),
  paste0("- max.cells.per.ident = ", ifelse(is.null(PARAMS$MARKER_MAX_CELLS_PER_IDENT), "NULL", PARAMS$MARKER_MAX_CELLS_PER_IDENT))
)
writeLines(method_lines, con = file.path(OUT_METHODS, "node1_methods_summary.txt"))

# ==============================================================================
# 16. SAVE FINAL OBJECT AND SESSION INFO
# ==============================================================================

write_progress(10, "Saving final Seurat object and session info")

DefaultAssay(scRNA) <- analysis_assay
qsave(scRNA, file.path(OUT_RDS, "scRNA_node1_final.qs"))
# saveRDS(scRNA, file.path(OUT_RDS, "test_scRNA2_final_processed.rds"))  # duplicate, removed per user
capture.output(sessionInfo(), file = file.path(OUT_LOGS, "sessionInfo_node01.txt"))

elapsed_min <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)

write_progress(11, "Pipeline completed")
log_msg("Final object : ", file.path(OUT_RDS, "scRNA_node1_final.qs"))
log_msg("Final object : ", file.path(OUT_RDS, "scRNA_node1_final.qs"))
log_msg("Final cells  : ", ncol(scRNA))
log_msg("Final genes  : ", nrow(scRNA))
log_msg("Default assay: ", DefaultAssay(scRNA))
log_msg("Clusters     : ", length(unique(scRNA$seurat_clusters)))
log_msg("Resolution   : ", PARAMS$CLUST_RES)
log_msg("Elapsed      : ", elapsed_min, " minutes")
log_msg("Output dir   : ", OUT)
log_msg("[RESUME] All done.")

################################################################################
# END
################################################################################

