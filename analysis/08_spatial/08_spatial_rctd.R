#!/usr/bin/env Rscript
################################################################################
# node08_spatial_v2.R - RCTD 14-subtype deconvolution + spatial QC + MISTy
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
  library(mistyR)
  library(Matrix)
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(ggpubr)
  library(pheatmap)
  library(stringr)
  library(scales)
})

set.seed(1234)

# ==============================================================================
# 0. GLOBAL CONFIGURATION
# ==============================================================================

TEST_MODE <- TRUE

OUT <- "/path/to/project/results/08_spatial"
SP  <- "/path/to/public_data/spatial/GSE220442/counts_and_images"

SC_ANNOTATED <- "/path/to/project/results/02_annotation/rds/scRNA_annotated.qs"
ENDO_QS <- "/path/to/project/results/04_subtypes/rds/subtype_Endothelial.qs"
ASTRO_QS <- "/path/to/project/results/04_subtypes/rds/subtype_Astrocytes.qs"
MICRO_QS <- "/path/to/project/results/04_subtypes/rds/subtype_Microglia.qs"

OUT_REF <- file.path(OUT, "reference_rctd"); OUT_RDS <- file.path(OUT, "rds_rctd")
OUT_FIG <- file.path(OUT, "figures_rctd"); OUT_TAB <- file.path(OUT, "tables_rctd")
OUT_LOG <- file.path(OUT, "logs"); OUT_MISTY <- file.path(OUT, "misty_rctd")
for (d in c(OUT_REF, OUT_RDS, OUT_FIG, OUT_TAB, OUT_LOG, OUT_MISTY)) dir.create(d, showWarnings=FALSE, recursive=TRUE)

SAMPLES <- c("1-1"="Control","18-64"="Control","2-5"="Control","2-3"="AD","2-8"="AD","T4857"="AD")
GROUP_COLORS <- c("Control"="#4DBBD5","AD"="#E64B35")
MAX_CELLS_PER_TYPE <- 10000; RCTD_MAX_CORES <- 1; RCTD_DOUBLET_MODE <- "full"

# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================

log_msg <- function(m) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), m))
save_csv <- function(x, name) write.csv(x, file.path(OUT_TAB, name), row.names=FALSE, quote=FALSE)

save_fig <- function(p, name, w=8, h=6) {
  ggsave(file.path(OUT_FIG, paste0(name,".pdf")), p, width=w, height=h, dpi=300, limitsize=FALSE, useDingbats=FALSE)
  ggsave(file.path(OUT_FIG, paste0(name,".png")), p, width=w, height=h, dpi=300, limitsize=FALSE)
  log_msg(paste("Saved:", name))
}

save_heatmap <- function(ph, name, w=8, h=7) {
  pdf(file.path(OUT_FIG, paste0(name,".pdf")), width=w, height=h, useDingbats=FALSE); print(ph); dev.off()
  png(file.path(OUT_FIG, paste0(name,".png")), width=w, height=h, units="in", res=300); print(ph); dev.off()
  log_msg(paste("Saved heatmap:", name))
}

theme_publication <- function(base=11) theme_classic(base_size=base) + theme(
  plot.title=element_text(face="bold",hjust=0.5,size=base+1), plot.subtitle=element_text(hjust=0.5),
  axis.title=element_text(face="bold",color="black"), axis.text=element_text(color="black"),
  legend.title=element_text(face="bold"), legend.text=element_text(color="black"),
  strip.background=element_rect(fill="grey95",color="grey70"), strip.text=element_text(face="bold",color="black"),
  panel.border=element_rect(color="black",fill=NA,linewidth=0.3)
)

get_counts <- function(obj) tryCatch(GetAssayData(obj, assay="RNA", layer="counts"), error=function(e) GetAssayData(obj, assay="RNA", slot="counts"))
get_data   <- function(obj) tryCatch(GetAssayData(obj, assay="RNA", layer="data"), error=function(e) GetAssayData(obj, assay="RNA", slot="data"))

as_matrix <- function(mat) {
  row_pos <- mat@i; col_pos <- findInterval(seq(mat@x)-1, mat@p[-1])
  tmp <- matrix(0, nrow=mat@Dim[1], ncol=mat@Dim[2])
  for (k in seq_along(mat@x)) tmp[row_pos[k]+1, col_pos[k]+1] <- mat@x[k]
  rownames(tmp) <- mat@Dimnames[[1]]; colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

sample_cells <- function(obj, group_col, max_n=10000) {
  grps <- as.character(obj[[group_col]][,1]); cells <- colnames(obj)
  keep <- unlist(lapply(unique(grps), function(g) {
    cg <- cells[grps==g]; if (length(cg)>max_n) sample(cg, max_n) else cg
  }))
  subset(obj, cells=keep)
}

# ==============================================================================
# 2. LOGGING
# ==============================================================================

sink(file.path(OUT_LOG, "session_v2.txt"), split=TRUE)
cat("node08 RCTD v2 - 14 subtypes + spatial QC + MISTy\n")
cat("Start:", format(Sys.time()), "\n")
cat("RCTD mode:", RCTD_DOUBLET_MODE, "| cores:", RCTD_MAX_CORES, "\n\n")

# ==============================================================================
# 3. BUILD 14-SUBTYPE RCTD REFERENCE
# ==============================================================================

log_msg("=== Step 1: Building 14-subtype RCTD Reference ===")

log_msg("Loading scRNA annotated (for non-NVU cells)...")
sc <- qread(SC_ANNOTATED)
gc()

# Strategy: extract counts + labels per compartment
# Combine into single matrix WITHOUT Seurat merge (avoids memory explosion)

combined_counts <- NULL
combined_types  <- c()
type_counts <- list()

# --- Non-NVU: major cell types ---
for (ct_name in c("Oligodendrocytes","OPCs","Inhibitory_neurons","Excitatory_neurons")) {
  log_msg(paste("  Extracting", ct_name, "..."))
  sub_obj <- subset(sc, subset=celltype==ct_name)
  cnt <- get_counts(sub_obj)
  n_cells <- ncol(cnt)
  type_counts[[ct_name]] <- n_cells
  if (is.null(combined_counts)) {
    combined_counts <- cnt
  } else {
    common_genes <- intersect(rownames(combined_counts), rownames(cnt))
    combined_counts <- cbind(combined_counts[common_genes,,drop=FALSE], cnt[common_genes,,drop=FALSE])
  }
  combined_types <- c(combined_types, rep(ct_name, n_cells))
  rm(sub_obj, cnt); gc()
}
rm(sc); gc()

# --- NVU: subtypes ---
log_msg("Loading NVU subtype objects...")

add_nvu_cells <- function(obj, subtype_col) {
  subtypes <- as.character(obj[[subtype_col]][,1])
  for (st in unique(subtypes)) {
    log_msg(paste("  Extracting", st, "..."))
    sub_obj <- subset(obj, subset=!!sym(subtype_col)==st)
    cnt <- get_counts(sub_obj)
    n_cells <- ncol(cnt)
    type_counts[[st]] <<- n_cells
    common_genes <- intersect(rownames(combined_counts), rownames(cnt))
    combined_counts <<- cbind(combined_counts[common_genes,,drop=FALSE], cnt[common_genes,,drop=FALSE])
    combined_types  <<- c(combined_types, rep(st, n_cells))
    rm(sub_obj, cnt); gc()
  }
}

endo <- qread(ENDO_QS); add_nvu_cells(endo, "ec_subtype"); rm(endo); gc()
astro <- qread(ASTRO_QS); add_nvu_cells(astro, "astro_subtype"); rm(astro); gc()
micro <- qread(MICRO_QS); add_nvu_cells(micro, "micro_subtype"); rm(micro); gc()

# Remove types with <25 cells
names(combined_types) <- colnames(combined_counts)
tbl <- table(combined_types)
log_msg(paste("Subtype counts:", paste(names(tbl), tbl, sep="=", collapse=", ")))

bad_types <- names(tbl)[tbl < 25]
if (length(bad_types) > 0) {
  log_msg(paste("Removing types with <25 cells:", paste(bad_types, collapse=", ")))
  keep_cells <- !combined_types %in% bad_types
  combined_counts <- combined_counts[, keep_cells, drop=FALSE]
  combined_types  <- combined_types[keep_cells]
}

# Downsample per type (max 10000) to keep RCTD happy
set.seed(1234)
keep_idx <- unlist(lapply(unique(combined_types), function(ct) {
  idx <- which(combined_types==ct)
  if (length(idx) > MAX_CELLS_PER_TYPE) sample(idx, MAX_CELLS_PER_TYPE) else idx
}))
combined_counts <- combined_counts[, keep_idx, drop=FALSE]
combined_types  <- combined_types[keep_idx]
combined_types   <- factor(combined_types)

ref_nUMI <- Matrix::colSums(combined_counts)

log_msg(sprintf("Reference: %d cells, %d types", ncol(combined_counts), length(levels(combined_types))))

reference <- Reference(combined_counts, combined_types, ref_nUMI)
log_msg("RCTD Reference built successfully")
qsave(reference, file.path(OUT_REF, "rctd_reference_14type.qs"))

# --- Reference QC ---
log_msg("=== Step 2: Reference QC ===")

ref_count_df <- data.frame(subtype=names(tbl), n=as.integer(tbl))
save_csv(ref_count_df, "TableS1_reference_subtype_counts.csv")

# Composition barplot
p_comp <- ggplot(ref_count_df, aes(x=reorder(subtype,-n), y=n, fill=subtype)) +
  geom_bar(stat="identity") + theme_publication() + coord_flip() +
  labs(title="Reference composition", x="", y="Number of cells") + theme(legend.position="none")
save_fig(p_comp, "FigS1_reference_composition", 8, 6)

# Marker validation: build minimal Seurat for DotPlot
log_msg("Building minimal object for marker QC...")
ref_seu <- CreateSeuratObject(counts=combined_counts, assay="RNA", verbose=FALSE)
ref_seu$subtype <- combined_types
ref_seu <- NormalizeData(ref_seu, verbose=FALSE)
Idents(ref_seu) <- "subtype"

ref_markers <- list(
  Endo_ActLow=c("CLDN5","PECAM1","VWF"), Endo_ActMid=c("ANGPT2","KDR","FLT1"),
  Endo_ActHigh=c("VCAM1","SELE","ICAM1"), Capillary=c("CLDN5","MFSD2A","SLC2A1"),
  SMC=c("ACTA2","MYH11","TAGLN"), Astro_Homeo=c("AQP4","ALDH1L1","SLC1A3"),
  Astro_Inter=c("GFAP","VIM"), Astro_React=c("CD44","SERPINA3","VIM"),
  Micro_Homeo=c("P2RY12","CX3CR1","TMEM119"), Micro_DAM=c("APOE","SPP1","TREM2"),
  Oligo=c("MBP","MOG","PLP1"), OPC=c("PDGFRA","CSPG4"),
  Inhibitory=c("GAD1","GAD2"), Excitatory=c("SLC17A7","SATB2")
)

p_dot <- DotPlot(ref_seu, features=unlist(ref_markers), group.by="subtype") +
  RotatedAxis() + theme_publication() +
  labs(title="Reference marker validation", x="Marker genes", y="Subtype")
save_fig(p_dot, "FigS2_reference_marker_dotplot", 16, 6)

rm(ref_seu); gc()
rm(combined_counts); gc()

# ==============================================================================
# 4. PER-SAMPLE SPATIAL QC + RCTD DECONVOLUTION
# ==============================================================================

log_msg("=== Step 3: Spatial QC + RCTD per sample ===")
mt_pattern <- "^MT-"; ribo_pattern <- "^RP[SL]|^MRP"

all_weights <- list()

for (sid in names(SAMPLES)) {
  log_msg(paste("--- Processing", sid, paste0("[", SAMPLES[sid], "]"), "---"))
  
  # Load spatial data
  h5_path <- file.path(SP, sid, "filtered_feature_bc_matrix.h5")
  if (!file.exists(h5_path)) {
    log_msg(paste("  WARNING: h5 not found at", h5_path, "- checking directory..."))
    h5_path <- file.path(SP, sid, "filtered_feature_bc_matrix.h5")
  }
  
  stRNA <- tryCatch({
    Load10X_Spatial(data.dir=file.path(SP, sid), filename="filtered_feature_bc_matrix.h5",
                    slice=sid, filter.matrix=TRUE)
  }, error=function(e) {
    log_msg(paste("  Load10X_Spatial failed:", e$message, "- trying alt approach"))
    # Try reading from directory
    counts <- Read10X_h5(h5_path)
    spatial_dir <- file.path(SP, sid, "spatial")
    # Manual object construction approach
    img <- Read10X_Image(spatial_dir, image.name="tissue_lowres_image.png")
    seu <- CreateSeuratObject(counts=counts, assay="Spatial")
    seu[[sid]] <- img
    return(seu)
  })
  
  # Filter on-tissue spots
  coords_tissue <- GetTissueCoordinates(stRNA, cols=c("row","col","tissue"), scale=NULL)
  stRNA$tissue <- coords_tissue[colnames(stRNA), "tissue"]
  stRNA$tissue <- ifelse(stRNA$tissue == 1, "on_tissue", "not_on_tissue")
  log_msg(paste("  Tissue spots:", sum(stRNA$tissue=="on_tissue"), "/", ncol(stRNA)))
  stRNA <- subset(stRNA, subset=tissue=="on_tissue")
  
  # Remove MT and ribosomal genes
  mt_genes  <- rownames(stRNA)[grepl(mt_pattern, rownames(stRNA))]
  rps_genes <- rownames(stRNA)[grepl("^RPS", rownames(stRNA))]
  mrp_genes <- rownames(stRNA)[grepl("^MRP", rownames(stRNA))]
  rpl_genes <- rownames(stRNA)[grepl("^RPL", rownames(stRNA))]
  rb_genes  <- c(rps_genes, mrp_genes, rpl_genes)
  stRNA <- stRNA[!rownames(stRNA) %in% c(rb_genes, mt_genes), ]
  
  # Filter low-quality spots (genes expressed in <10 spots, then spot-level filters)
  stRNA <- stRNA[rowSums(as.matrix(stRNA[["Spatial"]]$counts) > 0) > 10, ]
  stRNA$nFeature_Spatial_filt <- colSums(stRNA[["Spatial"]]$counts > 0)
  stRNA$nCount_Spatial_filt  <- colSums(stRNA[["Spatial"]]$counts)
  stRNA$percent.mt_filt <- PercentageFeatureSet(stRNA, pattern=mt_pattern)
  
  # Spot QC filtering
  stRNA <- subset(stRNA, subset=nFeature_Spatial_filt > 50 & nCount_Spatial_filt > 150)
  log_msg(paste("  After spot QC:", ncol(stRNA), "spots"))
  
  stRNA$orig.ident <- sid; stRNA$group <- SAMPLES[sid]
  
  # Prepare for RCTD
  coords <- GetTissueCoordinates(stRNA, cols=c("row","col"), scale=NULL)
  colnames(coords) <- c("xcoord","ycoord")
  
  sp_counts <- as_matrix(stRNA[["Spatial"]]$counts)
  
  # Align spots
  common_spots <- intersect(rownames(coords), colnames(sp_counts))
  coords <- coords[common_spots, , drop=FALSE]
  sp_counts <- sp_counts[, common_spots, drop=FALSE]
  sp_nUMI <- colSums(sp_counts)
  
  # Build puck and run RCTD
  puck <- SpatialRNA(coords, sp_counts, sp_nUMI)
  myRCTD <- create.RCTD(puck, reference, max_cores=RCTD_MAX_CORES)
  
  log_msg("  Running RCTD...")
  myRCTD <- run.RCTD(myRCTD, doublet_mode=RCTD_DOUBLET_MODE)
  
  # Save RCTD object
  saveRDS(myRCTD, file.path(OUT_RDS, paste0(sid, "_rctd.rds")))
  
  # Extract weights
  w <- as.data.frame(myRCTD@results$weights)
  all_weights[[sid]] <- w
  save_csv(w, paste0("TableS_", sid, "_RCTD_weights.csv"))
  
  log_msg(sprintf("  Done: %d spots, weight range [%.4f, %.4f]", 
                  nrow(w), min(w), max(w)))
}

# Save combined weights
saveRDS(all_weights, file.path(OUT_RDS, "all_weights_list.rds"))

# ==============================================================================
# 5. SPATIAL VISUALIZATION (TIERED)
# ==============================================================================

log_msg("=== Step 4: Tiered spatial visualization ===")
subtype_colors <- c(
  "Endothelial_ActLow"="#2166AC","Endothelial_ActMid"="#4393C3","Endothelial_ActHigh"="#92C5DE",
  "Capillary"="#053061","SMC"="#B2182B",
  "Astrocytes_Homeostatic"="#E69F00","Astrocytes_Intermediate"="#F0E442","Astrocytes_Reactive"="#D55E00",
  "Microglia_Homeostatic"="#1B9E77","Microglia_DAM"="#66C2A5",
  "Oligodendrocytes"="#7570B3","OPCs"="#A6CEE3",
  "Inhibitory_neurons"="#E78AC3","Excitatory_neurons"="#B3B3B3"
)

for (sid in names(SAMPLES)) {
  log_msg(paste("Visualizing", sid))
  
  # Load RCTD
  myRCTD <- readRDS(file.path(OUT_RDS, paste0(sid, "_rctd.rds")))
  w <- myRCTD@results$weights
  
  # Load and align spatial object
  stRNA <- tryCatch({
    Load10X_Spatial(data.dir=file.path(SP, sid), filename="filtered_feature_bc_matrix.h5",
                    slice=sid, filter.matrix=TRUE)
  }, error=function(e) NULL)
  
  if (is.null(stRNA)) next
  stRNA <- subset(stRNA, cells=intersect(colnames(stRNA), rownames(w)))
  
  # Attach weights
  common <- intersect(colnames(stRNA), rownames(w))
  w_sub <- w[common, , drop=FALSE]
  stRNA <- stRNA[, common]
  
  for (ct in colnames(w_sub)) {
    ct_clean <- gsub("[^A-Za-z0-9_]","_", ct)
    stRNA[[ct_clean]] <- w_sub[, ct]
  }
  
  # Tier 1: NVU subtype proportion maps (individual, continuous scale)
  nvu_subtypes <- intersect(colnames(w_sub), 
    c("Endothelial_ActLow","Endothelial_ActMid","Endothelial_ActHigh","Capillary","SMC",
      "Astrocytes_Homeostatic","Astrocytes_Intermediate","Astrocytes_Reactive",
      "Microglia_Homeostatic","Microglia_DAM"))
  
  for (ct in nvu_subtypes) {
    ct_clean <- gsub("[^A-Za-z0-9_]","_", ct)
    p <- SpatialFeaturePlot(stRNA, features=ct_clean, image.alpha=0.3, pt.size.factor=1.6,
                            alpha=c(0.1, 1)) +
      scale_fill_viridis_c(option="B", limits=c(0, quantile(stRNA[[ct_clean]][,1], 0.99, na.rm=TRUE)),
                           oob=scales::squish, name="Proportion") +
      ggtitle(paste(sid, "-", ct), subtitle=paste0(SAMPLES[sid], " | Visium MTG")) +
      theme(plot.title=element_text(size=10, face="bold"))
    save_fig(p, paste0("Fig_Spatial_", sid, "_", ct_clean), 6, 5)
  }
  
  # Tier 2: Combined dominant map (collapsed to 7 major for readability)
  w_collapsed <- data.frame(
    Endothelial  = rowSums(w_sub[, grep("Endothelial|Capillary|SMC", colnames(w_sub)), drop=FALSE]),
    Astrocytes   = rowSums(w_sub[, grep("Astrocytes", colnames(w_sub)), drop=FALSE]),
    Microglia    = rowSums(w_sub[, grep("Microglia", colnames(w_sub)), drop=FALSE]),
    Oligodendrocytes = w_sub[, "Oligodendrocytes"],
    OPCs         = w_sub[, "OPCs"],
    Inhibitory   = w_sub[, "Inhibitory_neurons"],
    Excitatory   = w_sub[, "Excitatory_neurons"]
  )
  dominant <- apply(w_collapsed, 1, which.max)
  dominant_label <- colnames(w_collapsed)[dominant]
  
  stRNA$dominant_7 <- factor(dominant_label, levels=colnames(w_collapsed))
  
  major_colors <- c("Endothelial"="#0072B2","Astrocytes"="#E69F00","Microglia"="#009E73",
                    "Oligodendrocytes"="#7570B3","OPCs"="#A6CEE3",
                    "Inhibitory"="#E78AC3","Excitatory"="#B3B3B3")
  
  p_dom <- SpatialDimPlot(stRNA, group.by="dominant_7", image.alpha=0.3, pt.size.factor=1.6,
                          cols=major_colors) +
    ggtitle(paste(sid, "- Dominant cell type (7 classes)"), subtitle=SAMPLES[sid])
  save_fig(p_dom, paste0("Fig_Dominant_7class_", sid), 7, 5.5)
  
  rm(stRNA); gc()
}

# ==============================================================================
# 6. SAMPLE-LEVEL DIFFERENTIAL ABUNDANCE
# ==============================================================================

log_msg("=== Step 5: Sample-level differential abundance ===")

# Aggregate per sample
prop_list <- list()
for (sid in names(SAMPLES)) {
  w <- all_weights[[sid]]
  prop_list[[sid]] <- data.frame(
    sample=sid, group=SAMPLES[sid],
    t(colMeans(w))
  )
}
prop_sample <- do.call(rbind, lapply(prop_list, function(x) x))
prop_long <- pivot_longer(prop_sample, cols=-c(sample, group), names_to="subtype", values_to="proportion")

# Grouped barplot
p_bar <- ggplot(prop_long, aes(x=subtype, y=proportion, fill=group)) +
  geom_boxplot(outlier.shape=NA, alpha=0.6) +
  geom_jitter(width=0.15, size=2, alpha=0.8) +
  scale_fill_manual(values=GROUP_COLORS) + theme_publication() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title="Sample-level cell type abundance (AD vs Control)", x="", y="Mean RCTD proportion")
save_fig(p_bar, "Fig_Subtype_Abundance_AD_vs_Control", 12, 5)

# Wilcoxon test per subtype
da_results <- data.frame()
for (ct in colnames(prop_sample)[-c(1,2)]) {
  ad_vals  <- prop_sample[prop_sample$group=="AD", ct]
  ctrl_vals <- prop_sample[prop_sample$group=="Control", ct]
  if (length(ad_vals)>=3 && length(ctrl_vals)>=3) {
    wt <- wilcox.test(ad_vals, ctrl_vals)
    da_results <- rbind(da_results, data.frame(
      subtype=ct, mean_AD=mean(ad_vals), mean_Control=mean(ctrl_vals),
      log2FC=log2(mean(ad_vals)/mean(ctrl_vals)),
      p_value=wt$p.value, W=wt$statistic
    ))
  }
}
da_results$FDR <- p.adjust(da_results$p_value, method="BH")
save_csv(da_results, "Table_SampleLevel_DifferentialAbundance.csv")

# Volcano plot
da_results$sig <- ifelse(da_results$FDR < 0.05, "FDR<0.05", "NS")
p_volcano <- ggplot(da_results, aes(x=log2FC, y=-log10(FDR), color=sig, label=subtype)) +
  geom_point(size=3) + geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50") +
  geom_text_repel(size=3, max.overlaps=20) + scale_color_manual(values=c("FDR<0.05"="red","NS"="grey50")) +
  theme_publication() + labs(title="Differential abundance (sample-level)", x="log2FC (AD/Control)")
save_fig(p_volcano, "Fig_DA_Volcano", 8, 6)

# ==============================================================================
# 7. MARKER SANITY CHECK
# ==============================================================================

log_msg("=== Step 6: Marker sanity check ===")

marker_modules <- list(
  Endo_ActLow=c("CLDN5","PECAM1"), Endo_ActMid=c("ANGPT2","KDR"),
  Endo_ActHigh=c("VCAM1","SELE"), Capillary=c("CLDN5","MFSD2A"),
  SMC=c("ACTA2","MYH11"), Astro_Homeo=c("AQP4","ALDH1L1"),
  Astro_Inter=c("GFAP","VIM"), Astro_React=c("CD44","SERPINA3"),
  Micro_Homeo=c("P2RY12","CX3CR1"), Micro_DAM=c("APOE","SPP1"),
  Oligo=c("MBP","MOG"), OPC=c("PDGFRA"),
  Inhibitory=c("GAD1","GAD2"), Excitatory=c("SLC17A7")
)

cor_results <- data.frame()

for (sid in names(SAMPLES)) {
  log_msg(paste("Marker check for", sid))
  myRCTD <- readRDS(file.path(OUT_RDS, paste0(sid, "_rctd.rds")))
  w <- myRCTD@results$weights
  
  stRNA <- tryCatch({
    Load10X_Spatial(data.dir=file.path(SP, sid), filename="filtered_feature_bc_matrix.h5",
                    slice=sid, filter.matrix=TRUE)
  }, error=function(e) NULL)
  if (is.null(stRNA)) next
  
  common <- intersect(colnames(stRNA), rownames(w))
  stRNA <- stRNA[, common]; w_sub <- w[common, ]
  
  # Compute marker scores
  stRNA <- NormalizeData(stRNA, assay="Spatial", verbose=FALSE)
  for (mod_name in names(marker_modules)) {
    genes <- marker_modules[[mod_name]]
    genes <- intersect(genes, rownames(stRNA))
    if (length(genes) > 0) {
      stRNA[[paste0(mod_name,"_score")]] <- colMeans(stRNA[["Spatial"]]$data[genes, , drop=FALSE])
    }
  }
  
  # Correlate RCTD proportion vs marker score
  for (mod_name in names(marker_modules)) {
    score_name <- paste0(mod_name, "_score")
    if (!score_name %in% colnames(stRNA@meta.data)) next
    
    # Find matching RCTD column
    matching_rctd <- grep(gsub("_","", mod_name), gsub("[^A-Za-z]","", colnames(w_sub)), value=TRUE, ignore.case=TRUE)
    if (length(matching_rctd)==0) next
    
    ct_col <- matching_rctd[1]
    sp <- cor.test(stRNA[[score_name]][,1], w_sub[, ct_col], method="spearman")
    
    cor_results <- rbind(cor_results, data.frame(
      sample=sid, marker_module=mod_name, RCTD_col=ct_col,
      spearman_rho=sp$estimate, p_value=sp$p.value
    ))
  }
  
  # Sanity check figure
  san_cols <- intersect(names(marker_modules), unique(gsub("_score","",
    grep("_score$", colnames(stRNA@meta.data), value=TRUE))))
  if (length(san_cols) >= 4) {
    san_cols <- san_cols[1:min(4, length(san_cols))]
    plot_list <- list()
    for (i in seq_along(san_cols)) {
      mn <- san_cols[i]; sc <- paste0(mn,"_score")
      matching_rctd <- grep(gsub("_","",mn), gsub("[^A-Za-z]","", colnames(w_sub)), value=TRUE, ignore.case=TRUE)
      if (length(matching_rctd)==0) next
      ct_col <- matching_rctd[1]
      df <- data.frame(proportion=w_sub[, ct_col], score=stRNA[[sc]][,1])
      plot_list[[i]] <- ggplot(df, aes(x=proportion, y=score)) +
        geom_point(size=0.5, alpha=0.5) + geom_smooth(method="lm", se=TRUE, color="red") +
        theme_publication(9) + labs(x=paste(mn,"RCTD proportion"), y=paste(mn,"marker score"))
    }
    if (length(plot_list) > 0) {
      p_san <- wrap_plots(plot_list) + plot_annotation(title=paste(sid,"- RCTD proportion vs marker score"))
      save_fig(p_san, paste0("Fig_Sanity_", sid), 10, 8)
    }
  }
  
  rm(stRNA, myRCTD); gc()
}

if (nrow(cor_results) > 0) {
  cor_results$FDR <- p.adjust(cor_results$p_value, method="BH")
  save_csv(cor_results, "Table_MarkerSanity_Correlation.csv")
}

# Summary heatmap
if (nrow(cor_results) > 0) {
  cor_mat <- reshape2::acast(cor_results, marker_module ~ sample, value.var="spearman_rho", fill=0)
  ph_cor <- pheatmap(cor_mat, color=colorRampPalette(c("#2166AC","white","#B2182B"))(100),
                     main="RCTD proportion vs marker score correlation (Spearman rho)")
  save_heatmap(ph_cor, "Fig_MarkerSanity_Summary", 8, 7)
}

# ==============================================================================
# 8. MISTy SPATIAL INTERACTION ANALYSIS
# ==============================================================================

log_msg("=== Step 7: MISTy spatial interaction ===")

run_misty_sample <- function(sid, weight_list, misty_out_dir) {
  w <- weight_list[[sid]]
  if (is.null(w)) return(NULL)
  
  # Load spatial object
  stRNA <- tryCatch({
    Load10X_Spatial(data.dir=file.path(SP, sid), filename="filtered_feature_bc_matrix.h5",
                    slice=sid, filter.matrix=TRUE)
  }, error=function(e) NULL)
  if (is.null(stRNA)) return(NULL)
  
  # Align
  common <- intersect(colnames(stRNA), rownames(w))
  stRNA <- stRNA[, common]; w_sub <- w[common, ]
  
  # Clean colnames for MISTy (replace special chars)
  colnames(w_sub) <- gsub("[^A-Za-z0-9_]", "_", colnames(w_sub))
  
  # Attach weights as assay
  stRNA[["RCTD_props"]] <- CreateAssayObject(t(w_sub))
  
  # Run MISTy
  assay <- "RCTD_props"
  DefaultAssay(stRNA) <- assay
  useful_features <- rownames(stRNA)
  
  misty_out <- paste0(misty_out_dir, "/", sid, "_misty_", assay)
  
  run_misty_seurat(
    visium.slide = stRNA,
    view.assays  = list("main"=assay, "juxta"=assay, "para"=assay),
    view.features = list("main"=useful_features, "juxta"=useful_features, "para"=useful_features),
    view.types   = list("main"="intra", "juxta"="juxta", "para"="para"),
    view.params  = list("main"=NULL, "juxta"=5, "para"=15),
    spot.ids     = NULL,
    out.alias    = misty_out
  )
  
  return(misty_out)
}

# Run MISTy per sample
misty_outputs <- list()
for (sid in names(SAMPLES)) {
  log_msg(paste("MISTy for", sid))
  misty_outputs[[sid]] <- tryCatch({
    run_misty_sample(sid, all_weights, OUT_MISTY)
  }, error=function(e) {
    log_msg(paste("  MISTy failed for", sid, ":", e$message))
    return(NULL)
  })
}

# Collect and visualize MISTy results
for (sid in names(misty_outputs)) {
  if (is.null(misty_outputs[[sid]])) next
  
  mout <- misty_outputs[[sid]]
  misty_res <- tryCatch(collect_results(mout), error=function(e) NULL)
  if (is.null(misty_res)) next
  
  # Improvement stats
  p_imp <- plot_improvement_stats(misty_res)
  save_fig(p_imp, paste0("Fig_MISTy_", sid, "_improvement"), 7, 5)
  
  # View contributions
  p_view <- plot_view_contributions(misty_res)
  save_fig(p_view, paste0("Fig_MISTy_", sid, "_views"), 8, 6)
  
  # Interaction heatmaps
  for (view_name in c("intra", "juxta_5", "para_15")) {
    p_int <- tryCatch(
      plot_interaction_heatmap(misty_res, view_name, cutoff=0),
      error=function(e) NULL
    )
    if (!is.null(p_int)) save_fig(p_int, paste0("Fig_MISTy_", sid, "_", view_name), 8, 7)
  }
}

# ==============================================================================
# 9. OUTPUT MANIFEST
# ==============================================================================

log_msg("=== Step 8: Output manifest ===")

fig_files <- list.files(OUT_FIG, pattern="\\.pdf$|\\.png$")
tab_files <- list.files(OUT_TAB, pattern="\\.csv$")
rds_files <- list.files(OUT_RDS, pattern="\\.rds$")

manifest <- data.frame(
  category=c(rep("figure", length(fig_files)), rep("table", length(tab_files)), rep("RDS", length(rds_files))),
  filename=c(fig_files, tab_files, rds_files),
  path=c(file.path(OUT_FIG, fig_files), file.path(OUT_TAB, tab_files), file.path(OUT_RDS, rds_files))
)
save_csv(manifest, "output_manifest.csv")

# Run summary
sink(file.path(OUT_LOG, "run_summary_v2.txt"))
cat("node08 RCTD v2 Summary\n======================\n")
cat("Completed:", format(Sys.time()), "\n")
cat("Samples processed:", length(names(SAMPLES)), "\n")
cat("Reference cell types:", length(levels(reference@cell_types)), "\n")
cat("Figures:", length(fig_files), "\n")
cat("Tables:", length(tab_files), "\n")
sink()

sink()
log_msg("=== PIPELINE COMPLETE ===")
