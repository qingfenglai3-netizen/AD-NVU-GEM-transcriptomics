#!/usr/bin/env Rscript
################################################################################
# _worker_subtype.R - isolated subprocess for one cell-type subtype analysis
# Runs ONE subtype (Cerebrovascular cells/Astrocytes/Microglia), saves qs, exits.
# Memory fully released to OS upon process exit.
################################################################################
set.seed(123456)
options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 200 * 1024^3)
paths_to_add <- "/path/to/R_libraries"
for (p in paths_to_add) {
  if (dir.exists(p) && !(p %in% .libPaths())) .libPaths(c(p, .libPaths()))
}

suppressPackageStartupMessages({
  library(Seurat); library(qs)
  library(ggplot2); library(patchwork); library(ggpubr)
  library(dplyr); library(tidyr); library(tibble)
  library(harmony); library(pheatmap); library(lme4)
  library(RColorBrewer); library(scales); library(viridis)
  library(speckle)
})

# ---- Config from env vars ----
CELLTYPE <- Sys.getenv("WORKER_CELLTYPE")
CKPT_QS  <- Sys.getenv("WORKER_CKPT_QS")
OUT      <- Sys.getenv("WORKER_OUT")
OUT_FIG  <- Sys.getenv("WORKER_OUT_FIG")
OUT_ED   <- Sys.getenv("WORKER_OUT_ED")
OUT_TAB  <- Sys.getenv("WORKER_OUT_TAB")
OUT_RDS  <- Sys.getenv("WORKER_OUT_RDS")
OUT_SRC  <- Sys.getenv("WORKER_OUT_SRC")
LOG_FILE <- Sys.getenv("WORKER_LOG_FILE")

for (d in c(OUT_FIG, OUT_ED, OUT_TAB, OUT_RDS, OUT_SRC))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) {
  msg <- sprintf("[%s][WORKER:%s] %s", format(Sys.time(), "%H:%M:%S"), CELLTYPE, paste0(...))
  cat(msg, "\n"); flush.console()
}

# ---- Helpers (exact copies from main script) ----
publication_theme <- function(base = 7) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      axis.text   = element_text(color = "black", size = base),
      axis.title  = element_text(color = "black", size = base + 0.5),
      axis.line   = element_line(color = "black", linewidth = 0.35),
      axis.ticks  = element_line(color = "black", linewidth = 0.35),
      legend.text = element_text(size = base - 0.5),
      legend.title= element_text(size = base),
      legend.key.size = unit(0.32, "cm"),
      strip.background = element_blank(),
      strip.text  = element_text(size = base, face = "bold"),
      plot.title  = element_text(size = base + 1, face = "bold", hjust = 0.5),
      panel.grid  = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
}

group_palette <- c("CN" = "#4DBBD5", "AD" = "#E64B35",
                   "NC" = "#4DBBD5", "UC" = "#E64B35",
                   "Ctrl" = "#4DBBD5",
                   "Control" = "#4DBBD5", "Disease" = "#E64B35", "ND" = "#4DBBD5")

standardize_group <- function(x) {
  v <- tolower(trimws(as.character(x)))
  v[v %in% c("cn", "nc", "nd", "ctrl", "control", "normal", "healthy")] <- "CN"
  v[v %in% c("ad", "uc", "disease", "case", "alz", "alzheimers")] <- "AD"
  v
}

save_fig <- function(p, name, dir, w = 6, h = 4) {
  pdf_f  <- file.path(dir, paste0(name, ".pdf"))
  png_f  <- file.path(dir, paste0(name, ".png"))
  svg_f  <- file.path(dir, paste0(name, ".svg"))
  tiff_f <- file.path(dir, paste0(name, ".tiff"))
  ggsave(pdf_f, p, width = w, height = h, units = "in", device = cairo_pdf)
  ggsave(png_f, p, width = w, height = h, units = "in", dpi = 300)
  tryCatch(ggsave(svg_f, p, width = w, height = h, units = "in"),
           error = function(e) log_msg(sprintf("SVG skipped (%s)", name)))
  tryCatch(ggsave(tiff_f, p, width = w, height = h, units = "in", dpi = 600,
                  compression = "lzw", device = "tiff"),
           error = function(e) log_msg(sprintf("TIFF skipped (%s)", name)))
  log_msg(sprintf("Saved: %s", name))
}

safe_write_csv <- function(x, f)
  tryCatch(write.csv(x, f, row.names = FALSE),
           error = function(e) log_msg(sprintf("CSV save failed: %s", f)))

save_source <- function(df, fig_name)
  safe_write_csv(df, file.path(OUT_SRC, paste0(fig_name, "_source_data.csv")))

ave_col <- function(mat, key) {
  key <- as.character(key)
  cn <- colnames(mat)
  if (key %in% cn) return(key)
  gk <- paste0("g", key)
  if (gk %in% cn) return(gk)
  if (key %in% as.character(as.numeric(cn))) return(key)
  return(NA_character_)
}

prepare_v5 <- function(obj, a = "RNA") {
  DefaultAssay(obj) <- a
  if (inherits(obj[[a]], "Assay5")) obj <- JoinLayers(obj, assay = a)
  obj
}

harmony_emb <- function(obj, gbv = "orig.ident", du = 1:15, th = 2) {
  pca_emb <- Embeddings(obj, "pca")[, du]
  harm_emb <- harmony::RunHarmony(pca_emb, obj@meta.data,
                                  vars_use = gbv, theta = th,
                                  verbose = FALSE)
  colnames(harm_emb) <- paste0("harmony_", 1:ncol(harm_emb))
  rownames(harm_emb) <- colnames(obj)
  obj[["harmony"]] <- CreateDimReducObject(harm_emb, key = "harmony_", assay = "RNA")
  obj
}

subcluster_one <- function(obj, ct_name, res = 0.5, dims_pca = 1:30, dims_h = 1:15) {
  obj <- prepare_v5(obj)
  obj <- NormalizeData(obj, verbose = FALSE) %>%
    FindVariableFeatures(nfeatures = 2000, verbose = FALSE)
  hvg <- setdiff(VariableFeatures(obj),
                 grep("^MT-|^RPS|^RPL|^MALAT1|^FOS|^JUN|^HSP",
                      rownames(obj), value = TRUE))
  VariableFeatures(obj) <- hvg
  npcs <- min(30, max(5, ncol(obj) - 1))
  obj <- ScaleData(obj, verbose = FALSE) %>% RunPCA(npcs = npcs, verbose = FALSE)
  if ("orig.ident" %in% colnames(obj[[]]) && length(unique(obj$orig.ident)) > 1) {
    dh <- intersect(seq_len(npcs), dims_h)
    obj <- harmony_emb(obj, "orig.ident", dh, 2)
    red <- "harmony"
  } else {
    red <- "pca"
  }
  use_dims <- intersect(seq_len(npcs), 1:10)
  obj <- FindNeighbors(obj, reduction = red, dims = use_dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = res, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = red, dims = use_dims, verbose = FALSE)
  obj
}

assign_subtype <- function(obj, marker_list, cluster_col = NULL) {
  if (is.null(cluster_col)) {
    snn_cols <- grep("_snn_res", colnames(obj[[]]), value = TRUE)
    cluster_col <- tail(snn_cols, 1)
  }
  Idents(obj) <- obj[[cluster_col]][[1]]
  avg <- AverageExpression(obj, assays = "RNA", layer = "data", verbose = FALSE)$RNA
  res <- list()
  for (cl in levels(Idents(obj))) {
    cn <- ave_col(avg, cl)
    if (is.na(cn)) { res[[cl]] <- list(top = "Unassigned", margin = 0, scores = rep(0, length(marker_list))); next }
    sc <- sapply(marker_list, function(g) {
      gx <- intersect(g, rownames(avg))
      if (length(gx) == 0) return(NA_real_)
      mean(avg[gx, cn], na.rm = TRUE)
    })
    sc[is.na(sc)] <- 0
    o <- order(sc, decreasing = TRUE)
    res[[cl]] <- list(top = names(sc)[o[1]],
                      margin = sc[o[1]] - sc[o[2]],
                      scores = sc)
  }
  res
}

run_subtype_analysis <- function(seu_in, ct_name, marker_list, sub_col, sub_pal, part_label) {
  log_msg(sprintf("================ %s: %s subtypes ================", part_label, ct_name))
  obj <- subset(seu_in, celltype == ct_name)
  if (ncol(obj) < 30) {
    log_msg(sprintf("Too few cells (%d) - skip", ncol(obj))); return(NULL)
  }
  obj <- subcluster_one(obj, ct_name)

  marker_present <- lapply(marker_list, function(g) intersect(g, rownames(obj)))
  marker_present <- marker_present[sapply(marker_present, length) >= 1]
  for (nm in names(marker_present))
    log_msg(sprintf("%s markers available: %d (%s)", nm, length(marker_present[[nm]]),
                    paste(marker_present[[nm]], collapse = ",")))

  ass <- assign_subtype(obj, marker_present)

  if (ct_name == "Astrocytes") {
    for (cl in names(ass)) {
      s <- ass[[cl]]$scores
      h <- s["Homeostatic"]; if(is.na(h)||h<=0) h <- 0.01
      r <- s["Reactive"]
      ratio <- r / h
      if (ratio > 0.5) ass[[cl]]$top <- "Reactive"
      else if (ratio > 0.25) ass[[cl]]$top <- "Intermediate"
      else ass[[cl]]$top <- "Homeostatic"
      s2 <- sort(s, decreasing=TRUE)
      ass[[cl]]$margin <- s2[1] - s2[2]
    }
    log_msg("Using ratio-based assignment for Astrocytes (H/I/R)")
  }

  scores_mat <- do.call(rbind, lapply(ass, function(x) x$scores))
  rownames(scores_mat) <- names(ass)
  margin_vec <- sapply(ass, function(x) x$margin)
  top_vec    <- sapply(ass, function(x) x$top)

  cl_log <- data.frame(
    cluster = names(ass), assigned = top_vec, margin = margin_vec,
    confidence = cut(margin_vec, breaks = c(-Inf, 0.1, 0.3, Inf),
                     labels = c("Low","Medium","High"))
  )
  cl_log <- cbind(cl_log, scores_mat)
  safe_write_csv(cl_log, file.path(OUT_TAB, sprintf("subtype_assignment_%s.csv", ct_name)))

  cluster_col <- grep("_snn_res", colnames(obj[[]]), value = TRUE)
  cluster_col <- tail(cluster_col, 1)
  cl_ids <- as.character(obj[[cluster_col]][[1]])
  obj[[sub_col]] <- unname(top_vec[cl_ids])
  obj[[sub_col]][is.na(obj[[sub_col]])] <- "Unassigned"

  used <- intersect(names(sub_pal), unique(obj[[sub_col]][[1]]))
  pal_use <- sub_pal[used]
  if (!"Unassigned" %in% names(pal_use) && "Unassigned" %in% obj[[sub_col]][[1]])
    pal_use <- c(pal_use, "Unassigned" = "grey75")

  # UMAP
  umap_df <- as.data.frame(Embeddings(obj, "umap"))
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  umap_df$subtype <- obj[[sub_col]][[1]]
  if ("Group" %in% colnames(obj[[]])) umap_df$group <- as.character(obj$Group)

  p_um <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = subtype)) +
    geom_point(size = 0.5, alpha = 0.85) +
    scale_color_manual(values = pal_use) +
    guides(color = guide_legend(override.aes = list(size = 2.5))) +
    labs(title = sprintf("%s subtypes", ct_name), x = "UMAP 1", y = "UMAP 2", color = NULL) +
    publication_theme(7) + theme(aspect.ratio = 1)
  save_fig(p_um, sprintf("Fig_%s_UMAP", ct_name), OUT_FIG, w = 5.2, h = 4.5)
  save_source(umap_df, sprintf("Fig_%s_UMAP", ct_name))

  # DotPlot
  present_subtypes <- unique(obj[[sub_col]][[1]])
  marker_for_dot <- marker_present[intersect(names(marker_present), present_subtypes)]
  marker_for_dot <- rev(marker_for_dot)
  if (length(marker_for_dot) > 0 && length(present_subtypes) > 1) {
    Idents(obj) <- obj[[sub_col]][[1]]
    obj <- SetIdent(obj, value = factor(obj[[sub_col]][[1]], levels = intersect(names(pal_use), present_subtypes)))
    p_dot <- DotPlot(obj, features = marker_for_dot, dot.scale = 4) +
      scale_color_gradientn(colors = c("grey92", "#FDAE61", "#D73027", "#67001F")) +
      scale_size(range = c(0.5, 5)) +
      labs(x = NULL, y = NULL, title = sprintf("%s markers", ct_name)) +
      publication_theme(6.5) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
    total_features <- length(unique(unlist(marker_for_dot)))
    save_fig(p_dot, sprintf("Fig_%s_DotPlot", ct_name), OUT_FIG,
             w = max(5, 0.25 * total_features + 2.5), h = 3.5)
    save_source(p_dot$data, sprintf("Fig_%s_DotPlot", ct_name))

    # Heatmap
    pg <- unique(unlist(marker_present))
    mat_avg <- as.matrix(AverageExpression(obj, assays = "RNA", layer = "data",
                                           features = pg, group.by = sub_col,
                                           verbose = FALSE)$RNA)
    if (nrow(mat_avg) >= 2 && ncol(mat_avg) >= 2) {
      ms <- t(scale(t(mat_avg)))
      ms[ms > 2] <- 2; ms[ms < -2] <- -2
      ms[is.na(ms)] <- 0
      pal_heat <- pal_use
      names(pal_heat) <- gsub("_", "-", names(pal_heat))
      ann_col <- data.frame(Subtype = colnames(ms), row.names = colnames(ms))
      pheatmap(ms, cluster_rows = FALSE, cluster_cols = FALSE,
               color = colorRampPalette(c("#053061","#2166AC","#F7F7F7","#B2182B","#67001F"))(100),
               annotation_col = ann_col,
               annotation_colors = list(Subtype = pal_heat),
               annotation_names_col = FALSE, border_color = NA,
               cellwidth = 40, cellheight = 14,
               filename = file.path(OUT_FIG, sprintf("Fig_%s_Heatmap.pdf", ct_name)),
               width = 6, height = max(4, 0.22 * length(pg) + 1.5))
      log_msg(sprintf("Heatmap saved: Fig_%s_Heatmap.pdf", ct_name))
      save_source(as.data.frame(ms) %>% rownames_to_column("gene"),
                  sprintf("Fig_%s_Heatmap", ct_name))
    }

    # Composition by group
    if ("Group" %in% colnames(obj[[]])) {
      comp <- as.data.frame(obj[[]]) %>%
        dplyr::count(Group, !!sym(sub_col)) %>%
        group_by(Group) %>%
        mutate(prop = n / sum(n)) %>% ungroup()
      colnames(comp)[2] <- "subtype"
      p_comp <- ggplot(comp, aes(Group, prop, fill = subtype)) +
        geom_col(color = "white", linewidth = 0.3) +
        scale_fill_manual(values = pal_use) +
        scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
        labs(x = NULL, y = "Proportion", fill = NULL,
             title = sprintf("%s composition by group", ct_name)) +
        publication_theme(7)
      save_fig(p_comp, sprintf("Fig_%s_Composition", ct_name), OUT_FIG, w = 4, h = 4)
      save_source(comp, sprintf("Fig_%s_Composition", ct_name))

      # AD vs CN per subtype (Propeller FDR: pseudobulk + limma eBayes)
      comp_by_cell <- as.data.frame(obj[[]])[, c("orig.ident", "Group", sub_col), drop=FALSE]
      colnames(comp_by_cell)[3] <- "subtype"
      comp_by_cell$Group <- standardize_group(comp_by_cell$Group)
      sample_prop <- comp_by_cell %>%
        dplyr::count(orig.ident, Group, subtype) %>%
        group_by(orig.ident) %>%
        mutate(prop = n / sum(n)) %>% ungroup()
      
      prop_res <- tryCatch({
        suppressMessages(propeller(clusters = comp_by_cell$subtype,
                                   sample = comp_by_cell$orig.ident,
                                   group = comp_by_cell$Group,
                                   transform = "logit", robust = TRUE, trend = FALSE))
      }, error = function(e) NULL)
      
      subtypes_present <- unique(sample_prop$subtype)
      anno_list <- lapply(subtypes_present, function(st) {
        dat <- sample_prop %>% filter(subtype == st)
        cn_mean <- mean(dat$prop[dat$Group == "CN"]) * 100
        ad_mean <- mean(dat$prop[dat$Group == "AD"]) * 100
        delta <- ad_mean - cn_mean
        if (!is.null(prop_res) && st %in% rownames(prop_res)) {
          label <- sprintf("%+.1f%%\nFDR=%.2f", delta, prop_res[st, "FDR"])
        } else {
          label <- sprintf("%+.1f%%", delta)
        }
        data.frame(subtype = st, label = label, stringsAsFactors = FALSE)
      })
      sig_df <- do.call(rbind, anno_list)
      
      p_sig <- ggplot(sample_prop, aes(Group, prop, fill = Group)) +
        geom_boxplot(width = 0.5, alpha = 0.7, outlier.size = 0.5) +
        geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
        geom_text(data = sig_df, aes(label = label, x = 1.5, y = Inf, vjust = 2),
                  inherit.aes = FALSE, size = 2.8, hjust = 0.5) +
        facet_wrap(~subtype, nrow = 1, scales = "free_y") +
        scale_fill_manual(values = group_palette) +
        labs(x = NULL, y = "Proportion (per sample)",
             title = paste0(ct_name, ": AD vs CN (Propeller FDR)")) +
        publication_theme(6) +
        theme(strip.text = element_text(size = 5, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 5))
      save_fig(p_sig, sprintf("Fig_%s_ADvsCN", ct_name), OUT_FIG, w = 2*length(subtypes_present)+1, h = 3.5)
    }

    # DEGs per subtype
    Idents(obj) <- obj[[sub_col]][[1]]
    if (length(unique(Idents(obj))) >= 2) {
      degs <- tryCatch(
        FindAllMarkers(obj, only.pos = TRUE, logfc.threshold = 0.25,
                       min.pct = 0.1, verbose = FALSE),
        error = function(e) NULL
      )
      if (!is.null(degs) && nrow(degs) > 0) {
        safe_write_csv(degs, file.path(OUT_TAB, sprintf("subtype_DEGs_%s.csv", ct_name)))
      }
    }
  }

  qs::qsave(obj, file.path(OUT_RDS, sprintf("scRNA_%s_annotated.qs", tolower(substr(ct_name,1,5)))),
            preset = "fast")
  log_msg(sprintf("%s subtype distribution:", ct_name))
  print(table(obj[[sub_col]][[1]]))
  obj
}

# ---- Marker definitions ----
ec_markers <- list(
  Pericyte  = c("PDGFRB","RGS5","NOTCH3"),
  SMC       = c("ACTA2","TAGLN","CNN1"),
  Arterial  = c("GJA5","EFNB2","BMX"),
  Capillary = c("MFSD2A","RGCC","CA4"),
  Venous    = c("NR2F2","ACKR1")
)
ec_pal <- c("Pericyte"="#8B4513","SMC"="#FF7F00","Arterial"="#E41A1C",
            "Capillary"="#4DAF4A","Venous"="#377EB8","Unassigned"="grey75")

astro_markers <- list(
  Homeostatic  = c("ALDH1L1","SLC1A2","GLUL"),
  Intermediate = c("VIM","CLU"),
  Reactive = c("GFAP","CD44")
)
astro_pal <- c("Homeostatic"="#66C2A5","Intermediate"="#FFD92F",
               "Reactive"="#FC8D62","Unassigned"="grey75")

micro_markers <- list(
  Homeostatic  = c("P2RY12","TMEM119","CX3CR1"),
  DAM          = c("APOE","LPL","CST7","SPP1"),
  Inflammatory = c("IL1B","CCL2","TNF")
)
micro_pal <- c("Homeostatic"="#66C2A5","DAM"="#FC8D62",
               "Inflammatory"="#E78AC3","Unassigned"="grey75")

# ---- Dispatch ----
log_msg(sprintf("Worker started for %s", CELLTYPE))
log_msg(sprintf("Loading checkpoint: %s", CKPT_QS))

seu_all <- qs::qread(CKPT_QS)
log_msg(sprintf("Loaded %d cells", ncol(seu_all)))

markers <- switch(CELLTYPE,
  "Cerebrovascular cells" = ec_markers,
  "Endothelial" = ec_markers,
  "Astrocytes"  = astro_markers,
  "Microglia"   = micro_markers,
  stop("Unknown celltype: ", CELLTYPE))

sub_col <- switch(CELLTYPE,
  "Cerebrovascular cells" = "vascular_subtype",
  "Endothelial" = "vascular_subtype",
  "Astrocytes"  = "astro_subtype",
  "Microglia"   = "micro_subtype")

sub_pal <- switch(CELLTYPE,
  "Cerebrovascular cells" = ec_pal,
  "Endothelial" = ec_pal,
  "Astrocytes"  = astro_pal,
  "Microglia"   = micro_pal)

log_msg("Running subtype analysis...")
res <- tryCatch(
  run_subtype_analysis(seu_all, CELLTYPE, markers, sub_col, sub_pal, "SUBTYPE"),
  error = function(e) {
    log_msg(sprintf("CRASH: %s", e$message))
    traceback()
    NULL
  })

if (!is.null(res)) {
  out_path <- file.path(OUT_RDS, sprintf("subtype_%s.qs", CELLTYPE))
  qs::qsave(res, out_path, preset = "fast")
  log_msg(sprintf("Saved: %s", out_path))
  log_msg("Worker SUCCESS. Exiting.")
  quit(save = "no", status = 0)
} else {
  log_msg("Worker FAILED. Exiting with error.")
  quit(save = "no", status = 1)
}
