#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
node07 Final: SCTour trajectory analysis for AD-NVU project
===========================================================

Three disease-relevant NVU trajectories:

1. Microglial DAM trajectory:
   Homeostatic microglia -> activated / DAM-like microglia

2. Astrocytic reactive trajectory:
   Homeostatic astrocytes -> intermediate -> reactive / DAA-like astrocytes

3. Capillary endothelial activation trajectory:
   Capillary EC homeostatic-like state -> activated / angiogenic / inflammatory EC state

This script is designed to generate publication-quality main and supplementary
materials suitable for publication-family manuscript preparation.

Author: analysis team
Updated: 2026-05-25
Environment: Python + scanpy + sctour
"""

# ============================================================
# 0. Imports
# ============================================================

import os
import sys
import json
import time
import platform
import warnings
import traceback
import subprocess
from datetime import datetime

import numpy as np
import pandas as pd

import scanpy as sc
import anndata as ad
import sctour as sct

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from scipy import sparse
from scipy.stats import spearmanr, mannwhitneyu, kruskal

warnings.filterwarnings("ignore")


# ============================================================
# 1. Global configuration
# ============================================================

BASE_DIR = r"/path/to/project/results/node07_sctour/sctour"

OUT_DIR = r"/path/to/project/results/node07_sctour/result/node07_SCTour_Final"

os.makedirs(OUT_DIR, exist_ok=True)


# ----------------------------
# Important computational params
# ----------------------------

RANDOM_SEED = 20260525

N_HVG = 2000
MIN_CELLS = 50

# SCTour epoch:
# None = SCTour default.
# For final manuscript, you may set 400-1000 depending on convergence.
SCTOUR_NEPOCH = None

# Dense conversion memory guard.
# SCTour often needs dense input. However, converting large sparse matrices
# blindly can crash the machine. This guard estimates the dense float32 size.
MAX_DENSE_GB = 20

# Whether to draw exploratory vector field on raw UMAP.
# IMPORTANT:
# Latent vector field should be the primary figure.
# Raw UMAP vector field can be visually attractive but may be mathematically
# less rigorous unless properly projected.
PLOT_RAW_UMAP_VECTOR_FIELD = False

# If True, save processed AnnData h5ad for each trajectory.
SAVE_H5AD = True

# Figure DPI
DPI_MAIN = 300
DPI_SUPP = 300

# UMAP / neighbors params for SCTour latent visualization
LATENT_N_NEIGHBORS = 15
LATENT_UMAP_MIN_DIST = 0.15

# Pseudotime bin number for heatmaps and composition analysis
N_PTIME_BINS = 20


# ============================================================
# 2. Matplotlib / seaborn style: publication-like
# ============================================================

def set_publication_style():
    """
    Set a clean, publication-like plotting style.
    """
    matplotlib.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "Arial",
        "font.size": 7,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 6,
        "figure.titlesize": 9,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "savefig.dpi": DPI_MAIN,
        "figure.dpi": DPI_MAIN,
    })

    sns.set_theme(
        context="paper",
        style="white",
        font="Arial",
        rc={
            "axes.linewidth": 0.6,
            "axes.edgecolor": "black",
            "xtick.bottom": True,
            "ytick.left": True,
        }
    )


set_publication_style()


# ============================================================
# 3. Color palettes
# ============================================================

GROUP_COLORS = {
    "CN": "#4C78A8",
    "NC": "#4C78A8",
    "Control": "#4C78A8",
    "CTRL": "#4C78A8",
    "AD": "#E45756",
    "UC": "#E45756",
    "Disease": "#E45756",
    "MCI": "#F2A541",
}

# Okabe-Ito / publication-friendly discrete palette
DEFAULT_CATEGORY_COLORS = [
    "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7",
    "#56B4E9", "#F0E442", "#000000", "#999999", "#882255",
    "#44AA99", "#AA4499", "#117733", "#332288", "#DDCC77",
]

STATE_COLORS = {
    # Microglia
    "Homeostatic": "#3B6EA8",
    "Homeostatic_Micro": "#3B6EA8",
    "DAM": "#C44E52",
    "Disease_Associated": "#C44E52",
    "Activated": "#D55E00",
    "Intermediate": "#E69F00",

    # Astrocyte
    "Reactive": "#C44E52",
    "Reactive_DAA": "#C44E52",
    "DAA": "#C44E52",
    "Astro_Homeostatic": "#3B6EA8",

    # Endo
    "Capillary": "#4C78A8",
    "Activated_Capillary": "#C44E52",
    "Inflammatory": "#D55E00",
    "Angiogenic": "#E69F00",
}


# ============================================================
# 4. Trajectory configuration
# ============================================================

CELL_TYPE_CONFIG = {
    "micro": {
        "subdir": "micro",
        "prefix": "Micro",
        "label": "Microglia",
        "subtype_col": "micro_subtype",
        "cap_only": False,
        "early_markers": ["P2RY12", "CX3CR1", "TMEM119", "SALL1", "GPR34"],
        "late_markers": ["APOE", "SPP1", "CST7", "LPL", "ITGAX", "TYROBP", "TREM2"],
        "title": "Microglial DAM trajectory",
    },

    "astro": {
        "subdir": "astro",
        "prefix": "Astro",
        "label": "Astrocyte",
        "subtype_col": "astro_subtype",
        "cap_only": False,
        "early_markers": ["AQP4", "SLC1A2", "SLC1A3", "ALDH1L1", "GJA1"],
        "late_markers": ["GFAP", "CD44", "VIM", "C3", "SERPINA3", "LCN2"],
        "title": "Astrocytic reactive trajectory",
    },

    "endo": {
        "subdir": "endo",
        "prefix": "Endo",
        "label": "Capillary endothelial cell",
        "subtype_col": "CapEC_activation_state",
        "endo_subtype_col": "endo_subtype",
        "cap_only": True,
        "capillary_label": "Capillary",
        "early_markers": ["CLDN5", "FLT1", "KDR", "PECAM1", "RAMP2"],
        "late_markers": ["ANGPT2", "VCAM1", "ICAM1", "SELE", "VWF", "VEGFA"],
        "title": "Capillary endothelial activation trajectory",
    },
}


# ============================================================
# 5. Utility functions
# ============================================================

def make_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def save_figure_all(fig, base_path, dpi=300):
    """
    Save figure as PNG, PDF, SVG, and TIFF for manuscript editing.
    """
    for ext in ["png", "pdf", "svg", "tiff"]:
        out = f"{base_path}.{ext}"
        try:
            kwargs = {"dpi": dpi, "bbox_inches": "tight"}
            if ext == "tiff":
                kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
            fig.savefig(out, **kwargs)
        except Exception as e:
            print(f"[WARN] Failed to save {out}: {e}")
    print(f"[FIG] {base_path}.png/pdf/svg/tiff")


def harmonize_group_labels(values):
    """
    Use final manuscript group names while retaining compatibility with earlier exports.
    """
    return pd.Series(values).astype(str).replace({
        "NC": "CN",
        "UC": "AD",
        "Control": "CN",
        "CTRL": "CN",
        "Disease": "AD",
    }).values


def save_table(df, path, index=False):
    df.to_csv(path, index=index)
    print(f"[TABLE] {path}")


def write_json(obj, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)
    print(f"[JSON] {path}")


def log_versions(out_dir):
    version_file = os.path.join(out_dir, "software_versions.txt")
    with open(version_file, "w", encoding="utf-8") as f:
        f.write(f"Run time: {now()}\n")
        f.write(f"Python: {sys.version}\n")
        f.write(f"Platform: {platform.platform()}\n")
        f.write(f"numpy: {np.__version__}\n")
        f.write(f"pandas: {pd.__version__}\n")
        f.write(f"scanpy: {sc.__version__}\n")
        f.write(f"anndata: {ad.__version__}\n")
        f.write(f"matplotlib: {matplotlib.__version__}\n")
        f.write(f"seaborn: {sns.__version__}\n")
        try:
            f.write(f"sctour: {sct.__version__}\n")
        except Exception:
            f.write("sctour: version unknown\n")
    print(f"[INFO] Software versions saved: {version_file}")


def safe_obs_category(adata, col):
    """
    Convert obs column to string category safely.
    """
    if col in adata.obs.columns:
        adata.obs[col] = adata.obs[col].astype(str).astype("category")


def find_existing_col(adata, candidates):
    """
    Return first existing column from candidate list.
    """
    for c in candidates:
        if c in adata.obs.columns:
            return c
    return None


def make_category_palette(values):
    """
    Generate stable category palette.
    """
    values = [str(v) for v in values]
    palette = {}
    for i, v in enumerate(values):
        if v in STATE_COLORS:
            palette[v] = STATE_COLORS[v]
        elif v in GROUP_COLORS:
            palette[v] = GROUP_COLORS[v]
        else:
            palette[v] = DEFAULT_CATEGORY_COLORS[i % len(DEFAULT_CATEGORY_COLORS)]
    return palette


def estimate_dense_gb(n_cells, n_genes, dtype=np.float32):
    bytes_per = np.dtype(dtype).itemsize
    return n_cells * n_genes * bytes_per / 1024**3


def get_expr_vector(adata, gene):
    """
    Return expression vector for one gene from adata.X.
    Assumes gene is in adata.var_names.
    """
    x = adata[:, gene].X
    if sparse.issparse(x):
        x = x.toarray()
    x = np.asarray(x).reshape(-1)
    return x


def minmax_scale(x):
    x = np.asarray(x, dtype=float)
    xmin = np.nanmin(x)
    xmax = np.nanmax(x)
    if not np.isfinite(xmin) or not np.isfinite(xmax) or xmax <= xmin:
        return np.zeros_like(x)
    return (x - xmin) / (xmax - xmin)


def sanitize_filename(x):
    x = str(x)
    bad = ["/", "\\", ":", "*", "?", "\"", "<", ">", "|", " "]
    for b in bad:
        x = x.replace(b, "_")
    return x


# ============================================================
# 6. Data loading
# ============================================================

def load_cell_data(subdir, prefix):
    """
    Load counts matrix + genes + barcodes + metadata + UMAP.

    Expected files:
      {prefix}_counts.mtx      genes x cells from R export
      {prefix}_genes.tsv
      {prefix}_barcodes.tsv
      {prefix}_metadata.csv
      {prefix}_umap.csv
    """
    cell_dir = os.path.join(BASE_DIR, subdir)

    mtx_path = os.path.join(cell_dir, f"{prefix}_counts.mtx")
    gene_path = os.path.join(cell_dir, f"{prefix}_genes.tsv")
    barcode_path = os.path.join(cell_dir, f"{prefix}_barcodes.tsv")
    meta_path = os.path.join(cell_dir, f"{prefix}_metadata.csv")
    umap_path = os.path.join(cell_dir, f"{prefix}_umap.csv")

    for p in [mtx_path, gene_path, barcode_path, meta_path, umap_path]:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing required file: {p}")

    print(f"[LOAD] {prefix}: reading matrix...")
    adata = sc.read_mtx(mtx_path).T

    genes = pd.read_csv(gene_path, header=None, sep="\t")
    barcodes = pd.read_csv(barcode_path, header=None, sep="\t")

    gene_names = genes.iloc[:, 0].astype(str).values
    barcode_names = barcodes.iloc[:, 0].astype(str).values

    adata.var_names = gene_names
    adata.obs_names = barcode_names
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    metadata = pd.read_csv(meta_path)

    # Robust metadata alignment
    possible_cell_cols = ["cell_id", "barcode", "cell", "Cell", "CellID", "orig_barcode"]
    cell_col = None
    for c in possible_cell_cols:
        if c in metadata.columns:
            cell_col = c
            break

    if cell_col is not None:
        metadata[cell_col] = metadata[cell_col].astype(str)
        metadata = metadata.set_index(cell_col)
        metadata = metadata.reindex(adata.obs_names)
        if metadata.isnull().all(axis=None):
            raise ValueError(
                f"Metadata alignment failed using column {cell_col}. "
                f"Please check barcode names."
            )
    else:
        if metadata.shape[0] != adata.n_obs:
            raise ValueError(
                f"Metadata rows ({metadata.shape[0]}) do not match cells ({adata.n_obs}), "
                f"and no cell_id/barcode column was found."
            )
        metadata.index = adata.obs_names

    adata.obs = metadata.copy()
    for group_col in ["Group", "group"]:
        if group_col in adata.obs.columns:
            adata.obs[group_col] = harmonize_group_labels(adata.obs[group_col])

    # UMAP alignment
    umap_df = pd.read_csv(umap_path)

    if {"umap_1", "umap_2"}.issubset(set(umap_df.columns)):
        umap_cols = ["umap_1", "umap_2"]
    elif {"UMAP_1", "UMAP_2"}.issubset(set(umap_df.columns)):
        umap_cols = ["UMAP_1", "UMAP_2"]
    else:
        raise ValueError(f"UMAP file must contain umap_1/umap_2 or UMAP_1/UMAP_2: {umap_path}")

    umap_cell_col = None
    for c in possible_cell_cols:
        if c in umap_df.columns:
            umap_cell_col = c
            break

    if umap_cell_col is not None:
        umap_df[umap_cell_col] = umap_df[umap_cell_col].astype(str)
        umap_df = umap_df.set_index(umap_cell_col)
        umap_df = umap_df.reindex(adata.obs_names)
    else:
        if umap_df.shape[0] != adata.n_obs:
            raise ValueError(
                f"UMAP rows ({umap_df.shape[0]}) do not match cells ({adata.n_obs}), "
                f"and no cell_id/barcode column found."
            )
        umap_df.index = adata.obs_names

    adata.obsm["X_umap_raw"] = umap_df[umap_cols].values.astype(np.float32)
    adata.obsm["X_umap"] = adata.obsm["X_umap_raw"].copy()

    # Store raw counts layer
    adata.layers["counts"] = adata.X.copy()

    print(f"[LOAD] {prefix}: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    return adata


# ============================================================
# 7. Preprocessing
# ============================================================

def preprocess_for_sctour(adata, label, out_dir, n_hvg=N_HVG):
    """
    Preprocess for SCTour:
    - QC metrics
    - HVG selection on count matrix
    - dense conversion with memory guard
    - float32 conversion
    """
    print(f"[PREPROCESS] {label}")

    if adata.n_obs < MIN_CELLS:
        raise ValueError(f"{label}: too few cells ({adata.n_obs}) for trajectory analysis.")

    # Basic QC
    try:
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    except Exception as e:
        print(f"[WARN] QC metrics failed: {e}")

    qc_cols = [c for c in ["n_genes_by_counts", "total_counts", "pct_counts_mt"] if c in adata.obs.columns]
    if len(qc_cols) > 0:
        qc_df = adata.obs[qc_cols].copy()
        qc_df.insert(0, "cell_id", adata.obs_names)
        save_table(qc_df, os.path.join(out_dir, f"{label}_qc_metrics.csv"))

    # HVG selection
    print(f"[PREPROCESS] {label}: HVG selection n_top_genes={n_hvg}")

    try:
        sc.pp.highly_variable_genes(
            adata,
            flavor="seurat_v3",
            n_top_genes=n_hvg,
            layer="counts",
            subset=False,
        )
    except Exception as e:
        print(f"[WARN] seurat_v3 HVG failed, fallback to cell_ranger. Reason: {e}")
        sc.pp.highly_variable_genes(
            adata,
            flavor="cell_ranger",
            n_top_genes=n_hvg,
            subset=False,
        )

    hvg_df = adata.var[["highly_variable"]].copy()
    if "highly_variable_rank" in adata.var.columns:
        hvg_df["highly_variable_rank"] = adata.var["highly_variable_rank"]
    hvg_df.insert(0, "gene", adata.var_names)
    save_table(hvg_df, os.path.join(out_dir, f"{label}_hvg_table.csv"))

    adata = adata[:, adata.var["highly_variable"].values].copy()

    # SCTour usually uses raw counts; preserve counts subset
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()

    dense_gb = estimate_dense_gb(adata.n_obs, adata.n_vars, np.float32)
    print(f"[PREPROCESS] {label}: dense float32 estimate = {dense_gb:.2f} GB")

    if dense_gb > MAX_DENSE_GB:
        raise MemoryError(
            f"{label}: dense conversion would require ~{dense_gb:.2f} GB, "
            f"exceeding MAX_DENSE_GB={MAX_DENSE_GB}. "
            f"Please reduce N_HVG, subset cells, or increase memory limit."
        )

    if sparse.issparse(adata.X):
        print(f"[PREPROCESS] {label}: sparse -> dense float32")
        adata.X = adata.X.toarray().astype(np.float32)
    else:
        adata.X = np.asarray(adata.X).astype(np.float32)

    return adata


# ============================================================
# 8. Pseudotime orientation
# ============================================================

def orient_pseudotime_by_markers(adata, label, early_markers, late_markers, out_dir):
    """
    Orient pseudotime using early and late marker signatures.

    Logic:
      - Compute average late marker expression signature.
      - Compute average early marker expression signature.
      - Disease progression score = late_score - early_score.
      - If Spearman(ptime_raw, progression_score) < 0, flip pseudotime.

    This is more robust than checking one or two markers.
    """

    p_raw = np.asarray(adata.obs["ptime_raw"].values, dtype=float)
    p_raw = minmax_scale(p_raw)

    early_avail = [g for g in early_markers if g in adata.var_names]
    late_avail = [g for g in late_markers if g in adata.var_names]

    rows = []

    for g in early_avail + late_avail:
        expr = get_expr_vector(adata, g)
        rho, pval = spearmanr(p_raw, expr)
        rows.append({
            "trajectory": label,
            "gene": g,
            "marker_class": "early" if g in early_avail else "late",
            "spearman_rho_with_ptime_raw": rho,
            "p_value": pval,
        })

    corr_df = pd.DataFrame(rows)
    save_table(corr_df, os.path.join(out_dir, f"{label}_marker_pseudotime_correlation_raw.csv"))

    if len(late_avail) == 0:
        print(f"[ORIENT] {label}: no late markers found. Keeping raw pseudotime.")
        adata.obs["ptime"] = p_raw
        orientation_info = {
            "trajectory": label,
            "flipped": False,
            "reason": "no_late_markers_found",
            "early_markers_used": early_avail,
            "late_markers_used": late_avail,
            "spearman_progression_score": None,
        }
        return adata, orientation_info, corr_df

    late_mat = np.vstack([get_expr_vector(adata, g) for g in late_avail]).T
    late_score = np.nanmean(late_mat, axis=1)

    if len(early_avail) > 0:
        early_mat = np.vstack([get_expr_vector(adata, g) for g in early_avail]).T
        early_score = np.nanmean(early_mat, axis=1)
    else:
        early_score = np.zeros_like(late_score)

    progression_score = late_score - early_score
    rho, pval = spearmanr(p_raw, progression_score)

    if np.isfinite(rho) and rho < 0:
        print(f"[ORIENT] {label}: flipping pseudotime. Spearman rho={rho:.3f}, p={pval:.2e}")
        adata.obs["ptime"] = 1.0 - p_raw
        flipped = True
    else:
        print(f"[ORIENT] {label}: keeping raw pseudotime. Spearman rho={rho:.3f}, p={pval:.2e}")
        adata.obs["ptime"] = p_raw
        flipped = False

    adata.obs["late_marker_score"] = minmax_scale(late_score)
    adata.obs["early_marker_score"] = minmax_scale(early_score)
    adata.obs["progression_marker_score"] = minmax_scale(progression_score)

    orientation_info = {
        "trajectory": label,
        "flipped": flipped,
        "reason": "marker_signature_orientation",
        "early_markers_used": early_avail,
        "late_markers_used": late_avail,
        "spearman_progression_score": float(rho) if np.isfinite(rho) else None,
        "spearman_p_value": float(pval) if np.isfinite(pval) else None,
    }

    write_json(orientation_info, os.path.join(out_dir, f"{label}_pseudotime_orientation.json"))

    return adata, orientation_info, corr_df


# ============================================================
# 9. SCTour training
# ============================================================

def train_sctour(adata, label, config, out_dir, nepoch=SCTOUR_NEPOCH):
    """
    Train SCTour and orient pseudotime.
    """
    print(f"[SCTOUR] Training {label}: {adata.n_obs:,} cells x {adata.n_vars:,} HVGs")

    np.random.seed(RANDOM_SEED)

    t0 = time.time()

    trainer = sct.train.Trainer(
        adata,
        loss_mode="nb",
        alpha_recon_lec=0.5,
        alpha_recon_lode=0.5,
        nepoch=nepoch,
    )

    trainer.train()

    elapsed_min = (time.time() - t0) / 60
    print(f"[SCTOUR] {label}: training finished in {elapsed_min:.1f} min")

    ptime_raw = trainer.get_time()
    ptime_raw = minmax_scale(ptime_raw)

    if np.any(~np.isfinite(ptime_raw)):
        raise ValueError(f"{label}: SCTour returned non-finite pseudotime.")

    adata.obs["ptime_raw"] = ptime_raw

    adata, orientation_info, corr_df = orient_pseudotime_by_markers(
        adata=adata,
        label=label,
        early_markers=config["early_markers"],
        late_markers=config["late_markers"],
        out_dir=out_dir,
    )

    adata.obs["ptime"] = minmax_scale(adata.obs["ptime"].values)

    train_info = {
        "trajectory": label,
        "n_cells": int(adata.n_obs),
        "n_hvg": int(adata.n_vars),
        "nepoch": nepoch,
        "elapsed_min": elapsed_min,
        "random_seed": RANDOM_SEED,
        "sctour_loss_mode": "nb",
        "alpha_recon_lec": 0.5,
        "alpha_recon_lode": 0.5,
        "orientation": orientation_info,
    }

    write_json(train_info, os.path.join(out_dir, f"{label}_sctour_training_info.json"))

    return adata, trainer


# ============================================================
# 10. Latent space and vector field
# ============================================================

def compute_latent_and_vector_field(adata, trainer, label, out_dir):
    """
    Compute SCTour latent space and vector field.
    """
    print(f"[LATENT] {label}")

    mix_zs, zs, pred_zs = trainer.get_latentsp(alpha_z=0.5, alpha_predz=0.5)

    adata.obsm["X_TNODE"] = np.asarray(mix_zs, dtype=np.float32)
    adata.obsm["X_TNODE_z"] = np.asarray(zs, dtype=np.float32)
    adata.obsm["X_TNODE_predz"] = np.asarray(pred_zs, dtype=np.float32)

    vf = trainer.get_vector_field(
        adata.obs["ptime"].values,
        adata.obsm["X_TNODE"]
    )

    adata.obsm["X_VF"] = np.asarray(vf, dtype=np.float32)

    # Rebuild UMAP from latent representation without reordering cells
    sc.pp.neighbors(adata, use_rep="X_TNODE", n_neighbors=LATENT_N_NEIGHBORS)
    sc.tl.umap(adata, min_dist=LATENT_UMAP_MIN_DIST, random_state=RANDOM_SEED)
    adata.obsm["X_umap_sctour"] = adata.obsm["X_umap"].copy()

    # Restore default X_umap to SCTour latent UMAP for scanpy plotting;
    # raw UMAP is preserved as X_umap_raw.
    adata.obsm["X_umap"] = adata.obsm["X_umap_sctour"].copy()

    # Save latent and vector field tables
    latent_df = pd.DataFrame(
        adata.obsm["X_TNODE"],
        index=adata.obs_names,
        columns=[f"latent_{i + 1}" for i in range(adata.obsm["X_TNODE"].shape[1])]
    )
    latent_df.insert(0, "cell_id", adata.obs_names)
    save_table(latent_df, os.path.join(out_dir, f"{label}_latent_coordinates.csv"))

    vf_df = pd.DataFrame(
        adata.obsm["X_VF"],
        index=adata.obs_names,
        columns=[f"vf_{i + 1}" for i in range(adata.obsm["X_VF"].shape[1])]
    )
    vf_df.insert(0, "cell_id", adata.obs_names)
    save_table(vf_df, os.path.join(out_dir, f"{label}_vector_field.csv"))

    return adata


# ============================================================
# 11. Save core outputs
# ============================================================

def save_core_outputs(adata, label, out_dir):
    """
    Save per-cell pseudotime and metadata.
    """
    pt_df = pd.DataFrame({
        "cell_id": adata.obs_names,
        "pseudotime": adata.obs["ptime"].values,
        "pseudotime_raw": adata.obs["ptime_raw"].values,
    })

    extra_cols = [
        "celltype2",
        "Group",
        "orig.ident",
        "sample_id",
        "donor",
        "donor_id",
        "late_marker_score",
        "early_marker_score",
        "progression_marker_score",
    ]

    for c in extra_cols:
        if c in adata.obs.columns:
            pt_df[c] = adata.obs[c].astype(str).values if str(adata.obs[c].dtype) == "category" else adata.obs[c].values

    save_table(pt_df, os.path.join(out_dir, f"{label}_cell_pseudotime.csv"))

    if SAVE_H5AD:
        h5ad_path = os.path.join(out_dir, f"{label}_sctour_final.h5ad")
        adata.write_h5ad(h5ad_path, compression="gzip")
        print(f"[H5AD] {h5ad_path}")


# ============================================================
# 12. Plot helpers
# ============================================================

def scanpy_umap_from_key(adata, key, color, ax, title, cmap=None, palette=None,
                         legend_loc="right margin", size=None, alpha=0.8):
    """
    Temporarily set adata.obsm['X_umap'] to a selected embedding key and call sc.pl.umap.
    """
    old_umap = adata.obsm["X_umap"].copy()
    adata.obsm["X_umap"] = adata.obsm[key].copy()

    sc.pl.umap(
        adata,
        color=color,
        ax=ax,
        show=False,
        frameon=False,
        title=title,
        cmap=cmap,
        palette=palette,
        legend_loc=legend_loc,
        size=size,
        alpha=alpha,
    )

    adata.obsm["X_umap"] = old_umap


def add_panel_label(ax, label):
    ax.text(
        -0.08, 1.05, label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
        ha="right"
    )


# ============================================================
# 13. Main figure panels
# ============================================================

def plot_raw_umap_panels(adata, label, config, out_dir):
    """
    Raw UMAP colored by subtype, group, and sample.
    """
    print(f"[PLOT] {label}: raw UMAP panels")

    group_col = find_existing_col(adata, ["Group", "group", "Diagnosis", "diagnosis"])
    sample_col = find_existing_col(adata, ["orig.ident", "sample_id", "Sample", "sample", "donor", "donor_id"])

    cols = ["celltype2"]
    titles = ["Cell state"]

    if group_col:
        cols.append(group_col)
        titles.append("Group")

    if sample_col:
        cols.append(sample_col)
        titles.append("Sample / donor")

    n = len(cols)
    fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 3.0))
    if n == 1:
        axes = [axes]

    for i, (c, title) in enumerate(zip(cols, titles)):
        if c == "celltype2":
            vals = adata.obs[c].astype(str).unique().tolist()
            palette = make_category_palette(vals)
        elif group_col and c == group_col:
            vals = adata.obs[c].astype(str).unique().tolist()
            palette = make_category_palette(vals)
        else:
            palette = None

        scanpy_umap_from_key(
            adata,
            key="X_umap_raw",
            color=c,
            ax=axes[i],
            title=title,
            palette=palette,
            legend_loc="right margin" if c != "celltype2" else "on data",
            size=8,
            alpha=0.75,
        )
        add_panel_label(axes[i], chr(65 + i))

    plt.tight_layout()
    save_figure_all(fig, os.path.join(out_dir, f"{label}_Fig_raw_umap_panels"), dpi=DPI_MAIN)
    plt.close(fig)


def plot_sctour_main_panels(adata, label, config, out_dir):
    """
    Main SCTour figure:
      A raw UMAP state
      B latent UMAP state
      C latent UMAP pseudotime
      D latent vector field
    """
    print(f"[PLOT] {label}: SCTour main panels")

    fig, axes = plt.subplots(1, 4, figsize=(12.0, 3.0))

    state_values = adata.obs["celltype2"].astype(str).unique().tolist()
    state_palette = make_category_palette(state_values)

    scanpy_umap_from_key(
        adata,
        key="X_umap_raw",
        color="celltype2",
        ax=axes[0],
        title="Original UMAP",
        palette=state_palette,
        legend_loc="on data",
        size=8,
        alpha=0.75,
    )
    add_panel_label(axes[0], "A")

    scanpy_umap_from_key(
        adata,
        key="X_umap_sctour",
        color="celltype2",
        ax=axes[1],
        title="SCTour latent UMAP",
        palette=state_palette,
        legend_loc="on data",
        size=8,
        alpha=0.75,
    )
    add_panel_label(axes[1], "B")

    scanpy_umap_from_key(
        adata,
        key="X_umap_sctour",
        color="ptime",
        ax=axes[2],
        title="Pseudotime",
        cmap="viridis",
        legend_loc="right margin",
        size=8,
        alpha=0.85,
    )
    add_panel_label(axes[2], "C")

    try:
        # Primary rigorous vector field visualization in SCTour latent space.
        sct.vf.plot_vector_field(
            adata,
            zs_key="X_TNODE",
            vf_key="X_VF",
            use_rep_neigh="X_TNODE",
            color="celltype2",
            show=False,
            ax=axes[3],
            legend_loc="none",
            frameon=False,
            size=20,
            alpha=0.25,
        )
        axes[3].set_title("Latent vector field")
    except Exception as e:
        axes[3].text(0.5, 0.5, f"Vector field plot failed:\n{e}",
                     ha="center", va="center", fontsize=6)
        axes[3].axis("off")

    add_panel_label(axes[3], "D")

    fig.suptitle(config["title"], y=1.04, fontsize=10)
    plt.tight_layout()
    save_figure_all(fig, os.path.join(out_dir, f"{label}_Fig_sctour_main_panels"), dpi=DPI_MAIN)
    plt.close(fig)


def plot_pseudotime_raw_and_latent(adata, label, out_dir):
    """
    Pseudotime on raw UMAP and SCTour UMAP.
    """
    print(f"[PLOT] {label}: pseudotime raw vs latent")

    fig, axes = plt.subplots(1, 2, figsize=(6.2, 3.0))

    scanpy_umap_from_key(
        adata,
        key="X_umap_raw",
        color="ptime",
        ax=axes[0],
        title="Original UMAP",
        cmap="viridis",
        size=8,
        alpha=0.85,
    )
    add_panel_label(axes[0], "A")

    scanpy_umap_from_key(
        adata,
        key="X_umap_sctour",
        color="ptime",
        ax=axes[1],
        title="SCTour latent UMAP",
        cmap="viridis",
        size=8,
        alpha=0.85,
    )
    add_panel_label(axes[1], "B")

    plt.tight_layout()
    save_figure_all(fig, os.path.join(out_dir, f"{label}_Supp_pseudotime_raw_vs_latent"), dpi=DPI_SUPP)
    plt.close(fig)


# ============================================================
# 14. Statistical and supplementary plots
# ============================================================

def plot_pseudotime_by_group(adata, label, out_dir):
    """
    Violin/boxplot and density of pseudotime by group.
    """
    group_col = find_existing_col(adata, ["Group", "group", "Diagnosis", "diagnosis"])
    if group_col is None:
        print(f"[SKIP] {label}: no group column found")
        return

    df = pd.DataFrame({
        "pseudotime": adata.obs["ptime"].values,
        "group": adata.obs[group_col].astype(str).values,
    })

    groups = sorted(df["group"].unique().tolist())
    if len(groups) < 2:
        print(f"[SKIP] {label}: fewer than two groups")
        return

    palette = make_category_palette(groups)

    # Violin + box
    fig, ax = plt.subplots(figsize=(3.0, 3.0))
    sns.violinplot(
        data=df,
        x="group",
        y="pseudotime",
        palette=palette,
        inner=None,
        cut=0,
        linewidth=0.6,
        ax=ax,
    )
    sns.boxplot(
        data=df,
        x="group",
        y="pseudotime",
        width=0.25,
        showcaps=True,
        boxprops={"facecolor": "white", "edgecolor": "black", "linewidth": 0.6},
        whiskerprops={"linewidth": 0.6},
        medianprops={"color": "black", "linewidth": 0.8},
        showfliers=False,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Pseudotime")
    ax.set_title("Pseudotime by group")
    sns.despine(ax=ax)
    plt.tight_layout()
    save_figure_all(fig, os.path.join(out_dir, f"{label}_Supp_pseudotime_by_group_violin"), dpi=DPI_SUPP)
    plt.close(fig)

    # Density
    fig, ax = plt.subplots(figsize=(3.2, 2.6))
    for g in groups:
        sns.kdeplot(
            data=df[df["group"] == g],
            x="pseudotime",
            fill=True,
            alpha=0.25,
            linewidth=1.0,
            color=palette[g],
            label=g,
            ax=ax,
        )
    ax.set_xlabel("Pseudotime")
    ax.set_ylabel("Density")
    ax.set_title("Pseudotime density")
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    plt.tight_layout()
    save_figure_all(fig, os.path.join(out_dir, f"{label}_Supp_pseudotime_by_group_density"), dpi=DPI_SUPP)
    plt.close(fig)

    # Statistics
    stat_rows = []
    if len(groups) == 2:
        x = df.loc[df["group"] == groups[0], "pseudotime"].values
        y = df.loc[df["group"] == groups[1], "pseudotime"].values
        stat, p = mannwhitneyu(x, y, alternative="two-sided")
        stat_rows.append({
            "test": "Mann-Whitney U",
            "group1": groups[0],
            "group2": groups[1],
            "statistic": stat,
            "p_value": p,
            "n_group1": len(x),
            "n_group2": len(y),
            "median_group1": np.median(x),
            "median_group2": np.median(y),
        })
    else:
        vals = [df.loc[df["group"] == g, "pseudotime"].values for g in groups]
        stat, p = kruskal(*vals)
        stat_rows.append({
            "test": "Kruskal-Wallis",
            "groups": ";".join(groups),
            "statistic": stat,
            "p_value": p,
        })

    save_table(pd.DataFrame(stat_rows), os.path.join(out_dir, f"{label}_pseudotime_by_group_statistics.csv"))


def plot_pseudotime_by_state(adata, label, out_dir):
    """
    Pseudotime distribution across annotated states/subtypes.
    """
    if "celltype2" not in adata.obs.columns:
        return

    df = pd.DataFrame({
        "pseudotime": adata.obs["ptime"].values,
        "state": adata.obs["celltype2"].astype(str).values,
    })

    order = df.groupby("state")["pseudotime"].median().sort_values().index.tolist()
    palette = make_category_palette(order)

    fig, ax = plt.subplots(figsize=(max(3.2, 0.35 * len(order)), 3.0))

    sns.boxplot(
        data=df,
        x="state",
        y="pseudotime",
        order=order,
        palette=palette,
        showfliers=False,
        linewidth=0.6,
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel("Pseudotime")
    ax.set_title("Pseudotime by cell state")
    ax.tick_params(axis="x", rotation=45)
    sns.despine(ax=ax)
    plt.tight_layout()
    save_figure_all(fig, os.path.join(out_dir, f"{label}_Supp_pseudotime_by_state"), dpi=DPI_SUPP)
    plt.close(fig)

    summary = df.groupby("state").agg(
        n_cells=("pseudotime", "size"),
        mean_pseudotime=("pseudotime", "mean"),
        median_pseudotime=("pseudotime", "median"),
        sd_pseudotime=("pseudotime", "std"),
    ).reset_index()

    save_table(summary, os.path.join(out_dir, f"{label}_state_pseudotime_summary.csv"))


def plot_marker_trends(adata, label, config, out_dir):
    """
    Plot marker expression trends along pseudotime.
    """
    markers = []
    for g in config["early_markers"] + config["late_markers"]:
        if g in adata.var_names and g not in markers:
            markers.append(g)

    if len(markers) == 0:
        print(f"[SKIP] {label}: no configured markers found")
        return

    print(f"[PLOT] {label}: marker trends, {len(markers)} markers")

    df_list = []
    ptime = adata.obs["ptime"].values

    for g in markers:
        expr = get_expr_vector(adata, g)
        df_list.append(pd.DataFrame({
            "pseudotime": ptime,
            "expression": expr,
            "gene": g,
            "marker_class": "early" if g in config["early_markers"] else "late",
        }))

    df = pd.concat(df_list, axis=0, ignore_index=True)

    ncol = 3
    nrow = int(np.ceil(len(markers) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.2 * nrow), sharex=True)
    axes = np.asarray(axes).reshape(-1)

    for i, g in enumerate(markers):
        ax = axes[i]
        sub = df[df["gene"] == g]
        color = "#3B6EA8" if g in config["early_markers"] else "#C44E52"

        # Sample for scatter to keep plots light
        if sub.shape[0] > 5000:
            sub_scatter = sub.sample(5000, random_state=RANDOM_SEED)
        else:
            sub_scatter = sub

        ax.scatter(
            sub_scatter["pseudotime"],
            sub_scatter["expression"],
            s=2,
            alpha=0.12,
            color=color,
            linewidths=0,
        )

        sns.regplot(
            data=sub,
            x="pseudotime",
            y="expression",
            scatter=False,
            lowess=True,
            color="black",
            line_kws={"linewidth": 1.0},
            ax=ax,
        )

        ax.set_title(g)
        ax.set_xlabel("Pseudotime")
        ax.set_ylabel("Expression")
        sns.despine(ax=ax)

    for j in range(len(markers), len(axes)):
        axes[j].axis("off")

    plt.tight_layout()
    save_figure_all(fig, os.path.join(out_dir, f"{label}_Supp_marker_trends_along_pseudotime"), dpi=DPI_SUPP)
    plt.close(fig)


def plot_marker_heatmap_binned(adata, label, config, out_dir):
    """
    Heatmap of binned marker expression along pseudotime.
    """
    markers = []
    for g in config["early_markers"] + config["late_markers"]:
        if g in adata.var_names and g not in markers:
            markers.append(g)

    if len(markers) == 0:
        return

    ptime = adata.obs["ptime"].values
    bins = pd.cut(ptime, bins=N_PTIME_BINS, labels=False, include_lowest=True)

    heat = []
    for g in markers:
        expr = get_expr_vector(adata, g)
        tmp = pd.DataFrame({"bin": bins, "expr": expr})
        y = tmp.groupby("bin")["expr"].mean().reindex(range(N_PTIME_BINS)).values
        y = minmax_scale(y)
        heat.append(y)

    heat = np.vstack(heat)

    fig, ax = plt.subplots(figsize=(4.2, max(2.2, 0.22 * len(markers))))

    sns.heatmap(
        heat,
        cmap="viridis",
        xticklabels=[str(i + 1) for i in range(N_PTIME_BINS)],
        yticklabels=markers,
        cbar_kws={"label": "Scaled expression"},
        ax=ax,
    )

    ax.set_xlabel("Pseudotime bin")
    ax.set_ylabel("")
    ax.set_title("Marker dynamics along pseudotime")
    plt.tight_layout()
    save_figure_all(fig, os.path.join(out_dir, f"{label}_Supp_marker_binned_heatmap"), dpi=DPI_SUPP)
    plt.close(fig)

    heat_df = pd.DataFrame(
        heat,
        index=markers,
        columns=[f"ptime_bin_{i + 1}" for i in range(N_PTIME_BINS)]
    )
    heat_df.insert(0, "gene", heat_df.index)
    save_table(heat_df, os.path.join(out_dir, f"{label}_marker_binned_expression.csv"), index=False)


def compute_pseudotime_bin_composition(adata, label, out_dir):
    """
    Composition of cell states across pseudotime bins.
    """
    if "celltype2" not in adata.obs.columns:
        return

    df = pd.DataFrame({
        "ptime": adata.obs["ptime"].values,
        "state": adata.obs["celltype2"].astype(str).values,
    })

    df["ptime_bin"] = pd.cut(df["ptime"], bins=N_PTIME_BINS, labels=False, include_lowest=True) + 1

    comp = df.groupby(["ptime_bin", "state"]).size().reset_index(name="n")
    comp["fraction"] = comp.groupby("ptime_bin")["n"].transform(lambda x: x / x.sum())

    save_table(comp, os.path.join(out_dir, f"{label}_pseudotime_bin_state_composition.csv"))

    states = comp["state"].unique().tolist()
    palette = make_category_palette(states)

    pivot = comp.pivot(index="ptime_bin", columns="state", values="fraction").fillna(0)
    pivot = pivot[states]

    fig, ax = plt.subplots(figsize=(4.5, 2.8))

    bottom = np.zeros(pivot.shape[0])
    x = pivot.index.values

    for state in states:
        ax.bar(
            x,
            pivot[state].values,
            bottom=bottom,
            color=palette[state],
            label=state,
            width=0.9,
            linewidth=0,
        )
        bottom += pivot[state].values

    ax.set_xlabel("Pseudotime bin")
    ax.set_ylabel("Fraction")
    ax.set_title("State composition along pseudotime")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    sns.despine(ax=ax)
    plt.tight_layout()
    save_figure_all(fig, os.path.join(out_dir, f"{label}_Supp_state_composition_along_pseudotime"), dpi=DPI_SUPP)
    plt.close(fig)


def compute_sample_level_summary(adata, label, out_dir):
    """
    Sample/donor-level pseudotime summary.
    This is important because cell-level distributions can inflate significance.
    """
    sample_col = find_existing_col(adata, ["orig.ident", "sample_id", "Sample", "sample", "donor", "donor_id"])
    group_col = find_existing_col(adata, ["Group", "group", "Diagnosis", "diagnosis"])

    if sample_col is None:
        print(f"[SKIP] {label}: no sample column found")
        return

    df = pd.DataFrame({
        "sample": adata.obs[sample_col].astype(str).values,
        "pseudotime": adata.obs["ptime"].values,
    })

    if group_col is not None:
        df["group"] = adata.obs[group_col].astype(str).values

    summary = df.groupby("sample").agg(
        n_cells=("pseudotime", "size"),
        mean_pseudotime=("pseudotime", "mean"),
        median_pseudotime=("pseudotime", "median"),
        sd_pseudotime=("pseudotime", "std"),
        q25_pseudotime=("pseudotime", lambda x: np.quantile(x, 0.25)),
        q75_pseudotime=("pseudotime", lambda x: np.quantile(x, 0.75)),
    ).reset_index()

    if group_col is not None:
        gmap = df.groupby("sample")["group"].first().to_dict()
        summary["group"] = summary["sample"].map(gmap)

    save_table(summary, os.path.join(out_dir, f"{label}_sample_level_pseudotime_summary.csv"))

    if group_col is not None and summary["group"].nunique() >= 2:
        groups = sorted(summary["group"].unique().tolist())
        palette = make_category_palette(groups)

        fig, ax = plt.subplots(figsize=(3.0, 3.0))
        sns.boxplot(
            data=summary,
            x="group",
            y="median_pseudotime",
            palette=palette,
            showfliers=False,
            linewidth=0.6,
            ax=ax,
        )
        sns.stripplot(
            data=summary,
            x="group",
            y="median_pseudotime",
            color="black",
            size=3,
            alpha=0.7,
            jitter=0.15,
            ax=ax,
        )
        ax.set_xlabel("")
        ax.set_ylabel("Sample median pseudotime")
        ax.set_title("Donor-level pseudotime")
        sns.despine(ax=ax)
        plt.tight_layout()
        save_figure_all(fig, os.path.join(out_dir, f"{label}_Supp_sample_level_pseudotime_by_group"), dpi=DPI_SUPP)
        plt.close(fig)


# ============================================================
# 15. Optional raw UMAP vector field
# ============================================================

def plot_optional_raw_umap_vector_field(adata, label, out_dir):
    """
    Exploratory raw UMAP vector field.

    Note:
    This is intentionally disabled by default because latent vector field
    is the rigorous SCTour representation. Direct raw UMAP arrows may be
    visually useful but should be interpreted cautiously.
    """
    if not PLOT_RAW_UMAP_VECTOR_FIELD:
        return

    print(f"[PLOT] {label}: exploratory raw UMAP vector field")

    try:
        # Approximate using first two VF dimensions.
        # This is not a true UMAP-projected vector field.
        adata.obsm["X_VF_2D_APPROX"] = adata.obsm["X_VF"][:, :2].copy()

        fig, ax = plt.subplots(figsize=(3.5, 3.2))

        sct.vf.plot_vector_field(
            adata,
            zs_key="X_umap_raw",
            vf_key="X_VF_2D_APPROX",
            use_rep_neigh="X_TNODE",
            color="celltype2",
            show=False,
            ax=ax,
            legend_loc="none",
            frameon=False,
            size=15,
            alpha=0.20,
        )
        ax.set_title("Exploratory VF on original UMAP")
        plt.tight_layout()
        save_figure_all(fig, os.path.join(out_dir, f"{label}_Exploratory_raw_umap_vector_field"), dpi=DPI_SUPP)
        plt.close(fig)

    except Exception as e:
        print(f"[WARN] {label}: raw UMAP vector field failed: {e}")


# ============================================================
# 16. Per-trajectory pipeline
# ============================================================

def run_trajectory(cell_key, config):
    """
    Complete pipeline for one trajectory.
    """
    label = cell_key
    label_pretty = config["label"]

    traj_dir = make_dir(os.path.join(OUT_DIR, label))
    fig_dir = make_dir(os.path.join(traj_dir, "figures"))
    table_dir = make_dir(os.path.join(traj_dir, "tables"))
    h5ad_dir = make_dir(os.path.join(traj_dir, "h5ad"))

    print("\n" + "=" * 80)
    print(f"[START] {label}: {config['title']}")
    print("=" * 80)

    t0 = time.time()

    try:
        # ----------------------------
        # Step 1. Load data
        # ----------------------------
        adata = load_cell_data(config["subdir"], config["prefix"])

        # ----------------------------
        # Step 2. Capillary-only filtering for endothelial analysis
        # ----------------------------
        if config.get("cap_only", False):
            endo_subtype_col = config.get("endo_subtype_col", "endo_subtype")
            capillary_label = config.get("capillary_label", "Capillary")

            if endo_subtype_col not in adata.obs.columns:
                raise ValueError(
                    f"{label}: cap_only=True but column {endo_subtype_col} not found."
                )

            n0 = adata.n_obs
            adata = adata[adata.obs[endo_subtype_col].astype(str) == capillary_label].copy()
            print(f"[FILTER] {label}: capillary only {n0:,} -> {adata.n_obs:,}")

            if adata.n_obs < MIN_CELLS:
                raise ValueError(f"{label}: too few capillary cells after filtering.")

        # ----------------------------
        # Step 3. Define celltype2
        # ----------------------------
        subtype_col = config["subtype_col"]
        if subtype_col not in adata.obs.columns:
            raise ValueError(f"{label}: subtype column not found: {subtype_col}")

        adata.obs["celltype2"] = adata.obs[subtype_col].astype(str).values

        # Clean important columns
        for c in ["celltype2", "Group", "group", "orig.ident", "sample_id", "donor", "donor_id"]:
            safe_obs_category(adata, c)

        # Subtype distribution table
        subtype_table = (
            adata.obs["celltype2"]
            .astype(str)
            .value_counts()
            .reset_index()
        )
        subtype_table.columns = ["celltype2", "n_cells"]
        save_table(subtype_table, os.path.join(table_dir, f"{label}_subtype_distribution.csv"))

        print("[INFO] Subtype distribution:")
        print(subtype_table.to_string(index=False))

        # ----------------------------
        # Step 4. Raw UMAP QC panels before HVG subsetting
        # ----------------------------
        plot_raw_umap_panels(adata, label, config, fig_dir)

        # ----------------------------
        # Step 5. Preprocess
        # ----------------------------
        adata = preprocess_for_sctour(
            adata=adata,
            label=label,
            out_dir=table_dir,
            n_hvg=N_HVG,
        )

        # ----------------------------
        # Step 6. SCTour training
        # ----------------------------
        adata, trainer = train_sctour(
            adata=adata,
            label=label,
            config=config,
            out_dir=table_dir,
            nepoch=SCTOUR_NEPOCH,
        )

        # ----------------------------
        # Step 7. Latent and vector field
        # ----------------------------
        adata = compute_latent_and_vector_field(
            adata=adata,
            trainer=trainer,
            label=label,
            out_dir=table_dir,
        )

        # ----------------------------
        # Step 8. Save core tables
        # ----------------------------
        # If h5ad should be placed in h5ad_dir:
        old_save_h5ad = SAVE_H5AD
        save_core_outputs(adata, label, table_dir)

        if SAVE_H5AD:
            h5ad_path = os.path.join(h5ad_dir, f"{label}_sctour_final.h5ad")
            adata.write_h5ad(h5ad_path, compression="gzip")
            print(f"[H5AD] {h5ad_path}")

        # ----------------------------
        # Step 9. Main and supplementary figures
        # ----------------------------
        plot_sctour_main_panels(adata, label, config, fig_dir)
        plot_pseudotime_raw_and_latent(adata, label, fig_dir)
        plot_pseudotime_by_group(adata, label, fig_dir)
        plot_pseudotime_by_state(adata, label, fig_dir)
        plot_marker_trends(adata, label, config, fig_dir)
        plot_marker_heatmap_binned(adata, label, config, fig_dir)
        compute_pseudotime_bin_composition(adata, label, fig_dir)
        compute_sample_level_summary(adata, label, fig_dir)
        plot_optional_raw_umap_vector_field(adata, label, fig_dir)

        # ----------------------------
        # Step 10. Run summary
        # ----------------------------
        elapsed = (time.time() - t0) / 60

        run_summary = {
            "trajectory": label,
            "label": label_pretty,
            "title": config["title"],
            "n_cells": int(adata.n_obs),
            "n_hvg": int(adata.n_vars),
            "elapsed_min": elapsed,
            "output_dir": traj_dir,
            "finished_time": now(),
        }

        write_json(run_summary, os.path.join(traj_dir, f"{label}_run_summary.json"))

        print(f"[DONE] {label}: elapsed {elapsed:.1f} min")
        return adata

    except Exception as e:
        print(f"[ERROR] {label} failed: {e}")
        traceback.print_exc()

        err_path = os.path.join(traj_dir, f"{label}_ERROR.txt")
        with open(err_path, "w", encoding="utf-8") as f:
            f.write(f"Time: {now()}\n")
            f.write(f"Trajectory: {label}\n")
            f.write(str(e) + "\n\n")
            f.write(traceback.format_exc())

        raise e


# ============================================================
# 17. Cross-trajectory summary
# ============================================================

def collect_cross_trajectory_summary():
    """
    Collect key output tables across trajectories.
    """
    print("[SUMMARY] collecting cross-trajectory summaries")

    rows = []

    for label in CELL_TYPE_CONFIG.keys():
        summary_path = os.path.join(OUT_DIR, label, f"{label}_run_summary.json")
        if os.path.exists(summary_path):
            with open(summary_path, "r", encoding="utf-8") as f:
                x = json.load(f)
            rows.append(x)

    if len(rows) > 0:
        df = pd.DataFrame(rows)
        save_table(df, os.path.join(OUT_DIR, "node07_cross_trajectory_run_summary.csv"))

    # Also generate a manuscript output index
    index_rows = []
    for root, dirs, files in os.walk(OUT_DIR):
        for fn in files:
            if fn.endswith((".csv", ".png", ".pdf", ".svg", ".json", ".h5ad", ".txt")):
                path = os.path.join(root, fn)
                index_rows.append({
                    "file": fn,
                    "path": path,
                    "size_MB": os.path.getsize(path) / 1024**2,
                })

    if len(index_rows) > 0:
        index_df = pd.DataFrame(index_rows)
        save_table(index_df, os.path.join(OUT_DIR, "node07_output_file_index.csv"))


# ============================================================
# 18. Save global parameters
# ============================================================

def save_global_parameters():
    params = {
        "BASE_DIR": BASE_DIR,
        "OUT_DIR": OUT_DIR,
        "RANDOM_SEED": RANDOM_SEED,
        "N_HVG": N_HVG,
        "MIN_CELLS": MIN_CELLS,
        "SCTOUR_NEPOCH": SCTOUR_NEPOCH,
        "MAX_DENSE_GB": MAX_DENSE_GB,
        "PLOT_RAW_UMAP_VECTOR_FIELD": PLOT_RAW_UMAP_VECTOR_FIELD,
        "SAVE_H5AD": SAVE_H5AD,
        "LATENT_N_NEIGHBORS": LATENT_N_NEIGHBORS,
        "LATENT_UMAP_MIN_DIST": LATENT_UMAP_MIN_DIST,
        "N_PTIME_BINS": N_PTIME_BINS,
        "CELL_TYPE_CONFIG": CELL_TYPE_CONFIG,
    }

    write_json(params, os.path.join(OUT_DIR, "node07_global_parameters.json"))


# ============================================================
# 19. Main
# ============================================================

def main():
    print("=" * 80)
    print("node07 Final: SCTour trajectory analysis")
    print(f"Start time: {now()}")
    print(f"Output: {OUT_DIR}")
    print("=" * 80)

    t0 = time.time()

    np.random.seed(RANDOM_SEED)

    save_global_parameters()
    log_versions(OUT_DIR)

    adata_results = {}

    # Run three trajectories
    for cell_key in ["micro", "astro", "endo"]:
        adata_results[cell_key] = run_trajectory(
            cell_key=cell_key,
            config=CELL_TYPE_CONFIG[cell_key],
        )

    collect_cross_trajectory_summary()

    publication_script = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "node07_publication_package.py"
    )
    if os.path.exists(publication_script):
        print(f"[PUBLICATION] Running final package generator: {publication_script}")
        subprocess.run([sys.executable, publication_script], check=True)
    else:
        print(f"[WARN] Publication package generator not found: {publication_script}")

    elapsed = (time.time() - t0) / 60

    print("\n" + "=" * 80)
    print("node07 Final complete!")
    print(f"Total elapsed time: {elapsed:.1f} min")
    print(f"End time: {now()}")
    print("=" * 80)


if __name__ == "__main__":
    main()
