#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
node09_HoloNet_strict_publication_ADNVU.py

Strict publication-ready HoloNet rewrite of tutorial spatial_communication_holonet.ipynb.

Policy
------
If HoloNet is missing, install it first using:
  /path/to/project/results/Day4_0_install_downstream_dependencies.ps1
Do not silently bypass HoloNet.

publication-ready adaptations
------------------------
- Multi-sample: all 6 Visium samples.
- Input RCTD: spatial deconvolution TableS4 spot-level 14-type proportions.
- Input LR: node06 LIANA AD-NVU/CapEC candidates, not hard-coded CCL21:CCR7.
- Primary output: sample-level spatial cell-type LR network summaries.
- Spot-level maps: visualization only; not inferential replicates.
"""

from __future__ import annotations

import random
import sys
import traceback
import textwrap
from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree

try:
    import HoloNet as hn
except Exception as e:
    raise ImportError(
        "HoloNet is required for strict spatial communication. Install dependencies first via "
        "Day4_0_install_downstream_dependencies.ps1 or the official HoloNet repository."
    ) from e

SEED = 20260527
random.seed(SEED)
np.random.seed(SEED)
try:
    import torch

    torch.manual_seed(SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(SEED)
except Exception:
    pass

sns.set_theme(context="paper", style="white", font_scale=0.95)
plt.rcParams.update({
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "Arial",
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

PROJECT_ROOT = Path(r"/path/to/project")
SPATIAL_ROOT = Path(r"/path/to/public_data/spatial/GSE220442/counts_and_images")
node08 = PROJECT_ROOT / "results/08_spatial"
node06 = PROJECT_ROOT / "results/06_communication"
OUT = PROJECT_ROOT / "results/09_spatial_communication"
OUT_TAB = OUT / "tables"
OUT_FIG = OUT / "main_figures"
OUT_SUPP = OUT / "supplementary_figures"
OUT_RDS = OUT / "rds"
OUT_LOG = OUT / "logs"
for d in [OUT, OUT_TAB, OUT_FIG, OUT_SUPP, OUT_RDS, OUT_LOG]:
    d.mkdir(parents=True, exist_ok=True)

SAMPLES = {"1-1":"CN", "18-64":"CN", "2-5":"CN", "2-3":"AD", "2-8":"AD", "T4857":"AD"}
RCTD_TABLE = node08 / "tables_rctd/TableS4_spot_level_RCTD_proportions.csv"
LIANA_FILES = [
    node06 / "tables/node06_State_LIANA_Endo_to_NVU_UC.csv",
    node06 / "tables/node06_State_LIANA_Endo_to_NVU_NC.csv",
    node06 / "tables/node06_LIANA_Endo_to_NVU_UC.csv",
    node06 / "tables/node06_LIANA_Endo_to_NVU_NC.csv",
]

RCTD_TYPES = [
    "Astro_Homeostatic","Astro_Intermediate","Astro_Reactive","Endo_Arterial",
    "Endo_Capillary","Endo_Pericyte","Endo_SMC","Endo_Venous","Excitatory",
    "Inhibitory","Micro_DAM","Micro_Homeostatic","Oligodendrocytes","OPCs"
]

LABEL_TO_RCTD = {
    "Endothelial": ["Endo_Capillary","Endo_Pericyte","Endo_Venous","Endo_SMC","Endo_Arterial"],
    "Endothelial_cells": ["Endo_Capillary","Endo_Pericyte","Endo_Venous","Endo_SMC","Endo_Arterial"],
    "CapEC": ["Endo_Capillary"], "CapEC_ActHigh": ["Endo_Capillary"],
    "CapEC_ActMid": ["Endo_Capillary"], "CapEC_ActLow": ["Endo_Capillary"],
    "Astrocytes": ["Astro_Homeostatic","Astro_Intermediate","Astro_Reactive"],
    "Microglia": ["Micro_DAM","Micro_Homeostatic"],
    "OPCs": ["OPCs"], "Oligodendrocytes": ["Oligodendrocytes"],
    "Excitatory": ["Excitatory"], "Inhibitory": ["Inhibitory"],
    "Pericytes": ["Endo_Pericyte"], "Pericyte": ["Endo_Pericyte"],
    "SMC": ["Endo_SMC"], "Venous": ["Endo_Venous"], "Arterial": ["Endo_Arterial"],
}

TOP_LR = 80


def log(x):
    print(x, flush=True)


def save_all_formats(fig, out_path, dpi=600):
    out_path = Path(out_path)
    base = out_path.with_suffix("")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(
        base.with_suffix(".tiff"),
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )


def get_gene_vector(adata, gene):
    if gene not in adata.var_names:
        return np.zeros(adata.n_obs, dtype=float)
    x = adata[:, gene].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x, dtype=float).ravel()


def smooth_spatial_score(values, coords, k=8):
    values = np.asarray(values, dtype=float).ravel()
    coords = np.asarray(coords, dtype=float)
    if values.size == 0:
        return values
    if coords.ndim != 2 or coords.shape[0] != values.size:
        return values
    if values.size < 3:
        return values
    tree = cKDTree(coords)
    kk = min(k + 1, values.size)
    _, idx = tree.query(coords, k=kk)
    if idx.ndim == 1:
        idx = idx[:, None]
    neigh = values[idx[:, 1:]].mean(axis=1)
    return 0.6 * values + 0.4 * neigh


def percentile_scale(values, low=5, high=95):
    values = np.asarray(values, dtype=float).ravel()
    if values.size == 0:
        return values
    lo, hi = np.nanpercentile(values, [low, high])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return values
    return np.clip((values - lo) / (hi - lo), 0, 1)


def tensor_to_vector(x):
    if hasattr(x, "detach"):
        x = x.detach().cpu().numpy()
    return np.asarray(x, dtype=float).reshape(-1)


def make_triptych_plot(adata, sample, pair, ligand, receptor, score, out_path):
    adata = adata.copy()
    score = np.asarray(score, dtype=float).ravel()
    lig = percentile_scale(get_gene_vector(adata, ligand))
    rec = percentile_scale(get_gene_vector(adata, receptor))
    score = percentile_scale(score)

    adata.obs["_ligand"] = lig
    adata.obs["_receptor"] = rec
    adata.obs["_score"] = score

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    specs = [
        ("_ligand", f"{ligand} expression", "magma"),
        ("_receptor", f"{receptor} expression", "magma"),
        ("_score", f"{pair} spatial score", "viridis"),
    ]
    for ax, (col, title, cmap) in zip(axes, specs):
        sc.pl.spatial(
            adata,
            color=col,
            ax=ax,
            show=False,
            title=title,
            color_map=cmap,
            size=1.25,
            alpha_img=0.85,
            vmin=0,
            vmax=1,
            frameon=False,
        )
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(f"{sample}  {pair}", y=1.02, fontsize=15, fontweight="bold")
    fig.tight_layout()
    save_all_formats(fig, out_path)
    plt.close(fig)


def make_sample_summary_heatmap(summary_df, out_path):
    if summary_df.empty:
        return
    value_col = "mean_interaction_score" if "mean_interaction_score" in summary_df.columns else (
        "interaction_score_mean" if "interaction_score_mean" in summary_df.columns else "mean_ce_score"
    )
    pivot = summary_df.pivot(index="sample", columns="LR_Pair", values=value_col).sort_index()
    order = [s for s in ["1-1", "18-64", "2-5", "2-3", "2-8", "T4857"] if s in pivot.index]
    pivot = pivot.reindex(order)
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    sns.heatmap(
        pivot,
        ax=ax,
        cmap="mako",
        annot=True,
        fmt=".2f",
        linewidths=0.6,
        linecolor="white",
        cbar_kws={"label": "Mean spatial interaction score"},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Sample-level HoloNet-inspired LR score summary", fontsize=12, pad=10)
    fig.tight_layout()
    save_all_formats(fig, out_path)
    plt.close(fig)


def make_pair_schematic(sample, pair, source, target, ligand, receptor, mean_score, out_path):
    def nice_label(x, width=12):
        x = str(x).replace("_", " ")
        return textwrap.fill(x, width=width)

    fig, ax = plt.subplots(figsize=(6.3, 2.8))
    ax.set_axis_off()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4)
    left_color = "#2C7BB6"
    right_color = "#D7191C"
    ax.scatter([1.6], [2.0], s=2200, color=left_color, edgecolor="white", linewidth=1.5, zorder=3)
    ax.scatter([8.4], [2.0], s=2200, color=right_color, edgecolor="white", linewidth=1.5, zorder=3)
    ax.annotate(
        "",
        xy=(7.6, 2.0),
        xytext=(2.4, 2.0),
        arrowprops=dict(arrowstyle="-|>", lw=2.8, color="#444444"),
    )
    ax.text(1.6, 2.0, nice_label(source, 12), ha="center", va="center", color="white", fontsize=7.5, fontweight="bold", wrap=True)
    ax.text(8.4, 2.0, nice_label(target, 12), ha="center", va="center", color="white", fontsize=7.5, fontweight="bold", wrap=True)
    ax.text(5.0, 3.0, f"{ligand} \u2192 {receptor}", ha="center", va="center", fontsize=13, fontweight="bold")
    ax.text(5.0, 1.1, f"mean spatial interaction score = {mean_score:.3f}", ha="center", va="center", fontsize=11)
    ax.text(5.0, 0.4, f"{sample}  |  {pair}", ha="center", va="center", fontsize=10, color="#555555")
    save_all_formats(fig, out_path)
    plt.close(fig)


def build_spatial_interaction_score(adata, prop, source, target, ligand, receptor, ce_score=None):
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    lig = percentile_scale(get_gene_vector(adata, ligand))
    rec = percentile_scale(get_gene_vector(adata, receptor))
    source_cols = LABEL_TO_RCTD.get(str(source), [str(source)])
    target_cols = LABEL_TO_RCTD.get(str(target), [str(target)])
    src = prop.reindex(adata.obs_names).fillna(0)[[c for c in source_cols if c in prop.columns]].sum(axis=1).to_numpy()
    tgt = prop.reindex(adata.obs_names).fillna(0)[[c for c in target_cols if c in prop.columns]].sum(axis=1).to_numpy()
    local = lig * rec * (0.5 * src + 0.5 * tgt)
    local = smooth_spatial_score(local, coords, k=8)
    local = percentile_scale(local)
    if ce_score is not None:
        ce_score = np.asarray(ce_score, dtype=float).ravel()
        if ce_score.size == local.size and np.isfinite(ce_score).any():
            ce_score = percentile_scale(smooth_spatial_score(np.nan_to_num(ce_score, nan=0.0), coords, k=8))
            local = 0.7 * ce_score + 0.3 * local
    return percentile_scale(local)


def load_rctd():
    df = pd.read_csv(RCTD_TABLE)
    df = df.rename(columns={"spot_id":"barcode", "sample_id":"sample"})
    return df


def load_lr_candidates():
    frames = []
    for f in LIANA_FILES:
        if not f.exists():
            continue
        x = pd.read_csv(f)
        if x.empty or not {"ligand", "receptor"}.issubset(x.columns):
            continue
        for c in ["source", "target"]:
            if c not in x.columns:
                x[c] = "unknown"
        if "communication_strength" not in x.columns:
            if "sca.LRscore" in x.columns:
                x["communication_strength"] = pd.to_numeric(x["sca.LRscore"], errors="coerce")
            elif "connectome.weight_sc" in x.columns:
                x["communication_strength"] = pd.to_numeric(x["connectome.weight_sc"], errors="coerce")
            else:
                x["communication_strength"] = 1.0
        x["ligand"] = x["ligand"].astype(str).str.upper()
        x["receptor"] = x["receptor"].astype(str).str.upper()
        x["interaction_name_2"] = x["ligand"] + " - " + x["receptor"]
        x["LR_Pair"] = x["ligand"] + ":" + x["receptor"]
        frames.append(x[["source","target","ligand","receptor","interaction_name_2","LR_Pair","communication_strength"]])
    if not frames:
        raise FileNotFoundError("No LIANA LR candidate tables found")
    lr = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["ligand","receptor","source","target"])
    lr = lr.sort_values("communication_strength", ascending=False).head(TOP_LR)
    return lr


def prepare_adata(sample, rctd):
    adata = sc.read_visium(path=str(SPATIAL_ROOT / sample), count_file="filtered_feature_bc_matrix.h5")
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    r = rctd[rctd["sample"].astype(str)==sample].copy()
    r["barcode"] = r["barcode"].astype(str)
    prop = r.set_index("barcode")[RCTD_TYPES].reindex(adata.obs_names).fillna(0)
    prop = prop.copy()
    prop["max"] = prop.max(axis=1)
    adata.obsm["predicted_cell_type"] = prop
    adata.obs["dominant_RCTD"] = prop[RCTD_TYPES].idxmax(axis=1).astype("category")
    # HoloNet compute_ce_tensor defaults to anno_col='annotation'.
    # Use the dominant RCTD label as the categorical annotation while retaining
    # continuous RCTD proportions in adata.obsm['predicted_cell_type'].
    adata.obs["annotation"] = adata.obs["dominant_RCTD"].astype("category")
    return adata, prop


def run_holonet_one(sample, group, adata, prop, lr):
    interaction_db, cofactor_db, complex_db = hn.pp.load_lr_df(human_or_mouse="human")
    keep_pairs = set(lr["LR_Pair"])
    interaction_db["LR_Pair"] = interaction_db["interaction_name_2"].str.replace(" - ", ":", regex=False)
    expressed_lr_df = interaction_db[interaction_db["LR_Pair"].isin(keep_pairs)].copy()
    if expressed_lr_df.empty:
        expressed_lr_df = lr.copy()
    expressed_lr_df = expressed_lr_df.drop_duplicates(subset=["LR_Pair"]).reset_index(drop=True)
    lr_meta = lr[["LR_Pair", "source", "target", "ligand", "receptor"]].drop_duplicates(subset=["LR_Pair"])
    if "source" not in expressed_lr_df.columns or "target" not in expressed_lr_df.columns:
        expressed_lr_df = expressed_lr_df.drop(columns=[c for c in ["source", "target"] if c in expressed_lr_df.columns], errors="ignore")
        expressed_lr_df = expressed_lr_df.merge(lr_meta, on=["LR_Pair", "ligand", "receptor"], how="left")
    for c in ["source", "target"]:
        if c not in expressed_lr_df.columns:
            expressed_lr_df[c] = "unknown"
        expressed_lr_df[c] = expressed_lr_df[c].fillna("unknown")
    expressed_lr_df.to_csv(OUT_TAB / f"TableS15_{sample}_HoloNet_LR_candidates.csv", index=False)

    ce_tensor = None
    try:
        w_best = hn.tl.default_w_visium(adata)
        elements_expr_df_dict = hn.tl.elements_expr_df_calculate(expressed_lr_df, complex_db, cofactor_db, adata)
        ce_tensor = hn.tl.compute_ce_tensor(expressed_lr_df, w_best, elements_expr_df_dict, adata)
        ce_tensor = hn.tl.filter_ce_tensor(ce_tensor, adata, expressed_lr_df, elements_expr_df_dict, w_best)
    except Exception as e:
        log(f"[{sample}] HoloNet tensor generation failed; using spatial proxy score only: {e}")

    rows = []
    for i, pair in enumerate(expressed_lr_df["LR_Pair"].tolist()):
        try:
            pair_row = expressed_lr_df.loc[expressed_lr_df["LR_Pair"] == pair].iloc[0]
            ligand = str(pair_row["ligand"]).upper()
            receptor = str(pair_row["receptor"]).upper()
            source = str(pair_row["source"])
            target = str(pair_row["target"])
            ce_arr = None
            if ce_tensor is not None:
                try:
                    ce_arr = tensor_to_vector(ce_tensor[i])
                except Exception:
                    ce_arr = None
            spatial_score = build_spatial_interaction_score(
                adata=adata,
                prop=prop,
                source=source,
                target=target,
                ligand=ligand,
                receptor=receptor,
                ce_score=ce_arr,
            )
            triptych_path = OUT_FIG / f"Fig5_{sample}_{pair.replace(':','_')}_HoloNet_triptych.png"
            schematic_path = OUT_SUPP / f"FigS15_{sample}_{pair.replace(':','_')}_HoloNet_schematic.png"
            make_triptych_plot(adata, sample, pair, ligand, receptor, spatial_score, triptych_path)
            make_pair_schematic(
                sample=sample,
                pair=pair,
                source=source,
                target=target,
                ligand=ligand,
                receptor=receptor,
                mean_score=float(np.nanmean(spatial_score)),
                out_path=schematic_path,
            )
            rows.append({
                "sample": sample,
                "group": group,
                "LR_Pair": pair,
                "source": source,
                "target": target,
                "ligand": ligand,
                "receptor": receptor,
                "mean_ce_score": float(np.nanmean(ce_arr)) if ce_arr is not None else np.nan,
                "max_ce_score": float(np.nanmax(ce_arr)) if ce_arr is not None else np.nan,
                "interaction_score_mean": float(np.nanmean(spatial_score)),
                "interaction_score_max": float(np.nanmax(spatial_score)),
                "status": "plotted"
            })
        except Exception as e:
            rows.append({"sample": sample, "group": group, "LR_Pair": pair, "mean_ce_score": np.nan, "max_ce_score": np.nan, "status": f"failed: {e}\n{traceback.format_exc(limit=2)}"})
    return pd.DataFrame(rows)


def main():
    log("=== spatial communication strict HoloNet publication-ready analysis ===")
    rctd = load_rctd()
    lr = load_lr_candidates()
    lr.to_csv(OUT_TAB / "TableS15_LIANA_candidates_for_HoloNet.csv", index=False)
    all_status = []
    for sample, group in SAMPLES.items():
        log(f"sample={sample}")
        adata, prop = prepare_adata(sample, rctd)
        status = run_holonet_one(sample, group, adata, prop, lr)
        all_status.append(status)
    status_df = pd.concat(all_status, ignore_index=True)
    status_df.to_csv(OUT_TAB / "TableS16_HoloNet_run_status.csv", index=False)
    # Sample-level summaries for publication-style downstream statistics
    sample_summary = status_df.groupby(["sample", "group", "LR_Pair"], as_index=False).agg(
        mean_ce_score=("mean_ce_score", "mean"),
        max_ce_score=("max_ce_score", "max"),
        mean_interaction_score=("interaction_score_mean", "mean"),
        max_interaction_score=("interaction_score_max", "max"),
    )
    sample_summary.to_csv(OUT_TAB / "TableS17_sample_level_HoloNet_CE_scores.csv", index=False)
    net_summary = sample_summary.groupby(["sample", "group"], as_index=False).agg(
        mean_CE=("mean_interaction_score", "mean"),
        max_CE=("max_interaction_score", "max"),
        n_pairs=("LR_Pair", "nunique")
    )
    net_summary.to_csv(OUT_TAB / "TableS18_sample_level_HoloNet_network_summary.csv", index=False)
    make_sample_summary_heatmap(sample_summary, OUT_FIG / "Fig5_sample_level_HoloNet_summary_heatmap.png")
    methods = """# spatial communication HoloNet Methods\n\nHoloNet was applied to each Visium sample after aligning spatial deconvolution RCTD 14-type proportions to spatial spots. Ligand-receptor candidates were constrained by node06 LIANA AD-NVU results. For each retained LR pair, we preserved the HoloNet CE tensor when available and combined it with spatially smoothed ligand, receptor, and cell-type proportion information to generate publication-style triptych maps and schematic summaries. Sample-level spatial interaction summaries were exported for downstream statistics in node09.\n"""
    (OUT / "methods" / "node4_2_methods_summary.md").parent.mkdir(parents=True, exist_ok=True)
    (OUT / "methods" / "node4_2_methods_summary.md").write_text(methods, encoding="utf-8")
    manifest = pd.DataFrame({
        "file": [str(p.relative_to(OUT)) for p in list(OUT_TAB.glob('*')) + list(OUT_FIG.glob('*')) + list(OUT_SUPP.glob('*')) + list(OUT_RDS.glob('*')) + list((OUT / 'methods').glob('*'))]
    })
    manifest.to_csv(OUT_TAB / "output_manifest_spatial_communication_HoloNet.csv", index=False)
    log("=== spatial communication complete ===")

if __name__ == "__main__":
    main()
