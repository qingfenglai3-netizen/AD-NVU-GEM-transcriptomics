#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
node09_HoloNet_full_NVU_landscape_ADNVU.py

Spatial communication analysis supporting Figure 7.
This script is deliberately separated from the existing final manuscript figures.
It evaluates the complete LIANA-prioritized NVU ligand-receptor landscape and
separates official HoloNet CE-supported axes from spatial LR support-only axes.

Outputs are written to:
  results/09_spatial_communication/full_NVU_landscape_extension
"""
from __future__ import annotations

import json
import random
import importlib.util
from pathlib import Path
import os

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

try:
    import torch
except Exception:
    torch = None

try:
    import HoloNet as hn
except Exception as exc:
    raise ImportError("HoloNet is required for the full NVU landscape extension.") from exc

SEED = 20260527
random.seed(SEED)
np.random.seed(SEED)
if torch is not None:
    torch.manual_seed(SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(SEED)

PROJECT_ROOT = Path(os.getenv("AD_NVU_PROJECT_ROOT", "/path/to/project"))
REPO_ROOT = Path(os.getenv("AD_NVU_REPO_ROOT", ".")).resolve()
BASE_SCRIPT = REPO_ROOT / "analysis/09_holonet/09a_holonet_base.py"
OUT = PROJECT_ROOT / "results/09_spatial_communication"
OUT_TAB = OUT / "tables"
OUT_FIG = OUT / "figures"
OUT_TENSOR = OUT / "tensors"
OUT_AUDIT = OUT / "audit"
for d in (OUT, OUT_TAB, OUT_FIG, OUT_TENSOR, OUT_AUDIT):
    d.mkdir(parents=True, exist_ok=True)

spec = importlib.util.spec_from_file_location("spatial_communication_base", BASE_SCRIPT)
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

SAMPLES = base.SAMPLES
RCTD_TYPES = base.RCTD_TYPES
NVU_CELLTYPES = [
    "Astro_Homeostatic", "Astro_Intermediate", "Astro_Reactive",
    "Endo_Arterial", "Endo_Capillary", "Endo_Pericyte", "Endo_SMC", "Endo_Venous",
    "Micro_DAM", "Micro_Homeostatic",
]
PRIORITY_KEYWORDS = [
    "ANGPT", "TEK", "TIE1", "DLL4", "NOTCH3", "VEGF", "FLT1", "KDR",
    "VWF", "ITGB1", "LRP1", "FN1", "CD44", "ITGA", "SPARC", "ENG",
    "EDN1", "ADGRL4", "TGM2", "SDC4", "PDGF", "CXCL", "DPP4", "SEMA",
]
MAX_OFFICIAL_CE_PAIRS = 24
MAX_ROWS_FIG = 32

matplotlib.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "Arial", "font.size": 6, "font.weight": "bold",
    "axes.labelsize": 6, "axes.titlesize": 7, "xtick.labelsize": 5.5, "ytick.labelsize": 5.5,
    "legend.fontsize": 5.5, "axes.linewidth": 0.55,
})
sns.set_theme(context="paper", style="white", font="Arial")


def save_all(fig, stem: Path, dpi=600):
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=dpi, bbox_inches="tight", facecolor="white", pil_kwargs={"compression":"tiff_lzw"})


def as_numpy(x):
    if hasattr(x, "detach"):
        return x.detach().cpu().numpy()
    return np.asarray(x)


def gene_vec(adata, gene):
    gene = str(gene).upper()
    if gene not in adata.var_names:
        return np.zeros(adata.n_obs, dtype=float)
    x = adata[:, gene].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x, dtype=float).ravel()


def expr_summary(adata, gene):
    v = gene_vec(adata, gene)
    return float(np.mean(v > 0)), float(np.nanmean(v)), float(np.nanmax(v) if v.size else 0)


def nvu_priority(pair, source, target):
    text = f"{pair}|{source}|{target}".upper()
    score = sum(1 for k in PRIORITY_KEYWORDS if k in text)
    if any(ct.upper() in text for ct in ["CAPEC", "ENDO", "PERICY", "SMC", "ASTRO", "MICRO"]):
        score += 2
    return score


def get_holonet_db():
    interaction_db, cofactor_db, complex_db = hn.pp.load_lr_df(human_or_mouse="human")
    interaction_db = interaction_db.copy()
    interaction_db["LR_Pair"] = interaction_db["interaction_name_2"].str.replace(" - ", ":", regex=False).str.upper()
    interaction_db["ligand"] = interaction_db["ligand"].astype(str).str.upper()
    interaction_db["receptor"] = interaction_db["receptor"].astype(str).str.upper()
    return interaction_db, cofactor_db, complex_db


def classify_candidates(lr, interaction_db):
    lr = lr.copy()
    lr["ligand"] = lr["ligand"].astype(str).str.upper()
    lr["receptor"] = lr["receptor"].astype(str).str.upper()
    lr["LR_Pair"] = lr["ligand"] + ":" + lr["receptor"]
    db_pairs = set(interaction_db["LR_Pair"])
    lr["holonet_db_exact"] = lr["LR_Pair"].isin(db_pairs)
    lr["nvu_priority_score"] = [nvu_priority(p, s, t) for p, s, t in zip(lr["LR_Pair"], lr["source"], lr["target"])]
    lr["candidate_tier"] = np.where(lr["holonet_db_exact"], "official_CE_candidate", "spatial_LR_support_only")
    return lr.sort_values(["holonet_db_exact", "nvu_priority_score", "communication_strength"], ascending=[False, False, False])


def official_expressed_df(lr_sample, interaction_db):
    exact_pairs = set(lr_sample.loc[lr_sample["holonet_db_exact"], "LR_Pair"])
    expressed = interaction_db[interaction_db["LR_Pair"].isin(exact_pairs)].copy()
    if expressed.empty:
        return expressed
    expressed = expressed.drop_duplicates(subset=["LR_Pair"]).reset_index(drop=True)
    meta = lr_sample[["LR_Pair", "source", "target", "communication_strength", "nvu_priority_score"]].drop_duplicates("LR_Pair")
    expressed = expressed.drop(columns=[c for c in ["source", "target", "communication_strength", "nvu_priority_score"] if c in expressed.columns], errors="ignore")
    expressed = expressed.merge(meta, on="LR_Pair", how="left")
    expressed["source"] = expressed["source"].fillna("unknown")
    expressed["target"] = expressed["target"].fillna("unknown")
    return expressed.sort_values(["nvu_priority_score", "communication_strength"], ascending=[False, False]).head(MAX_OFFICIAL_CE_PAIRS).reset_index(drop=True)


def network_from_ce(ce_one, prop):
    prop_arr = prop[RCTD_TYPES].to_numpy(dtype=float)
    prop_arr = np.nan_to_num(prop_arr, nan=0.0)
    row_sum = prop_arr.sum(axis=1, keepdims=True)
    prop_norm = np.divide(prop_arr, row_sum, out=np.zeros_like(prop_arr), where=row_sum > 0)
    arr = np.asarray(ce_one, dtype=float)
    if arr.ndim == 2 and arr.shape[0] == prop_norm.shape[0] and arr.shape[1] == prop_norm.shape[0]:
        denom = np.maximum(prop_norm.sum(axis=0)[:, None] * prop_norm.sum(axis=0)[None, :], 1e-9)
        return (prop_norm.T @ np.nan_to_num(arr, nan=0.0) @ prop_norm) / denom
    flat = arr.reshape(-1)
    if flat.size == prop_norm.shape[0]:
        weighted = prop_norm * flat[:, None]
        return np.diag(weighted.sum(axis=0) / np.maximum(prop_norm.sum(axis=0), 1e-9))
    return np.zeros((len(RCTD_TYPES), len(RCTD_TYPES)), dtype=float)


def hotspot_from_ce(ce_one, n):
    arr = np.asarray(ce_one, dtype=float)
    if arr.ndim == 2 and arr.shape[0] == n and arr.shape[1] == n:
        score = np.nan_to_num(arr, nan=0.0).sum(axis=0) + np.nan_to_num(arr, nan=0.0).sum(axis=1)
    else:
        flat = arr.reshape(-1)
        score = flat if flat.size == n else np.zeros(n)
    lo, hi = np.nanpercentile(score, [5, 95]) if score.size else (0, 0)
    if hi > lo:
        score = np.clip((score - lo) / (hi - lo), 0, 1)
    return score


def compute_official_ce(sample, adata, expressed, cofactor_db, complex_db):
    w_best = hn.tl.default_w_visium(adata)
    elements = hn.tl.elements_expr_df_calculate(expressed, complex_db, cofactor_db, adata)
    ce_tensor = hn.tl.compute_ce_tensor(expressed, w_best, elements, adata)
    ce_tensor = hn.tl.filter_ce_tensor(ce_tensor, adata, expressed, elements, w_best)
    ce = as_numpy(ce_tensor)
    np.savez_compressed(OUT_TENSOR / f"{sample}_official_CE_tensor.npz", ce_tensor=ce, lr_pairs=expressed["LR_Pair"].astype(str).to_numpy())
    return ce, float(w_best)


def plot_landscape_heatmap(sample_scores):
    d = sample_scores.copy()
    axis_rank = (d.groupby("LR_Pair", as_index=False)
                   .agg(mean_score=("spatial_score_mean", "mean"), max_strength=("communication_strength", "max"), holonet_db_exact=("holonet_db_exact", "max"), nvu_priority_score=("nvu_priority_score", "max"))
                   .sort_values(["holonet_db_exact", "nvu_priority_score", "mean_score", "max_strength"], ascending=[False, False, False, False]))
    keep = axis_rank.head(MAX_ROWS_FIG)["LR_Pair"].tolist()
    piv = d[d["LR_Pair"].isin(keep)].pivot_table(index="LR_Pair", columns="sample", values="spatial_score_mean", aggfunc="mean").reindex(keep)
    cols = [s for s in SAMPLES.keys() if s in piv.columns]
    piv = piv[cols]
    fig_h = max(4.0, 0.16 * len(piv) + 1.2)
    fig, ax = plt.subplots(figsize=(6.69, fig_h))
    sns.heatmap(piv, cmap="rocket_r", ax=ax, cbar_kws={"label":"Mean spatial LR score", "shrink":0.6}, linewidths=0.15, linecolor="white")
    ax.set_title("Full NVU spatial ligand-receptor landscape")
    ax.set_xlabel("Visium sample")
    ax.set_ylabel("Candidate LR axis")
    save_all(fig, OUT_FIG / "FULL01_full_NVU_spatial_LR_landscape_heatmap")
    plt.close(fig)


def plot_delta_lollipop(axis_summary):
    d = axis_summary.sort_values(["holonet_db_exact", "abs_AD_minus_CN", "mean_score"], ascending=[False, False, False]).head(24).iloc[::-1]
    colors = np.where(d["holonet_db_exact"], "#B23A48", "#3E7CB1")
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    ax.axvline(0, color="#555555", lw=0.6)
    ax.hlines(d["LR_Pair"], 0, d["AD_minus_CN_mean"], color="#999999", lw=0.7)
    ax.scatter(d["AD_minus_CN_mean"], d["LR_Pair"], s=np.clip(d["mean_score"]*280, 18, 160), c=colors, edgecolor="black", linewidth=0.25)
    ax.set_xlabel("AD - CN mean spatial LR score")
    ax.set_ylabel("")
    ax.set_title("Candidate axes with AD/CN directional contrast")
    sns.despine(ax=ax)
    save_all(fig, OUT_FIG / "FULL02_AD_minus_CN_candidate_axis_lollipop")
    plt.close(fig)


def plot_ce_contribution(contrib):
    if contrib.empty:
        return
    d = contrib.copy()
    axis_order = (d.groupby("LR_Pair")["contribution"].sum().sort_values(ascending=False).index.tolist())
    ct_order = (d.groupby("celltype")["contribution"].mean().sort_values(ascending=False).head(10).index.tolist())
    piv = d[d["celltype"].isin(ct_order)].pivot_table(index="LR_Pair", columns="celltype", values="contribution", aggfunc="mean").reindex(axis_order).fillna(0)
    fig, ax = plt.subplots(figsize=(6.69, max(2.5, 0.25*len(piv)+1)))
    sns.heatmap(piv, cmap="viridis", ax=ax, cbar_kws={"label":"Mean CE contribution", "shrink":0.65}, linewidths=0.15, linecolor="white")
    ax.set_title("Official HoloNet CE cell-state contribution")
    ax.set_xlabel("NVU / brain cell state")
    ax.set_ylabel("Official CE-supported LR")
    save_all(fig, OUT_FIG / "FULL03_official_CE_celltype_contribution_heatmap")
    plt.close(fig)


def plot_evidence_bubble(axis_summary):
    d = axis_summary.sort_values(["holonet_db_exact", "nvu_priority_score", "mean_score"], ascending=[False, False, False]).head(32)
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    colors = np.where(d["holonet_db_exact"], "#C75146", "#2F6F73")
    ax.scatter(d["communication_strength_max"], d["mean_score"], s=np.clip(d["n_source_target_edges"]*22, 24, 220), c=colors, edgecolor="black", linewidth=0.25, alpha=0.85)
    for _, r in d.head(14).iterrows():
        ax.text(r["communication_strength_max"]+0.015, r["mean_score"], r["LR_Pair"], fontsize=4.8, va="center")
    ax.set_xlabel("LIANA max communication strength")
    ax.set_ylabel("Mean spatial LR score")
    ax.set_title("Evidence map for NVU candidate axes")
    sns.despine(ax=ax)
    save_all(fig, OUT_FIG / "FULL04_NVU_axis_evidence_bubble")
    plt.close(fig)


def main():
    params = {
        "seed": SEED,
        "base_script": str(BASE_SCRIPT),
        "extension_type": "full_NVU_landscape_with_official_CE_and_spatial_support_layers",
        "max_official_ce_pairs": MAX_OFFICIAL_CE_PAIRS,
        "interpretation_boundary": "Official HoloNet CE only for exact HoloNet DB matches; remaining candidates are full NVU spatial LR support, not official CE.",
    }
    (OUT_AUDIT / "run_parameters.json").write_text(json.dumps(params, indent=2), encoding="utf-8")
    rctd = base.load_rctd()
    lr = base.load_lr_candidates()
    interaction_db, cofactor_db, complex_db = get_holonet_db()
    candidates = classify_candidates(lr, interaction_db)
    candidates.to_csv(OUT_TAB / "Table_FULL00_all_NVU_candidates_HoloNet_DB_classification.csv", index=False)

    sample_rows = []
    ce_rows = []
    contrib_rows = []
    failure_rows = []

    for sample, group in SAMPLES.items():
        print(f"[FULL NVU] sample={sample}", flush=True)
        adata, prop = base.prepare_adata(sample, rctd)
        varset = set(map(str.upper, adata.var_names))
        lr_sample = candidates.copy()
        lr_sample["ligand_detected_frac"] = [expr_summary(adata, g)[0] for g in lr_sample["ligand"]]
        lr_sample["receptor_detected_frac"] = [expr_summary(adata, g)[0] for g in lr_sample["receptor"]]
        lr_sample["ligand_in_var"] = lr_sample["ligand"].isin(varset)
        lr_sample["receptor_in_var"] = lr_sample["receptor"].isin(varset)

        official = official_expressed_df(lr_sample, interaction_db)
        official.to_csv(OUT_TAB / f"Table_FULL01_{sample}_official_HoloNet_CE_candidates.csv", index=False)
        ce = None
        if not official.empty:
            try:
                ce, w_best = compute_official_ce(sample, adata, official, cofactor_db, complex_db)
                ce_rows.append({"sample": sample, "group": group, "status": "official_ce_computed", "n_official_pairs": int(ce.shape[0]), "ce_shape": str(tuple(ce.shape)), "w_best": w_best})
            except Exception as e:
                failure_rows.append({"sample": sample, "stage": "official_ce", "message": str(e)})
                ce_rows.append({"sample": sample, "group": group, "status": "official_ce_failed", "n_official_pairs": int(len(official)), "message": str(e)})

        ce_pair_to_index = {p: i for i, p in enumerate(official["LR_Pair"].tolist())} if ce is not None else {}
        for _, row in lr_sample.iterrows():
            pair = row["LR_Pair"]
            lig = row["ligand"]
            rec = row["receptor"]
            src = row["source"]
            tgt = row["target"]
            ce_mean = np.nan
            ce_hotspot_mean = np.nan
            ce_hotspot_max = np.nan
            top_ce_contributor = "NA"
            if ce is not None and pair in ce_pair_to_index:
                idx = ce_pair_to_index[pair]
                score = hotspot_from_ce(ce[idx], adata.n_obs)
                ce_hotspot_mean = float(np.nanmean(score))
                ce_hotspot_max = float(np.nanmax(score))
                ce_mean = float(np.nanmean(ce[idx]))
                net = network_from_ce(ce[idx], prop)
                net_df = pd.DataFrame(net, index=RCTD_TYPES, columns=RCTD_TYPES)
                net_df.to_csv(OUT_TAB / f"Table_FULL02_{sample}_{pair.replace(':','_')}_official_CE_celltype_network.csv")
                contrib = np.nansum(np.abs(net), axis=0) + np.nansum(np.abs(net), axis=1)
                denom = float(np.nansum(contrib)) or 1.0
                top_ce_contributor = RCTD_TYPES[int(np.nanargmax(contrib))] if np.isfinite(contrib).any() else "NA"
                for ct, val in zip(RCTD_TYPES, contrib / denom):
                    contrib_rows.append({"sample": sample, "group": group, "LR_Pair": pair, "celltype": ct, "contribution": float(val)})
            try:
                spatial_score = base.build_spatial_interaction_score(
                    adata=adata, prop=prop, source=str(src), target=str(tgt), ligand=str(lig), receptor=str(rec),
                    ce_score=None,
                )
                spatial_mean = float(np.nanmean(spatial_score))
                spatial_max = float(np.nanmax(spatial_score))
            except Exception as e:
                spatial_mean = np.nan
                spatial_max = np.nan
                failure_rows.append({"sample": sample, "LR_Pair": pair, "stage": "spatial_score", "message": str(e)})
            support_status = "official_CE" if row["holonet_db_exact"] and pair in ce_pair_to_index else "spatial_LR_support"
            if not row["holonet_db_exact"]:
                reason = "not_in_HoloNet_human_LR_database_exact_pair"
            elif pair not in ce_pair_to_index:
                reason = "official_DB_match_but_CE_not_computed_or_filtered"
            else:
                reason = "official_CE_computed"
            sample_rows.append({
                "sample": sample, "group": group, "LR_Pair": pair, "source": src, "target": tgt,
                "ligand": lig, "receptor": rec, "communication_strength": row["communication_strength"],
                "nvu_priority_score": row["nvu_priority_score"], "n_source_target_edges": row.get("n_source_target_edges", np.nan),
                "holonet_db_exact": bool(row["holonet_db_exact"]), "support_status": support_status, "support_reason": reason,
                "ligand_in_var": bool(row["ligand_in_var"]), "receptor_in_var": bool(row["receptor_in_var"]),
                "ligand_detected_frac": row["ligand_detected_frac"], "receptor_detected_frac": row["receptor_detected_frac"],
                "spatial_score_mean": spatial_mean, "spatial_score_max": spatial_max,
                "official_ce_mean": ce_mean, "official_ce_hotspot_mean": ce_hotspot_mean, "official_ce_hotspot_max": ce_hotspot_max,
                "top_ce_contributor": top_ce_contributor,
            })

    sample_scores = pd.DataFrame(sample_rows)
    sample_scores.to_csv(OUT_TAB / "Table_FULL03_sample_level_full_NVU_spatial_LR_scores.csv", index=False)
    pd.DataFrame(ce_rows).to_csv(OUT_TAB / "Table_FULL04_official_CE_run_summary.csv", index=False)
    pd.DataFrame(contrib_rows).to_csv(OUT_TAB / "Table_FULL05_official_CE_celltype_contributions.csv", index=False)
    pd.DataFrame(failure_rows).to_csv(OUT_TAB / "Table_FULL06_failure_and_filter_reasons.csv", index=False)

    axis_summary = (sample_scores.groupby("LR_Pair", as_index=False)
                    .agg(communication_strength_max=("communication_strength", "max"),
                         mean_score=("spatial_score_mean", "mean"),
                         max_score=("spatial_score_max", "max"),
                         holonet_db_exact=("holonet_db_exact", "max"),
                         nvu_priority_score=("nvu_priority_score", "max"),
                         n_source_target_edges=("source", "nunique")))
    deltas = []
    for pair, d in sample_scores.groupby("LR_Pair"):
        cn = d.loc[d["group"].eq("CN"), "spatial_score_mean"].dropna()
        ad = d.loc[d["group"].eq("AD"), "spatial_score_mean"].dropna()
        deltas.append({"LR_Pair": pair, "CN_mean": float(cn.mean()) if len(cn) else np.nan,
                       "AD_mean": float(ad.mean()) if len(ad) else np.nan,
                       "AD_minus_CN_mean": float(ad.mean()-cn.mean()) if len(cn) and len(ad) else np.nan})
    delta_df = pd.DataFrame(deltas)
    axis_summary = axis_summary.merge(delta_df, on="LR_Pair", how="left")
    axis_summary["abs_AD_minus_CN"] = axis_summary["AD_minus_CN_mean"].abs()
    axis_summary = axis_summary.sort_values(["holonet_db_exact", "nvu_priority_score", "mean_score"], ascending=[False, False, False])
    axis_summary.to_csv(OUT_TAB / "Table_FULL07_axis_level_full_NVU_summary_and_highlights.csv", index=False)

    plot_landscape_heatmap(sample_scores)
    plot_delta_lollipop(axis_summary)
    plot_ce_contribution(pd.DataFrame(contrib_rows))
    plot_evidence_bubble(axis_summary)

    audit = f"""# Full NVU HoloNet/spatial LR landscape audit\n\nSeed: {SEED}\n\nThis extension evaluates all LIANA-prioritized NVU ligand-receptor candidates instead of only retained plotted axes. Candidate axes are separated into official HoloNet CE-compatible exact database matches and spatial LR support-only axes.\n\nKey boundary: official HoloNet CE is only claimed for exact HoloNet human LR database matches that successfully computed CE tensors. Other axes remain biologically important LIANA/RCTD/spatial LR candidates and should not be described as official HoloNet CE/FCE results.\n\nThe full MGC/FCE/target Delta-E layer is not claimed in this script. If Figure 7/8 needs a full HoloNet target-response claim, a separate MGC target-gene prediction extension must be run with predefined target genes and saved trained models.\n"""
    (OUT_AUDIT / "Full_NVU_landscape_audit.md").write_text(audit, encoding="utf-8")
    print(OUT, flush=True)


if __name__ == "__main__":
    main()
