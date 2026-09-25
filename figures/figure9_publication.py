#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 9: drug-target mapping, drug and pathway enrichment, and module classification."""

from __future__ import annotations
import sys
from pathlib import Path


import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch


PROJECT_ROOT = Path(r"<ANALYSIS_ROOT>")
SCREEN = PROJECT_ROOT / "fuxian_test" / "day6_GEM_method_enhancement_screen"
TRANS = PROJECT_ROOT / "fuxian_test" / "day6_GEM_translational"
JTM = PROJECT_ROOT / "fuxian_test" / "JTM_manuscript"

OUT = Path(r"<WORKDIR>/figures/Fig9")
DIRS = {
    "fig": OUT / "figures",
    "table": OUT / "tables",
    "main": OUT / "main",
    "supp": OUT / "supp",
}
for p in DIRS.values():
    p.mkdir(parents=True, exist_ok=True)

COL = {
    "gem": "#E64B35",
    "vascular": "#4DBBD5",
    "ion": "#00A087",
    "phospho": "#3C5488",
    "drug": "#F39B7F",
    "axis": "#8491B4",
    "grey": "#D8DCE2",
    "dark": "#2B2B2B",
    "line": "#7E8AA2",
}


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 6,
            "axes.labelsize": 6,
            "axes.titlesize": 6,
            "xtick.labelsize": 5.6,
            "ytick.labelsize": 5.6,
            "legend.fontsize": 5.6,
            "axes.linewidth": 0.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    sns.set_theme(context="paper", style="white", font="Arial")


def save_all(fig: plt.Figure, path: Path) -> None:
    for ext in (".png", ".pdf", ".svg", ".tiff"):
        kwargs = {"facecolor": "white", "transparent": False}
        if ext in (".png", ".tiff"):
            kwargs["dpi"] = 600
        if ext == ".tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        if ext == '.tiff':
            import io
            from PIL import Image
            buffer = io.BytesIO()
            fig.savefig(buffer, format='png', dpi=600, facecolor='white', transparent=False)
            buffer.seek(0)
            with Image.open(buffer) as rendered:
                rendered.convert('RGB').save(path.with_suffix(ext), compression='tiff_lzw', dpi=(600,600))
        else:
            fig.savefig(path.with_suffix(ext), **kwargs)
    plt.close(fig)


def clean_term(term: str) -> str:
    raw = str(term)
    low = raw.lower()
    if "g protein-coupled glutamate" in low:
        return "GPCR glutamate receptor signaling"
    if "glutamate receptor signaling" in low:
        return "Glutamate receptor signaling"
    if "synaptic transmission" in low:
        return "Synaptic transmission regulation"
    if "class c/3" in low or "metabotropic glutamate" in low:
        return "Metabotropic glutamate receptors"
    if "gpcr downstream" in low:
        return "GPCR downstream signaling"
    if "signaling by gpcr" in low:
        return "Signaling by GPCR"
    if "pdgf" in low:
        return "PDGF signaling"
    if "inositol phosphate" in low:
        return "Inositol phosphate metabolism"
    if "diacylglycerol" in low:
        return "Diacylglycerol metabolism"
    if "adenylate cyclase" in low:
        return "Adenylate cyclase regulation"
    if "electrical coupling" in low:
        return "Electrical coupling"
    if "calcium ion transport" in low or "calcium channel" in low:
        return "Calcium-channel target set"
    if "quisqualate" in low:
        return "Quisqualate target set"
    if "cediranib" in low:
        return "Cediranib perturbation"
    if "nintedanib" in low:
        return "Nintedanib target set"
    if "nimodipine" in low:
        return "Nimodipine target set"
    if "nisoldipine" in low:
        return "Nisoldipine target set"
    if "felodipine" in low:
        return "Felodipine target set"
    if "nilotinib" in low:
        return "Nilotinib target set"
    if "vatalanib" in low:
        return "Vatalanib target set"
    out = raw.replace(" R-HSA-", " ").replace("_", " ").strip()
    return out if len(out) <= 38 else out[:35] + "..."


def short_drug(name: str) -> str:
    x = str(name).upper()
    repl = {
        "VASCULAR WALL PDGF-NOTCH AXIS": "PDGF-NOTCH axis",
        "PHOSPHOINOSITIDE-CGMP AXIS": "PI-cGMP axis",
        "ION-GPCR AXIS": "Ion-GPCR axis",
        "FLAVOXATE HYDROCHLORIDE": "FLAVOXATE",
        "DIPYRIDAMOLE": "DIPYRIDAMOLE",
        "PENTOXIFYLLINE": "PENTOXIFYLLINE",
        "PARATHYROID HORMONE": "PTH",
        "TERIPARATIDE ACETATE": "TERIPARATIDE",
    }
    return repl.get(x, x.title() if len(x) > 12 else x)


def z01(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce").fillna(0)
    lo, hi = x.min(), x.max()
    if math.isclose(float(lo), float(hi)):
        return pd.Series(0.0, index=x.index)
    return (x - lo) / (hi - lo)


def load_inputs() -> dict[str, pd.DataFrame]:
    data = {
        "edges": pd.read_csv(TRANS / "tables" / "TableS2_day6_GEM_drug_target_edges.csv"),
        "drug_priority": pd.read_csv(TRANS / "tables" / "TableS1_day6_GEM_drug2cell_priority.csv"),
        "gene_weights": pd.read_csv(TRANS / "tables" / "TableS3_day6_GEM_gene_druggability_weights.csv"),
        "ml": pd.read_csv(Path(r"<WORKDIR>/results/ML_labelblind_OOF_summary.csv")),
        "enrichr": pd.read_csv(SCREEN / "tables" / "GEM_Enrichr_drug_perturbagen_pathway_results.csv"),
    }
    edges=data["edges"].copy()
    data["dock"]=edges.groupby("gene",as_index=False).agg(n_drug_sets=("drug_name","nunique"),example_drugs=("drug_name",lambda x:";".join(dict.fromkeys(x))),max_edge_weight=("edge_weight","max"))
    data["dock"].to_csv(OUT/"tables/Figure9A_mapping_counts.csv",index=False)
    return data


def draw_panel_mapping(ax, dock, vina=None):
    fixed_display=["CACNA1C","PDGFRB","PTH1R","DGKB","PDE7B","NOTCH3","GRM3","GRM8"]
    keep = dock[dock.gene.isin(fixed_display)].sort_values(["n_drug_sets", "gene"], ascending=[False,True]).iloc[::-1].reset_index(drop=True)
    y=np.arange(len(keep))
    ax.barh(y,keep.n_drug_sets,color=COL["vascular"],height=.58)
    for i,r in keep.iterrows():
        examples={"CACNA1C":"verapamil / nimodipine","PDGFRB":"imatinib / nintedanib","PTH1R":"teriparatide","DGKB":"PI-cGMP axis","PDE7B":"dipyridamole","NOTCH3":"PDGF-NOTCH axis","GRM3":"Ion-GPCR axis","GRM8":"Ion-GPCR axis"}[r.gene]
        ax.text(float(r.n_drug_sets)+1.35,i,examples,va="center",fontsize=6.2,fontweight="bold",color=COL["dark"])
    ax.set_yticks(y); ax.set_yticklabels(keep.gene,fontsize=7.6,fontweight="bold")
    ax.set_xlabel("Mapped drug or mechanism sets",fontsize=7.4)
    ax.tick_params(axis="x",labelsize=7.2)
    ax.set_xlim(0,max(keep.n_drug_sets.max()+28,85))
    ax.grid(axis="x",color="#EAEAEA",lw=.5)


def select_enrichment(enrichr: pd.DataFrame) -> pd.DataFrame:
    sig = enrichr[enrichr["Adjusted P-value"] < 0.05].copy()
    wanted = sig[
        sig["library"].isin(["GO_Biological_Process_2023", "Reactome_2022", "DGIdb_Drug_Targets_2024", "DSigDB", "Drug_Perturbations_from_GEO_down"])
    ].copy()
    preferred_mech = [
        "glutamate",
        "gpcr",
        "pdgf",
        "inositol",
        "diacylglycerol",
        "adenylate cyclase",
        "electrical coupling",
    ]
    preferred_drug = [
        "quisqualate",
        "cediranib",
        "nintedanib",
        "nimodipine",
        "nisoldipine",
        "felodipine",
        "nilotinib",
        "vatalanib",
    ]
    mech_all = wanted[wanted["library"].isin(["GO_Biological_Process_2023", "Reactome_2022"])].copy()
    drug_all = wanted[wanted["library"].isin(["DGIdb_Drug_Targets_2024", "DSigDB", "Drug_Perturbations_from_GEO_down"])].copy()
    mech_all["pref"] = mech_all["Term"].str.lower().apply(lambda x: min([i for i, k in enumerate(preferred_mech) if k in x], default=99))
    drug_all["pref"] = drug_all["Term"].str.lower().apply(lambda x: min([i for i, k in enumerate(preferred_drug) if k in x], default=99))
    mech = mech_all[mech_all["pref"] < 99].sort_values(["pref", "Adjusted P-value"]).head(7)
    drug = drug_all[drug_all["pref"] < 99].sort_values(["pref", "Adjusted P-value"]).head(7)
    out = pd.concat([mech, drug], ignore_index=True)
    out["display"] = out["Term"].map(clean_term)
    out["neglog10_fdr"] = -np.log10(out["Adjusted P-value"].clip(lower=1e-300))
    out["category"] = np.where(out["library"].isin(["GO_Biological_Process_2023", "Reactome_2022"]), "Pathway", "Drug target")
    out.to_csv(DIRS["table"] / "Figure9C_GEM_drug_pathway_enrichment.csv", index=False)
    return out


def draw_panel_enrich(ax: plt.Axes, enr: pd.DataFrame) -> None:
    enr = enr.sort_values(["category", "Adjusted P-value"], ascending=[True, False]).reset_index(drop=True)
    y = np.arange(len(enr))
    colors = np.where(enr["category"].eq("Pathway"), COL["ion"], COL["drug"])
    ax.barh(y, enr["neglog10_fdr"], color=colors, height=0.62)
    labels = []
    for _, r in enr.iterrows():
        labels.append(r["display"])
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.4, fontweight="bold")
    ax.set_xlabel("-log10(FDR)", fontsize=7.4, fontweight="bold")
    ax.tick_params(axis="x", labelsize=7.2)
    for label in ax.get_xticklabels():
        label.set_fontweight("bold")
    ax.set_xlim(0, max(float(enr["neglog10_fdr"].max()) + 1.65, 3.65))
    ax.grid(axis="x", color="#EAEAEA", lw=0.5)
    legend = [
        Line2D([0], [0], marker="s", color="none", markerfacecolor=COL["ion"], markeredgecolor="none", markersize=5, label="Pathway"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=COL["drug"], markeredgecolor="none", markersize=5, label="Drug target"),
    ]
    ax.legend(
        handles=legend,
        loc="lower right",
        bbox_to_anchor=(1.02, 1.03),
        ncol=2,
        frameon=False,
        handlelength=0.8,
        columnspacing=0.7,
        labelspacing=0.45,
        borderaxespad=0,
        prop={"family": "Arial", "size": 7.2, "weight": "bold"},
    )


def ml_pivot(ml: pd.DataFrame) -> pd.DataFrame:
    tab = ml.pivot_table(index="pipeline", columns="cv", values="auc", aggfunc="max")
    if "10fold" in tab.columns:
        tab["mean_AUC"] = tab["10fold"]
    else:
        tab["mean_AUC"] = tab.mean(axis=1, skipna=True)
    tab = tab.sort_values("mean_AUC", ascending=False)
    return tab


def draw_panel_ml_heatmap(ax: plt.Axes, ml: pd.DataFrame, top_n: int = 36) -> pd.DataFrame:
    selector_order = ["all", "top3", "top5", "top8", "top12", "top16"]
    model_order = [
        "Logistic_L2",
        "Logistic_L1",
        "ElasticNet_logistic",
        "Linear_SVM",
        "RBF_SVM",
        "RandomForest",
        "ExtraTrees",
        "GradientBoosting",
        "DecisionTree",
        "KNN_3",
        "KNN_5",
        "GaussianNB",
    ]
    show = ml.pivot_table(index="model", columns="selector", values="auc", aggfunc="max")
    show = show.reindex(index=[m for m in model_order if m in show.index], columns=[s for s in selector_order if s in show.columns])
    show.columns = ["All genes", "Top 3", "Top 5", "Top 8", "Top 12", "Top 16"][: len(show.columns)]
    sns.heatmap(
        show,
        ax=ax,
        cmap="RdYlBu_r",
        vmin=0.55,
        vmax=0.92,
        cbar=False,
        linewidths=0.2,
        linecolor="white",
        annot=True,
        fmt=".3f",
        annot_kws={"fontsize": 5.4, "fontfamily": "Arial", "fontweight": "bold"},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=35, labelsize=7.6)
    ax.tick_params(axis="y", labelsize=7.6)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
    return show


def make_main_figure(data: dict[str, pd.DataFrame]) -> None:
    set_style()
    enr = select_enrichment(data["enrichr"])
    fig = plt.figure(figsize=(7.5, 7.2))
    ax_a = fig.add_axes([0.16, 0.70, 0.80, 0.26])
    ax_b = fig.add_axes([0.29, 0.13, 0.17, 0.43])
    ax_c = fig.add_axes([0.70, 0.18, 0.28, 0.38])

    draw_panel_mapping(ax_a, data["dock"], data.get("vina"))
    draw_panel_enrich(ax_b, enr)
    full_ml = draw_panel_ml_heatmap(ax_c, data["ml"], top_n=36)

    cax = fig.add_axes([0.74, 0.060, 0.20, 0.015])
    cb = fig.colorbar(ax_c.collections[0], cax=cax, orientation="horizontal")
    cb.set_label("AUC", labelpad=1)
    cb.set_ticks([0.6,0.7,0.8,0.9])
    from matplotlib.text import Text
    for t in fig.findobj(match=Text):
        t.set_fontsize(8)
    for letter, x, y in [('A', 0.015, 0.99), ('B', 0.015, 0.59), ('C', 0.51, 0.59)]:
        fig.text(x, y, letter, fontsize=12, fontweight='bold', ha='left', va='top')

    save_all(fig, DIRS["fig"] / "Figure9_GEM_drugtarget_first_main_ABCD")
    for ext in (".png", ".pdf", ".svg", ".tiff"):
        src = DIRS["fig"] / f"Figure9_GEM_drugtarget_first_main_ABCD{ext}"
        dst = DIRS["main"] / f"Figure9_GEM_drugtarget_first_main_ABCD{ext}"
        dst.write_bytes(src.read_bytes())

    full_ml.to_csv(DIRS["table"] / "Figure9D_all_model_pipeline_auc_matrix.csv")


def make_supp_ml_heatmap(ml: pd.DataFrame) -> None:
    set_style()
    selector_order = ["all", "top3", "top5", "top8", "top12", "top16"]
    model_order = [
        "Logistic_L2",
        "Logistic_L1",
        "ElasticNet_logistic",
        "Linear_SVM",
        "RBF_SVM",
        "RandomForest",
        "ExtraTrees",
        "GradientBoosting",
        "DecisionTree",
        "KNN_3",
        "KNN_5",
        "GaussianNB",
    ]
    tab = ml.pivot_table(index="model", columns="selector", values="auc", aggfunc="max")
    tab = tab.reindex(index=model_order, columns=selector_order)
    tab.columns = ["All genes", "Top 3", "Top 5", "Top 8", "Top 12", "Top 16"]
    fig, ax = plt.subplots(figsize=(4.9, 3.5))
    fig.subplots_adjust(left=0.28, right=0.88, top=0.92, bottom=0.18)
    sns.heatmap(
        tab,
        ax=ax,
        cmap="RdYlBu_r",
        vmin=0.55,
        vmax=0.92,
        cbar_kws={"label": "AUC", "shrink": 0.38, "pad": 0.015},
        linewidths=0.12,
        linecolor="white",
        annot=True,
        fmt=".3f",
        annot_kws={"fontsize": 5.0, "fontfamily": "Arial"},
        yticklabels=True,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=35, labelsize=5.8)
    ax.tick_params(axis="y", labelsize=5.6)
    ax.collections[0].colorbar.ax.tick_params(labelsize=5.2)
    ax.collections[0].colorbar.set_label("AUC", fontsize=5.6)
    fig.text(0.045, 0.985, "A", fontsize=12, fontweight="bold", ha="left", va="top")
    save_all(fig, DIRS["fig"] / "Supplementary_Figure9_all_GEM_model_AUC_heatmap")
    for ext in (".png", ".pdf", ".svg", ".tiff"):
        src = DIRS["fig"] / f"Supplementary_Figure9_all_GEM_model_AUC_heatmap{ext}"
        dst = DIRS["supp"] / f"Supplementary_Figure9_all_GEM_model_AUC_heatmap{ext}"
        dst.write_bytes(src.read_bytes())


def main():
    data=load_inputs()
    data["vina"]=None
    make_main_figure(data)
if __name__=="__main__": main()
