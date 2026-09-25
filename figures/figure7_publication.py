#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rebuild Figure 7 v3 as a polished, data-redrawn publication/MANUSCRIPT-style figure.
This version minimizes screenshot panels and redraws scTour, RCTD and NVU LR
summary panels from source data where possible.
"""
from __future__ import annotations

from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

ROOT = Path(r"/path/to/project/results")
MANUSCRIPT = ROOT / "publication_outputs"
OUT = MANUSCRIPT / "main_figures"
OUT.mkdir(parents=True, exist_ok=True)
SCTOUR = ROOT / "node07_sctour/outputs/source_data"
RCTD = ROOT / "08_spatial/tables_rctd/TableS4_spot_level_RCTD_proportions.csv"
SPATIAL_ROOT = Path(r"/path/to/public_data/spatial/GSE220442/counts_and_images")
FULL = ROOT / "09_spatial_communication/full_NVU_landscape_extension"
SCTOUR_FINAL = ROOT / "node07_sctour/result/node07_SCTour_Final"
SCTOUR_FULL_VASC = ROOT / "node07_sctour/result/node07_SCTour_FullVascular"
SCTOUR_MAIN_FIGS = ROOT / "node07_sctour/main_figures"

FIG_W = 17 / 2.54
FIG_H = 24 / 2.54
SEED = 20260527
rng = np.random.default_rng(SEED)

matplotlib.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "font.family": "Arial",
    "font.size": 6,
    "font.weight": "bold",
    "axes.labelsize": 6,
    "axes.titlesize": 6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "axes.linewidth": 0.5,
})
sns.set_theme(context="paper", style="white", font="Arial")
matplotlib.rcParams.update({
    "font.family": "Arial",
    "font.size": 6,
    "font.weight": "bold",
    "axes.labelsize": 6,
    "axes.titlesize": 6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
})

BG = "#FFFFFF"
CARD = "#FFFFFF"
INK = "#25313A"
MUTED = "#66737D"
CN = "#4DBBD5"
AD = "#E64B35"
TEAL = "#25706D"
GOLD = "#D6A642"
PURPLE = "#6D4C7D"
ORANGE = "#CF7B3A"
GREEN = "#4E9A63"

# Manuscript-wide subtype colors aligned with Figure 4.
SCTOUR_PANEL_COLORS = {
    "micro": {"Homeostatic": "#E6AB6F", "DAM": "#D55E00"},
    "astro": {"Homeostatic": "#D8A7B1", "Intermediate": "#B07AA1", "Reactive": "#CC79A7"},
    "endo": {
        "CapEC_ActLow": "#4C78A8", "CapEC_ActMid": "#59A14F",
        "CapEC_ActHigh": "#E15759"
    },
    "vascular": {
        "Arterial": "#CC6677", "Capillary": "#B79F00", "CapEC_ActLow": "#D8C85A",
        "CapEC_ActMid": "#B79F00", "CapEC_ActHigh": "#8F7D00",
        "Pericyte": "#44AA99", "SMC": "#4477AA", "Venous": "#AA4499"
    },
}


def label(ax, s, dx=0.018, dy=0.012):
    """Place panel letters in figure coordinates above-left of the panel."""
    pos = ax.get_position()
    ax.figure.text(pos.x0 - dx, pos.y1 + dy, s, fontsize=12, fontweight="bold",
                   ha="left", va="bottom", color="black", clip_on=False, zorder=20)


def clean(ax):
    sns.despine(ax=ax)
    ax.tick_params(length=2, width=0.45, colors=INK, pad=2.2)
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)


def card_title(ax, title, subtitle=None):
    # MANUSCRIPT-facing figures in this manuscript use minimal in-panel titles.
    return


def enforce_text_style(fig):
    """Final manuscript-wide typography pass: panel letters 12 pt, all else 6 pt, bold."""
    for txt in fig.findobj(matplotlib.text.Text):
        txt.set_fontfamily("Arial")
        txt.set_fontweight("bold")
        if txt.get_text() in {"A", "B", "C", "D", "E", "F"}:
            txt.set_fontsize(12)
        else:
            txt.set_fontsize(6)


def sample_scatter(df, n=4500):
    if len(df) <= n:
        return df.copy()
    return df.sample(n=n, random_state=SEED)


def plot_density_trajectory(ax, d, title):
    d = sample_scatter(d, 5000)
    d = d.sort_values("pseudotime")
    x = d["latent_umap_1"].to_numpy(float)
    y = d["latent_umap_2"].to_numpy(float)
    pt = d["pseudotime"].to_numpy(float)
    sc = ax.scatter(x, y, c=pt, s=1.35, cmap="viridis", linewidths=0, alpha=0.86,
                    vmin=0, vmax=1, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.text(0.03, 0.96, title, transform=ax.transAxes, ha="left", va="top",
            fontsize=6, fontweight="bold", color=INK,
            bbox=dict(boxstyle="round,pad=0.18", facecolor=(1, 1, 1, 0.86), edgecolor="none"))
    return sc


def plot_sctour_vector_panel(ax, key, title):
    """Draw the already computed SCTour latent vector field for Figure 7A.

    Preferred path: replot from saved h5ad through SCTour, matching the node07
    final script's rigorous latent-vector-field representation.
    Fallback: use the previously generated SCTour vector-field PNG if SCTour is
    unavailable in the current Python environment.
    """
    h5ad = SCTOUR_FINAL / key / "h5ad" / f"{key}_sctour_final.h5ad"
    png = SCTOUR_MAIN_FIGS / f"{key}_sctour_vf.png"
    try:
        import scanpy as sc
        import sctour as sct
        adata = sc.read_h5ad(h5ad)
        color = "celltype2"
        if color in adata.obs.columns:
            cats_present = [c for c in SCTOUR_PANEL_COLORS.get(key, {}) if c in set(adata.obs[color].astype(str))]
            if cats_present:
                adata.obs[color] = pd.Categorical(adata.obs[color].astype(str), categories=cats_present, ordered=True)
                adata.uns[f"{color}_colors"] = [SCTOUR_PANEL_COLORS[key][c] for c in cats_present]
        sct.vf.plot_vector_field(
            adata,
            zs_key="X_TNODE",
            vf_key="X_VF",
            use_rep_neigh="X_TNODE",
            color=color,
            ax=ax,
            show=False,
            legend_loc=None,
            size=7,
            alpha=0.58,
            stream_density=0.72,
            stream_color="#555555",
            stream_linewidth=0.24,
            stream_arrowsize=0.45,
        )
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()
        ax.collections[-1].set_rasterized(True)
        # Literature-style on-data subtype labels, similar to boxed labels used
        # in common single-cell/spatial plotting workflows.
        label_col = "celltype2" if "celltype2" in adata.obs.columns else None
        if label_col is not None and "X_umap" in adata.obsm:
            # scTour plots the vector field on the displayed UMAP-like
            # embedding, so labels should stay in X_umap coordinates. The small
            # arterial state is manually nudged below to avoid the venous label.
            emb = np.asarray(adata.obsm["X_umap"])
            cats = adata.obs[label_col].astype(str)
            label_map = {
                "Homeostatic": "Homeostatic",
                "DAM": "DAM",
                "Intermediate": "Intermediate",
                "Reactive": "Reactive",
                "Capillary": "Capillary",
                "Pericyte": "Pericyte",
                "Venous": "Venous",
                "SMC": "SMC",
                "Arterial": "Arterial",
                "CapEC_ActHigh": "High",
                "CapEC_ActMid": "Mid",
                "CapEC_ActLow": "Low",
            }
            for cat in sorted(cats.unique()):
                idx = np.where(cats.values == cat)[0]
                if len(idx) < 20:
                    continue
                x, y = np.nanmedian(emb[idx, 0]), np.nanmedian(emb[idx, 1])
                dx, dy = {
                    "CapEC_ActLow": (-0.42, -0.26),
                    "CapEC_ActMid": (0.04, 0.03),
                    "CapEC_ActHigh": (0.42, 0.28),
                }.get(cat, (0.0, 0.0))
                box_fc = "white"
                if key == "endo":
                    box_fc = SCTOUR_PANEL_COLORS[key].get(cat, "white")
                ax.text(x + dx, y + dy, label_map.get(cat, cat), ha="center", va="center",
                        fontsize=4.8, fontweight="bold",
                        color=("white" if key == "endo" and cat != "CapEC_ActMid" else INK),
                        zorder=30,
                        transform=ax.transData,
                        bbox=dict(boxstyle="round,pad=0.11", fc=box_fc, ec="#6A6A6A",
                                  lw=0.30, alpha=0.78))
    except Exception:
        img = plt.imread(png)
        # Crop only non-informative title/legend margins from the archived
        # source output; do not crop plotted data.
        h, w = img.shape[:2]
        img = img[int(h*0.05):int(h*0.90), int(w*0.02):int(w*0.78)]
        ax.imshow(img)
    ax.set_title("", loc="center")
    ax.set_title("", loc="right")
    ax.set_title(title, loc="left", pad=2, fontsize=6, fontweight="bold", color=INK)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def build_A(fig, spec):
    order = [("micro", "Microglia"), ("astro", "Astrocytes"), ("endo", "CapEC activation")]
    sg = spec.subgridspec(1, 3, width_ratios=[1,1,1], wspace=0.08)
    axes=[]
    for i,(key,title) in enumerate(order):
        ax=fig.add_subplot(sg[0,i])
        plot_sctour_vector_panel(ax, key, title)
        axes.append(ax)
    label(axes[0], "A", dx=0.020, dy=0.014)


def build_B(ax):
    df=pd.read_csv(SCTOUR / "Figure_9D_source_data.csv")
    df["trajectory_label"] = df["trajectory_label"].replace({"Capillaries":"CapEC activation"})
    df = df[["trajectory_label", "Group", "median_pseudotime"]]
    order=["Microglia","Astrocytes","CapEC activation"]
    xpos={k:i for i,k in enumerate(order)}
    pvals={}
    for key in order:
        ad_vals=df[(df["trajectory_label"].eq(key)) & (df["Group"].eq("AD"))]["median_pseudotime"].dropna().to_numpy(float)
        cn_vals=df[(df["trajectory_label"].eq(key)) & (df["Group"].eq("CN"))]["median_pseudotime"].dropna().to_numpy(float)
        try:
            from scipy.stats import mannwhitneyu
            pvals[key]=mannwhitneyu(ad_vals, cn_vals, alternative="two-sided").pvalue
        except Exception:
            pvals[key]=np.nan
        for group, color, off in [("CN", CN, -0.15), ("AD", AD, 0.15)]:
            vals=df[(df["trajectory_label"].eq(key)) & (df["Group"].eq(group))]["median_pseudotime"].dropna().to_numpy(float)
            x=np.full(len(vals), xpos[key]+off)+rng.normal(0,0.028,len(vals))
            ax.scatter(x, vals, s=8.5, c=color, alpha=0.72, edgecolor="white", linewidth=0.2, zorder=2)
            if len(vals): ax.plot([xpos[key]+off-0.09,xpos[key]+off+0.09],[np.median(vals),np.median(vals)], color="black", lw=0.8, zorder=3)
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order)
    ax.set_ylabel("Median pseudotime")
    ax.set_ylim(-0.02,1.05)
    for key in order:
        p=pvals.get(key, np.nan)
        txt = f"P={p:.3f}" if np.isfinite(p) and p >= 0.001 else (f"P={p:.1e}" if np.isfinite(p) else "P=NA")
        ax.text(xpos[key], 1.015, txt, ha="center", va="bottom", fontsize=6, fontweight="bold", color=INK)
    ax.grid(axis="y", lw=0.35, color="#D7D3C8", alpha=0.8)
    card_title(ax, "Sample-level disease-state shift")
    ax.legend(handles=[plt.Line2D([0],[0], marker="o", color="w", markerfacecolor=CN, markersize=4, label="CN"), plt.Line2D([0],[0], marker="o", color="w", markerfacecolor=AD, markersize=4, label="AD")], frameon=False, loc="upper left", bbox_to_anchor=(1.01,1.00), borderaxespad=0, handletextpad=0.25, labelspacing=0.18)
    clean(ax); label(ax,"B", dx=0.026, dy=0.014)


def read_spatial_coords(sample="2-3"):
    # Directly read 10x spatial tissue positions without requiring scanpy.
    spatial_dir = SPATIAL_ROOT / sample / "spatial"
    candidates = [spatial_dir / "tissue_positions.csv", spatial_dir / "tissue_positions_list.csv"]
    p = next((x for x in candidates if x.exists()), None)
    if p is None:
        raise FileNotFoundError(f"No tissue positions found for {sample}")
    raw = pd.read_csv(p, header=None)
    if raw.shape[1] >= 6:
        raw.columns = ["barcode","in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres"] + list(raw.columns[6:])
    coords = raw[["barcode","pxl_col_in_fullres","pxl_row_in_fullres"]].copy()
    return coords


def plot_tissue_mask_enriched_rctd(ax, d, col, title, extent, he_img=None, high_q=0.92):
    """Plot literature-style continuous RCTD proportions over pale H&E.

    This is a display-only change. RCTD values are not smoothed, rescaled, or
    interpolated; every measured tissue spot is colored by its continuous RCTD
    proportion and each panel carries an independent colorbar, matching common
    spatial transcriptomics display practice.
    """
    vals = d[col].to_numpy(float)
    ax.set_facecolor("#FFFFFF")
    if he_img is not None:
        he = he_img[..., :3] if he_img.ndim == 3 else he_img
        if he.ndim == 3:
            he = np.dot(he[..., :3], [0.299, 0.587, 0.114])
        he = he.astype(float)
        if he.max() > 1:
            he = he / 255.0
        he = 1.0 - (1.0 - he) * 0.38
        ax.imshow(he, extent=extent, cmap="gray", vmin=0, vmax=1, alpha=0.60, zorder=0)
    vmax = np.nanpercentile(vals, 98)
    sc = ax.scatter(d["pxl_col_in_fullres"], d["pxl_row_in_fullres"],
                    c=vals, s=1.35, cmap="Spectral_r", vmin=0, vmax=vmax,
                    linewidths=0, alpha=0.92, rasterized=True, zorder=2)
    # Extend the vertical plotting range so the horizontal colorbar and its
    # endpoint labels remain inside each spatial subpanel instead of spilling
    # outside the panel. This lengthens the spatial subpanel itself rather than
    # increasing row spacing.
    y_extra = (extent[2] - extent[3]) * 0.18
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2] + y_extra, extent[3])
    ax.set_aspect("equal")
    ax.set_anchor("N")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.text(0.03, 0.95, title, transform=ax.transAxes, ha="left", va="top",
            fontsize=6, fontweight="bold", color=INK,
            bbox=dict(boxstyle="round,pad=0.10", facecolor=(1, 1, 1, 0.76),
                      edgecolor="none"))
    cax = ax.inset_axes([0.19, 0.060, 0.62, 0.050])
    cb = ax.figure.colorbar(sc, cax=cax, orientation="horizontal")
    cb.set_ticks([0, vmax])
    cb.set_ticklabels([f"{0:.2g}", f"{vmax:.2g}"])
    cb.ax.tick_params(labelsize=5.0, length=1.1, width=0.35, pad=0.7)
    cb.outline.set_linewidth(0.35)


def build_C(fig, spec):
    sample="2-3"
    rctd=pd.read_csv(RCTD)
    coords=read_spatial_coords(sample)
    d=rctd[rctd["sample_id"].astype(str).eq(sample)].merge(coords, left_on="spot_id", right_on="barcode", how="inner")
    spatial_dir = SPATIAL_ROOT / sample / "spatial"
    scales = json.loads((spatial_dir / "scalefactors_json.json").read_text())
    img_path = spatial_dir / "tissue_hires_image.png"
    scale_key = "tissue_hires_scalef"
    if not img_path.exists():
        img_path = spatial_dir / "tissue_lowres_image.png"
        scale_key = "tissue_lowres_scalef"
    img = plt.imread(img_path)
    sf = scales[scale_key]
    img_extent = [0, img.shape[1] / sf, img.shape[0] / sf, 0]
    features=[
        ("Endo_Capillary","Capillary EC"),
        ("Endo_Pericyte","Pericyte"),
        ("Endo_SMC","SMC"),
        ("Endo_Venous","Venous"),
        ("Astro_Reactive","Reactive astro"),
        ("Micro_DAM","DAM microglia"),
    ]
    sg=spec.subgridspec(3,2,width_ratios=[1,1],height_ratios=[1,1,1],wspace=0.08,hspace=0.055)
    axes=[]
    # Original-paper-aligned rendering: continuous all-spot RCTD proportions
    # over pale grayscale H&E.
    high_q = 0.92
    positions=[(0,0),(0,1),(1,0),(1,1),(2,0),(2,1)]
    for i,(col,title) in enumerate(features):
        ax=fig.add_subplot(sg[positions[i][0],positions[i][1]])
        plot_tissue_mask_enriched_rctd(ax, d, col, title, img_extent, he_img=img, high_q=high_q)
        axes.append(ax)
    label(axes[0],"C",dx=0.018,dy=0.014)


def build_D(ax, letter="D"):
    scores=pd.read_csv(FULL / "tables/Table_FULL03_sample_level_full_NVU_spatial_LR_scores.csv")
    axis=pd.read_csv(FULL / "tables/Table_FULL07_axis_level_full_NVU_summary_and_highlights.csv")
    order=["VWF:LRP1","SPARC:ENG","SPARC:FGFR1","VWF:ITGB1","TGM2:SDC4","TGM2:ITGB1","PDGFC:FLT1","FN1:ITGA6","FN1:MAG","ANGPT2:TEK","DLL4:NOTCH3","VEGFC:FLT1","ANGPT2:TIE1"]
    order=[x for x in order if x in set(axis["LR_Pair"])]
    samples=["1-1","18-64","2-5","2-3","2-8","T4857"]
    sample_labels=["CN-1","CN-2","CN-3","AD-1","AD-2","AD-3"]
    piv=scores[scores["LR_Pair"].isin(order)].pivot_table(index="LR_Pair",columns="sample",values="spatial_score_mean",aggfunc="mean").reindex(order)[samples]
    im=ax.imshow(piv.values, aspect="auto", cmap="rocket_r", vmin=0, vmax=np.nanpercentile(piv.values,96))
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order)
    ax.set_xticks(range(len(samples))); ax.set_xticklabels(sample_labels, rotation=35, ha="right", rotation_mode="anchor")
    ax.tick_params(length=0, pad=3)
    ax.axvline(2.5,color="white",lw=1.3)
    ax.text(1,-1.25,"CN",ha="center",va="center",color=CN,fontsize=6,fontweight="bold")
    ax.text(4,-1.25,"AD",ha="center",va="center",color=AD,fontsize=6,fontweight="bold")
    card_title(ax,"Full NVU spatial ligand-receptor landscape")
    cb=plt.colorbar(im, ax=ax, fraction=0.018, pad=0.018); cb.set_label("LR score", labelpad=1)
    label(ax, letter, dx=0.026, dy=0.014)


def build_E(ax, letter="E"):
    axis=pd.read_csv(FULL / "tables/Table_FULL07_axis_level_full_NVU_summary_and_highlights.csv")
    keep=["VWF:LRP1","ANGPT2:TIE1","VEGFC:FLT1","VWF:ITGB1","SPARC:ENG","SPARC:FGFR1","PDGFC:FLT1","TGM2:ITGB1","ANXA2:TLR2","ANGPT2:TEK","DLL4:NOTCH3"]
    d=axis[axis["LR_Pair"].isin(keep)].copy().sort_values("AD_minus_CN_mean")
    y=np.arange(len(d)); vals=d["AD_minus_CN_mean"].to_numpy(float)
    colors=np.where(vals>=0,AD,CN)
    ax.axvline(0,color=INK,lw=0.65)
    ax.hlines(y,0,vals,color="#C7C2B8",lw=0.8)
    ax.scatter(vals,y,s=np.clip(d["mean_score"].to_numpy(float)*360,20,140),c=colors,edgecolor="black",linewidth=0.25,zorder=3)
    official=d["holonet_db_exact"].astype(bool).to_numpy()
    ax.scatter(vals[official],y[official],s=np.clip(d.loc[official,"mean_score"].to_numpy(float)*520,36,175),facecolors="none",edgecolor=GOLD,linewidth=1.0,zorder=4)
    ax.set_yticks(y); ax.set_yticklabels(d["LR_Pair"])
    ax.set_xlabel("AD - CN mean spatial LR score")
    card_title(ax,"Directional NVU communication remodeling")
    ax.grid(axis="x", color="#D8D2C4", lw=0.35, alpha=0.8)
    handles=[
        plt.Line2D([0],[0], marker="o", color="none", markerfacecolor=AD,
                   markeredgecolor="black", markersize=4, label="AD-enriched"),
        plt.Line2D([0],[0], marker="o", color="none", markerfacecolor=CN,
                   markeredgecolor="black", markersize=4, label="CN-enriched"),
    ]
    ax.legend(handles=handles, frameon=False, loc="upper left",
              bbox_to_anchor=(0.02,0.98), borderaxespad=0, handletextpad=0.20,
              labelspacing=0.12, markerscale=0.9)
    clean(ax); label(ax, letter, dx=0.026, dy=0.014)


def build_F(ax, letter="F"):
    contrib=pd.read_csv(FULL / "tables/Table_FULL05_official_CE_celltype_contributions.csv")
    keep_ct=["Endo_Capillary","Endo_Pericyte","Endo_Venous","Endo_SMC","Astro_Reactive","Astro_Homeostatic","Micro_DAM","OPCs"]
    label_map={"Endo_Capillary":"Capillary","Endo_Pericyte":"Pericyte","Endo_Venous":"Venous","Endo_SMC":"SMC","Astro_Reactive":"Reactive astro","Astro_Homeostatic":"Homeostatic astro","Micro_DAM":"DAM micro","OPCs":"OPC"}
    d=contrib[contrib["celltype"].isin(keep_ct)].groupby(["LR_Pair","celltype"],as_index=False)["contribution"].mean()
    xmap={ct:i for i,ct in enumerate(keep_ct)}; yorder=["ANGPT2:TEK","DLL4:NOTCH3"]; ymap={p:i for i,p in enumerate(yorder)}
    sc_last = None
    for _,r in d.iterrows():
        if r["LR_Pair"] not in ymap: continue
        sc_last = ax.scatter(xmap[r["celltype"]], ymap[r["LR_Pair"]], s=900*r["contribution"], c=r["contribution"], cmap="viridis", vmin=d["contribution"].min(), vmax=d["contribution"].max(), edgecolor="black", linewidth=0.25)
    ax.set_xticks(range(len(keep_ct))); ax.set_xticklabels([label_map[x] for x in keep_ct], rotation=35, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(yorder))); ax.set_yticklabels(yorder)
    ax.set_xlim(-0.65,len(keep_ct)-0.35); ax.set_ylim(-0.55,2.15)
    card_title(ax,"Official HoloNet CE contribution", "Bubble size/color encode mean contribution")
    ax.grid(color="#E0DACE", lw=0.35)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.tick_params(length=0)
    if sc_last is not None:
        cax = inset_axes(ax, width="34%", height="5.2%", loc="upper left", borderpad=0.10,
                         bbox_to_anchor=(0.02, 0.00, 1, 1), bbox_transform=ax.transAxes)
        cb = plt.colorbar(sc_last, cax=cax, orientation="horizontal")
        cb.set_ticks([d["contribution"].min(), d["contribution"].max()])
        cb.set_ticklabels(["Low", "High"])
        cb.set_label("")
        cax.set_title("Contribution", fontsize=6, fontweight="bold", pad=1)
        # A small size key next to the colorbar; colorbar and size key are parallel.
        ax.scatter([3.72, 4.13], [2.04, 2.04], s=[35, 110], c="#6CBFA5",
                   edgecolor="black", linewidth=0.25, clip_on=False, zorder=5)
        ax.text(4.42, 2.04, "Size", ha="left", va="center", fontsize=6,
                fontweight="bold", color=INK, clip_on=False)
    label(ax, letter, dx=0.018, dy=0.014)


def main():
    fig=plt.figure(figsize=(FIG_W,FIG_H), dpi=600, facecolor=BG)
    gs=GridSpec(4,2,figure=fig,height_ratios=[0.82,1.30,0.96,1.02],width_ratios=[1.04,0.96],hspace=0.44,wspace=0.48)
    fig.subplots_adjust(left=0.155, right=0.935, top=0.965, bottom=0.075)
    build_A(fig, gs[0,:])
    axb=fig.add_subplot(gs[1,0], facecolor=CARD); build_B(axb)
    build_C(fig, gs[1:3,1])
    axe=fig.add_subplot(gs[2,0], facecolor=CARD); build_E(axe, letter="D")
    axf=fig.add_subplot(gs[3,:], facecolor=CARD); build_F(axf, letter="E")
    enforce_text_style(fig)
    base=OUT/"Figure_7_state_trajectory_spatial_evidence"
    for stem in (base,):
        fig.savefig(stem.with_suffix(".png"),dpi=600,facecolor=BG)
        fig.savefig(stem.with_suffix(".pdf"),facecolor=BG)
        fig.savefig(stem.with_suffix(".svg"),facecolor=BG)
        fig.savefig(stem.with_suffix(".tiff"),dpi=600,facecolor=BG,pil_kwargs={"compression":"tiff_lzw"})
    plt.close(fig)
    supp=plt.figure(figsize=(8.5/2.54, 8.0/2.54), dpi=600, facecolor=BG)
    axd=supp.add_subplot(111, facecolor=CARD)
    build_D(axd, letter="E")
    enforce_text_style(supp)
    supp_base=MANUSCRIPT/"supplementary_figures"/"Supplementary_Figure_11E_spatial_LR_heatmap_for_merge"
    supp.savefig(supp_base.with_suffix(".png"), dpi=600, facecolor=BG)
    supp.savefig(supp_base.with_suffix(".pdf"), facecolor=BG)
    supp.savefig(supp_base.with_suffix(".svg"), facecolor=BG)
    supp.savefig(supp_base.with_suffix(".tiff"), dpi=600, facecolor=BG, pil_kwargs={"compression":"tiff_lzw"})
    plt.close(supp)
    print(base)

if __name__=="__main__":
    main()
