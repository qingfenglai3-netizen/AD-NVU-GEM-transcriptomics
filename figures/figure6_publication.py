from __future__ import annotations

from pathlib import Path
import math
import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyArrowPatch


BASE = Path(r"/path/to/project/results")
DAY = BASE / "06_communication" / "outputs"
MANUSCRIPT = BASE / "publication_outputs"
SRC_OUT = MANUSCRIPT / "source_data"
OUT_BASE = MANUSCRIPT / "main_figures" / "Figure_6_NVU_cell_cell_communication"

ARTERIAL = "Cerebrovascular_Arterial"
CIRCLE_MIN_CELLS = 30
CIRCLE_ORDER = [
    "Excitatory",
    "Inhibitory",
    "Oligodendrocytes",
    "OPCs",
    "Astro_Homeostatic",
    "Astro_Intermediate",
    "Astro_Reactive",
    "Micro_Homeostatic",
    "Micro_DAM",
    "CapEC_ActLow",
    "CapEC_ActMid",
    "CapEC_ActHigh",
    "Cerebrovascular_Pericyte",
    "Cerebrovascular_SMC",
    "Cerebrovascular_Venous",
]
CIRCLE_COLORS = {
    # Manuscript-wide family colors aligned with Figures 2-4.
    "Neuron": "#8A8A8A",
    "Oligo/OPC": "#0072B2",
    "Astrocyte": "#CC79A7",
    "Microglia": "#D55E00",
    "Vascular": "#B79F00",
}

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.1,
        "ytick.labelsize": 6.1,
        "legend.fontsize": 6,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
    }
)


def clean_label(x: str) -> str:
    repl = {
        "CapEC_ActHigh": "Cap high",
        "CapEC_ActMid": "Cap mid",
        "CapEC_ActLow": "Cap low",
        "Cerebrovascular_Pericyte": "Pericyte",
        "Cerebrovascular_SMC": "SMC",
        "Cerebrovascular_Venous": "Venous EC",
        "Astro_Reactive": "Reactive astro",
        "Astro_Intermediate": "Intermediate astro",
        "Astro_Homeostatic": "Homeostatic astro",
        "Micro_DAM": "DAM micro",
        "Micro_Homeostatic": "Homeostatic micro",
        "Oligodendrocytes": "Oligodendrocytes",
        "Excitatory": "Excitatory",
        "Inhibitory": "Inhibitory",
        "OPCs": "OPCs",
    }
    return repl.get(x, x.replace("Cerebrovascular_", "").replace("_", " "))


def family(x: str) -> str:
    if x in {"Excitatory", "Inhibitory"}:
        return "Neuron"
    if x in {"OPCs", "Oligodendrocytes"}:
        return "Oligo/OPC"
    if x.startswith("Astro"):
        return "Astrocyte"
    if x.startswith("Micro"):
        return "Microglia"
    if x.startswith("CapEC") or x.startswith("Cerebrovascular"):
        return "Vascular"
    return "Other"


def short_lr(x: str) -> str:
    return x.replace("_", "/")


def counts(group: str) -> dict[str, int]:
    df = pd.read_csv(DAY / "tables" / f"node06_State_LIANA_cellcounts_after_downsampling_{group}.csv")
    return dict(zip(df["label"], df["n_cells"]))


def circle_counts(group: str) -> dict[str, int]:
    df = pd.read_csv(DAY / "tables" / f"node06_State_LIANA_cellcounts_after_downsampling_{group}.csv")
    df = df[(df["label"] != ARTERIAL) & (df["n_cells"] >= CIRCLE_MIN_CELLS)].copy()
    return dict(zip(df["label"], df["n_cells"]))


def circle_aggregate(group: str) -> pd.DataFrame:
    df = pd.read_csv(DAY / "tables" / f"node06_State_LIANA_aggregate_all_{group}.csv")
    df = df[(df["source"] != ARTERIAL) & (df["target"] != ARTERIAL)].copy()
    ct = circle_counts(group)
    for side in ["source", "target"]:
        df[f"{side}_n"] = df[side].map(ct).fillna(0).astype(int)
    df = df[(df["source_n"] >= CIRCLE_MIN_CELLS) & (df["target_n"] >= CIRCLE_MIN_CELLS)].copy()
    df["source_family"] = df["source"].map(family)
    df.to_csv(SRC_OUT / f"Figure_6A_v4_circle_network_{group}_source_data.csv", index=False)
    return df


def draw_circle_self_loop(ax, p: np.ndarray, color: str, width: float, angle: float) -> None:
    outward = p / (np.linalg.norm(p) + 1e-9)
    tangent = np.array([-outward[1], outward[0]])
    center = p + outward * 0.17
    start = center - tangent * 0.085 - outward * 0.005
    end = center + tangent * 0.085 - outward * 0.005
    rad = 2.15 if angle > 0 else -2.15
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            connectionstyle=f"arc3,rad={rad}",
            arrowstyle="-|>",
            mutation_scale=6.8,
            lw=width,
            color=color,
            alpha=0.58,
            zorder=4,
        )
    )


def draw_circle_panel(ax, df: pd.DataFrame, ct: dict[str, int], title: str) -> None:
    nodes = [n for n in CIRCLE_ORDER if n in set(df["source"]).union(df["target"])]
    angles = np.linspace(math.pi / 2, math.pi / 2 - 2 * math.pi, len(nodes), endpoint=False)
    pos = {n: np.array([math.cos(a), math.sin(a)]) for n, a in zip(nodes, angles)}

    self_edges = df[df["source"] == df["target"]].nlargest(15, "mean_strength")
    cross_edges = df[df["source"] != df["target"]].nlargest(52, "mean_strength")
    edges = pd.concat([self_edges, cross_edges], ignore_index=True)
    wmin, wmax = edges["mean_strength"].min(), edges["mean_strength"].max()

    for r in cross_edges.itertuples():
        p1, p2 = pos[r.source], pos[r.target]
        i, j = nodes.index(r.source), nodes.index(r.target)
        rad = 0.16 + 0.015 * min(abs(i - j), 7)
        if i > j:
            rad = -rad
        width = 0.25 + 2.6 * ((r.mean_strength - wmin) / (wmax - wmin + 1e-9))
        ax.add_patch(
            FancyArrowPatch(
                p1,
                p2,
                connectionstyle=f"arc3,rad={rad}",
                arrowstyle="-|>",
                mutation_scale=5.2,
                lw=width,
                color=CIRCLE_COLORS.get(r.source_family, "0.5"),
                alpha=0.27,
                shrinkA=13,
                shrinkB=13,
                zorder=1,
            )
        )

    for r in self_edges.itertuples():
        width = 0.85 + 3.5 * ((r.mean_strength - wmin) / (wmax - wmin + 1e-9))
        draw_circle_self_loop(ax, pos[r.source], CIRCLE_COLORS.get(r.source_family, "0.5"), width, angles[nodes.index(r.source)])

    vals = np.array([ct.get(n, 0) for n in nodes], dtype=float)
    node_sizes = 80 + 170 * np.sqrt(vals / vals.max())
    for n, size in zip(nodes, node_sizes):
        p = pos[n]
        fam = family(n)
        ax.scatter(p[0], p[1], s=size, color=CIRCLE_COLORS.get(fam, "0.6"), edgecolor="white", linewidth=0.7, zorder=5)
        ha = "left" if p[0] >= 0 else "right"
        va = "bottom" if p[1] >= 0.15 else ("top" if p[1] <= -0.15 else "center")
        ax.text(p[0] * 1.18, p[1] * 1.18, clean_label(n), ha=ha, va=va, fontsize=5.6)

    ax.set_title(title, fontweight="bold", pad=3, fontsize=8)
    ax.set_xlim(-1.42, 1.42)
    ax.set_ylim(-1.48, 1.62)
    ax.set_aspect("equal")
    ax.axis("off")


def robust_edges_from_aggregate() -> pd.DataFrame:
    ad = pd.read_csv(DAY / "tables" / "node06_State_LIANA_aggregate_all_AD.csv")
    cn = pd.read_csv(DAY / "tables" / "node06_State_LIANA_aggregate_all_CN.csv")
    ad_counts, cn_counts = counts("AD"), counts("CN")
    key = ["source", "target"]
    df = ad[key + ["mean_strength", "n_lr"]].merge(
        cn[key + ["mean_strength", "n_lr"]],
        on=key,
        how="outer",
        suffixes=("_AD", "_CN"),
    ).fillna(0)
    df = df[(df["source"] != ARTERIAL) & (df["target"] != ARTERIAL)].copy()
    df["delta_AD_vs_CN"] = df["mean_strength_AD"] - df["mean_strength_CN"]
    for side in ["source", "target"]:
        df[f"{side}_n_AD"] = df[side].map(ad_counts).fillna(0).astype(int)
        df[f"{side}_n_CN"] = df[side].map(cn_counts).fillna(0).astype(int)
    df["min_group_cell_count"] = df[["source_n_AD", "source_n_CN", "target_n_AD", "target_n_CN"]].min(axis=1)
    df = df[df["min_group_cell_count"] >= 100].copy()
    df["source_family"] = df["source"].map(family)
    df["target_family"] = df["target"].map(family)
    df.to_csv(SRC_OUT / "Figure_6A_v3_full_subtype_communication_delta_source_data.csv", index=False)
    return df


def robust_raw_lr() -> pd.DataFrame:
    ad = pd.read_csv(
        DAY / "tables" / "node06_State_LIANA_raw_AD.csv",
        usecols=["source", "target", "ligand", "receptor", "communication_strength"],
    )
    cn = pd.read_csv(
        DAY / "tables" / "node06_State_LIANA_raw_CN.csv",
        usecols=["source", "target", "ligand", "receptor", "communication_strength"],
    )
    ad_counts, cn_counts = counts("AD"), counts("CN")
    key = ["source", "target", "ligand", "receptor"]
    df = ad.merge(cn, on=key, how="outer", suffixes=("_AD", "_CN")).fillna(0)
    df = df[(df["source"] != ARTERIAL) & (df["target"] != ARTERIAL)].copy()
    df["delta_AD_vs_CN"] = df["communication_strength_AD"] - df["communication_strength_CN"]
    for side in ["source", "target"]:
        df[f"{side}_n_AD"] = df[side].map(ad_counts).fillna(0).astype(int)
        df[f"{side}_n_CN"] = df[side].map(cn_counts).fillna(0).astype(int)
    df["min_group_cell_count"] = df[["source_n_AD", "source_n_CN", "target_n_AD", "target_n_CN"]].min(axis=1)
    df = df[df["min_group_cell_count"] >= 100].copy()
    df["source_family"] = df["source"].map(family)
    df["target_family"] = df["target"].map(family)
    df["lr_pair"] = df["ligand"] + " - " + df["receptor"].map(short_lr)
    return df


def build_panel_a() -> pd.DataFrame:
    df = robust_edges_from_aggregate()
    order = [
        "Excitatory",
        "Inhibitory",
        "Oligodendrocytes",
        "OPCs",
        "Astro_Homeostatic",
        "Astro_Intermediate",
        "Astro_Reactive",
        "Micro_Homeostatic",
        "Micro_DAM",
        "CapEC_ActLow",
        "CapEC_ActMid",
        "CapEC_ActHigh",
        "Cerebrovascular_Pericyte",
        "Cerebrovascular_Venous",
    ]
    mat = df.pivot_table(index="source", columns="target", values="delta_AD_vs_CN", fill_value=0)
    mat = mat.reindex(index=order, columns=order).fillna(0)
    mat.to_csv(SRC_OUT / "Figure_6A_v3_full_subtype_heatmap_matrix.csv")
    return mat


def build_panel_b(raw: pd.DataFrame) -> pd.DataFrame:
    df = raw.copy()
    df["is_cross_family"] = df["source_family"] != df["target_family"]
    df = df[df["is_cross_family"]].copy()
    ad_top = df.nlargest(7, "delta_AD_vs_CN")
    cn_top = df.nsmallest(5, "delta_AD_vs_CN")
    out = pd.concat([cn_top, ad_top], ignore_index=True)
    out["direction"] = np.where(out["delta_AD_vs_CN"] > 0, "AD-enriched", "CN-enriched")
    out["label"] = out.apply(lambda r: f"{clean_label(r.source)}>{clean_label(r.target)} | {r.lr_pair}", axis=1)
    out.to_csv(SRC_OUT / "Figure_6B_v3_global_cross_compartment_top_LR.csv", index=False)
    return out


def build_panel_c(raw: pd.DataFrame) -> pd.DataFrame:
    focus = raw[
        ((raw["source_family"].isin(["Vascular", "Astrocyte", "Microglia"])) |
         (raw["target_family"].isin(["Vascular", "Astrocyte", "Microglia"])))
        & (raw["source_family"] != raw["target_family"])
    ].copy()
    ad_top = focus.nlargest(10, "delta_AD_vs_CN")
    ad_top["label"] = ad_top.apply(lambda r: f"{clean_label(r.source)}>{clean_label(r.target)} | {r.lr_pair}", axis=1)
    ad_top.to_csv(SRC_OUT / "Figure_6C_v3_NVU_BBB_focused_AD_enriched_LR.csv", index=False)
    return ad_top


def build_panel_d() -> pd.DataFrame:
    df = pd.read_csv(DAY / "tables" / "node06_CapEC_ActHigh_vs_ActLow_LIANA.csv")
    low_count_boundary = {ARTERIAL, "Cerebrovascular_SMC"}
    df = df[(~df["target"].isin(low_count_boundary)) & (df["group"] == "AD") & (df["direction"] == "ActHigh_up")].copy()
    df["target_clean"] = df["target"].map(clean_label)
    df["lr_pair"] = df["ligand"] + " - " + df["receptor"].map(short_lr)
    df["label"] = df["target_clean"] + " | " + df["lr_pair"]
    df = df.sort_values("delta_high_vs_low", ascending=False).head(10)
    df.to_csv(SRC_OUT / "Figure_6D_v3_activated_capillary_LR.csv", index=False)
    return df


def shorten_term(x: str) -> str:
    x = x.replace("extracellular matrix", "ECM")
    x = x.replace("organization", "org.")
    x = x.replace("regulation of ", "reg. ")
    x = x.replace("positive ", "pos. ")
    x = x.replace("cellular response to", "response to")
    if len(x) > 40:
        return x[:37] + "..."
    return x


def build_panel_e() -> pd.DataFrame:
    df = pd.read_csv(DAY / "tables" / "Table_S31_node06_GO_Reactome_enrichment_all.csv")
    df = df[df["p.adjust"] < 0.05].copy()
    keep = {
        "Cerebrovascular": ["Extracellular matrix organization", "endothelial cell migration", "sprouting angiogenesis", "leukocyte migration"],
        "Astro": ["Extracellular matrix organization", "Collagen formation", "regulation of angiogenesis", "positive regulation of synapse assembly"],
        "Micro": ["regulation of immune effector process", "myeloid leukocyte activation", "positive regulation of cytokine-mediated signaling pathway", "regulation of inflammatory response"],
    }
    rows = []
    for fam, terms in keep.items():
        part = df[df["family"].eq(fam)]
        for term in terms:
            hit = part[part["Description"].str.lower().eq(term.lower())]
            if hit.empty:
                hit = part[part["Description"].str.contains(term, case=False, regex=False, na=False)]
            if not hit.empty:
                rows.append(hit.sort_values(["p.adjust", "Count"], ascending=[True, False]).head(1))
    out = pd.concat(rows, ignore_index=True).drop_duplicates(["family", "Description"])
    out["neglog10_fdr"] = -np.log10(out["p.adjust"].clip(lower=1e-300))
    out["term_short"] = out["Description"].map(shorten_term)
    out.to_csv(SRC_OUT / "Figure_6E_v3_curated_pathway_support.csv", index=False)
    return out


def panel_label(ax, label: str, x: float = -0.12, y: float = 1.08) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")


def wrap_labels(labels, width=29):
    return ["\n".join(textwrap.wrap(str(x), width=width, break_long_words=False)) for x in labels]


def set_wrapped_yticklabels(ax, labels, width=29, fontsize=5.1, linespacing=1.25):
    texts = ax.set_yticklabels(wrap_labels(labels, width), fontsize=fontsize)
    for text in texts:
        text.set_linespacing(linespacing)
    return texts


def main() -> None:
    SRC_OUT.mkdir(parents=True, exist_ok=True)
    mat = build_panel_a()
    circle_cn = circle_aggregate("CN")
    circle_ad = circle_aggregate("AD")
    circle_cn_counts = circle_counts("CN")
    circle_ad_counts = circle_counts("AD")
    raw = robust_raw_lr()
    raw.to_csv(SRC_OUT / "Figure_6_v3_full_filtered_LIANA_LR_source_data.csv", index=False)
    top_lr = build_panel_b(raw)
    nvu_lr = build_panel_c(raw)
    cap_lr = build_panel_d()
    pathways = build_panel_e()

    fig = plt.figure(figsize=(6.69, 7.55))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.05, 1.16, 1.0], width_ratios=[1.06, 1.0], hspace=0.42, wspace=0.70)
    fig.subplots_adjust(left=0.22, right=0.965, top=0.965, bottom=0.08)

    ax_a = fig.add_subplot(gs[0, :])
    ax_a.axis("off")
    ax_cn = ax_a.inset_axes([-0.01, 0.03, 0.48, 0.91])
    ax_ad = ax_a.inset_axes([0.50, 0.03, 0.48, 0.91])
    draw_circle_panel(ax_cn, circle_cn, circle_cn_counts, "CN")
    draw_circle_panel(ax_ad, circle_ad, circle_ad_counts, "AD")
    handles = [
        mpl.lines.Line2D([0], [0], marker="o", color="w", label=k, markerfacecolor=v, markersize=4.8)
        for k, v in CIRCLE_COLORS.items()
    ]
    ax_a.legend(handles=handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.10), fontsize=5.7)
    panel_label(ax_a, "A", x=-0.08, y=1.06)

    ax_b = fig.add_subplot(gs[1, 0])
    plot = top_lr.iloc[::-1].copy()
    colors = plot["direction"].map({"AD-enriched": "#d95f5f", "CN-enriched": "#4f83b1"})
    ax_b.barh(range(len(plot)), plot["delta_AD_vs_CN"], color=colors, edgecolor="0.25", linewidth=0.35)
    ax_b.axvline(0, color="0.25", lw=0.7)
    ax_b.set_yticks(range(len(plot)))
    set_wrapped_yticklabels(ax_b, plot["label"], width=30, fontsize=5.0, linespacing=1.28)
    ax_b.set_xlabel("Delta LR strength (AD-CN)")
    ax_b.grid(axis="x", color="0.90", lw=0.5)
    panel_label(ax_b, "B")

    ax_c = fig.add_subplot(gs[1, 1])
    plot = nvu_lr.iloc[::-1].copy()
    ax_c.barh(range(len(plot)), plot["delta_AD_vs_CN"], color="#d95f5f", edgecolor="0.25", linewidth=0.35)
    ax_c.set_yticks(range(len(plot)))
    set_wrapped_yticklabels(ax_c, plot["label"], width=29, fontsize=4.9, linespacing=1.28)
    ax_c.set_xlabel("Delta LR strength (AD-CN)")
    ax_c.grid(axis="x", color="0.90", lw=0.5)
    panel_label(ax_c, "C")

    ax_d = fig.add_subplot(gs[2, 0])
    plot = cap_lr.iloc[::-1].copy()
    ax_d.barh(range(len(plot)), plot["delta_high_vs_low"], color="#8ab6d6", edgecolor="0.25", linewidth=0.35)
    ax_d.set_yticks(range(len(plot)))
    set_wrapped_yticklabels(ax_d, plot["label"], width=30, fontsize=5.1, linespacing=1.18)
    ax_d.set_xlabel("ActHigh minus ActLow LR strength")
    ax_d.grid(axis="x", color="0.90", lw=0.5)
    panel_label(ax_d, "D")

    ax_e = fig.add_subplot(gs[2, 1])
    fam_order = ["Cerebrovascular", "Astro", "Micro"]
    terms = list(dict.fromkeys(pathways["term_short"]))
    pathways["x"] = pathways["family"].map({f: i for i, f in enumerate(fam_order)})
    pathways["y"] = pathways["term_short"].map({t: i for i, t in enumerate(reversed(terms))})
    sc = ax_e.scatter(
        pathways["x"],
        pathways["y"],
        s=18 + pathways["Count"] * 2.2,
        c=pathways["neglog10_fdr"],
        cmap="plasma",
        edgecolor="0.25",
        linewidth=0.3,
    )
    ax_e.set_xticks(range(len(fam_order)))
    ax_e.set_xticklabels(["Vascular", "Astrocyte", "Microglia"], rotation=15, ha="right")
    ax_e.set_yticks(range(len(terms)))
    set_wrapped_yticklabels(ax_e, list(reversed(terms)), width=24, fontsize=5.0, linespacing=1.12)
    ax_e.grid(color="0.90", lw=0.5)
    cb2 = fig.colorbar(sc, ax=ax_e, fraction=0.045, pad=0.02)
    cb2.set_label("-log10 FDR")
    panel_label(ax_e, "E")

    for ext in ["png", "pdf", "svg", "tiff"]:
        kwargs = {"facecolor": "white"}
        if ext in {"png", "tiff"}:
            kwargs["dpi"] = 600
        fig.savefig(f"{OUT_BASE}.{ext}", **kwargs)
    print(OUT_BASE)


if __name__ == "__main__":
    main()
