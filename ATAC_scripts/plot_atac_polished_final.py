#!/usr/bin/env python3
from pathlib import Path
import textwrap

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


ROOT = Path("/mnt/d/lxk/project/lishan-20260613")
ATAC = ROOT / "analysis/ATAC"
FIG = ROOT / "fig" / "ATAC_final_polished"
FIG.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.linewidth": 0.55,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
})

COLORS = {
    "Promoter": "#4C78A8",
    "UTR": "#59A14F",
    "Exon": "#E15759",
    "Intron": "#F28E2B",
    "Downstream": "#B07AA1",
    "Intergenic": "#EDC948",
    "Other": "#BAB0AC",
    "consensus": "#6B7280",
    "PE_open": "#B22222",
    "PE_closed": "#2B6CB0",
    "ATAC": "#B22222",
    "scRNA": "#4C78A8",
}


def wrap(s, width=45):
    return "\n".join(textwrap.wrap(str(s).replace("_", " "), width=width))


def simplify_annotation(x):
    x = str(x)
    if x.startswith("Promoter"):
        return "Promoter"
    if "UTR" in x:
        return "UTR"
    if "Exon" in x:
        return "Exon"
    if "Intron" in x:
        return "Intron"
    if x.startswith("Downstream"):
        return "Downstream"
    if "Intergenic" in x:
        return "Intergenic"
    return "Other"


def parse_ratio(x):
    if isinstance(x, str) and "/" in x:
        a, b = x.split("/")[:2]
        return float(a) / float(b)
    return float(x)


anno = pd.read_csv(ROOT / "fig/ATAC_windowsR/ChIPseeker_all_peak_annotations.csv")
anno["simple"] = anno["annotation"].map(simplify_annotation)
anno["set"] = pd.Categorical(anno["set"], ["consensus", "PE_open", "PE_closed"], ordered=True)

linked = pd.read_csv(ATAC / "annotation/atac_scrna_linked_gene_summary.tsv", sep="\t")
linked_genes = set(linked["gene"])
linked_cells = linked["celltype"].drop_duplicates().tolist()

kegg = pd.read_csv(ROOT / "fig/ATAC_windowsR/clusterProfiler_KEGG_enrichment.csv")
kegg["GeneRatioValue"] = kegg["GeneRatio"].map(parse_ratio)
kegg["mlog10padj"] = -np.log10(kegg["p.adjust"].clip(lower=1e-300))
kegg_top = kegg.sort_values("p.adjust").head(14).copy()

scrna_rows = []
for cell in linked_cells:
    f = ROOT / "DEG" / cell / "GSEA.csv"
    if not f.exists():
        continue
    df = pd.read_csv(f, encoding="utf-8-sig")
    if not {"Description", "NES", "p.adjust"}.issubset(df.columns):
        continue
    df["celltype"] = cell
    df["mlog10padj"] = -np.log10(pd.to_numeric(df["p.adjust"], errors="coerce").clip(lower=1e-300))
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    df = df.dropna(subset=["NES", "mlog10padj"])
    key = "INFLAM|TNFA|NFKB|HYPOXIA|ANGIO|T_CELL|LEUKOCYTE|MIGRATION|ADHESION|DEVELOPMENT|EMT|VASCULAR|CYTOKINE|INTERFERON|COMPLEMENT|METABOL"
    rel = df[df["Description"].str.contains(key, case=False, regex=True, na=False)]
    if rel.empty:
        rel = df
    scrna_rows.append(rel.sort_values("p.adjust").head(4))
scrna = pd.concat(scrna_rows, ignore_index=True) if scrna_rows else pd.DataFrame()
scrna = scrna.sort_values(["celltype", "p.adjust"]).head(22).copy()


def save(fig, name, pages):
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    pages.savefig(fig, bbox_inches="tight")
    plt.close(fig)


with PdfPages(FIG / "ATAC_final_polished_all.pdf") as pages:
    # 1 Polished annotation composition
    order = ["Promoter", "UTR", "Exon", "Intron", "Downstream", "Intergenic", "Other"]
    count_tab = anno.groupby(["set", "simple"], observed=False).size().unstack(fill_value=0)
    count_tab = count_tab.reindex(index=["consensus", "PE_open", "PE_closed"], columns=order, fill_value=0)
    pct_tab = count_tab.div(count_tab.sum(axis=1), axis=0) * 100
    fig, ax = plt.subplots(figsize=(7.2, 3.25), constrained_layout=True)
    left = np.zeros(3)
    sets = ["consensus", "PE_open", "PE_closed"]
    y = np.arange(len(sets))
    for cat in order:
        vals = pct_tab[cat].reindex(sets).fillna(0).values
        ax.barh(y, vals, left=left, color=COLORS[cat], edgecolor="white", linewidth=0.4, label=cat)
        for yi, lft, val in zip(y, left, vals):
            if val >= 8:
                ax.text(lft + val / 2, yi, f"{val:.0f}%", ha="center", va="center", color="white", fontsize=7)
        left += np.array(vals)
    ax.set_yticks(y, ["Consensus peaks", "PE-open candidates", "PE-closed candidates"])
    ax.set_xlabel("Peak fraction (%)")
    ax.set_xlim(0, 100)
    ax.legend(ncol=7, loc="upper center", bbox_to_anchor=(0.5, -0.22), frameon=False, columnspacing=0.9)
    ax.set_title("Genomic annotation of ATAC peak sets", fontweight="bold")
    save(fig, "01_peak_annotation_composition_polished", pages)

    # 2 Donut grid with clean legend
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.8), constrained_layout=True)
    for ax, st in zip(axes, sets):
        vals = [len(anno[(anno["set"] == st) & (anno["simple"] == cat)]) for cat in order]
        wedges, _ = ax.pie(vals, colors=[COLORS[c] for c in order], startangle=90,
                           wedgeprops={"width": 0.42, "edgecolor": "white", "linewidth": 0.5})
        ax.text(0, 0, st.replace("_", "\n"), ha="center", va="center", fontweight="bold", fontsize=8)
        ax.set_title(f"n={sum(vals):,}", fontsize=8)
    fig.legend(wedges, order, ncol=7, loc="lower center", bbox_to_anchor=(0.5, -0.02), frameon=False)
    fig.suptitle("Peak annotation overview", fontweight="bold")
    save(fig, "02_peak_annotation_donut_grid_polished", pages)

    # 3 Distance to TSS clean distribution
    fig, ax = plt.subplots(figsize=(6.6, 3.0), constrained_layout=True)
    bins = [-1e9, -100000, -10000, -3000, -1000, 0, 1000, 3000, 10000, 100000, 1e9]
    labels = [">100 kb upstream", "10-100 kb up", "3-10 kb up", "1-3 kb up", "0-1 kb up",
              "0-1 kb down", "1-3 kb down", "3-10 kb down", "10-100 kb down", ">100 kb downstream"]
    dist_tab = []
    for st in sets:
        sub = anno[anno["set"] == st].copy()
        sub["bin"] = pd.cut(sub["distanceToTSS"], bins=bins, labels=labels)
        fr = sub["bin"].value_counts(normalize=True).reindex(labels).fillna(0) * 100
        for lab, val in fr.items():
            dist_tab.append({"set": st, "bin": lab, "pct": val})
    dt = pd.DataFrame(dist_tab)
    x = np.arange(len(labels))
    width = 0.24
    for i, st in enumerate(sets):
        vals = dt[dt["set"] == st]["pct"].values
        ax.bar(x + (i - 1) * width, vals, width=width, color=COLORS[st], label=st)
    ax.axvline(4.5, color="#111827", lw=0.6)
    ax.set_xticks(x, [wrap(l, 11) for l in labels], rotation=35, ha="right")
    ax.set_ylabel("Peak fraction (%)")
    ax.legend(frameon=False, ncol=3)
    ax.set_title("Distance from peaks to nearest TSS", fontweight="bold")
    save(fig, "03_distance_to_TSS_polished", pages)

    # 4 KEGG: wide, wrapped, no overlap
    kt = kegg_top.sort_values("GeneRatioValue")
    fig, ax = plt.subplots(figsize=(7.0, 4.8), constrained_layout=True)
    y = np.arange(len(kt))
    sc = ax.scatter(kt["GeneRatioValue"], y, s=np.clip(kt["Count"] * 6, 35, 220),
                    c=kt["mlog10padj"], cmap="Reds", edgecolor="#374151", linewidth=0.35)
    ax.set_yticks(y, [wrap(x, 44) for x in kt["Description"]])
    ax.set_xlabel("GeneRatio")
    ax.set_title("ATAC-linked genes: KEGG enrichment", fontweight="bold")
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("-log10(adj. P)")
    for size in [10, 20, 30]:
        ax.scatter([], [], s=size * 6, c="white", edgecolor="#374151", label=f"{size} genes")
    ax.legend(title="Count", loc="lower right", frameon=False)
    save(fig, "04_ATAC_linked_KEGG_dotplot_polished", pages)

    # 5 scRNA GSEA linked cell types
    if not scrna.empty:
        scrna = scrna.copy()
        scrna["label"] = scrna["celltype"] + " | " + scrna["Description"]
        scrna = scrna.sort_values("NES")
        fig, ax = plt.subplots(figsize=(7.4, 5.6), constrained_layout=True)
        y = np.arange(len(scrna))
        sc = ax.scatter(scrna["NES"], y, s=np.clip(scrna["mlog10padj"] * 18, 35, 180),
                        c=scrna["celltype"].astype("category").cat.codes, cmap="tab10",
                        edgecolor="#374151", linewidth=0.35)
        ax.axvline(0, color="#111827", lw=0.65)
        ax.set_yticks(y, [wrap(x, 58) for x in scrna["label"]])
        ax.set_xlabel("scRNA GSEA NES")
        ax.set_title("scRNA pathways in ATAC-linked cell types", fontweight="bold")
        save(fig, "05_scRNA_linked_GSEA_dotplot_polished", pages)

    # 6 Integrated panel
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.2), constrained_layout=True)
    kt2 = kegg_top.sort_values("GeneRatioValue").tail(10)
    y1 = np.arange(len(kt2))
    axes[0].barh(y1, kt2["GeneRatioValue"], color="#B22222", alpha=0.82)
    axes[0].set_yticks(y1, [wrap(x, 34) for x in kt2["Description"]])
    axes[0].set_xlabel("KEGG GeneRatio")
    axes[0].set_title("ATAC-linked genes", fontweight="bold")
    if not scrna.empty:
        sc2 = scrna.sort_values("mlog10padj", ascending=False).head(10).sort_values("NES")
        y2 = np.arange(len(sc2))
        colors = ["#B22222" if v > 0 else "#2B6CB0" for v in sc2["NES"]]
        axes[1].barh(y2, sc2["NES"], color=colors, alpha=0.82)
        axes[1].axvline(0, color="#111827", lw=0.65)
        axes[1].set_yticks(y2, [wrap(c + " | " + d, 34) for c, d in zip(sc2["celltype"], sc2["Description"])])
        axes[1].set_xlabel("scRNA GSEA NES")
        axes[1].set_title("Matched scRNA cell states", fontweight="bold")
    fig.suptitle("Integrated ATAC-scRNA pathway view", fontweight="bold", fontsize=11)
    save(fig, "06_integrated_ATAC_scRNA_pathway_panel", pages)

print(f"Wrote polished figures to {FIG}")
