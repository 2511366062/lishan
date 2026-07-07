#!/usr/bin/env python3
from pathlib import Path
import textwrap

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap


ROOT = Path("/mnt/d/lxk/project/lishan-20260613")
ATAC = ROOT / "analysis/ATAC"
FIG_ROOT = ROOT / "fig" / "ATAC_linked_gene_sets"
FIG_ROOT.mkdir(parents=True, exist_ok=True)

SAMPLES = pd.read_csv(ATAC / "metadata/samples.tsv", sep="\t")
SAMPLE_COLS = SAMPLES["sample"].tolist()
COUNTS = pd.read_csv(ATAC / "counts/consensus_peak_counts.tsv", sep="\t")
ANNOT = pd.read_csv(ATAC / "annotation/consensus_peaks_nearest_gene.tsv", sep="\t")
RELAXED = pd.read_csv(ATAC / "annotation/atac_scrna_linked_genes_relaxed_all_peaks.tsv", sep="\t")
STRICT = pd.read_csv(ATAC / "annotation/atac_scrna_linked_gene_summary.tsv", sep="\t")

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

CMAP = LinearSegmentedColormap.from_list("blue_white_red", ["#2166AC", "#F7F7F7", "#B2182B"])
COL_OPEN = "#B22222"
COL_CLOSED = "#2B6CB0"


def wrap(x, width=42):
    return "\n".join(textwrap.wrap(str(x), width=width))


def cpm_matrix():
    mat = COUNTS.set_index(["chrom", "start", "end"])[SAMPLE_COLS].astype(float)
    return mat.div(mat.sum(axis=0), axis=1) * 1e6


CPM = cpm_matrix()
GENE_MAP = (
    ANNOT.sort_values("distance", key=lambda s: s.abs())
    .drop_duplicates(["chrom", "start", "end", "gene"])
    .set_index(["chrom", "start", "end"])["gene"].astype(str)
)


def gene_accessibility_matrix(genes):
    rows, names = [], []
    for gene in genes:
        idx = [x for x in GENE_MAP[GENE_MAP == gene].index if x in CPM.index]
        if not idx:
            continue
        rows.append(np.log2(CPM.loc[idx].sum(axis=0) + 1).values)
        names.append(gene)
    if not rows:
        return pd.DataFrame(columns=SAMPLE_COLS)
    return pd.DataFrame(rows, index=names, columns=SAMPLE_COLS)


def zscore_rows(df):
    arr = df.values.astype(float)
    arr = (arr - arr.mean(axis=1, keepdims=True)) / (arr.std(axis=1, keepdims=True) + 1e-9)
    return arr


def save(fig, path, pages=None):
    fig.savefig(path, bbox_inches="tight")
    if pages is not None:
        pages.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def make_set(name, table, gene_col="gene"):
    out = FIG_ROOT / name
    out.mkdir(parents=True, exist_ok=True)
    table = table.copy()
    genes = table[gene_col].drop_duplicates().tolist()
    table.to_csv(out / f"{name}_linked_gene_table.tsv", sep="\t", index=False)

    with PdfPages(out / f"{name}_linked_gene_figures_all.pdf") as pages:
        # Effect bar
        bar = table.sort_values("score_relaxed" if "score_relaxed" in table.columns else "best_score", ascending=True)
        if "score_relaxed" not in bar.columns:
            bar["score_relaxed"] = bar["best_score"]
        fig_h = max(3.2, 0.28 * len(bar) + 1.2)
        fig, ax = plt.subplots(figsize=(7.2, fig_h), constrained_layout=True)
        y = np.arange(len(bar))
        colors = [COL_OPEN if v > 0 else COL_CLOSED for v in bar["atac_log2fc"]]
        labels = bar["gene"] + " (" + bar["celltype"] + ")"
        ax.barh(y, bar["atac_log2fc"], color=colors, height=0.66)
        ax.axvline(0, color="#111827", lw=0.65)
        ax.set_yticks(y, [wrap(x, 38) for x in labels])
        ax.set_xlabel("ATAC log2FC (PE1 vs NP mean)")
        ax.set_title(f"{name}: ATAC-linked scRNA genes", fontweight="bold")
        save(fig, out / f"{name}_01_ATAC_effect_bar.pdf", pages)

        # ATAC-scRNA scatter
        fig, ax = plt.subplots(figsize=(4.4, 3.8), constrained_layout=True)
        size_col = "score_relaxed" if "score_relaxed" in table.columns else "best_score"
        ax.scatter(table["scrna_log2fc"], table["atac_log2fc"],
                   s=np.clip(table[size_col] * 18, 28, 180),
                   c=[COL_OPEN if v > 0 else COL_CLOSED for v in table["atac_log2fc"]],
                   edgecolor="#111827", linewidth=0.35, alpha=0.88)
        for _, r in table.sort_values(size_col, ascending=False).head(18).iterrows():
            ax.text(r["scrna_log2fc"], r["atac_log2fc"], " " + r["gene"], fontsize=7, va="center")
        ax.axhline(0, color="#9CA3AF", lw=0.6)
        ax.axvline(0, color="#9CA3AF", lw=0.6)
        ax.set_xlabel("scRNA log2FC")
        ax.set_ylabel("ATAC log2FC")
        ax.set_title(f"{name}: ATAC-scRNA concordance", fontweight="bold")
        save(fig, out / f"{name}_02_ATAC_scRNA_scatter.pdf", pages)

        # Heatmap
        mat = gene_accessibility_matrix(genes)
        if not mat.empty:
            spread = mat.max(axis=1) - mat.min(axis=1)
            mat = mat.loc[spread.sort_values(ascending=False).index]
            fig_h = max(3.2, 0.22 * mat.shape[0] + 1.2)
            fig, ax = plt.subplots(figsize=(3.8, fig_h), constrained_layout=True)
            im = ax.imshow(zscore_rows(mat), aspect="auto", cmap=CMAP, vmin=-2, vmax=2)
            ax.set_xticks(np.arange(len(SAMPLE_COLS)), SAMPLE_COLS)
            ax.set_yticks(np.arange(mat.shape[0]), mat.index)
            cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
            cb.set_label("row z-score")
            ax.set_title(f"{name}: accessibility near linked genes", fontweight="bold")
            save(fig, out / f"{name}_03_gene_accessibility_heatmap.pdf", pages)

        # Candidate composition
        if "candidate" in table.columns:
            comp = table["candidate"].value_counts()
            fig, ax = plt.subplots(figsize=(3.8, 2.9), constrained_layout=True)
            ax.bar(comp.index, comp.values, color=[COL_OPEN if x == "PE_open" else COL_CLOSED if x == "PE_closed" else "#A0A0A0" for x in comp.index])
            ax.set_ylabel("Gene-celltype links")
            ax.set_title(f"{name}: link classes", fontweight="bold")
            ax.tick_params(axis="x", rotation=25)
            save(fig, out / f"{name}_04_link_class_counts.pdf", pages)

    return out


strict = RELAXED[(RELAXED["candidate"] != "stable") & (RELAXED["concordant"])].copy()
strict = strict.sort_values("score_relaxed", ascending=False)
medium = RELAXED[(RELAXED["abs_atac_log2fc"] >= 0.5) & (RELAXED["concordant"])].copy()
medium = medium.sort_values("score_relaxed", ascending=False)

strict_out = make_set("strict", strict)
medium_out = make_set("medium_absATAC0.5", medium)

print(f"Wrote strict figures to {strict_out}")
print(f"Wrote medium figures to {medium_out}")
