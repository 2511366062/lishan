#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


ROOT = Path("/mnt/d/lxk/project/lishan-20260613")
ATAC = ROOT / "analysis/ATAC"
FIG = ATAC / "figures"
FIG.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
})

COLORS = {"NP": "#4C78A8", "PE": "#D55E00"}
CMAP = LinearSegmentedColormap.from_list("nature_atac", ["#2B6CB0", "#F7F7F7", "#B22222"])


def read_table(path, **kwargs):
    if not Path(path).exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", **kwargs)


samples = pd.read_csv(ATAC / "metadata/samples.tsv", sep="\t")
qc = read_table(ATAC / "qc/sample_qc_metrics.tsv")
counts = read_table(ATAC / "counts/consensus_peak_counts.tsv")
annot = read_table(ATAC / "annotation/consensus_peaks_nearest_gene.tsv")
linked = read_table(ATAC / "annotation/atac_scrna_linked_genes.tsv")

pdf_path = FIG / "ATAC_PE_NP_scRNA_linked_summary.pdf"
with PdfPages(pdf_path) as pdf:
    fig = plt.figure(figsize=(7.2, 6.8), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig)
    fig.suptitle("Bulk ATAC-seq overview: PE-linked chromatin accessibility", fontsize=11, fontweight="bold")

    ax = fig.add_subplot(gs[0, 0])
    if not qc.empty:
        q = qc.copy()
        if "group" not in q.columns:
            q = q.merge(samples[["sample", "group"]], on="sample", how="left")
        x = np.arange(len(q))
        vals = q["final_mapped_fragments"].astype(float) / 1e6
        ax.bar(x, vals, color=[COLORS.get(g, "0.5") for g in q["group"]], width=0.72)
        ax.set_xticks(x, q["sample"], rotation=0)
        ax.set_ylabel("Deduplicated fragments (M)")
        ax.set_title("Library depth")
    else:
        ax.text(0.5, 0.5, "QC metrics not available yet", ha="center", va="center")
        ax.set_axis_off()

    ax = fig.add_subplot(gs[0, 1])
    if not counts.empty:
        mat = counts.set_index(["chrom", "start", "end"])[samples["sample"].tolist()].astype(float)
        cpm = mat.div(mat.sum(axis=0), axis=1) * 1e6
        x = StandardScaler().fit_transform(np.log2(cpm + 1).T)
        pca = PCA(n_components=2).fit(x)
        coords = pca.transform(x)
        for i, row in samples.iterrows():
            ax.scatter(coords[i, 0], coords[i, 1], s=55, color=COLORS[row.group], edgecolor="black", linewidth=0.5)
            ax.text(coords[i, 0], coords[i, 1], " " + row["sample"], va="center", fontsize=8)
        ax.axhline(0, color="0.85", lw=0.5)
        ax.axvline(0, color="0.85", lw=0.5)
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        ax.set_title("Consensus peak accessibility")
    else:
        ax.text(0.5, 0.5, "Peak count matrix not available yet", ha="center", va="center")
        ax.set_axis_off()

    ax = fig.add_subplot(gs[1, 0])
    if not linked.empty:
        sub = linked.dropna(subset=["atac_log2fc"]).copy()
        sub = sub.sort_values("score", ascending=False).head(16)
        y = np.arange(len(sub))[::-1]
        ax.barh(y, sub["atac_log2fc"], color=[("#B22222" if v > 0 else "#2B6CB0") for v in sub["atac_log2fc"]], height=0.65)
        labels = sub["gene"].astype(str) + " (" + sub["celltype"].astype(str) + ")"
        ax.set_yticks(y, labels)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_xlabel("ATAC log2FC (PE1 vs NP mean)")
        ax.set_title("scRNA DEG genes with nearby ATAC change")
    else:
        ax.text(0.5, 0.5, "scRNA-linked peak annotation not available yet", ha="center", va="center")
        ax.set_axis_off()

    ax = fig.add_subplot(gs[1, 1])
    if not counts.empty and not annot.empty:
        mat = counts.set_index(["chrom", "start", "end"])[samples["sample"].tolist()].astype(float)
        cpm = mat.div(mat.sum(axis=0), axis=1) * 1e6
        gene_map = annot.set_index(["chrom", "start", "end"])["gene"].astype(str)
        top_genes = linked.sort_values("score", ascending=False)["gene"].drop_duplicates().head(20).tolist() if not linked.empty else []
        rows = []
        names = []
        for gene in top_genes:
            idx = gene_map[gene_map == gene].index
            idx = [i for i in idx if i in cpm.index]
            if not idx:
                continue
            v = np.log2(cpm.loc[idx].sum(axis=0) + 1)
            rows.append(v.values)
            names.append(gene)
        if rows:
            arr = np.vstack(rows)
            arr = (arr - arr.mean(axis=1, keepdims=True)) / (arr.std(axis=1, keepdims=True) + 1e-9)
            im = ax.imshow(arr, aspect="auto", cmap=CMAP, vmin=-2, vmax=2)
            ax.set_xticks(np.arange(len(samples)), samples["sample"])
            ax.set_yticks(np.arange(len(names)), names)
            ax.set_title("Accessibility near linked genes")
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
            cb.set_label("row z-score")
        else:
            ax.text(0.5, 0.5, "No linked genes in count matrix", ha="center", va="center")
            ax.set_axis_off()
    else:
        ax.text(0.5, 0.5, "Gene-linked heatmap not available yet", ha="center", va="center")
        ax.set_axis_off()

    pdf.savefig(fig)
    plt.close(fig)

print(pdf_path)
