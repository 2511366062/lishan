#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


ROOT = Path("/mnt/d/lxk/project/lishan-20260613")
ATAC = ROOT / "analysis/ATAC"
OUT = ATAC / "figures" / "figure_set"
OUT.mkdir(parents=True, exist_ok=True)

SAMPLES = pd.read_csv(ATAC / "metadata/samples.tsv", sep="\t")
SAMPLE_COLS = SAMPLES["sample"].tolist()
GROUP = dict(zip(SAMPLES["sample"], SAMPLES["group"]))

COUNTS = pd.read_csv(ATAC / "counts/consensus_peak_counts.tsv", sep="\t")
DESC = pd.read_csv(ATAC / "counts/consensus_peak_descriptive_PE_vs_NP.tsv", sep="\t")
ANNOT = pd.read_csv(ATAC / "annotation/consensus_peaks_nearest_gene.tsv", sep="\t")
LINKED = pd.read_csv(ATAC / "annotation/atac_scrna_linked_genes.tsv", sep="\t")
LINKED_SUM = pd.read_csv(ATAC / "annotation/atac_scrna_linked_gene_summary.tsv", sep="\t")
QC = pd.read_csv(ATAC / "qc/sample_qc_metrics.tsv", sep="\t")
SCRNA = pd.read_csv(ATAC / "annotation/scrna_sig_genes_all.tsv", sep="\t")
MANUAL = pd.read_csv(ATAC / "annotation/manual_marker_gene_panel.tsv", sep="\t")

mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "DejaVu Sans",
    "font.size": 7.5,
    "axes.linewidth": 0.55,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "legend.frameon": False,
})

COL = {
    "NP": "#4C78A8",
    "PE": "#D55E00",
    "PE_open": "#B22222",
    "PE_closed": "#2B6CB0",
    "stable": "#C9CDD3",
    "dark": "#202124",
}
CMAP = LinearSegmentedColormap.from_list("blue_white_red", ["#2166AC", "#F7F7F7", "#B2182B"])


def figsave(fig, name, pages):
    fig.suptitle(name.replace("_", " "), fontsize=10, fontweight="bold", y=0.995)
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    pages.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def add_panel_label(ax, label):
    ax.text(-0.12, 1.08, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")


def cpm_matrix():
    mat = COUNTS.set_index(["chrom", "start", "end"])[SAMPLE_COLS].astype(float)
    return mat.div(mat.sum(axis=0), axis=1) * 1e6


def zscore_rows(df):
    arr = df.values.astype(float)
    return (arr - arr.mean(axis=1, keepdims=True)) / (arr.std(axis=1, keepdims=True) + 1e-9)


def style_ax(ax):
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.45)
    ax.set_axisbelow(True)


with PdfPages(OUT / "ATAC_PE_NP_figure_set_all.pdf") as pages:
    # 1 QC overview
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.45), constrained_layout=True)
    q = QC.copy()
    metrics = [
        ("final_mapped_fragments", "Deduplicated fragments (M)", 1e6),
        ("peaks", "MACS3 peaks (k)", 1e3),
        ("frip", "FRiP", 1),
    ]
    for ax, (col, ylabel, scale) in zip(axs, metrics):
        vals = q[col].astype(float) / scale
        ax.bar(q["sample"], vals, color=[COL[GROUP[s]] for s in q["sample"]], width=0.68)
        ax.set_ylabel(ylabel)
        style_ax(ax)
    add_panel_label(axs[0], "A")
    figsave(fig, "01_QC_library_peak_FRiP", pages)

    # 2 PCA
    cpm = cpm_matrix()
    logcpm = np.log2(cpm + 1)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(StandardScaler().fit_transform(logcpm.T))
    fig, ax = plt.subplots(figsize=(3.7, 3.2), constrained_layout=True)
    for i, s in enumerate(SAMPLE_COLS):
        ax.scatter(coords[i, 0], coords[i, 1], s=72, c=COL[GROUP[s]], edgecolor="black", lw=0.5, zorder=3)
        ax.text(coords[i, 0], coords[i, 1], f" {s}", va="center")
    ax.axhline(0, color="#D1D5DB", lw=0.6)
    ax.axvline(0, color="#D1D5DB", lw=0.6)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    add_panel_label(ax, "A")
    figsave(fig, "02_consensus_peak_PCA", pages)

    # 3 Sample correlation
    corr = logcpm.corr(method="pearson")
    fig, ax = plt.subplots(figsize=(3.35, 3.15), constrained_layout=True)
    im = ax.imshow(corr.values, cmap=CMAP, vmin=0.75, vmax=1.0)
    ax.set_xticks(range(len(SAMPLE_COLS)), SAMPLE_COLS)
    ax.set_yticks(range(len(SAMPLE_COLS)), SAMPLE_COLS)
    for i in range(len(SAMPLE_COLS)):
        for j in range(len(SAMPLE_COLS)):
            ax.text(j, i, f"{corr.iloc[i, j]:.2f}", ha="center", va="center", fontsize=7)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("Pearson r")
    add_panel_label(ax, "A")
    figsave(fig, "03_sample_correlation_heatmap", pages)

    # 4 Accessibility distribution
    fig, ax = plt.subplots(figsize=(4.3, 3.0), constrained_layout=True)
    data = [logcpm[s].sample(min(6000, logcpm.shape[0]), random_state=4).values for s in SAMPLE_COLS]
    parts = ax.violinplot(data, showmeans=False, showextrema=False, widths=0.75)
    for body, s in zip(parts["bodies"], SAMPLE_COLS):
        body.set_facecolor(COL[GROUP[s]])
        body.set_edgecolor("none")
        body.set_alpha(0.75)
    ax.boxplot(data, widths=0.18, showfliers=False, patch_artist=True,
               boxprops={"facecolor": "white", "edgecolor": "black", "linewidth": 0.5},
               medianprops={"color": "black", "linewidth": 0.7},
               whiskerprops={"linewidth": 0.5}, capprops={"linewidth": 0.5})
    ax.set_xticks(range(1, len(SAMPLE_COLS) + 1), SAMPLE_COLS)
    ax.set_ylabel("log2(CPM + 1)")
    style_ax(ax)
    add_panel_label(ax, "A")
    figsave(fig, "04_peak_accessibility_distribution", pages)

    # 5 MA scatter
    d = DESC.copy()
    fig, ax = plt.subplots(figsize=(4.2, 3.2), constrained_layout=True)
    for cat in ["stable", "PE_closed", "PE_open"]:
        sub = d[d["candidate"] == cat]
        ax.scatter(np.log2(sub["mean_cpm"] + 1), sub["atac_log2fc"], s=5 if cat == "stable" else 10,
                   c=COL[cat], alpha=0.25 if cat == "stable" else 0.75, linewidths=0)
    ax.axhline(0, color="black", lw=0.55)
    ax.axhline(1, color="#9CA3AF", lw=0.45, ls="--")
    ax.axhline(-1, color="#9CA3AF", lw=0.45, ls="--")
    ax.set_xlabel("mean accessibility log2(CPM + 1)")
    ax.set_ylabel("ATAC log2FC (PE1 vs NP mean)")
    add_panel_label(ax, "A")
    figsave(fig, "05_PE_vs_NP_MA_candidate_peaks", pages)

    # 6 z-score volcano-like plot
    fig, ax = plt.subplots(figsize=(4.2, 3.2), constrained_layout=True)
    for cat in ["stable", "PE_closed", "PE_open"]:
        sub = d[d["candidate"] == cat]
        ax.scatter(sub["atac_log2fc"], sub["np_z"], s=5 if cat == "stable" else 10,
                   c=COL[cat], alpha=0.22 if cat == "stable" else 0.78, linewidths=0, label=cat.replace("_", " "))
    ax.axvline(0, color="black", lw=0.55)
    ax.axvline(1, color="#9CA3AF", lw=0.45, ls="--")
    ax.axvline(-1, color="#9CA3AF", lw=0.45, ls="--")
    ax.axhline(1.5, color="#9CA3AF", lw=0.45, ls="--")
    ax.axhline(-1.5, color="#9CA3AF", lw=0.45, ls="--")
    ax.set_xlabel("ATAC log2FC")
    ax.set_ylabel("NP-referenced z-score")
    ax.legend(loc="upper left", markerscale=2)
    add_panel_label(ax, "A")
    figsave(fig, "06_candidate_peak_zscore_landscape", pages)

    # 7 Candidate counts and chromosome distribution
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.8), constrained_layout=True)
    cand_counts = d["candidate"].value_counts().reindex(["PE_open", "PE_closed", "stable"]).fillna(0)
    axs[0].bar(cand_counts.index.str.replace("_", " "), cand_counts.values,
               color=[COL[x] for x in cand_counts.index])
    axs[0].set_ylabel("Peak count")
    axs[0].tick_params(axis="x", rotation=25)
    style_ax(axs[0])
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    cd = d[d["candidate"] != "stable"].groupby(["chrom", "candidate"]).size().unstack(fill_value=0).reindex(chrom_order).fillna(0)
    bottom = np.zeros(len(cd))
    for cat in ["PE_open", "PE_closed"]:
        vals = cd[cat].values if cat in cd else np.zeros(len(cd))
        axs[1].bar(range(len(cd)), vals, bottom=bottom, color=COL[cat], width=0.8, label=cat.replace("_", " "))
        bottom += vals
    axs[1].set_xticks(range(len(cd)), [c.replace("chr", "") for c in cd.index], fontsize=6)
    axs[1].set_ylabel("Candidate peaks")
    axs[1].legend(loc="upper right")
    style_ax(axs[1])
    add_panel_label(axs[0], "A")
    add_panel_label(axs[1], "B")
    figsave(fig, "07_candidate_peak_counts_by_chromosome", pages)

    # 8 Distance to nearest gene
    nearest = ANNOT.sort_values("distance", key=lambda s: s.abs()).drop_duplicates("peak_id")
    nearest["abs_distance_kb"] = nearest["distance"].abs() / 1000
    fig, ax = plt.subplots(figsize=(4.25, 3.1), constrained_layout=True)
    bins = np.linspace(0, min(250, nearest["abs_distance_kb"].quantile(0.99)), 45)
    for cat in ["stable", "PE_open", "PE_closed"]:
        sub = nearest[nearest["candidate"] == cat]
        ax.hist(sub["abs_distance_kb"], bins=bins, density=True, histtype="stepfilled",
                alpha=0.22 if cat == "stable" else 0.38, color=COL[cat], label=cat.replace("_", " "))
    ax.set_xlabel("Distance to nearest gene (kb)")
    ax.set_ylabel("Density")
    ax.legend(loc="upper right")
    add_panel_label(ax, "A")
    figsave(fig, "08_peak_to_gene_distance_distribution", pages)

    # 9 scRNA-linked gene concordance bar
    lsum = LINKED_SUM.sort_values("best_score", ascending=True)
    fig, ax = plt.subplots(figsize=(4.65, 3.25), constrained_layout=True)
    y = np.arange(len(lsum))
    colors = [COL["PE_open"] if x > 0 else COL["PE_closed"] for x in lsum["atac_log2fc"]]
    labels = lsum["gene"] + " (" + lsum["celltype"] + ")"
    ax.barh(y, lsum["atac_log2fc"], color=colors, height=0.68)
    ax.set_yticks(y, labels)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("ATAC log2FC (PE1 vs NP mean)")
    add_panel_label(ax, "A")
    figsave(fig, "09_scRNA_linked_gene_ATAC_effects", pages)

    # 10 linked gene heatmap
    gene_map = nearest.set_index(["chrom", "start", "end"])["gene"].astype(str)
    gene_rows = []
    names = []
    for gene in LINKED_SUM.sort_values("best_score", ascending=False)["gene"].drop_duplicates():
        idx = [x for x in gene_map[gene_map == gene].index if x in cpm.index]
        if not idx:
            continue
        gene_rows.append(np.log2(cpm.loc[idx].sum(axis=0) + 1).values)
        names.append(gene)
    fig, ax = plt.subplots(figsize=(3.7, 3.2), constrained_layout=True)
    arr = zscore_rows(pd.DataFrame(gene_rows, columns=SAMPLE_COLS))
    im = ax.imshow(arr, aspect="auto", cmap=CMAP, vmin=-2, vmax=2)
    ax.set_xticks(range(len(SAMPLE_COLS)), SAMPLE_COLS)
    ax.set_yticks(range(len(names)), names)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("row z-score")
    add_panel_label(ax, "A")
    figsave(fig, "10_scRNA_linked_gene_accessibility_heatmap", pages)

    # 11 ATAC vs scRNA effect concordance
    fig, ax = plt.subplots(figsize=(3.65, 3.25), constrained_layout=True)
    ax.scatter(LINKED["scrna_log2fc"], LINKED["atac_log2fc"],
               c=[COL["PE_open"] if x > 0 else COL["PE_closed"] for x in LINKED["atac_log2fc"]],
               s=np.clip(LINKED["score"] * 18, 25, 160), edgecolor="black", lw=0.4, alpha=0.88)
    for _, r in LINKED_SUM.iterrows():
        ax.text(r["scrna_log2fc"], r["atac_log2fc"], " " + r["gene"], fontsize=7, va="center")
    ax.axhline(0, color="#9CA3AF", lw=0.55)
    ax.axvline(0, color="#9CA3AF", lw=0.55)
    ax.set_xlabel("scRNA log2FC")
    ax.set_ylabel("ATAC log2FC")
    add_panel_label(ax, "A")
    figsave(fig, "11_ATAC_scRNA_concordance_scatter", pages)

    # 12 marker panel accessibility
    marker_genes = [g for g in MANUAL["gene"].tolist() if g in set(nearest["gene"])]
    marker_rows, marker_names = [], []
    for gene in marker_genes:
        idx = [x for x in gene_map[gene_map == gene].index if x in cpm.index]
        if not idx:
            continue
        marker_rows.append(np.log2(cpm.loc[idx].sum(axis=0) + 1).values)
        marker_names.append(gene)
    if marker_rows:
        marker_df = pd.DataFrame(marker_rows, index=marker_names, columns=SAMPLE_COLS)
        spread = marker_df.max(axis=1) - marker_df.min(axis=1)
        marker_df = marker_df.loc[spread.sort_values(ascending=False).head(24).index]
        fig, ax = plt.subplots(figsize=(3.8, 5.2), constrained_layout=True)
        arr = zscore_rows(marker_df)
        im = ax.imshow(arr, aspect="auto", cmap=CMAP, vmin=-2, vmax=2)
        ax.set_xticks(range(len(SAMPLE_COLS)), SAMPLE_COLS)
        ax.set_yticks(range(marker_df.shape[0]), marker_df.index)
        cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
        cb.set_label("row z-score")
        add_panel_label(ax, "A")
        figsave(fig, "12_marker_gene_accessibility_heatmap", pages)

    # 13 overlap-like binary peak presence summary
    present = (COUNTS[SAMPLE_COLS] > 0).astype(int)
    patterns = present.astype(str).agg("".join, axis=1).value_counts().head(12)
    fig, axs = plt.subplots(2, 1, figsize=(5.4, 3.7), gridspec_kw={"height_ratios": [2.2, 1]}, constrained_layout=True)
    axs[0].bar(range(len(patterns)), patterns.values, color="#444444", width=0.72)
    axs[0].set_ylabel("Consensus peaks")
    axs[0].set_xticks([])
    style_ax(axs[0])
    mat = np.array([[int(c) for c in p] for p in patterns.index]).T
    axs[1].imshow(mat, aspect="auto", cmap=LinearSegmentedColormap.from_list("presence", ["#E5E7EB", "#111827"]), vmin=0, vmax=1)
    axs[1].set_yticks(range(len(SAMPLE_COLS)), SAMPLE_COLS)
    axs[1].set_xticks(range(len(patterns)), range(1, len(patterns) + 1))
    axs[1].set_xlabel("Peak presence pattern")
    add_panel_label(axs[0], "A")
    figsave(fig, "13_peak_presence_overlap_patterns", pages)

print(f"Wrote figures to {OUT}")
