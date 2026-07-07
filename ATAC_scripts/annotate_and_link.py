#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/mnt/d/lxk/project/lishan-20260613")
ATAC = ROOT / "analysis/ATAC"
samples = pd.read_csv(ATAC / "metadata/samples.tsv", sep="\t")
count_path = ATAC / "counts/consensus_peak_counts.tsv"
nearest_path = ATAC / "annotation/consensus_peaks_nearest_gene.raw.tsv"
scrna_path = ATAC / "annotation/scrna_sig_genes_all.tsv"

counts = pd.read_csv(count_path, sep="\t")
peak_id = counts[["chrom", "start", "end"]].astype(str).agg(":".join, axis=1)
sample_cols = samples["sample"].tolist()
mat = counts[sample_cols].astype(float)
cpm = mat.div(mat.sum(axis=0), axis=1) * 1e6
log_cpm = np.log2(cpm + 1)
np_cols = samples.loc[samples["group"] == "NP", "sample"].tolist()
pe_cols = samples.loc[samples["group"] == "PE", "sample"].tolist()
if len(pe_cols) != 1:
    raise SystemExit("This descriptive linker expects one PE sample.")
pe = pe_cols[0]

stats = counts[["chrom", "start", "end"]].copy()
stats["peak_id"] = peak_id
stats["np_mean_log2cpm"] = log_cpm[np_cols].mean(axis=1)
stats["np_sd_log2cpm"] = log_cpm[np_cols].std(axis=1)
stats["pe_log2cpm"] = log_cpm[pe]
stats["atac_log2fc"] = stats["pe_log2cpm"] - stats["np_mean_log2cpm"]
stats["np_z"] = stats["atac_log2fc"] / (stats["np_sd_log2cpm"] + 0.25)
stats["mean_cpm"] = cpm.mean(axis=1)
stats["candidate"] = np.select(
    [
        (stats["mean_cpm"] >= 0.5) & (stats["atac_log2fc"] >= 1.0) & (stats["np_z"] >= 1.5),
        (stats["mean_cpm"] >= 0.5) & (stats["atac_log2fc"] <= -1.0) & (stats["np_z"] <= -1.5),
    ],
    ["PE_open", "PE_closed"],
    default="stable",
)
stats.to_csv(ATAC / "counts/consensus_peak_descriptive_PE_vs_NP.tsv", sep="\t", index=False)

nearest_cols = [
    "chrom", "start", "end", "peak_name",
    "gene_chrom", "gene_start", "gene_end", "gene", "gene_id", "strand", "distance",
]
nearest = pd.read_csv(nearest_path, sep="\t", header=None, names=nearest_cols)
nearest = nearest[["chrom", "start", "end", "gene", "gene_id", "strand", "distance"]]
annot = nearest.merge(stats, on=["chrom", "start", "end"], how="left")
annot.to_csv(ATAC / "annotation/consensus_peaks_nearest_gene.tsv", sep="\t", index=False)

scrna = pd.read_csv(scrna_path, sep="\t")
scrna = scrna.rename(columns={"log2FoldChange": "scrna_log2fc", "padj": "scrna_padj"})
linked = annot.merge(scrna, on="gene", how="inner")
linked = linked[linked["candidate"] != "stable"].copy()
linked["direction_concordant"] = np.sign(linked["atac_log2fc"]) == np.sign(linked["scrna_log2fc"])
linked["score"] = (
    linked["atac_log2fc"].abs()
    * np.minimum(-np.log10(linked["scrna_padj"].clip(lower=1e-300)), 50)
    / (1 + linked["distance"].abs() / 100000)
)
linked = linked.sort_values(["direction_concordant", "score"], ascending=[False, False])
linked.to_csv(ATAC / "annotation/atac_scrna_linked_genes.tsv", sep="\t", index=False)

gene_summary = (
    linked.groupby(["gene", "celltype"], as_index=False)
    .agg(
        best_score=("score", "max"),
        atac_log2fc=("atac_log2fc", lambda x: x.iloc[x.abs().argmax()]),
        scrna_log2fc=("scrna_log2fc", "first"),
        scrna_padj=("scrna_padj", "first"),
        nearest_peak_distance=("distance", lambda x: x.abs().min()),
        n_candidate_peaks=("peak_id", "nunique"),
        concordant=("direction_concordant", "max"),
    )
    .sort_values(["concordant", "best_score"], ascending=[False, False])
)
gene_summary.to_csv(ATAC / "annotation/atac_scrna_linked_gene_summary.tsv", sep="\t", index=False)

print(f"candidate peaks: {(stats['candidate'] != 'stable').sum()}")
print(f"linked rows: {len(linked)}")
