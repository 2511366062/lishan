#!/usr/bin/env python3
from pathlib import Path
import gzip
import re
from collections import defaultdict

import numpy as np
import pandas as pd
import pyBigWig
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


ROOT = Path("/mnt/d/lxk/project/lishan-20260613")
ATAC = ROOT / "analysis/ATAC"
STRICT_OUT = ROOT / "fig/ATAC_linked_gene_sets/strict/browser_tracks"
MEDIUM_OUT = ROOT / "fig/ATAC_linked_gene_sets/medium_absATAC0.5/browser_tracks"
STRICT_OUT.mkdir(parents=True, exist_ok=True)
MEDIUM_OUT.mkdir(parents=True, exist_ok=True)

SAMPLES = pd.read_csv(ATAC / "metadata/samples.tsv", sep="\t")
SAMPLE_COLS = SAMPLES["sample"].tolist()
GROUP = dict(zip(SAMPLES["sample"], SAMPLES["group"]))
BW = {s: str(ATAC / "bigwig" / f"{s}.CPM.bw") for s in SAMPLE_COLS}
DESC = pd.read_csv(ATAC / "counts/consensus_peak_descriptive_PE_vs_NP.tsv", sep="\t")
STRICT_GENES = ["PLA2G2A", "CXCL13", "DAPP1", "LYVE1", "TNFSF10"]
relaxed = pd.read_csv(ATAC / "annotation/atac_scrna_linked_genes_relaxed_all_peaks.tsv", sep="\t")
MEDIUM_GENES = (
    relaxed[(relaxed["abs_atac_log2fc"] >= 0.5) & (relaxed["concordant"])]
    .loc[lambda d: ~d["gene"].astype(str).str.startswith("ENSG")]
    .sort_values("score_relaxed", ascending=False)["gene"]
    .drop_duplicates()
    .head(20)
    .tolist()
)
GENES = list(dict.fromkeys(STRICT_GENES + MEDIUM_GENES))
GTF = "/home/vapor/database/source/human/genes.filtered.gtf"
COL = {"NP": "#4C78A8", "PE": "#D55E00", "PE_open": "#B22222", "PE_closed": "#2B6CB0", "stable": "#BDBDBD"}

mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "DejaVu Sans", "font.size": 7.5})


def parse_attrs(text):
    return dict(re.findall(r'([A-Za-z0-9_]+) "([^"]+)"', text))


def load_gene_models(genes):
    wanted = set(genes)
    models = defaultdict(lambda: {"exon": []})
    opener = gzip.open if GTF.endswith(".gz") else open
    with opener(GTF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in {"gene", "exon"}:
                continue
            attrs = parse_attrs(f[8])
            gene = attrs.get("gene_name", attrs.get("gene_id", ""))
            if gene not in wanted:
                continue
            chrom, feature, start, end, strand = f[0], f[2], int(f[3]) - 1, int(f[4]), f[6]
            models[gene]["chrom"] = chrom
            models[gene]["strand"] = strand
            if feature == "gene":
                models[gene]["gene"] = (start, end)
            else:
                models[gene]["exon"].append((start, end))
    return models


def fetch_bw(path, chrom, start, end, bins):
    with pyBigWig.open(path) as bw:
        vals = bw.stats(chrom, max(0, start), end, nBins=bins, type="mean")
    return np.array([0 if v is None or np.isnan(v) else v for v in vals], dtype=float)


def plot_gene(gene, out_dir, title_prefix):
    m = models.get(gene, {})
    if "gene" not in m:
        return
    chrom = m["chrom"]
    gstart, gend = m["gene"]
    width = max(8000, gend - gstart)
    start = max(0, gstart - width // 2)
    end = gend + width // 2
    bins = 430
    x = np.linspace(start, end, bins) / 1e6
    fig = plt.figure(figsize=(7.2, 3.6), constrained_layout=True)
    gs = fig.add_gridspec(6, 1, height_ratios=[1, 1, 1, 1, 0.45, 0.55])
    tracks = {s: fetch_bw(BW[s], chrom, start, end, bins) for s in SAMPLE_COLS}
    ymax = max(float(v.max()) for v in tracks.values()) * 1.08 + 1e-6
    for i, s in enumerate(SAMPLE_COLS):
        ax = fig.add_subplot(gs[i, 0])
        ax.fill_between(x, tracks[s], step="mid", color=COL[GROUP[s]], alpha=0.9)
        ax.set_ylim(0, ymax)
        ax.set_ylabel(s, rotation=0, ha="right", va="center")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axp = fig.add_subplot(gs[4, 0])
    local = DESC[(DESC["chrom"] == chrom) & (DESC["start"] < end) & (DESC["end"] > start)]
    for _, r in local.iterrows():
        axp.add_patch(Rectangle((r["start"] / 1e6, 0.18), max((r["end"] - r["start"]) / 1e6, 0.00004), 0.64,
                                color=COL.get(r["candidate"], "#BDBDBD"), alpha=0.9))
    axp.set_xlim(start / 1e6, end / 1e6)
    axp.set_yticks([])
    axp.set_ylabel("peaks", rotation=0, ha="right", va="center")
    axp.set_xticks([])
    axg = fig.add_subplot(gs[5, 0])
    axg.hlines(0.5, gstart / 1e6, gend / 1e6, color="#111827", lw=0.8)
    for es, ee in m["exon"]:
        if ee >= start and es <= end:
            axg.add_patch(Rectangle((max(es, start) / 1e6, 0.32), (min(ee, end) - max(es, start)) / 1e6, 0.36, color="#111827"))
    axg.text((gstart + gend) / 2 / 1e6, 0.88, gene, ha="center", va="bottom", fontweight="bold")
    axg.set_xlim(start / 1e6, end / 1e6)
    axg.set_yticks([])
    axg.set_xlabel(f"{chrom} position (Mb)")
    fig.suptitle(f"{title_prefix}: {gene}", fontweight="bold")
    fig.savefig(out_dir / f"browser_{gene}.pdf", bbox_inches="tight")
    plt.close(fig)


models = load_gene_models(GENES)
for gene in STRICT_GENES:
    plot_gene(gene, STRICT_OUT, "Strict ATAC-scRNA linked gene")
for gene in MEDIUM_GENES:
    plot_gene(gene, MEDIUM_OUT, "Medium ATAC-scRNA linked gene")

print(f"Wrote strict browser tracks to {STRICT_OUT}")
print(f"Wrote medium browser tracks to {MEDIUM_OUT}")
