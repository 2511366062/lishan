#!/usr/bin/env python3
import os
from pathlib import Path

import pandas as pd


ROOT = Path("/mnt/d/lxk/project/lishan-20260613")
DEG_DIR = ROOT / "DEG"
OUT_DIR = ROOT / "analysis/ATAC/annotation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

priority_cells = [
    "Trophoblast",
    "EVT",
    "Endo",
    "Vascular_Endo",
    "Myeloid",
    "HLAII_Mye",
    "Infla_Mono",
    "Decidualized_DSC",
    "Contractile_DSC",
    "Fibroblast",
    "T",
]

rows = []
for sig in sorted(DEG_DIR.glob("*/sig.csv")):
    cell = sig.parent.name
    df = pd.read_csv(sig)
    if df.empty:
        continue
    df = df[df["gene"].astype(str).str.match(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")]
    df = df.dropna(subset=["log2FoldChange", "padj"])
    if df.empty:
        continue
    df["celltype"] = cell
    df["abs_lfc"] = df["log2FoldChange"].abs()
    rows.append(df)

if not rows:
    raise SystemExit("No DEG sig.csv rows found.")

all_sig = pd.concat(rows, ignore_index=True)
all_sig["priority"] = all_sig["celltype"].apply(
    lambda x: priority_cells.index(x) if x in priority_cells else 999
)
all_sig = all_sig.sort_values(["priority", "padj", "abs_lfc"], ascending=[True, True, False])
all_sig.to_csv(OUT_DIR / "scrna_sig_genes_all.tsv", sep="\t", index=False)

top_by_cell = []
for cell, sub in all_sig.groupby("celltype", sort=False):
    top_by_cell.append(
        sub.sort_values(["padj", "abs_lfc"], ascending=[True, False]).head(25)
    )
top = pd.concat(top_by_cell, ignore_index=True)
top.to_csv(OUT_DIR / "scrna_sig_genes_top_by_cell.tsv", sep="\t", index=False)

panel_genes = [
    "FLT1", "ENG", "PAPPA2", "FN1", "CCL21", "HLA-G", "KRT7", "MMP11",
    "IGFBP1", "PRL", "ACTA2", "TAGLN", "VWF", "PECAM1", "CDH5", "ACKR1",
    "LYVE1", "PROX1", "APOE", "FOLR2", "C1QA", "C1QB", "C1QC", "HLA-DRA",
    "HLA-DPA1", "HLA-DQA1", "LYZ", "S100A8", "S100A9", "TREM2", "NKG7",
    "GNLY", "TRAC", "CXCL13", "XCL1",
]
observed = set(all_sig["gene"])
panel = pd.DataFrame({"gene": panel_genes})
panel["in_scrna_sig"] = panel["gene"].isin(observed)
panel.to_csv(OUT_DIR / "manual_marker_gene_panel.tsv", sep="\t", index=False)

with open(OUT_DIR / "scrna_gene_panel.txt", "w") as handle:
    for gene in pd.concat([panel["gene"], top["gene"]]).drop_duplicates():
        handle.write(f"{gene}\n")

print(f"Wrote {OUT_DIR / 'scrna_sig_genes_all.tsv'}")
print(f"Wrote {OUT_DIR / 'scrna_gene_panel.txt'}")
