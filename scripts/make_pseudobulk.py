import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

try:
    import config
except ImportError:
    config = None


def parse_args():
    parser = argparse.ArgumentParser(description="Make pseudo-bulk count matrix for one cell type.")
    parser.add_argument("--base-dir", default=None)
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--celltype", required=True)
    parser.add_argument("--out-prefix", default=None)
    parser.add_argument("--celltype-col", default="celltype")
    parser.add_argument("--sample-col", default="sample")
    parser.add_argument("--counts-layer", default="counts")
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--min-counts", type=int, default=1000)
    return parser.parse_args()


def safe_name(x):
    return re.sub(r"[^A-Za-z0-9_]+", "_", x).strip("_")


def get_base_dir(args):
    if args.base_dir:
        return args.base_dir
    if config is not None and hasattr(config, "PY_BASE_DIR"):
        return str(config.PY_BASE_DIR)
    if config is not None and hasattr(config, "BASE_DIR"):
        return str(config.BASE_DIR)
    return "."


def get_h5ad_path(args):
    if args.h5ad:
        return args.h5ad
    if config is not None and hasattr(config, "H5AD_DIR"):
        return str(config.H5AD_DIR / "rna_final.h5ad")
    return "./h5ad/rna_final.h5ad"


def get_group_colors(adata):
    if "group_colors" not in adata.uns:
        return {}

    if pd.api.types.is_categorical_dtype(adata.obs["group"]):
        groups = list(adata.obs["group"].cat.categories)
    elif config is not None and hasattr(config, "GROUP_ORDER"):
        groups = list(config.GROUP_ORDER)
    else:
        groups = sorted(adata.obs["group"].astype(str).unique())

    colors = list(adata.uns["group_colors"])
    return dict(zip(groups, colors))


def main():
    args = parse_args()
    base_dir = get_base_dir(args)
    h5ad_path = get_h5ad_path(args)
    os.chdir(base_dir)

    celltype = args.celltype
    prefix = safe_name(args.out_prefix if args.out_prefix else celltype)
    out_dir = Path("./DEG") / prefix
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(h5ad_path)
    adata = adata[adata.obs[args.celltype_col].astype(str).isin([celltype])].copy()

    adata.obs[args.sample_col] = adata.obs[args.sample_col].astype(str)
    adata.obs["group"] = adata.obs[args.sample_col].str.replace(r"\d+", "", regex=True)
    group_colors = get_group_colors(adata)

    if args.counts_layer in adata.layers:
        count_mat = adata.layers[args.counts_layer]
    else:
        count_mat = adata.X

    if not sparse.issparse(count_mat):
        count_mat = sparse.csr_matrix(count_mat)

    pb_counts = {}
    metadata = []

    for sample in sorted(adata.obs[args.sample_col].unique()):
        cell_mask = (adata.obs[args.sample_col] == sample).values
        n_cells = int(cell_mask.sum())
        if n_cells == 0:
            continue

        group = re.sub(r"\d+", "", sample)
        pb_sample = f"{sample}_{prefix}"
        summed_counts = np.asarray(count_mat[cell_mask].sum(axis=0)).ravel()

        pb_counts[pb_sample] = summed_counts
        metadata.append(
            {
                "pb_sample": pb_sample,
                "sample": sample,
                "group": group,
                "group_color": group_colors.get(group, ""),
                "celltype": celltype,
                "n_cells": n_cells,
                "total_counts": int(summed_counts.sum()),
            }
        )

    count_df = pd.DataFrame(pb_counts, index=adata.var_names)
    count_df.index.name = "gene"

    meta_df = pd.DataFrame(metadata)
    meta_df["keep"] = (
        (meta_df["n_cells"] >= args.min_cells)
        & (meta_df["total_counts"] >= args.min_counts)
    )

    keep_samples = meta_df.loc[meta_df["keep"], "pb_sample"].tolist()
    count_df_filtered = count_df[keep_samples]
    meta_df_filtered = meta_df[meta_df["keep"]].copy()

    count_df_filtered.to_csv(out_dir / "counts.csv")
    meta_df_filtered.to_csv(out_dir / "metadata.csv", index=False)

    print("celltype:", celltype)
    print("base_dir:", base_dir)
    print("h5ad:", h5ad_path)
    print("all count matrix:", count_df.shape)
    print("filtered count matrix:", count_df_filtered.shape)
    print(meta_df)


if __name__ == "__main__":
    main()
