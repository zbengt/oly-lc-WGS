#!/usr/bin/env python3
"""
Plot PCA (PC1 vs PC2) and IBS heatmap from bam summary tables.

Usage:
    python plot_bam_relatedness.py \
        --pca output/02_bam_summary/tables/bam_pca_components.tsv \
        --ibs output/02_bam_summary/tables/ibs_matrix.tsv \
        --out output/02_bam_summary/figures/bam_relatedness_demo.png

Saves a two-panel PNG: PCA scatter (left) colored by location and IBS heatmap (right),
with hierarchical clustering ordering on the heatmap to highlight related sample groups.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy as sch
from scipy.spatial.distance import squareform


def read_pca(path):
    df = pd.read_csv(path, sep="\t", dtype={"sample_id": str})
    # remove blanks or placeholder rows if present
    df = df[~df["sample_id"].str.startswith("Blank")]
    return df.set_index("sample_id")


def read_ibs(path):
    # IBS matrix: first column is index, first row is header (empty first cell)
    df = pd.read_csv(path, sep="\t", index_col=0)
    # ensure order and numeric
    df = df.astype(float)
    return df


def plot(pca_df, ibs_df, out_path, annotate_thresh=1.0):
    # align samples present in both
    samples = [s for s in pca_df.index if s in ibs_df.index]
    pca_df = pca_df.loc[samples]
    ibs_df = ibs_df.loc[samples, samples]

    # prepare colors by location
    locations = pca_df["location"].astype(str)
    unique_locs = list(locations.unique())
    palette = dict(zip(unique_locs, sns.color_palette("tab20", n_colors=len(unique_locs))))
    colors = locations.map(palette)

    # clustering for heatmap ordering (use 1 - IBS as dist)
    # enforce numeric and symmetry
    mat = ibs_df.values
    # convert to distance matrix
    dist = 1.0 - mat
    # ensure diagonal zeros
    np.fill_diagonal(dist, 0.0)
    # condensed form
    try:
        condensed = squareform(dist)
        linkage = sch.linkage(condensed, method="average")
        dendro = sch.dendrogram(linkage, no_plot=True)
        order = dendro["leaves"]
    except Exception:
        # fallback: keep original order
        order = list(range(len(samples)))

    ordered_samples = [samples[i] for i in order]
    ordered_mat = ibs_df.loc[ordered_samples, ordered_samples]

    # make figure with two panels
    fig = plt.figure(figsize=(14, 6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.0], wspace=0.3)

    # PCA scatter
    ax0 = fig.add_subplot(gs[0, 0])
    sc = ax0.scatter(pca_df["PC1"], pca_df["PC2"], c=[palette[l] for l in pca_df["location"]], s=60, edgecolor='k', linewidth=0.3)
    ax0.axhline(0, color="0.8", linestyle="--", linewidth=0.6)
    ax0.axvline(0, color="0.8", linestyle="--", linewidth=0.6)
    ax0.set_xlabel("PC1")
    ax0.set_ylabel("PC2")
    ax0.set_title("PCA (PC1 vs PC2)")

    # add legend for locations
    handles = []
    for loc in unique_locs:
        handles.append(plt.Line2D([0], [0], marker='o', color='w', label=loc,
                                  markerfacecolor=palette[loc], markersize=8, markeredgecolor='k', markeredgewidth=0.3))
    ax0.legend(handles=handles, bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0., title='location')

    # annotate outliers (large absolute PC values)
    for sid, row in pca_df.iterrows():
        if abs(row["PC1"]) >= annotate_thresh or abs(row["PC2"]) >= annotate_thresh:
            ax0.annotate(sid, (row["PC1"], row["PC2"]), fontsize=8, alpha=0.9)

    # IBS heatmap
    ax1 = fig.add_subplot(gs[0, 1])
    sns.heatmap(ordered_mat, ax=ax1, cmap="viridis", vmin=ordered_mat.values.min(), vmax=ordered_mat.values.max(), cbar_kws={"label": "IBS"})
    ax1.set_xticks(np.arange(len(ordered_samples)) + 0.5)
    ax1.set_yticks(np.arange(len(ordered_samples)) + 0.5)
    # label sparsely to avoid crowding
    if len(ordered_samples) <= 40:
        ax1.set_xticklabels(ordered_samples, rotation=90, fontsize=6)
        ax1.set_yticklabels(ordered_samples, rotation=0, fontsize=6)
    else:
        # show every nth label
        step = max(1, len(ordered_samples)//40)
        ax1.set_xticklabels([s if (i%step==0) else '' for i,s in enumerate(ordered_samples)], rotation=90, fontsize=6)
        ax1.set_yticklabels([s if (i%step==0) else '' for i,s in enumerate(ordered_samples)], rotation=0, fontsize=6)
    ax1.set_title("IBS heatmap (clustered)")

    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Saved figure to {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pca", required=True, help="Path to bam_pca_components.tsv")
    parser.add_argument("--ibs", required=True, help="Path to ibs_matrix.tsv")
    parser.add_argument("--out", required=True, help="Output PNG path")
    args = parser.parse_args()

    pca_df = read_pca(args.pca)
    ibs_df = read_ibs(args.ibs)
    plot(pca_df, ibs_df, args.out)


if __name__ == '__main__':
    main()
