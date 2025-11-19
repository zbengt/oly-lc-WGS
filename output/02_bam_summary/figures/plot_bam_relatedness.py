#!/usr/bin/env python3
"""
Plot BAM-based genetic relatedness summary figure.

This script produces a two-panel figure summarizing genetic relatedness from
BAM-derived PCA and IBS (identity-by-state) data. The left panel shows a PCA
scatter plot (PC1 vs PC2) colored by sampling location, with extreme samples
annotated. The right panel displays a hierarchically-clustered heatmap of the
IBS similarity matrix.

Usage:
    python3 output/02_bam_summary/figures/plot_bam_relatedness.py \
        --pca output/02_bam_summary/tables/bam_pca_components.tsv \
        --ibs output/02_bam_summary/tables/ibs_matrix.tsv \
        --out output/02_bam_summary/figures/bam_relatedness_demo.png

Requirements:
    - pandas
    - numpy
    - matplotlib
    - seaborn
    - scipy

Behavior:
    - Filters out placeholder "Blank" samples from PCA input
    - Robustly parses IBS matrix, dropping unnamed trailing columns if present
    - Performs hierarchical clustering using distance = 1 - IBS with average linkage
    - Falls back to original ordering if clustering fails
    - Annotates samples with abs(PC1) >= 1.0 or abs(PC2) >= 1.0 by default
    - Saves figure at user-specified path and prints confirmation message
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot BAM-based genetic relatedness summary (PCA + IBS heatmap)."
    )
    parser.add_argument(
        "--pca",
        type=Path,
        required=True,
        help="Path to PCA components TSV file (sample_id, location, PC1, PC2, PC3).",
    )
    parser.add_argument(
        "--ibs",
        type=Path,
        required=True,
        help="Path to IBS matrix TSV file (symmetric matrix with sample IDs).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Path to save output figure (PNG format recommended).",
    )
    parser.add_argument(
        "--annotate_thresh",
        type=float,
        default=1.0,
        help="Threshold for annotating extreme samples (abs(PC1) or abs(PC2) >= thresh). Default: 1.0",
    )
    return parser.parse_args()


def load_pca_data(pca_path):
    """
    Load PCA components and filter out 'Blank' samples.
    
    Args:
        pca_path: Path to PCA components TSV file
        
    Returns:
        DataFrame with columns: sample_id, location, PC1, PC2, PC3
    """
    pca_df = pd.read_csv(pca_path, sep="\t")
    
    # Filter out Blank samples
    pca_df = pca_df[~pca_df["location"].str.contains("Blank", case=False, na=False)]
    
    return pca_df


def load_ibs_matrix(ibs_path):
    """
    Load and parse IBS matrix, handling potential unnamed trailing columns.
    
    Args:
        ibs_path: Path to IBS matrix TSV file
        
    Returns:
        DataFrame with sample IDs as both index and columns
    """
    # Read the file to inspect structure
    ibs_df = pd.read_csv(ibs_path, sep="\t", index_col=0)
    
    # Drop any unnamed or empty columns that may result from trailing tabs
    unnamed_cols = [col for col in ibs_df.columns if str(col).startswith("Unnamed")]
    if unnamed_cols:
        ibs_df = ibs_df.drop(columns=unnamed_cols)
    
    # Ensure the matrix is symmetric by keeping only samples present in both index and columns
    common_samples = ibs_df.index.intersection(ibs_df.columns)
    ibs_df = ibs_df.loc[common_samples, common_samples]
    
    # Filter out Blank samples
    non_blank_samples = [s for s in ibs_df.index if "Blank" not in str(s)]
    ibs_df = ibs_df.loc[non_blank_samples, non_blank_samples]
    
    return ibs_df


def cluster_ibs_matrix(ibs_df):
    """
    Perform hierarchical clustering on IBS matrix using distance = 1 - IBS.
    
    Args:
        ibs_df: IBS similarity matrix (samples x samples)
        
    Returns:
        Tuple of (reordered_ibs_df, clustering_success)
    """
    try:
        # Convert IBS similarity to distance
        distance_matrix = 1.0 - ibs_df.values
        
        # Ensure distance matrix is symmetric and valid
        distance_matrix = np.maximum(distance_matrix, 0)  # Clip negative values
        distance_matrix = (distance_matrix + distance_matrix.T) / 2  # Ensure symmetry
        
        # Convert to condensed distance matrix for linkage
        condensed_dist = squareform(distance_matrix, checks=False)
        
        # Perform hierarchical clustering with average linkage
        linkage_matrix = linkage(condensed_dist, method="average")
        
        # Get dendrogram to extract optimal leaf ordering
        dend = dendrogram(linkage_matrix, no_plot=True)
        cluster_order = dend["leaves"]
        
        # Reorder the IBS matrix
        new_order = ibs_df.index[cluster_order]
        reordered_ibs = ibs_df.loc[new_order, new_order]
        
        return reordered_ibs, True
        
    except Exception as e:
        print(f"Warning: Clustering failed ({e}). Using original order.", file=sys.stderr)
        return ibs_df, False


def plot_relatedness_figure(pca_df, ibs_df, output_path, annotate_thresh=1.0):
    """
    Create two-panel figure with PCA scatter and IBS heatmap.
    
    Args:
        pca_df: PCA components DataFrame
        ibs_df: IBS matrix DataFrame (should be pre-clustered if desired)
        output_path: Path to save the figure
        annotate_thresh: Threshold for annotating extreme samples
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Panel 1: PCA scatter plot
    unique_locations = pca_df["location"].unique()
    colors = sns.color_palette("tab10", n_colors=len(unique_locations))
    location_colors = dict(zip(unique_locations, colors))
    
    for location in unique_locations:
        loc_data = pca_df[pca_df["location"] == location]
        ax1.scatter(
            loc_data["PC1"],
            loc_data["PC2"],
            c=[location_colors[location]],
            label=location,
            alpha=0.7,
            s=50,
            edgecolors="black",
            linewidths=0.5,
        )
    
    # Annotate extreme samples
    extreme_samples = pca_df[
        (pca_df["PC1"].abs() >= annotate_thresh) | (pca_df["PC2"].abs() >= annotate_thresh)
    ]
    
    for _, row in extreme_samples.iterrows():
        ax1.annotate(
            row["sample_id"],
            (row["PC1"], row["PC2"]),
            fontsize=7,
            alpha=0.7,
            xytext=(5, 5),
            textcoords="offset points",
        )
    
    ax1.set_xlabel("PC1", fontsize=12)
    ax1.set_ylabel("PC2", fontsize=12)
    ax1.set_title("PCA of BAM-derived Variation", fontsize=14, fontweight="bold")
    ax1.legend(loc="best", fontsize=8, framealpha=0.9)
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0, color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
    ax1.axvline(x=0, color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
    
    # Panel 2: IBS heatmap
    sns.heatmap(
        ibs_df,
        ax=ax2,
        cmap="RdYlBu_r",
        vmin=0.8,
        vmax=1.0,
        cbar_kws={"label": "IBS Similarity"},
        xticklabels=False,
        yticklabels=False,
        square=True,
    )
    ax2.set_title("IBS Similarity Matrix (Clustered)", fontsize=14, fontweight="bold")
    ax2.set_xlabel("Samples", fontsize=12)
    ax2.set_ylabel("Samples", fontsize=12)
    
    plt.tight_layout()
    
    # Save figure
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    
    print(f"Figure saved to: {output_path}")


def main():
    """Main execution function."""
    args = parse_args()
    
    # Validate input files exist
    if not args.pca.exists():
        print(f"Error: PCA file not found: {args.pca}", file=sys.stderr)
        sys.exit(1)
    
    if not args.ibs.exists():
        print(f"Error: IBS file not found: {args.ibs}", file=sys.stderr)
        sys.exit(1)
    
    # Load data
    print(f"Loading PCA data from {args.pca}...")
    pca_df = load_pca_data(args.pca)
    print(f"  Loaded {len(pca_df)} samples after filtering blanks")
    
    print(f"Loading IBS matrix from {args.ibs}...")
    ibs_df = load_ibs_matrix(args.ibs)
    print(f"  Loaded {len(ibs_df)} x {len(ibs_df)} IBS matrix after filtering blanks")
    
    # Cluster IBS matrix
    print("Performing hierarchical clustering on IBS matrix...")
    ibs_clustered, clustering_success = cluster_ibs_matrix(ibs_df)
    if clustering_success:
        print("  Clustering successful")
    
    # Create figure
    print(f"Creating figure with annotate_thresh={args.annotate_thresh}...")
    plot_relatedness_figure(pca_df, ibs_clustered, args.out, args.annotate_thresh)
    
    print("Done!")


if __name__ == "__main__":
    main()
