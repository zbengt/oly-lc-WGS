#!/usr/bin/env python3
"""
Variant summary and quality assessment pipeline (VCF-driven).

This optional step complements `01_align_and_visualize.py` by generating a suite
of variant-level and per-sample metrics to assess genetic variation and data
quality. It consumes the joint variant calls produced in the previous step and
outputs summary statistics, tables, and visualizations to
`output/03_variant_summary/`.

Analyses performed:
1. `bcftools stats` over the filtered multi-sample VCF.
2. PLINK-based per-sample heterozygosity (`--het`) and missingness (`--missing`).
3. PLINK allele frequency stratification by sampling location (`--freq --within`).
4. Aggregated tables and plots (heterozygosity and missingness by location,
   allele-frequency distributions) to quickly assess variation and data quality.

All paths are resolved relative to the repository root. Inputs within `data/`
and `output/01_align_and_visualize/` are treated as read-only. Outputs, logs,
and metadata are saved under `output/03_variant_summary/` in accordance with the
conventions in `INSTRUCTIONS.md`.
"""

import argparse
import json
import logging
import multiprocessing
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description="Summarize variant-level diversity and per-sample quality metrics."
    )
    parser.add_argument(
        "--input-vcf",
        type=Path,
        default=Path("output/01_align_and_visualize/variants/all_samples.filtered.vcf.gz"),
        help="Multi-sample filtered VCF produced by step 01.",
    )
    parser.add_argument(
        "--plink-prefix",
        type=Path,
        default=Path("output/01_align_and_visualize/variants/plink_dataset"),
        help="Prefix of PLINK binary files (bed/bim/fam) produced by step 01.",
    )
    parser.add_argument(
        "--sample-metadata",
        type=Path,
        default=Path("output/01_align_and_visualize/metrics/sample_metadata.tsv"),
        help="Sample metadata table with sample_id and location columns.",
    )
    parser.add_argument(
        "--reference-fasta",
        type=Path,
        default=Path("data/genome/Olurida_v081.fa"),
        help="Reference FASTA used for variant calling (required for bcftools stats Ts/Tv).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/03_variant_summary"),
        help="Destination directory for reports, tables, and figures.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=min(50, multiprocessing.cpu_count()),
        help="Threads allocated to bcftools/plink (max 50).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run analyses even if outputs already exist.",
    )
    return parser.parse_args()


def configure_logging(log_path: Path) -> None:
    """Set up logging to file and stdout."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    handlers = [logging.FileHandler(log_path, mode="w"), logging.StreamHandler(sys.stdout)]
    logging.basicConfig(level=logging.INFO, format=log_format, handlers=handlers)


def check_dependencies(dependencies: Iterable[str]) -> None:
    """Ensure required externals are available."""
    missing = [exe for exe in dependencies if shutil.which(exe) is None]
    if missing:
        raise RuntimeError(
            "Missing required executables: "
            + ", ".join(missing)
            + ". Install them before rerunning."
        )


def run_command(
    command: Iterable[str] | str,
    *,
    cwd: Path | None = None,
    env: Dict[str, str] | None = None,
) -> None:
    """Execute command with logging and error propagation."""
    if isinstance(command, str):
        log_cmd = command
    else:
        log_cmd = " ".join(command)
    logging.info("Running command: %s", log_cmd)
    merged_env = os.environ.copy()
    if env:
        merged_env.update({k: str(v) for k, v in env.items()})
    subprocess.run(
        command,
        cwd=str(cwd) if cwd else None,
        env=merged_env,
        check=True,
    )


def to_relative_path(path: Path, base: Path) -> str:
    """Return `path` as a string relative to `base` when possible."""
    try:
        return str(path.relative_to(base))
    except ValueError:
        return str(path)


def load_sample_metadata(sample_metadata_path: Path) -> pd.DataFrame:
    """Load sample metadata TSV containing sample_id and location."""
    if not sample_metadata_path.exists():
        raise FileNotFoundError(f"Sample metadata not found: {sample_metadata_path}")
    df = pd.read_csv(sample_metadata_path, sep="\t")
    required_cols = {"sample_id", "location"}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise ValueError(
            f"Sample metadata missing required columns: {', '.join(sorted(missing_cols))}"
        )
    return df


def write_table(df: pd.DataFrame, path: Path) -> None:
    """Write a dataframe to TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def run_bcftools_stats(
    vcf_path: Path,
    reference_fasta: Path,
    stats_path: Path,
    threads: int,
    force: bool,
) -> Path:
    """Execute bcftools stats and return the output path."""
    if stats_path.exists() and not force:
        logging.info("bcftools stats output exists; skipping.")
        return stats_path

    stats_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "bcftools",
        "stats",
        "--threads",
        str(max(1, threads)),
        "-F",
        str(reference_fasta),
        "-s",
        "-",
        str(vcf_path),
    ]
    with stats_path.open("w") as handle:
        logging.info("Running bcftools stats with output -> %s", stats_path)
        subprocess.run(cmd, check=True, stdout=handle)
    return stats_path


def create_within_file(metadata: pd.DataFrame, cluster_path: Path) -> Path:
    """Create PLINK --within cluster file mapping samples to locations."""
    cluster_path.parent.mkdir(parents=True, exist_ok=True)
    with cluster_path.open("w") as handle:
        handle.write("FID IID CLUSTER\n")
        for _, row in metadata.iterrows():
            sample_id = str(row["sample_id"])
            location = str(row["location"])
            handle.write(f"{sample_id} {sample_id} {location}\n")
    return cluster_path


def run_plink_stats(
    plink_prefix: Path,
    cluster_file: Path,
    plink_dir: Path,
    threads: int,
    force: bool,
) -> Tuple[Path, Path, Path]:
    """Run PLINK commands for heterozygosity, missingness, and stratified frequencies."""
    plink_dir.mkdir(parents=True, exist_ok=True)
    heter_prefix = plink_dir / "heterozygosity"
    missing_prefix = plink_dir / "missingness"
    freq_prefix = plink_dir / "frequency"

    heter_file = heter_prefix.with_suffix(".het")
    missing_file = missing_prefix.with_suffix(".imiss")
    freq_strat_file = freq_prefix.with_suffix(".frq.strat")

    base_args = [
        "plink",
        "--bfile",
        str(plink_prefix),
        "--allow-extra-chr",
        "--threads",
        str(max(1, threads)),
        "--double-id",
    ]

    if (heter_file.exists() and missing_file.exists() and freq_strat_file.exists()) and not force:
        logging.info("PLINK outputs exist; skipping recomputation.")
        return heter_file, missing_file, freq_strat_file

    run_command(
        base_args + ["--het", "--out", str(heter_prefix)]
    )

    run_command(
        base_args + ["--missing", "--out", str(missing_prefix)]
    )

    run_command(
        base_args
        + [
            "--freq",
            "--within",
            str(cluster_file),
            "--family",
            "--out",
            str(freq_prefix),
        ]
    )

    return heter_file, missing_file, freq_strat_file


def summarise_heterozygosity(
    heter_path: Path,
    metadata: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Compute heterozygosity rates per sample and per location."""
    het_df = pd.read_csv(heter_path, delim_whitespace=True)
    het_df.rename(columns={"IID": "sample_id"}, inplace=True)
    het_df = het_df.merge(metadata, on="sample_id", how="left")
    het_df["heterozygosity_rate"] = 1.0 - (het_df["O(HOM)"] / het_df["N(NM)"])
    het_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    het_df["heterozygosity_rate"].fillna(0.0, inplace=True)

    loc_df = (
        het_df.groupby("location")
        .agg(
            samples=("sample_id", "count"),
            mean_heterozygosity=("heterozygosity_rate", "mean"),
            std_heterozygosity=("heterozygosity_rate", "std"),
        )
        .reset_index()
    )
    loc_df["std_heterozygosity"].fillna(0.0, inplace=True)
    return het_df, loc_df


def summarise_missingness(
    missing_path: Path,
    metadata: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Compute missingness per sample and per location."""
    miss_df = pd.read_csv(missing_path, delim_whitespace=True)
    miss_df.rename(columns={"IID": "sample_id"}, inplace=True)
    miss_df = miss_df.merge(metadata, on="sample_id", how="left")
    miss_df["call_rate"] = 1.0 - miss_df["F_MISS"]

    loc_df = (
        miss_df.groupby("location")
        .agg(
            samples=("sample_id", "count"),
            mean_call_rate=("call_rate", "mean"),
            std_call_rate=("call_rate", "std"),
        )
        .reset_index()
    )
    loc_df["std_call_rate"].fillna(0.0, inplace=True)
    return miss_df, loc_df


def summarise_frequency(freq_path: Path) -> pd.DataFrame:
    """Aggregate stratified allele frequencies to per-location summaries."""
    freq_df = pd.read_csv(freq_path, delim_whitespace=True)
    required_cols = {"CLST", "MAF", "NCHROBS"}
    missing = required_cols - set(freq_df.columns)
    if missing:
        raise ValueError(
            f"Expected columns missing from PLINK frequency output: {', '.join(sorted(missing))}"
        )
    freq_df.rename(columns={"CLST": "location"}, inplace=True)
    summary = (
        freq_df.groupby("location")
        .agg(
            variants=("SNP", "nunique"),
            mean_maf=("MAF", "mean"),
            median_maf=("MAF", "median"),
            sd_maf=("MAF", "std"),
        )
        .reset_index()
    )
    summary["sd_maf"].fillna(0.0, inplace=True)
    return summary


def plot_location_metric(
    df: pd.DataFrame,
    location_col: str,
    metric_col: str,
    ylabel: str,
    title: str,
    figure_path: Path,
) -> Path:
    """Generate a simple bar plot of location-level metrics."""
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib.pyplot as plt
    except ImportError as err:
        raise ImportError("matplotlib is required for visualization.") from err

    df_sorted = df.sort_values(metric_col, ascending=False)
    plt.figure(figsize=(10, 6))
    plt.bar(df_sorted[location_col], df_sorted[metric_col], color="#2878B5")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(figure_path, dpi=300)
    plt.close()
    return figure_path


def assemble_metadata(
    script_start: float,
    params: Dict[str, str],
    outputs: Dict[str, List[Path]],
    metadata_path: Path,
    repo_root: Path,
) -> None:
    """Write metadata JSON capturing run details."""
    metadata = {
        "script": Path(__file__).name,
        "date": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": time.time() - script_start,
        "parameters": params,
        "outputs": {
            key: [to_relative_path(path, repo_root) for path in paths]
            for key, paths in outputs.items()
        },
    }
    metadata_path.write_text(json.dumps(metadata, indent=2))


def main() -> None:
    args = parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    os.chdir(repo_root)

    output_dir = args.output_dir if args.output_dir.is_absolute() else repo_root / args.output_dir
    stats_dir = output_dir / "stats"
    tables_dir = output_dir / "tables"
    figures_dir = output_dir / "figures"
    logs_dir = output_dir / "logs"
    tmp_dir = output_dir / "tmp"

    output_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    configure_logging(logs_dir / "pipeline.log")
    script_start = time.time()

    try:
        check_dependencies(["bcftools", "plink"])

        vcf_path = args.input_vcf if args.input_vcf.is_absolute() else repo_root / args.input_vcf
        plink_prefix = args.plink_prefix if args.plink_prefix.is_absolute() else repo_root / args.plink_prefix
        sample_metadata_path = (
            args.sample_metadata if args.sample_metadata.is_absolute() else repo_root / args.sample_metadata
        )
        reference_fasta = (
            args.reference_fasta if args.reference_fasta.is_absolute() else repo_root / args.reference_fasta
        )

        for path in [
            vcf_path,
            plink_prefix.with_suffix(".bed"),
            plink_prefix.with_suffix(".bim"),
            plink_prefix.with_suffix(".fam"),
        ]:
            if not path.exists():
                raise FileNotFoundError(f"Required input not found: {path}")

        metadata_df = load_sample_metadata(sample_metadata_path)

        cluster_file = create_within_file(metadata_df, tmp_dir / "sample_clusters.txt")

        stats_path = run_bcftools_stats(
            vcf_path,
            reference_fasta,
            stats_dir / "bcftools_stats.txt",
            threads=max(1, min(args.threads, 50)),
            force=args.force,
        )

        heter_path, missing_path, freq_path = run_plink_stats(
            plink_prefix,
            cluster_file,
            stats_dir / "plink",
            threads=max(1, min(args.threads, 50)),
            force=args.force,
        )

        het_sample_df, het_location_df = summarise_heterozygosity(heter_path, metadata_df)
        miss_sample_df, miss_location_df = summarise_missingness(missing_path, metadata_df)
        freq_location_df = summarise_frequency(freq_path)

        write_table(het_sample_df, tables_dir / "heterozygosity_per_sample.tsv")
        write_table(het_location_df, tables_dir / "heterozygosity_by_location.tsv")
        write_table(miss_sample_df, tables_dir / "missingness_per_sample.tsv")
        write_table(miss_location_df, tables_dir / "missingness_by_location.tsv")
        write_table(freq_location_df, tables_dir / "allele_frequency_by_location.tsv")

        heter_fig = plot_location_metric(
            het_location_df,
            "location",
            "mean_heterozygosity",
            "Mean heterozygosity",
            "Per-location heterozygosity",
            figures_dir / "heterozygosity_by_location.png",
        )

        miss_fig = plot_location_metric(
            miss_location_df,
            "location",
            "mean_call_rate",
            "Mean call rate (1 - missingness)",
            "Per-location call rate",
            figures_dir / "call_rate_by_location.png",
        )

        assemble_metadata(
            script_start,
            params={
                "input_vcf": to_relative_path(vcf_path, repo_root),
                "plink_prefix": to_relative_path(plink_prefix, repo_root),
                "sample_metadata": to_relative_path(sample_metadata_path, repo_root),
                "reference_fasta": to_relative_path(reference_fasta, repo_root),
                "threads": max(1, min(args.threads, 50)),
                "force": args.force,
            },
            outputs={
                "stats": [stats_path, heter_path, missing_path, freq_path],
                "tables": [
                    tables_dir / "heterozygosity_per_sample.tsv",
                    tables_dir / "heterozygosity_by_location.tsv",
                    tables_dir / "missingness_per_sample.tsv",
                    tables_dir / "missingness_by_location.tsv",
                    tables_dir / "allele_frequency_by_location.tsv",
                ],
                "figures": [heter_fig, miss_fig],
                "logs": [logs_dir / "pipeline.log"],
            },
            metadata_path=output_dir / "metadata.json",
            repo_root=repo_root,
        )

        logging.info("Variant summary pipeline complete. Outputs written to %s", output_dir)

    except Exception as exc:  # pylint: disable=broad-except
        logging.exception("Variant summary pipeline failed: %s", exc)
        raise


if __name__ == "__main__":
    main()

