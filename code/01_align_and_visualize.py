#!/usr/bin/env python3
"""
Low-coverage WGS alignment and connectedness visualization pipeline.

This script discovers all paired-end FASTQ files in `data/raw/`, aligns them
against the Olympia oyster reference genome found in `data/genome/`, performs
joint variant calling, and generates a summary visualization of genetic
variation and connectedness across sampling locations inferred from file
prefixes. All outputs, logs, and metadata are written to
`output/01_align_and_visualize/` following the repository execution
conventions defined in `INSTRUCTIONS.md`.

Steps
-----
1. Dependency validation (`bwa`, `samtools`, `bcftools`, `plink`, Python libs).
2. Reference preparation inside the output directory (relative symlink +
   indices).
3. Sample discovery and metadata table generation from FASTQ filenames.
4. Sequential alignment (BWA-MEM ➜ SAMtools sort/index) with conservative
   memory usage suitable for a 200 MB constraint.
5. Coverage metrics collection and aggregation.
6. Joint variant calling (bcftools mpileup/call/filter) with depth limits.
7. PCA and pairwise-IBS estimation via PLINK followed by matplotlib figures.
8. Metadata logging (`metadata.json`) describing inputs, parameters, runtimes.

Usage
-----
    python code/01_align_and_visualize.py \
        --threads 32 \
        --threads-per-sample 4 \
        --memory-mb 200

All paths are relative to the repository root; do not modify files in `data/`.
"""

import argparse
import io
import json
import logging
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Sample:
    """Container for paired FASTQ paths and inferred metadata."""

    sample_id: str
    location: str
    r1: Path
    r2: Path


def to_relative_path(path: Path, base: Path) -> str:
    """Return `path` as a string relative to `base` when possible."""
    try:
        return str(path.relative_to(base))
    except ValueError:
        return str(path)


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Align low-coverage WGS FASTQs to the Olympia oyster reference "
            "and generate summary connectedness visualizations."
        )
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=Path("data/raw"),
        help="Directory containing raw FASTQ files (read-only).",
    )
    parser.add_argument(
        "--genome-fasta",
        type=Path,
        default=Path("data/genome/Olurida_v081.fa"),
        help="Reference genome FASTA (read-only source).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/01_align_and_visualize"),
        help="Destination directory for pipeline outputs.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=min(50, multiprocessing.cpu_count()),
        help="Total threads available for multi-threaded steps (max 50).",
    )
    parser.add_argument(
        "--threads-per-sample",
        type=int,
        default=2,
        help="Threads allocated to per-sample alignment/sorting steps.",
    )
    parser.add_argument(
        "--memory-mb",
        type=int,
        default=200,
        help="Total memory budget in megabytes (default reflects constraint).",
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=20,
        help="Minimum mapping quality for coverage statistics.",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=200,
        help="Maximum depth per position passed to bcftools mpileup.",
    )
    parser.add_argument(
        "--min-total-depth",
        type=int,
        default=10,
        help="Minimum INFO/DP required to keep a variant site.",
    )
    parser.add_argument(
        "--min-qual",
        type=float,
        default=30.0,
        help="Minimum QUAL score required to retain a variant site.",
    )
    parser.add_argument(
        "--plink-min-contig-length",
        type=int,
        default=100000,
        help="Minimum contig length to include in PLINK conversion (filters excessive small scaffolds).",
    )
    parser.add_argument(
        "--plink-max-contigs",
        type=int,
        default=200,
        help="Maximum number of contigs to include for PLINK after applying length filter (longest retained).",
    )
    parser.add_argument(
        "--plink2",
        action="store_true",
        help="Use plink2 instead of plink (if installed) for PCA/IBS steps.",
    )
    parser.add_argument(
        "--skip-blanks",
        action="store_true",
        help="Exclude samples whose prefixes start with 'Blank'.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run steps even if outputs already exist.",
    )
    return parser.parse_args()


def configure_logging(log_path: Path) -> None:
    """Configure logging to file and stdout."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    handlers = [logging.FileHandler(log_path, mode="w"), logging.StreamHandler(sys.stdout)]
    logging.basicConfig(level=logging.INFO, format=log_format, handlers=handlers)


def check_dependencies(dependencies: Iterable[str]) -> None:
    """Ensure required external binaries are available."""
    missing = [exe for exe in dependencies if shutil.which(exe) is None]
    if missing:
        raise RuntimeError(
            "Missing required executables: "
            + ", ".join(missing)
            + ". Please install them and re-run."
        )


def safe_symlink(source: Path, target: Path) -> None:
    """
    Create a relative symlink `target` -> `source` without touching the original data.
    If the link already exists and points to the same file, it is reused.
    """
    if target.exists():
        if target.is_symlink() and target.resolve() == source.resolve():
            return
        raise FileExistsError(f"Cannot create symlink; {target} already exists.")
    rel_source = os.path.relpath(source, start=target.parent)
    target.symlink_to(rel_source)


def prepare_reference(
    genome_fasta: Path,
    reference_dir: Path,
    threads: int,
    force: bool,
) -> Path:
    """
    Symlink the reference FASTA into the output tree and build required indices.

    Returns the path to the symlinked FASTA to be used by downstream tools.
    """
    reference_dir.mkdir(parents=True, exist_ok=True)
    reference_fasta = reference_dir / genome_fasta.name
    if not reference_fasta.exists():
        logging.info("Symlinking reference genome into %s", reference_fasta.parent)
        safe_symlink(genome_fasta, reference_fasta)

    bwa_indices = [
        reference_fasta.with_name(reference_fasta.name + ext)
        for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]
    ]
    if force or not all(idx.exists() for idx in bwa_indices):
        logging.info("Building BWA index for %s", reference_fasta)
        run_command(
            ["bwa", "index", str(reference_fasta)],
            env={"OMP_NUM_THREADS": str(max(1, threads))},
        )

    fai_path = reference_fasta.with_name(reference_fasta.name + ".fai")
    if force or not fai_path.exists():
        logging.info("Creating FASTA index with samtools faidx.")
        run_command(["samtools", "faidx", str(reference_fasta)])

    dict_path = reference_fasta.with_suffix(".dict")
    if force or not dict_path.exists():
        logging.info("Creating sequence dictionary with samtools dict.")
        run_command(["samtools", "dict", str(reference_fasta), "-o", str(dict_path)])

    return reference_fasta


def discover_samples(raw_dir: Path, skip_blanks: bool) -> List[Sample]:
    """Identify paired FASTQ samples and derive location metadata."""
    r1_files = list(raw_dir.glob("*_R1_*.fastq.gz"))
    if not r1_files:
        raise FileNotFoundError(f"No R1 FASTQs found in {raw_dir}.")

    sample_map: Dict[str, Dict[str, Path]] = {}
    pattern = re.compile(r"_S\d+_")
    for r1 in r1_files:
        sample_id = pattern.split(r1.name)[0]
        r2 = raw_dir / r1.name.replace("_R1_", "_R2_")
        if not r2.exists():
            logging.warning("Missing mate pair for %s; skipping sample.", r1.name)
            continue
        location = sample_id.rsplit("_", 1)[0] if "_" in sample_id else sample_id
        if skip_blanks and location.lower().startswith("blank"):
            logging.info("Skipping blank/control sample %s", sample_id)
            continue
        sample_map[sample_id] = {"location": location, "R1": r1, "R2": r2}

    samples = [
        Sample(sample_id=sid, location=data["location"], r1=data["R1"], r2=data["R2"])
        for sid, data in sorted(sample_map.items())
    ]
    if not samples:
        raise RuntimeError("No valid paired samples discovered.")
    logging.info("Discovered %d samples across %d locations.", len(samples), len({s.location for s in samples}))
    return samples


def run_command(
    command: Iterable[str | Path] | str,
    *,
    cwd: Optional[Path] = None,
    env: Optional[Dict[str, str]] = None,
    capture_output: bool = False,
    text: bool = True,
) -> subprocess.CompletedProcess:
    """Wrapper around subprocess.run with logging and error handling.

    Accepts Path objects in the command iterable and coerces them to strings.
    """
    if isinstance(command, str):
        log_cmd = command
        cmd_list: list[str] = [command]
    else:
        # Coerce all parts to str to avoid TypeErrors when joining/logging
        cmd_list = [str(part) for part in command]
        log_cmd = " ".join(cmd_list)
    logging.info("Running command: %s", log_cmd)
    merged_env = os.environ.copy()
    if env:
        merged_env.update({k: str(v) for k, v in env.items()})
    result = subprocess.run(
        cmd_list if not isinstance(command, str) else command,
        cwd=str(cwd) if cwd else None,
        env=merged_env,
        check=True,
        capture_output=capture_output,
        text=text,
    )
    return result


def align_sample(
    sample: Sample,
    reference_fasta: Path,
    align_dir: Path,
    threads_per_sample: int,
    memory_mb: int,
    force: bool,
) -> Path:
    """Align paired FASTQs with BWA-MEM and return sorted BAM path."""
    align_dir.mkdir(parents=True, exist_ok=True)
    bam_path = align_dir / f"{sample.sample_id}.sorted.bam"
    if bam_path.exists() and not force:
        logging.info("Alignment already exists for %s; skipping.", sample.sample_id)
        return bam_path

    threads = max(1, threads_per_sample)
    per_thread_mem = max(32, memory_mb // max(1, threads))
    pipeline = (
        "set -euo pipefail\n"
        f"bwa mem -t {threads} {reference_fasta} {sample.r1} {sample.r2} | "
        f"samtools sort -@ {threads} -m {per_thread_mem}M -o {bam_path} -"
    )
    run_command(["bash", "-lc", pipeline], env={"OMP_NUM_THREADS": str(threads)})
    run_command(["samtools", "index", "-@", str(threads), str(bam_path)])
    return bam_path


def collect_coverage(
    bam_path: Path,
    min_mapq: int,
) -> Dict[str, float]:
    """Gather coverage metrics using samtools coverage."""
    result = run_command(
        [
            "samtools",
            "coverage",
            "-q",
            str(min_mapq),
            str(bam_path),
        ],
        capture_output=True,
        text=True,
    )
    df = pd.read_csv(
        io.StringIO(result.stdout),
        sep="\t",
    )
    df.rename(columns=lambda col: str(col).lstrip("#"), inplace=True)
    if df.empty:
        raise RuntimeError(f"No coverage data produced for {bam_path.name}.")

    df["length"] = df["endpos"] - df["startpos"] + 1
    genome_length = df["length"].sum()
    if genome_length == 0:
        raise RuntimeError(f"Coverage computation yielded zero length for {bam_path.name}.")

    weighted = lambda column: float(np.average(df[column], weights=df["length"]))

    metrics = {
        "bam": str(bam_path),
        "numreads": float(df["numreads"].sum()),
        "covbases": float(df["covbases"].sum()),
        "coverage": float(df["covbases"].sum() / genome_length),
        "meandepth": weighted("meandepth"),
        "meanbaseq": weighted("meanbaseq"),
        "meanmapq": weighted("meanmapq"),
        "genome_length": float(genome_length),
    }
    return metrics


def write_table(rows: List[Dict[str, object]], output_path: Path) -> None:
    """Write a list of dictionaries to TSV using pandas."""
    if not rows:
        return
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep="\t", index=False)


def write_sample_sheet(samples: List[Sample], output_path: Path, repo_root: Path) -> None:
    """Persist sample metadata for downstream reuse."""
    rows = []
    for sample in samples:
        rows.append(
            {
                "sample_id": sample.sample_id,
                "location": sample.location,
                "r1": to_relative_path(sample.r1, repo_root),
                "r2": to_relative_path(sample.r2, repo_root),
            }
        )
    write_table(rows, output_path)


def call_variants(
    bam_paths: List[Path],
    reference_fasta: Path,
    variant_dir: Path,
    threads: int,
    max_depth: int,
    min_total_depth: int,
    min_qual: float,
    force: bool,
) -> Tuple[Path, Path]:
    """Run bcftools mpileup/call/filter on the provided BAM list."""
    variant_dir.mkdir(parents=True, exist_ok=True)
    bam_list_path = variant_dir / "bam_list.txt"
    bam_list_path.write_text("\n".join(str(bam) for bam in bam_paths) + "\n")

    raw_vcf = variant_dir / "all_samples.raw.vcf.gz"
    filtered_vcf = variant_dir / "all_samples.filtered.vcf.gz"

    if (raw_vcf.exists() and filtered_vcf.exists()) and not force:
        logging.info("Variant calling outputs present; skipping mpileup/call.")
        return raw_vcf, filtered_vcf

    # If raw VCF already exists (indexed or not) but filtered does not, skip mpileup and proceed
    if raw_vcf.exists() and not filtered_vcf.exists() and not force:
        logging.info("Raw VCF present; skipping mpileup/call and proceeding to filtering.")
        # Ensure raw VCF is indexed before filtering
        csi = raw_vcf.with_suffix(raw_vcf.suffix + ".csi")
        tbi = raw_vcf.with_suffix(raw_vcf.suffix + ".tbi")
        if not (csi.exists() or tbi.exists()):
            run_command(["bcftools", "index", "--threads", str(threads), raw_vcf])
        filter_cmd = [
            "bcftools",
            "filter",
            "--threads",
            str(threads),
            "-i",
            f"QUAL>{min_qual} && INFO/DP>{min_total_depth}",
            "-Oz",
            "-o",
            str(filtered_vcf),
            str(raw_vcf),
        ]
        run_command(filter_cmd)
        run_command(["bcftools", "index", "--threads", str(threads), filtered_vcf])
        return raw_vcf, filtered_vcf

    mpileup_cmd = (
        "set -euo pipefail\n"
        f"bcftools mpileup --threads {threads} "
        f"-Ou -f {reference_fasta} "
        f"--annotate FORMAT/AD,FORMAT/DP "
        f"--max-depth {max_depth} "
        f"-b {bam_list_path} "
        "| "
        f"bcftools call --threads {threads} -mv -Oz -o {raw_vcf}"
    )
    run_command(["bash", "-lc", mpileup_cmd])
    run_command(["bcftools", "index", "--threads", str(threads), raw_vcf])

    filter_cmd = [
        "bcftools",
        "filter",
        "--threads",
        str(threads),
        "-i",
        f"QUAL>{min_qual} && INFO/DP>{min_total_depth}",
        "-Oz",
        "-o",
        str(filtered_vcf),
        str(raw_vcf),
    ]
    run_command(filter_cmd)
    run_command(["bcftools", "index", "--threads", str(threads), filtered_vcf])
    return raw_vcf, filtered_vcf


def run_plink_pca(
    filtered_vcf: Path,
    variant_dir: Path,
    threads: int,
    force: bool,
    plink_min_contig_length: int,
    plink_max_contigs: int,
    use_plink2: bool,
) -> Tuple[Path, Path]:
    """Convert VCF to PLINK format and compute PCA + IBS matrices."""
    prefix = variant_dir / "plink_dataset"
    pca_prefix = variant_dir / "pca"
    ibs_prefix = variant_dir / "ibs"

    bed_file = prefix.with_suffix(".bed")
    eigenvec = pca_prefix.with_suffix(".eigenvec")
    ibs_matrix = ibs_prefix.with_suffix(".mibs")

    if not (bed_file.exists() and eigenvec.exists() and ibs_matrix.exists()) or force:
        # Subset VCF to long contigs to satisfy PLINK limits on distinct chromosome names.
        subset_vcf = variant_dir / "subset_for_plink.vcf.gz"
        subset_index_csi = subset_vcf.with_suffix(".csi")
        if not subset_vcf.exists() or force:
            logging.info(
                "Creating PLINK subset VCF retaining contigs length >= %d", plink_min_contig_length
            )
            # Extract contig lengths from header, select those meeting threshold.
            header_proc = run_command([
                "bcftools",
                "view",
                "-h",
                str(filtered_vcf),
            ], capture_output=True, text=True)
            contigs_with_len: list[tuple[str,int]] = []
            for line in header_proc.stdout.splitlines():
                if line.startswith("##contig="):
                    # Format: ##contig=<ID=Contig0,length=116746>
                    try:
                        inside = line.split("<", 1)[1].rsplit(">", 1)[0]
                        parts = dict(
                            kv.split("=") for kv in inside.split(",") if "=" in kv
                        )
                        cid = parts.get("ID")
                        length_val = int(parts.get("length", "0"))
                        if cid and length_val >= plink_min_contig_length:
                            contigs_with_len.append((cid,length_val))
                    except Exception:
                        continue
            if not contigs_with_len:
                raise RuntimeError(
                    "No contigs meet length threshold for PLINK conversion; reduce --plink-min-contig-length."
                )
            # Sort by length desc and retain top N
            contigs_with_len.sort(key=lambda x: x[1], reverse=True)
            selected = [cid for cid,_ in contigs_with_len[:plink_max_contigs]]
            logging.info(
                "Selected %d contigs (threshold=%d, max=%d). Shortest retained length=%d.",
                len(selected),
                plink_min_contig_length,
                plink_max_contigs,
                contigs_with_len[min(len(contigs_with_len), plink_max_contigs)-1][1],
            )
            regions_file = variant_dir / "plink_contigs.txt"
            regions_file.write_text("\n".join(selected) + "\n")
            view_cmd = [
                "bcftools",
                "view",
                "-Oz",
                "-o",
                str(subset_vcf),
                "-r",
                ",".join(selected),
                str(filtered_vcf),
            ]
            run_command(view_cmd)
            run_command([
                "bcftools",
                "index",
                "--threads",
                str(threads),
                str(subset_vcf),
            ])
        plink_input_vcf = subset_vcf if subset_vcf.exists() else filtered_vcf
        plink_bin = "plink2" if use_plink2 and shutil.which("plink2") else "plink"
        if use_plink2 and plink_bin != "plink2":
            logging.warning("--plink2 requested but 'plink2' not found; falling back to 'plink'.")
        convert_cmd = [
            plink_bin,
            "--vcf",
            str(plink_input_vcf),
            "--allow-extra-chr",
            "--double-id",
            "--set-missing-var-ids",
            "@:#:\\$1:\\$2",
            "--make-bed",
            "--out",
            str(prefix),
        ]
        run_command(convert_cmd, cwd=variant_dir)

        pca_cmd = [
            plink_bin,
            "--bfile",
            str(prefix),
            "--allow-extra-chr",
            "--geno",
            "0.5",
            "--mind",
            "0.9",
            "--pca",
            "10",
            "header",
            "--threads",
            str(threads),
            "--out",
            str(pca_prefix),
        ]
        run_command(pca_cmd, cwd=variant_dir)

        ibs_cmd = [
            plink_bin,
            "--bfile",
            str(prefix),
            "--allow-extra-chr",
            "--geno",
            "0.5",
            "--mind",
            "0.9",
            "--distance",
            "square",
            "ibs",
            "--threads",
            str(threads),
            "--out",
            str(ibs_prefix),
        ]
        run_command(ibs_cmd, cwd=variant_dir)

    return eigenvec, ibs_matrix


def generate_visualization(
    eigenvec_path: Path,
    ibs_matrix_path: Path,
    samples: List[Sample],
    figure_dir: Path,
) -> Path:
    """Create PCA scatter plot and IBS heatmap summarizing connectedness."""
    figure_dir.mkdir(parents=True, exist_ok=True)
    sample_lookup = {sample.sample_id: sample.location for sample in samples}

    eigenvec = pd.read_csv(
        eigenvec_path,
        delim_whitespace=True,
        header=0,
    )
    eigenvec.rename(columns={"#FID": "fid", "FID": "fid", "IID": "iid"}, inplace=True)
    eigenvec["location"] = eigenvec["iid"].map(lambda sid: sample_lookup.get(sid, "unknown"))
    if "PC1" not in eigenvec.columns or "PC2" not in eigenvec.columns:
        raise RuntimeError("PCA results missing PC1/PC2 columns.")

    ibs_ids = []
    with (ibs_matrix_path.with_suffix(".mibs.id")).open() as handle:
        for line in handle:
            parts = line.strip().split()
            if parts:
                ibs_ids.append(parts[-1])
    ibs_matrix = np.atleast_2d(np.loadtxt(ibs_matrix_path))
    if ibs_matrix.shape[0] != len(ibs_ids):
        raise RuntimeError("Mismatch between IBS matrix and ID file.")
    # Convert IBS distance to similarity score (1 - distance).
    ibs_similarity = 1.0 - ibs_matrix

    locations = sorted({sample.location for sample in samples})

    try:
        import matplotlib.cm as cm
        import matplotlib.pyplot as plt
    except ImportError as err:
        raise ImportError("matplotlib is required for visualization.") from err

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    cmap = cm.get_cmap("tab20", len(locations))
    color_assignments = {loc: cmap(i) for i, loc in enumerate(locations)}

    axes[0].set_title("Genomic PCA (PC1 vs PC2)")
    for location, group in eigenvec.groupby("location"):
        axes[0].scatter(
            group["PC1"],
            group["PC2"],
            label=location,
            s=30,
            alpha=0.8,
            color=color_assignments.get(location, "grey"),
        )
    axes[0].set_xlabel("PC1")
    axes[0].set_ylabel("PC2")
    axes[0].legend(loc="best", fontsize="small", frameon=False)

    im = axes[1].imshow(ibs_similarity, cmap="viridis", vmin=0, vmax=1)
    axes[1].set_title("Pairwise IBS Similarity")
    axes[1].set_xticks(range(len(ibs_ids)))
    axes[1].set_yticks(range(len(ibs_ids)))
    axes[1].set_xticklabels(ibs_ids, rotation=90, fontsize=6)
    axes[1].set_yticklabels(ibs_ids, fontsize=6)
    fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04, label="IBS similarity (1 - distance)")
    fig.tight_layout()

    figure_path = figure_dir / "genetic_connectedness.png"
    fig.savefig(figure_path, dpi=300)
    plt.close(fig)
    return figure_path


def assemble_metadata(
    script_start: float,
    params: Dict[str, object],
    outputs: Dict[str, List[Path]],
    metadata_path: Path,
    repo_root: Path,
) -> None:
    """Write metadata.json capturing run context."""
    relative_outputs: Dict[str, List[str]] = {}
    for key, paths in outputs.items():
        relative_outputs[key] = [to_relative_path(path, repo_root) for path in paths]

    metadata = {
        "script": Path(__file__).name,
        "date": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": time.time() - script_start,
        "parameters": params,
        "outputs": relative_outputs,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2))


def main() -> None:
    args = parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    os.chdir(repo_root)

    output_dir = args.output_dir if args.output_dir.is_absolute() else repo_root / args.output_dir
    align_dir = output_dir / "alignments"
    variant_dir = output_dir / "variants"
    metrics_dir = output_dir / "metrics"
    figure_dir = output_dir / "figures"
    reference_dir = output_dir / "reference"
    log_dir = output_dir / "logs"
    tmp_dir = output_dir / "tmp"
    output_dir.mkdir(parents=True, exist_ok=True)
    align_dir.mkdir(parents=True, exist_ok=True)
    variant_dir.mkdir(parents=True, exist_ok=True)
    metrics_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    configure_logging(log_dir / "pipeline.log")
    script_start = time.time()

    try:
        check_dependencies(["bwa", "samtools", "bcftools", "plink"])

        raw_dir = args.raw_dir if args.raw_dir.is_absolute() else repo_root / args.raw_dir
        genome_fasta = args.genome_fasta if args.genome_fasta.is_absolute() else repo_root / args.genome_fasta

        samples = discover_samples(raw_dir, skip_blanks=args.skip_blanks)
        write_sample_sheet(samples, metrics_dir / "sample_metadata.tsv", repo_root)

        total_threads = max(1, min(args.threads, 50))
        threads_per_sample = max(1, min(args.threads_per_sample, total_threads))

        reference_fasta = prepare_reference(
            genome_fasta,
            reference_dir,
            threads=total_threads,
            force=args.force,
        )

        bam_paths: List[Path] = []
        coverage_records: List[Dict[str, object]] = []
        coverage_summary_path = metrics_dir / "coverage_summary.tsv"
        
        # Check if we can skip coverage calculation
        if coverage_summary_path.exists() and not args.force:
            logging.info("Coverage summary exists; skipping coverage calculation.")
            # Still need to collect BAM paths for variant calling
            for sample in samples:
                bam_path = align_dir / f"{sample.sample_id}.sorted.bam"
                bam_paths.append(bam_path)
        else:
            for sample in samples:
                logging.info("Processing sample %s (%s).", sample.sample_id, sample.location)
                bam_path = align_sample(
                    sample,
                    reference_fasta,
                    align_dir,
                    threads_per_sample=threads_per_sample,
                    memory_mb=args.memory_mb,
                    force=args.force,
                )
                bam_paths.append(bam_path)
                coverage_metrics = collect_coverage(bam_path, min_mapq=args.min_mapq)
                coverage_metrics.update({"sample_id": sample.sample_id, "location": sample.location})
                coverage_records.append(coverage_metrics)

            write_table(coverage_records, coverage_summary_path)

        raw_vcf, filtered_vcf = call_variants(
            bam_paths,
            reference_fasta,
            variant_dir,
            threads=total_threads,
            max_depth=args.max_depth,
            min_total_depth=args.min_total_depth,
            min_qual=args.min_qual,
            force=args.force,
        )

        eigenvec_path, ibs_matrix_path = run_plink_pca(
            filtered_vcf,
            variant_dir,
            threads=total_threads,
            force=args.force,
            plink_min_contig_length=args.plink_min_contig_length,
            plink_max_contigs=args.plink_max_contigs,
            use_plink2=args.plink2,
        )

        figure_path = generate_visualization(
            eigenvec_path,
            ibs_matrix_path,
            samples,
            figure_dir,
        )

        assemble_metadata(
            script_start,
            params={
                "raw_dir": to_relative_path(raw_dir, repo_root),
                "genome_fasta": to_relative_path(genome_fasta, repo_root),
                "threads": total_threads,
                "threads_per_sample": threads_per_sample,
                "memory_mb": args.memory_mb,
                "min_mapq": args.min_mapq,
                "max_depth": args.max_depth,
                "min_total_depth": args.min_total_depth,
                "min_qual": args.min_qual,
                "skip_blanks": args.skip_blanks,
                "force": args.force,
            },
            outputs={
                "alignments": bam_paths,
                "variants": [raw_vcf, filtered_vcf],
                "figures": [figure_path],
                "metrics": [
                    metrics_dir / "sample_metadata.tsv",
                    metrics_dir / "coverage_summary.tsv",
                ],
            },
            metadata_path=output_dir / "metadata.json",
            repo_root=repo_root,
        )

        logging.info("Pipeline complete. Outputs written to %s", output_dir)

    except Exception as exc:  # pylint: disable=broad-except
        logging.exception("Pipeline failed: %s", exc)
        raise


if __name__ == "__main__":
    main()

