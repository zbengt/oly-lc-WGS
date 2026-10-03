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

from __future__ import annotations

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
from concurrent.futures import ThreadPoolExecutor
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
        "--min-mean-depth",
        type=float,
        default=1.0,
        help=(
            "Exclude samples whose genome-wide mean depth in metrics/coverage_summary.tsv "
            "is below this from the PLINK PCA/IBS and figure (0 disables)."
        ),
    )
    parser.add_argument(
        "--fix-read-groups",
        action="store_true",
        help=(
            "Add an @RG line (ID and SM set to the sample ID) in place to existing BAMs "
            "that lack one, using `samtools addreplacerg`, before variant calling."
        ),
    )
    parser.add_argument(
        "--force-variants",
        action="store_true",
        help=(
            "Re-run bcftools mpileup/call/filter and everything downstream even if the "
            "VCFs exist. Alignments and reference indices are kept."
        ),
    )
    parser.add_argument(
        "--force-plink",
        action="store_true",
        help="Re-run the PLINK subset, conversion, PCA, IBS, and figure from the existing filtered VCF.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run every step, including alignment and reference indexing, even if outputs exist.",
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


def load_samples_from_sheet(sheet_path: Path, repo_root: Path, skip_blanks: bool) -> List[Sample]:
    """Rebuild the sample list from a previous run's `sample_metadata.tsv`.

    Used when `data/raw/` holds no FASTQs (for example on a machine restored from
    the archived outputs). The FASTQ paths are carried along but need not exist;
    alignment cannot be re-run in this mode.
    """
    sheet = pd.read_csv(sheet_path, sep="\t")
    required = {"sample_id", "location", "r1", "r2"}
    missing_cols = required - set(sheet.columns)
    if missing_cols:
        raise ValueError(f"{sheet_path} lacks columns: {', '.join(sorted(missing_cols))}")
    samples: List[Sample] = []
    for row in sheet.itertuples(index=False):
        if skip_blanks and str(row.location).lower().startswith("blank"):
            logging.info("Skipping blank/control sample %s", row.sample_id)
            continue
        samples.append(
            Sample(
                sample_id=str(row.sample_id),
                location=str(row.location),
                r1=repo_root / str(row.r1),
                r2=repo_root / str(row.r2),
            )
        )
    if not samples:
        raise RuntimeError(f"No samples loaded from {sheet_path}.")
    logging.info(
        "Loaded %d samples across %d locations from %s.",
        len(samples), len({s.location for s in samples}), sheet_path,
    )
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


def bam_has_read_group(bam_path: Path, sample_id: str) -> bool:
    """Return True if the BAM header carries an @RG line with SM:<sample_id>."""
    try:
        result = subprocess.run(
            ["samtools", "view", "-H", str(bam_path)],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError:
        return False
    for line in result.stdout.splitlines():
        if line.startswith("@RG") and f"\tSM:{sample_id}" in line:
            return True
    return False


def add_read_group(bam_path: Path, sample_id: str, threads: int) -> None:
    """Rewrite `bam_path` in place with an @RG line naming the sample, then re-index."""
    tmp_bam = bam_path.with_name(bam_path.name + ".rg.tmp")
    read_group = f"@RG\tID:{sample_id}\tSM:{sample_id}\tPL:ILLUMINA"
    run_command(
        [
            "samtools", "addreplacerg",
            "-@", str(max(1, threads)),
            "-r", read_group,
            "-m", "overwrite_all",
            "-O", "BAM",
            "-o", str(tmp_bam),
            str(bam_path),
        ]
    )
    tmp_bam.replace(bam_path)
    run_command(["samtools", "index", "-@", str(max(1, threads)), str(bam_path)])


def ensure_read_groups(
    samples: List[Sample],
    bam_paths: List[Path],
    fix: bool,
    workers: int,
    threads_per_bam: int,
) -> None:
    """Check every BAM for an @RG line; add one (in parallel) when `fix` is set."""
    missing = [
        (sample, bam)
        for sample, bam in zip(samples, bam_paths)
        if not bam_has_read_group(bam, sample.sample_id)
    ]
    if not missing:
        logging.info("All %d BAMs carry an @RG line with the sample ID.", len(bam_paths))
        return
    if not fix:
        logging.warning(
            "%d of %d BAMs have no @RG line; bcftools names such samples by file path. "
            "Re-run with --fix-read-groups to add them in place.",
            len(missing),
            len(bam_paths),
        )
        return
    logging.info("Adding read groups to %d BAMs with %d parallel workers.", len(missing), workers)
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        list(
            pool.map(
                lambda pair: add_read_group(pair[1], pair[0].sample_id, threads_per_bam),
                missing,
            )
        )


def low_depth_samples(coverage_summary_path: Path, min_mean_depth: float) -> List[str]:
    """Sample IDs whose mean depth in the coverage summary falls below the cutoff."""
    if min_mean_depth <= 0 or not coverage_summary_path.exists():
        return []
    cov = pd.read_csv(coverage_summary_path, sep="\t")
    if not {"sample_id", "meandepth"} <= set(cov.columns):
        return []
    low = cov.loc[cov["meandepth"] < min_mean_depth, ["sample_id", "meandepth"]]
    for row in low.itertuples(index=False):
        logging.warning(
            "Excluding %s from PCA/IBS: mean depth %.2fx is below --min-mean-depth %.2fx.",
            row.sample_id, row.meandepth, min_mean_depth,
        )
    return [str(sid) for sid in low["sample_id"]]


def vcf_sample_names(vcf_path: Path) -> List[str]:
    """Return the sample names in a VCF header, in order."""
    result = run_command(["bcftools", "query", "-l", str(vcf_path)], capture_output=True)
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def map_vcf_sample_names(names: Iterable[str], known_ids: Iterable[str]) -> Dict[str, str]:
    """Map path-style VCF sample names (from BAMs without read groups) to sample IDs.

    Only names that are not already known IDs and whose basename, minus
    `.sorted.bam`/`.bam`, is a known ID are included.
    """
    known = set(known_ids)
    mapping: Dict[str, str] = {}
    for name in names:
        if name in known:
            continue
        base = Path(name).name
        for suffix in (".sorted.bam", ".bam"):
            if base.endswith(suffix):
                base = base[: -len(suffix)]
                break
        if base in known:
            mapping[name] = base
    return mapping


def rename_vcf_samples(vcf_path: Path, samples: List[Sample], threads: int) -> bool:
    """Rename path-style sample names in `vcf_path` to sample IDs, in place.

    Returns True if the file was rewritten (so downstream PLINK outputs, which
    carry the old names, must be regenerated).
    """
    names = vcf_sample_names(vcf_path)
    mapping = map_vcf_sample_names(names, (s.sample_id for s in samples))
    if not mapping:
        return False
    known = {s.sample_id for s in samples}
    unmapped = [n for n in names if n not in known and n not in mapping]
    if unmapped:
        logging.warning("%d VCF sample names could not be mapped to sample IDs: %s", len(unmapped), unmapped[:5])
    map_file = vcf_path.with_name(vcf_path.name + ".sample_rename.txt")
    map_file.write_text("".join(f"{old} {new}\n" for old, new in mapping.items()))
    renamed = vcf_path.with_name(vcf_path.name.replace(".vcf.gz", ".renamed.vcf.gz"))
    logging.info("Renaming %d path-style sample names in %s to sample IDs.", len(mapping), vcf_path.name)
    run_command(["bcftools", "reheader", "--threads", str(threads), "-s", str(map_file), "-o", str(renamed), str(vcf_path)])
    renamed.replace(vcf_path)
    for ext in (".csi", ".tbi"):
        idx = vcf_path.with_name(vcf_path.name + ext)
        if idx.exists():
            idx.unlink()
    run_command(["bcftools", "index", "--threads", str(threads), str(vcf_path)])
    return True


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
    # Tag reads with a read group so downstream tools (bcftools, PLINK) name the
    # sample by its ID rather than by the BAM file path. bwa expands the literal
    # "\t"; the single quotes keep bash from touching it.
    read_group = f"@RG\\tID:{sample.sample_id}\\tSM:{sample.sample_id}\\tPL:ILLUMINA"
    pipeline = (
        "set -euo pipefail\n"
        f"bwa mem -t {threads} -R '{read_group}' {reference_fasta} {sample.r1} {sample.r2} | "
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
    keep_ids: Optional[List[str]] = None,
) -> Tuple[Path, Path]:
    """Convert VCF to PLINK format and compute PCA + IBS matrices.

    `keep_ids`, when given, restricts PCA/IBS to those sample IDs via `--keep`
    (used to drop blanks and other undiscovered samples from an existing VCF).
    """
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

        keep_args: List[str] = []
        if keep_ids:
            keep_file = variant_dir / "plink_keep.txt"
            keep_file.write_text("".join(f"{sid} {sid}\n" for sid in keep_ids))
            keep_args = ["--keep", str(keep_file)]
            logging.info("Restricting PLINK PCA/IBS to %d samples via --keep.", len(keep_ids))

        pca_cmd = [
            plink_bin,
            "--bfile",
            str(prefix),
            *keep_args,
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
            *keep_args,
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

    def normalize_sample_id(raw_id: str) -> str:
        """Map a VCF/PLINK sample name back to the sample ID.

        BAMs aligned without a read group are named by bcftools after their file
        path (e.g. ``.../CS18_22_Wild_plate1_A7.sorted.bam``); strip that back to
        the bare sample ID so location lookups and axis labels still work.
        """
        name = Path(str(raw_id)).name
        for suffix in (".sorted.bam", ".bam"):
            if name.endswith(suffix):
                name = name[: -len(suffix)]
        return name

    eigenvec = pd.read_csv(
        eigenvec_path,
        sep=r"\s+",
        header=0,
    )
    eigenvec.rename(columns={"#FID": "fid", "FID": "fid", "IID": "iid"}, inplace=True)
    eigenvec["iid"] = eigenvec["iid"].map(normalize_sample_id)
    eigenvec["location"] = eigenvec["iid"].map(lambda sid: sample_lookup.get(sid, "unknown"))
    unknown = int((eigenvec["location"] == "unknown").sum())
    if unknown:
        logging.warning(
            "%d of %d PCA samples could not be matched to a location.", unknown, len(eigenvec)
        )
    if "PC1" not in eigenvec.columns or "PC2" not in eigenvec.columns:
        raise RuntimeError("PCA results missing PC1/PC2 columns.")

    ibs_ids = []
    with (ibs_matrix_path.with_suffix(".mibs.id")).open() as handle:
        for line in handle:
            parts = line.strip().split()
            if parts:
                ibs_ids.append(normalize_sample_id(parts[-1]))
    ibs_matrix = np.atleast_2d(np.loadtxt(ibs_matrix_path))
    if ibs_matrix.shape[0] != len(ibs_ids):
        raise RuntimeError("Mismatch between IBS matrix and ID file.")
    # PLINK `--distance square ibs` writes a .mibs file that already holds IBS
    # *similarity* (proportion of alleles shared, 1.0 on the diagonal), so it is
    # used as-is. Subtracting from 1 would turn it into a distance.
    ibs_similarity = ibs_matrix

    locations = sorted({sample.location for sample in samples})

    try:
        import matplotlib.pyplot as plt
    except ImportError as err:
        raise ImportError("matplotlib is required for visualization.") from err

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    cmap = plt.get_cmap("tab20", len(locations))
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
    fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04, label="IBS similarity (proportion of alleles shared)")
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

        sample_sheet_path = metrics_dir / "sample_metadata.tsv"
        try:
            samples = discover_samples(raw_dir, skip_blanks=args.skip_blanks)
            write_sample_sheet(samples, sample_sheet_path, repo_root)
        except FileNotFoundError:
            if not sample_sheet_path.exists():
                raise
            logging.warning(
                "No FASTQs found in %s; loading samples from %s instead. "
                "Alignment cannot be re-run without the raw reads.",
                raw_dir, sample_sheet_path,
            )
            samples = load_samples_from_sheet(sample_sheet_path, repo_root, args.skip_blanks)

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
            logging.info("Coverage summary exists; skipping alignment and coverage calculation.")
            for sample in samples:
                bam_path = align_dir / f"{sample.sample_id}.sorted.bam"
                if not bam_path.exists():
                    raise FileNotFoundError(
                        f"{bam_path} not found although coverage summary exists; "
                        "delete metrics/coverage_summary.tsv or re-run with --force to align."
                    )
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

        ensure_read_groups(
            samples,
            bam_paths,
            fix=args.fix_read_groups,
            workers=max(1, total_threads // threads_per_sample),
            threads_per_bam=threads_per_sample,
        )

        force_variants = args.force or args.force_variants
        raw_vcf, filtered_vcf = call_variants(
            bam_paths,
            reference_fasta,
            variant_dir,
            threads=total_threads,
            max_depth=args.max_depth,
            min_total_depth=args.min_total_depth,
            min_qual=args.min_qual,
            force=force_variants,
        )

        # VCFs built from BAMs without read groups carry file paths as sample
        # names; rewrite the header so PLINK and the figure use sample IDs.
        renamed = rename_vcf_samples(filtered_vcf, samples, threads=total_threads)
        vcf_names = set(vcf_sample_names(filtered_vcf))
        excluded = set(low_depth_samples(coverage_summary_path, args.min_mean_depth))
        keep_ids = [
            s.sample_id for s in samples if s.sample_id in vcf_names and s.sample_id not in excluded
        ]
        if not keep_ids:
            logging.warning("No VCF sample names match discovered sample IDs; figure labels may be wrong.")
            keep_ids = None
        elif len(keep_ids) == len(vcf_names):
            keep_ids = None  # nothing to drop

        eigenvec_path, ibs_matrix_path = run_plink_pca(
            filtered_vcf,
            variant_dir,
            threads=total_threads,
            force=force_variants or args.force_plink or renamed,
            plink_min_contig_length=args.plink_min_contig_length,
            plink_max_contigs=args.plink_max_contigs,
            use_plink2=args.plink2,
            keep_ids=keep_ids,
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
                "min_mean_depth": args.min_mean_depth,
                "excluded_low_depth_samples": sorted(excluded),
                "fix_read_groups": args.fix_read_groups,
                "force_variants": args.force_variants,
                "force_plink": args.force_plink,
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

