#!/usr/bin/env python3
"""Genotype side of the step 09 RDA: allele frequencies from genotype likelihoods.

Called from ``code/09_rda.Rmd`` (as SLURM jobs) with one of two subcommands.

``reduce``
    Streams the per-window Beagle files from step 06 in parallel. For every
    population set in ``output/09_rda/inputs/population-sets.tsv`` it estimates
    the minor-allele frequency of each unit (a location, or several pooled) by
    EM on the genotype likelihoods, keeps SNPs with no missing unit, a mean
    frequency inside ``[--min-maf, 1 - --min-maf]`` and non-zero variance, and
    standardises each SNP across units. Only the Gram matrix (units x units) is
    accumulated, so its principal coordinates reproduce an RDA on the full
    units x SNPs matrix exactly: eigenvalues, R2, F and permutation p depend on
    the response only through its row space, and total inertia equals the
    number of SNPs. The same is done for individuals (posterior-mean dosages
    under a Hardy-Weinberg prior at the all-sample frequency) for the
    individual-level tests. Frequencies for the sets named by ``--save-sets``
    are kept in ``work/`` for the loading scan.

``scan``
    Projects every retained SNP of the ``env`` set onto the constrained axes of
    the fitted environment-only RDA (``loading = standardised SNP . lc site
    score``) and runs the Capblancq and Forester (2021) outlier test: Mahalanobis
    distance under a robust (MCD) covariance, rescaled by its median against a
    chi-square reference, Benjamini-Hochberg q-values.

Needs only numpy and pandas (both in the step 06 ``angsd`` conda env).

Usage:
    python code/09_rda_genotypes.py reduce [--threads 16]
    python code/09_rda_genotypes.py scan [--threads 16]
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_NAME = "09_rda_genotypes.py"
GL_DIR = Path("output/06_angsd_structure/gl")
SAMPLES = Path("output/06_angsd_structure/samples/samples.tsv")
OUTPUT_DIR = Path("output/09_rda")
INPUTS_DIR = OUTPUT_DIR / "inputs"
WORK_DIR = OUTPUT_DIR / "work"
TABLES_DIR = OUTPUT_DIR / "tables"
LOGS_DIR = OUTPUT_DIR / "logs"
POPULATION_SETS = INPUTS_DIR / "population-sets.tsv"
WINDOW_GLOB = "[0-9][0-9][0-9]_*.beagle.gz"   # per-window files; excludes all.beagle.gz

log = logging.getLogger(SCRIPT_NAME)


def setup_logging(name: str) -> None:
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    log.setLevel(logging.INFO)
    for handler in (logging.FileHandler(LOGS_DIR / f"{name}.log", mode="w"),
                    logging.StreamHandler(sys.stdout)):
        handler.setFormatter(formatter)
        log.addHandler(handler)


# --------------------------------------------------------------------------- likelihood maths

def read_beagle(path: Path, n_ind: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Chromosome, position and (n_snp, n_ind, 3) likelihoods from an ANGSD Beagle file."""
    frame = pd.read_csv(path, sep="\t", dtype={"marker": str})
    if frame.shape[1] != 3 + 3 * n_ind:
        raise ValueError(f"{path}: {frame.shape[1] - 3} likelihood columns, expected {3 * n_ind}")
    marker = frame["marker"].str.rsplit("_", n=1, expand=True)
    gl = frame.iloc[:, 3:].to_numpy(np.float32).reshape(len(frame), n_ind, 3)
    return marker[0].to_numpy(), marker[1].astype(np.int64).to_numpy(), gl


def informative(gl: np.ndarray) -> np.ndarray:
    """False where an individual has no reads (all three likelihoods equal)."""
    return ~((np.abs(gl[..., 0] - gl[..., 1]) < 1e-6) & (np.abs(gl[..., 1] - gl[..., 2]) < 1e-6))


def em_frequency(gl: np.ndarray, iterations: int = 50, tolerance: float = 1e-5) -> np.ndarray:
    """Minor-allele frequency per SNP from the likelihoods of a group of individuals.

    EM under Hardy-Weinberg (Kim et al. 2011, as in ANGSD ``-doMaf``), using only
    individuals with data at the SNP. NaN where no individual has data.
    """
    weight = informative(gl).astype(np.float64)
    n_info = weight.sum(axis=1)
    gl = gl.astype(np.float64)
    f = np.full(gl.shape[0], 0.25)
    for _ in range(iterations):
        prior = np.stack([(1 - f) ** 2, 2 * f * (1 - f), f ** 2], axis=1)[:, None, :]
        post = gl * prior
        total = post.sum(axis=2)
        total[total == 0] = 1
        expected = (post[..., 1] + 2 * post[..., 2]) / total
        updated = (expected * weight).sum(axis=1) / np.maximum(2 * n_info, 1)
        done = np.nanmax(np.abs(updated - f)) < tolerance
        f = updated
        if done:
            break
    f[n_info == 0] = np.nan
    return f


def posterior_dosage(gl: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Posterior-mean minor-allele dosage (0-2) per individual, HWE prior at ``f``."""
    prior = np.stack([(1 - f) ** 2, 2 * f * (1 - f), f ** 2], axis=1)[:, None, :]
    post = gl.astype(np.float64) * prior
    total = post.sum(axis=2)
    total[total == 0] = 1
    return ((post[..., 1] + 2 * post[..., 2]) / total).T          # (n_ind, n_snp)


def standardise(matrix: np.ndarray, min_maf: float | None) -> tuple[np.ndarray, np.ndarray]:
    """Columns (SNPs) passing the filters, centred and scaled across rows (ddof 1).

    A column is kept when it has no missing value, non-zero variance and, unless
    ``min_maf`` is None, a mean inside ``[min_maf, 1 - min_maf]``.
    """
    sd = matrix.std(axis=0, ddof=1)
    keep = ~np.isnan(matrix).any(axis=0) & (sd > 1e-8)
    if min_maf is not None:
        mean = matrix.mean(axis=0)
        keep &= (mean >= min_maf) & (mean <= 1 - min_maf)
    sub = matrix[:, keep]
    return (sub - sub.mean(axis=0)) / sd[keep], keep


def principal_coordinates(gram: np.ndarray, n_snp: int) -> np.ndarray:
    """Principal coordinates of a Gram matrix; total variance is asserted to equal ``n_snp``."""
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    keep = values > 1e-8 * max(values.max(), 1)
    coords = vectors[:, keep] * np.sqrt(values[keep])
    inertia = float(coords.var(axis=0, ddof=1).sum())
    if abs(inertia - n_snp) > max(1e-3, 1e-6 * n_snp):
        raise AssertionError(f"inertia {inertia} != SNP count {n_snp}")
    return coords


# --------------------------------------------------------------------------- population sets

def load_sets(samples: pd.DataFrame) -> dict:
    """Population sets -> ordered units -> indices of their individuals in the Beagle files."""
    if not POPULATION_SETS.exists():
        raise FileNotFoundError(f"{POPULATION_SETS} not found; run the population-sets chunk of code/09_rda.Rmd")
    table = pd.read_csv(POPULATION_SETS, sep="\t")
    known = set(samples["location"])
    unknown = sorted(set(table["location"]) - known)
    if unknown:
        raise ValueError(f"population-sets.tsv names locations with no samples: {', '.join(unknown)}")
    sets = {}
    for name, rows in table.groupby("set", sort=False):
        units = {}
        for unit, members in rows.groupby("unit", sort=False):
            idx = np.where(samples["location"].isin(members["location"]))[0]
            units[unit] = idx
        sets[name] = units
    return sets


# --------------------------------------------------------------------------- reduce

def reduce_window(task: dict) -> dict:
    """Allele frequencies, filters and Gram increments for one Beagle window."""
    chrom, pos, gl = read_beagle(task["path"], task["n_ind"])
    min_maf = task["min_maf"]
    result = {"window": task["path"].name, "n_input": int(len(pos)), "sets": {}}

    f_all = em_frequency(gl)
    for name, units in task["sets"].items():
        freq = np.vstack([em_frequency(gl[:, idx, :]) for idx in units.values()])
        z, keep = standardise(freq, min_maf)
        result["sets"][name] = {"gram": z @ z.T, "n_snp": int(keep.sum())}
        if name in task["save_sets"]:
            out = WORK_DIR / name / task["path"].name.replace(".beagle.gz", ".npz")
            np.savez(out, freq=freq[:, keep].astype(np.float32), chrom=chrom[keep], pos=pos[keep])

    # individuals: SNPs filtered on the all-sample frequency, not on mean dosage
    common = (f_all >= min_maf) & (f_all <= 1 - min_maf)
    dosage = posterior_dosage(gl[common][:, task["individuals"], :], f_all[common])
    z, keep = standardise(dosage, None)
    result["individual"] = {"gram": z @ z.T, "n_snp": int(keep.sum())}
    return result


def run_reduce(args: argparse.Namespace) -> int:
    setup_logging("reduce")
    started = time.time()
    if not SAMPLES.exists():
        raise FileNotFoundError(f"{SAMPLES} not found; run code/06_angsd_structure.Rmd first")
    samples = pd.read_csv(SAMPLES, sep="\t")
    windows = sorted(GL_DIR.glob(WINDOW_GLOB))
    if not windows:
        raise FileNotFoundError(f"no per-window Beagle files ({WINDOW_GLOB}) in {GL_DIR}; run step 06 first")
    sets = load_sets(samples)
    save_sets = [s for s in args.save_sets.split(",") if s]
    for name in save_sets:
        if name not in sets:
            raise ValueError(f"--save-sets names unknown set {name}")
        (WORK_DIR / name).mkdir(parents=True, exist_ok=True)
    ind_set = sets[args.individual_set]
    individuals = np.sort(np.concatenate(list(ind_set.values())))
    log.info("%d windows, %d samples, sets: %s; individual set %s (%d individuals)", len(windows),
             len(samples), ", ".join(f"{k} ({len(v)} units)" for k, v in sets.items()),
             args.individual_set, len(individuals))

    tasks = [{"path": p, "n_ind": len(samples), "min_maf": args.min_maf, "sets": sets,
              "save_sets": save_sets, "individuals": individuals} for p in windows]
    grams = {name: np.zeros((len(units), len(units))) for name, units in sets.items()}
    counts = {name: 0 for name in sets}
    ind_gram = np.zeros((len(individuals), len(individuals)))
    ind_count, n_input, per_window = 0, 0, []
    with ProcessPoolExecutor(max_workers=args.threads) as pool:
        for result in pool.map(reduce_window, tasks):
            n_input += result["n_input"]
            row = {"window": result["window"], "n_input": result["n_input"]}
            for name, part in result["sets"].items():
                grams[name] += part["gram"]
                counts[name] += part["n_snp"]
                row[f"n_{name}"] = part["n_snp"]
            ind_gram += result["individual"]["gram"]
            ind_count += result["individual"]["n_snp"]
            row["n_individual"] = result["individual"]["n_snp"]
            per_window.append(row)
            log.info("  %s: %d SNPs in, %s", result["window"], result["n_input"],
                     ", ".join(f"{k} {v['n_snp']}" for k, v in result["sets"].items()))

    INPUTS_DIR.mkdir(parents=True, exist_ok=True)
    summary = {}
    for name, units in sets.items():
        coords = principal_coordinates(grams[name], counts[name])
        frame = pd.DataFrame(coords, index=list(units), columns=[f"PCo{i + 1}" for i in range(coords.shape[1])])
        frame.index.name = "unit"
        frame.to_csv(INPUTS_DIR / f"{name}-pcoa.tsv", sep="\t")
        summary[name] = {"n_units": len(units), "n_snp": counts[name], "rank": int(coords.shape[1])}
        log.info("set %s: %d units, %d SNPs, rank %d", name, len(units), counts[name], coords.shape[1])
    coords = principal_coordinates(ind_gram, ind_count)
    frame = pd.DataFrame(coords, columns=[f"PCo{i + 1}" for i in range(coords.shape[1])])
    frame.insert(0, "location", samples["location"].to_numpy()[individuals])
    frame.insert(0, "sample_id", samples["sample_id"].to_numpy()[individuals])
    frame.to_csv(INPUTS_DIR / "individual-pcoa.tsv", sep="\t", index=False)
    summary["individual"] = {"set": args.individual_set, "n_individuals": int(len(individuals)),
                             "n_snp": ind_count, "rank": int(coords.shape[1])}
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(per_window).to_csv(TABLES_DIR / "snps-per-window.tsv", sep="\t", index=False)

    with open(INPUTS_DIR / "reduce-metadata.json", "w") as handle:
        json.dump({"script": f"{SCRIPT_NAME} reduce", "date": datetime.now(timezone.utc).isoformat(),
                   "runtime_seconds": round(time.time() - started, 1),
                   "parameters": {"min_maf": args.min_maf, "threads": args.threads,
                                  "save_sets": save_sets, "individual_set": args.individual_set},
                   "inputs": [SAMPLES.as_posix(), POPULATION_SETS.as_posix(), f"{GL_DIR.as_posix()}/{WINDOW_GLOB}"],
                   "n_windows": len(windows), "n_snp_input": n_input, "sets": summary,
                   "software": {"python": platform.python_version(), "numpy": np.__version__,
                                "pandas": pd.__version__}}, handle, indent=2)
    log.info("done in %.0f s", time.time() - started)
    return 0


# --------------------------------------------------------------------------- scan

def chi2_sf(x: np.ndarray, df: int) -> np.ndarray:
    """Upper tail of the chi-square distribution for integer ``df`` (closed form)."""
    half = np.asarray(x, dtype=np.float64) / 2
    if df % 2 == 0:
        term = np.ones_like(half)
        total = np.ones_like(half)
        for i in range(1, df // 2):
            term = term * half / i
            total += term
        return np.exp(-half) * total
    erfc = np.frompyfunc(math.erfc, 1, 1)
    total = erfc(np.sqrt(half)).astype(np.float64)
    term = np.sqrt(half) / math.gamma(1.5)
    for i in range(1, (df - 1) // 2 + 1):
        if i > 1:
            term = term * half / (i + 0.5 - 1)
        total += np.exp(-half) * term
    return total


def chi2_median(df: int) -> float:
    low, high = 0.0, 10.0 * df + 10
    for _ in range(200):
        mid = (low + high) / 2
        if chi2_sf(np.array([mid]), df)[0] > 0.5:
            low = mid
        else:
            high = mid
    return (low + high) / 2


def mahalanobis(x: np.ndarray, centre: np.ndarray, cov: np.ndarray) -> np.ndarray:
    diff = x - centre
    return np.einsum("ij,jk,ik->i", diff, np.linalg.inv(cov), diff)


def minimum_covariance_determinant(x: np.ndarray, rng: np.random.Generator,
                                   starts: int = 50, steps: int = 30) -> tuple[np.ndarray, np.ndarray]:
    """FastMCD C-steps (Rousseeuw and Van Driessen 1999) from random (p+1)-subsets."""
    n, p = x.shape
    h = (n + p + 1) // 2
    best = (np.inf, None, None)
    for _ in range(starts):
        subset = rng.choice(n, p + 1, replace=False)
        centre, cov = x[subset].mean(axis=0), np.cov(x[subset], rowvar=False).reshape(p, p)
        if np.linalg.matrix_rank(cov) < p:
            cov = cov + np.eye(p) * 1e-6
        for _ in range(steps):
            subset = np.argsort(mahalanobis(x, centre, cov))[:h]
            new_centre, new_cov = x[subset].mean(axis=0), np.cov(x[subset], rowvar=False).reshape(p, p)
            converged = np.allclose(new_cov, cov) and np.allclose(new_centre, centre)
            centre, cov = new_centre, new_cov
            if converged:
                break
        det = np.linalg.slogdet(cov)[1]
        if det < best[0]:
            best = (det, centre, cov)
    return best[1], best[2]


def bh_qvalues(p: np.ndarray) -> np.ndarray:
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.empty_like(p)
    q[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    return np.clip(q, 0, 1)


def scan_window(task: dict) -> tuple[pd.DataFrame, np.ndarray]:
    data = np.load(task["path"], allow_pickle=True)
    z, keep = standardise(data["freq"].astype(np.float64), task["min_maf"])
    loadings = z.T @ task["scores"]
    return pd.DataFrame({"chrom": data["chrom"][keep], "pos": data["pos"][keep]}), loadings


def run_scan(args: argparse.Namespace) -> int:
    setup_logging("scan")
    started = time.time()
    rng = np.random.default_rng(args.seed)
    scores_path = TABLES_DIR / "axis-unit-scores.tsv"
    if not scores_path.exists():
        raise FileNotFoundError(f"{scores_path} not found; run the rda-models chunk of code/09_rda.Rmd")
    scores = pd.read_csv(scores_path, sep="\t", index_col=0)
    units = pd.read_csv(INPUTS_DIR / f"{args.set}-pcoa.tsv", sep="\t", index_col=0).index
    if list(scores.index) != list(units):
        raise ValueError(f"unit order in {scores_path} does not match {args.set}-pcoa.tsv")
    files = sorted((WORK_DIR / args.set).glob("*.npz"))
    if not files:
        raise FileNotFoundError(f"no frequency files in {WORK_DIR / args.set}; rerun reduce with --save-sets {args.set}")
    k = scores.shape[1]
    tasks = [{"path": f, "scores": scores.to_numpy(np.float64), "min_maf": args.min_maf} for f in files]
    positions, parts = [], []
    with ProcessPoolExecutor(max_workers=args.threads) as pool:
        for pos, loadings in pool.map(scan_window, tasks):
            positions.append(pos)
            parts.append(loadings)
    loci = pd.concat(positions, ignore_index=True)
    loadings = np.vstack(parts)
    log.info("%d SNPs projected onto %d axes", len(loci), k)

    sample = loadings[rng.choice(len(loadings), size=min(args.mcd_subsample, len(loadings)), replace=False)]
    centre, cov = minimum_covariance_determinant(sample, rng)
    distance = mahalanobis(loadings, centre, cov)
    inflation = float(np.median(distance) / chi2_median(k))
    p = chi2_sf(distance / inflation, k)
    q = bh_qvalues(p)
    for j in range(k):
        loci[scores.columns[j]] = loadings[:, j].round(5)
    loci["mahalanobis"] = distance.round(4)
    loci["p_value"] = p
    loci["q_value"] = q

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    loci.to_csv(WORK_DIR / "loadings-all.tsv.gz", sep="\t", index=False)
    top = loci.nsmallest(args.top, "p_value")
    top.to_csv(TABLES_DIR / "outlier-top-loci.tsv", sep="\t", index=False)
    thin = loci[(loci["p_value"] < 1e-3)
                | (rng.random(len(loci)) < args.manhattan_sample / max(len(loci), 1))]
    thin[["chrom", "pos", "p_value", "q_value"]].to_csv(TABLES_DIR / "manhattan-thinned.tsv.gz", sep="\t", index=False)
    summary = {"set": args.set, "n_axes": k, "n_snp": int(len(loci)),
               "n_q_lt_0.05": int((q < 0.05).sum()), "n_q_lt_0.10": int((q < 0.10).sum()),
               "n_q_lt_0.20": int((q < 0.20).sum()), "min_q": float(q.min()),
               "median_p": float(np.median(p)), "inflation": round(inflation, 4),
               "mcd_subsample": int(len(sample))}
    pd.DataFrame([summary]).to_csv(TABLES_DIR / "outlier-summary.tsv", sep="\t", index=False)
    with open(OUTPUT_DIR / "logs" / "scan-metadata.json", "w") as handle:
        json.dump({"script": f"{SCRIPT_NAME} scan", "date": datetime.now(timezone.utc).isoformat(),
                   "runtime_seconds": round(time.time() - started, 1),
                   "parameters": {k2: v for k2, v in vars(args).items() if k2 != "func"},
                   "results": summary}, handle, indent=2)
    log.info("outliers: %s", summary)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    reduce = sub.add_parser("reduce", help="allele frequencies and Gram-matrix principal coordinates")
    reduce.add_argument("--threads", type=int, default=16, help="parallel windows (default: 16, max 50)")
    reduce.add_argument("--min-maf", type=float, default=0.05,
                        help="minimum mean minor-allele frequency across units (default: 0.05)")
    reduce.add_argument("--save-sets", default="env",
                        help="comma-separated sets whose frequencies are kept for the scan (default: env)")
    reduce.add_argument("--individual-set", default="env",
                        help="population set whose individuals enter the individual-level Gram (default: env)")
    reduce.set_defaults(func=run_reduce)
    scan = sub.add_parser("scan", help="loading projection and robust Mahalanobis outlier test")
    scan.add_argument("--threads", type=int, default=16)
    scan.add_argument("--set", default="env", help="population set the RDA was fitted on (default: env)")
    scan.add_argument("--min-maf", type=float, default=0.05, help="must match reduce (default: 0.05)")
    scan.add_argument("--mcd-subsample", type=int, default=50000,
                      help="SNPs used to estimate the robust covariance (default: 50000)")
    scan.add_argument("--top", type=int, default=1000, help="loci written to outlier-top-loci.tsv (default: 1000)")
    scan.add_argument("--manhattan-sample", type=int, default=200000,
                      help="approximate SNPs kept for the Manhattan plot besides p < 1e-3 (default: 200000)")
    scan.add_argument("--seed", type=int, default=42)
    scan.set_defaults(func=run_scan)
    args = parser.parse_args()
    args.threads = min(args.threads, 50)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
