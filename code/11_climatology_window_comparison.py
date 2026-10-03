#!/usr/bin/env python3
"""Do the site predictors depend on the climatology window?

Steps 07 and 08 describe every site with one 2015-2018 climatology, but the
non-`18` collections were sampled in 2024. This script compares the step 08
site-env matrix built for the reference window with one built for another
window (by default 2021-2024, the years before the 2024 collections) and asks
whether the sites keep their values and, more importantly for the RDA, their
order. For every summary predictor it reports Pearson r, Spearman rho, the mean
shift and the largest rank change across the step 09 environmental-model sites,
and it lists sites whose Ecology station changed between windows.

The alternative window is produced by the existing steps with their output
redirected, e.g.:

    python code/07_orca_ecology_data.py --start-year 2021 --end-year 2024 \\
        --output-dir output/07_orca_ecology_data_2021-2024
    python code/08_environmental_predictors.py \\
        --step07-dir output/07_orca_ecology_data_2021-2024 \\
        --output-dir output/08_environmental_predictors_2021-2024
    python code/11_climatology_window_comparison.py

Usage:
    python code/11_climatology_window_comparison.py [--reference DIR] [--alternative DIR]
"""

from __future__ import annotations

import argparse
import json
import logging
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SCRIPT_NAME = "11_climatology_window_comparison.py"
OUTPUT_DIR = Path("output/11_climatology_window_comparison")
TABLES_DIR = OUTPUT_DIR / "tables"
FIGURES_DIR = OUTPUT_DIR / "figures"
LOGS_DIR = OUTPUT_DIR / "logs"
COLLECTION_YEAR = Path("output/09_rda/inputs/collection-year.tsv")

FAMILIES = ["temp", "sal", "do", "chl"]
STATS = ["annual_mean", "summer_mean", "winter_mean", "month_max", "month_min", "range"]

log = logging.getLogger(SCRIPT_NAME)


def setup_logging() -> None:
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    log.setLevel(logging.INFO)
    for handler in (logging.FileHandler(LOGS_DIR / "pipeline.log", mode="w"),
                    logging.StreamHandler(sys.stdout)):
        handler.setFormatter(formatter)
        log.addHandler(handler)


def window_label(step08_dir: Path) -> str:
    years = json.loads((step08_dir / "metadata.json").read_text())["parameters"]["years"]
    return f"{min(years)}-{max(years)}"


def agreement(ref: pd.Series, alt: pd.Series) -> dict:
    """Agreement between two site vectors over the sites present in both."""
    ok = ref.notna() & alt.notna()
    r, a = ref[ok], alt[ok]
    if len(r) < 3 or r.std() == 0 or a.std() == 0:
        return {"n_sites": int(len(r))}
    rank_r, rank_a = r.rank(), a.rank()
    return {"n_sites": int(len(r)),
            "pearson_r": float(np.corrcoef(r, a)[0, 1]),
            "spearman_rho": float(np.corrcoef(rank_r, rank_a)[0, 1]),
            "mean_reference": float(r.mean()), "mean_alternative": float(a.mean()),
            "mean_shift": float((a - r).mean()),
            "sd_reference": float(r.std()), "sd_alternative": float(a.std()),
            "max_abs_rank_change": int((rank_a - rank_r).abs().max())}


def plot_selected(values: pd.DataFrame, selected: list[str], labels: tuple[str, str], path: Path) -> None:
    fig, axes = plt.subplots(1, len(selected), figsize=(5 * len(selected), 5), squeeze=False)
    colours = {2018: "tab:blue", 2024: "tab:orange"}
    for ax, name in zip(axes[0], selected):
        d = values[(values["predictor"] == name) & values["in_env_model"]].dropna(subset=["reference", "alternative"])
        for year, g in d.groupby("collection_year"):
            ax.scatter(g["reference"], g["alternative"], color=colours.get(year, "grey"), label=f"collected {year}")
        for _, row in d.iterrows():
            ax.annotate(row["location"], (row["reference"], row["alternative"]), fontsize=6,
                        xytext=(3, 3), textcoords="offset points")
        lo = min(d["reference"].min(), d["alternative"].min())
        hi = max(d["reference"].max(), d["alternative"].max())
        ax.plot([lo, hi], [lo, hi], color="grey", lw=0.8, ls="--")
        rho = np.corrcoef(d["reference"].rank(), d["alternative"].rank())[0, 1]
        ax.set_title(f"{name}  (Spearman rho = {rho:.2f}, n = {len(d)})", fontsize=9)
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
    axes[0][0].legend(fontsize=7, frameon=False)
    fig.suptitle("Site predictors under two climatology windows (environmental-model sites)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", default="output/08_environmental_predictors",
                        help="step 08 output for the window used in the analysis")
    parser.add_argument("--alternative", default="output/08_environmental_predictors_2021-2024",
                        help="step 08 output for the window to compare against")
    args = parser.parse_args()
    ref_dir, alt_dir = Path(args.reference), Path(args.alternative)
    for d in (ref_dir, alt_dir):
        for f in ("tables/site-env-matrix.tsv", "metadata.json"):
            if not (d / f).exists():
                print(f"ERROR: {d / f} not found; run code/08_environmental_predictors.py for that window first.",
                      file=sys.stderr)
                return 1

    started = time.time()
    for directory in (TABLES_DIR, FIGURES_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    setup_logging()
    ref_label, alt_label = window_label(ref_dir), window_label(alt_dir)
    selected = json.loads((ref_dir / "metadata.json").read_text())["parameters"]["predictors"]
    ref = pd.read_csv(ref_dir / "tables/site-env-matrix.tsv", sep="\t").set_index("location")
    alt = pd.read_csv(alt_dir / "tables/site-env-matrix.tsv", sep="\t").set_index("location")
    if set(ref.index) != set(alt.index):
        log.error("site sets differ: %s", sorted(set(ref.index) ^ set(alt.index)))
        return 1
    alt = alt.loc[ref.index]
    env_sites = ref.index[ref["in_env_model"].astype(bool)]
    log.info("reference %s, alternative %s; %d sites, %d in the environmental model; selected %s",
             ref_label, alt_label, len(ref), len(env_sites), ", ".join(selected))

    if COLLECTION_YEAR.exists():
        year = pd.read_csv(COLLECTION_YEAR, sep="\t").set_index("location")["year"]
    else:   # same rule as step 09: labels containing "18" are the 2018 collections
        year = pd.Series(np.where(ref.index.str.contains("18"), 2018, 2024), index=ref.index)

    # station assignment per window
    stations = pd.DataFrame({
        "collection_year": year.reindex(ref.index).values,
        "in_env_model_reference": ref["in_env_model"].astype(bool).values,
        "in_env_model_alternative": alt["in_env_model"].astype(bool).values,
        "station_reference": ref["env_station"].fillna("").values,
        "station_alternative": alt["env_station"].fillna("").values,
        "n_profiles_reference": ref["ecology_n_profiles"].values,
        "n_profiles_alternative": alt["ecology_n_profiles"].values}, index=ref.index)
    stations["station_changed"] = stations["station_reference"] != stations["station_alternative"]
    stations.rename_axis("location").reset_index().to_csv(TABLES_DIR / "station-assignment.tsv",
                                                         sep="\t", index=False)
    for loc, row in stations[stations["station_changed"]].iterrows():
        log.info("%s: station %s (%s) -> %s (%s)", loc, row["station_reference"] or "none", ref_label,
                 row["station_alternative"] or "none", alt_label)

    # per-site values, long format
    predictors = [f"{f}_{s}" for f in FAMILIES for s in STATS if f"{f}_{s}" in ref and f"{f}_{s}" in alt]
    values = pd.concat([pd.DataFrame({
        "location": ref.index, "collection_year": year.reindex(ref.index).values,
        "in_env_model": ref.index.isin(env_sites), "predictor": p,
        "reference": ref[p].values, "alternative": alt[p].values}) for p in predictors], ignore_index=True)
    values["difference"] = values["alternative"] - values["reference"]
    values.to_csv(TABLES_DIR / "site-values.tsv", sep="\t", index=False, float_format="%.6g")

    # agreement per predictor: environmental-model sites, and those sites by collection year
    rows = []
    for p in predictors:
        for subset, sites in [("env model sites", env_sites),
                              ("env model sites, collected 2018", env_sites[year.reindex(env_sites).values == 2018]),
                              ("env model sites, collected 2024", env_sites[year.reindex(env_sites).values == 2024])]:
            rows.append({"predictor": p, "selected": p in selected, "sites": subset,
                         **agreement(ref.loc[sites, p], alt.loc[sites, p])})
    agree = pd.DataFrame(rows)
    agree.to_csv(TABLES_DIR / "predictor-agreement.tsv", sep="\t", index=False, float_format="%.4g")
    main_rows = agree[agree["sites"] == "env model sites"]
    for _, row in main_rows[main_rows["selected"]].iterrows():
        log.info("%s: pearson %.3f, spearman %.3f, mean shift %+.3f, max rank change %s",
                 row["predictor"], row["pearson_r"], row["spearman_rho"], row["mean_shift"],
                 row["max_abs_rank_change"])

    figure = FIGURES_DIR / "selected-predictors.png"
    plot_selected(values, selected, (ref_label, alt_label), figure)

    with open(OUTPUT_DIR / "metadata.json", "w") as handle:
        json.dump({
            "script": SCRIPT_NAME,
            "inputs": [(ref_dir / "tables/site-env-matrix.tsv").as_posix(),
                       (alt_dir / "tables/site-env-matrix.tsv").as_posix(),
                       COLLECTION_YEAR.as_posix() if COLLECTION_YEAR.exists() else "label rule ('18' = 2018)"],
            "parameters": {"reference_window": ref_label, "alternative_window": alt_label,
                           "selected_predictors": selected, "agreement_sites": "reference in_env_model"},
            "results": {
                "selected": {r["predictor"]: {k: r[k] for k in ("pearson_r", "spearman_rho", "mean_shift",
                                                                 "max_abs_rank_change")}
                             for _, r in main_rows[main_rows["selected"]].iterrows()},
                "stations_changed": stations.index[stations["station_changed"]].tolist(),
                "min_spearman_all_predictors": float(main_rows["spearman_rho"].min()),
                "median_spearman_all_predictors": float(main_rows["spearman_rho"].median()),
            },
            "produced_files": [(TABLES_DIR / f).as_posix() for f in
                               ("station-assignment.tsv", "site-values.tsv", "predictor-agreement.tsv")]
                              + [figure.as_posix()],
            "date": datetime.now(timezone.utc).isoformat(),
            "runtime_seconds": round(time.time() - started, 1),
            "software": {"python": platform.python_version(), "pandas": pd.__version__,
                         "numpy": np.__version__, "matplotlib": matplotlib.__version__},
        }, handle, indent=2, default=str)
    log.info("Wrote %s/", OUTPUT_DIR)
    return 0


if __name__ == "__main__":
    sys.exit(main())
