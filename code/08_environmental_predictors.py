#!/usr/bin/env python3
"""Site-by-predictor matrix for the step 09 RDA.

Takes the per-site Ecology climatologies from step 07 and turns them into the
explanatory side of a redundancy analysis:

1. fills temperature for sites with no Ecology station (Coos Bay) from a NOAA
   CO-OPS shore station over the same years, summarised with step 07's own
   climatology functions so the values are directly comparable,
2. flags which sites carry every selected predictor (``in_env_model``) and
   which share an Ecology station (identical predictor values),
3. builds distance-based Moran's eigenvector maps (dbMEM) from the site
   coordinates, once for the sites in the environmental model and once for all
   sites, to use as geographic conditioning terms,
4. screens the candidate predictors for collinearity: within-family and
   between-predictor correlations, backward VIF elimination, a comparison of
   named candidate sets, and the eigenvalue spectrum of the environment.

The selected set defaults to summer temperature, salinity seasonal range and
summer chlorophyll, the set chosen in the earlier Olurida_v081 analysis; the
VIF table shows whether it still holds on the current data.

Usage:
    python code/08_environmental_predictors.py [--predictors a,b,c] [--force]
"""

from __future__ import annotations

import argparse
import importlib.util
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

SCRIPT_NAME = "08_environmental_predictors.py"
SITE_COORDINATES = Path("output/04_environmental_data/site-coordinates.tsv")
STEP07_DIR = Path("output/07_orca_ecology_data")
STEP07_SCRIPT = Path("code/07_orca_ecology_data.py")
OUTPUT_DIR = Path("output/08_environmental_predictors")
RAW_DIR = OUTPUT_DIR / "raw"
TABLES_DIR = OUTPUT_DIR / "tables"
FIGURES_DIR = OUTPUT_DIR / "figures"
LOGS_DIR = OUTPUT_DIR / "logs"

COOPS_URL = ("https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=water_temperature"
             "&station={station}&begin_date={year}0101&end_date={year}1231&interval=h"
             "&units=metric&time_zone=gmt&format=csv&application=oly-lc-WGS")

# Shore stations used for sites that step 07 could not match to an Ecology
# station. Charleston OR is 0.4 km from the Coos Bay coordinate.
DEFAULT_COOPS = "Coos_Bay=9432780"

DEFAULT_PREDICTORS = "temp_summer_mean,sal_range,chl_summer_mean"
# One or two statistics per variable, chosen for distinct biological meaning:
# mean condition in the growing season, seasonality, and the oxygen minimum.
DEFAULT_CANDIDATES = ("temp_summer_mean,temp_range,sal_annual_mean,sal_range,"
                      "do_month_min,chl_summer_mean")
FAMILY_STATS = ["annual_mean", "summer_mean", "winter_mean", "month_max", "month_min", "range"]
FAMILIES = ["temp", "sal", "do", "chl"]

# Named sets compared in tables/candidate-sets.tsv in addition to the VIF survivors.
CANDIDATE_SETS = {
    "selected": None,  # filled from --predictors
    "thermal_hypoxia_chl": ["temp_summer_mean", "do_month_min", "chl_summer_mean"],
    "salinity_anchored": ["sal_annual_mean", "sal_range", "chl_summer_mean"],
    "thermal_salinity_chl_do": ["temp_summer_mean", "sal_range", "chl_summer_mean", "do_month_min"],
}

log = logging.getLogger(SCRIPT_NAME)


def setup_logging() -> None:
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    log.setLevel(logging.INFO)
    for handler in (logging.FileHandler(LOGS_DIR / "pipeline.log", mode="w"),
                    logging.StreamHandler(sys.stdout)):
        handler.setFormatter(formatter)
        log.addHandler(handler)


def load_step07():
    """Import step 07 as a module so the climatology definitions are shared, not copied."""
    spec = importlib.util.spec_from_file_location("step07", STEP07_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def haversine_matrix(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Pairwise great-circle distances in kilometres."""
    radius = 6371.0088
    p = np.radians(lat)
    dp = p[:, None] - p[None, :]
    dl = np.radians(lon[:, None] - lon[None, :])
    a = np.sin(dp / 2) ** 2 + np.cos(p[:, None]) * np.cos(p[None, :]) * np.sin(dl / 2) ** 2
    return 2 * radius * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


# --------------------------------------------------------------------------- CO-OPS temperature

def coops_daily(station: str, years: list[int], force: bool, step07) -> pd.DataFrame:
    """Daily mean water temperature from the CO-OPS API, one cached CSV per year."""
    frames = []
    for year in years:
        path = RAW_DIR / f"coops_{station}_{year}.csv"
        if force or not path.exists():
            step07.download(COOPS_URL.format(station=station, year=year), path)
        frame = pd.read_csv(path, skipinitialspace=True)
        frame.columns = [c.strip() for c in frame.columns]
        if "Water Temperature" not in frame.columns:
            log.warning("  CO-OPS %s %s: no temperature column (%s)", station, year,
                        path.read_text(errors="replace")[:120].strip())
            continue
        frame = frame.rename(columns={"Date Time": "time", "Water Temperature": "temperature"})
        frame["time"] = pd.to_datetime(frame["time"])
        frame["temperature"] = pd.to_numeric(frame["temperature"], errors="coerce")
        frames.append(frame[["time", "temperature"]])
        log.info("  CO-OPS %s %s: %d hourly values", station, year, frame["temperature"].notna().sum())
    if not frames:
        return pd.DataFrame(columns=["date", "temperature"])
    hourly = pd.concat(frames, ignore_index=True).dropna()
    hourly["date"] = hourly["time"].dt.normalize()
    return hourly.groupby("date", as_index=False)["temperature"].mean()


# --------------------------------------------------------------------------- dbMEM

def minimum_spanning_tree_max_edge(dist: np.ndarray) -> float:
    """Longest edge of the minimum spanning tree (Prim's algorithm)."""
    n = len(dist)
    in_tree = np.zeros(n, dtype=bool)
    in_tree[0] = True
    best = dist[0].copy()
    longest = 0.0
    for _ in range(n - 1):
        candidates = np.where(in_tree, np.inf, best)
        j = int(np.argmin(candidates))
        longest = max(longest, float(candidates[j]))
        in_tree[j] = True
        best = np.minimum(best, dist[j])
    return longest


def dbmem(lat: np.ndarray, lon: np.ndarray, labels: list[str], prefix: str) -> tuple[pd.DataFrame, list[dict]]:
    """Distance-based Moran's eigenvector maps (Dray, Legendre and Peres-Neto 2006).

    Distances above the longest minimum-spanning-tree edge are set to four times
    that threshold, the truncated matrix is Gower-centred, and the positive-
    eigenvalue axes with positive spatial autocorrelation (Moran's I above its
    expectation) are kept, broad scales first.
    """
    n = len(lat)
    dist = haversine_matrix(lat, lon)
    threshold = minimum_spanning_tree_max_edge(dist)
    truncated = np.where(dist <= threshold, dist, 4 * threshold)
    np.fill_diagonal(truncated, 0)
    centring = np.eye(n) - np.ones((n, n)) / n
    gower = -0.5 * centring @ (truncated ** 2) @ centring
    values, vectors = np.linalg.eigh(gower)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    positive = values > 1e-8 * np.abs(values).max()

    weights = ((dist <= threshold) & ~np.eye(n, dtype=bool)).astype(float)
    expected = -1 / (n - 1)
    columns, info = {}, []
    for j in np.where(positive)[0]:
        z = vectors[:, j] - vectors[:, j].mean()
        moran = (n / weights.sum()) * (z @ weights @ z) / (z @ z)
        keep = moran > expected
        name = f"{prefix}{len(columns) + 1}" if keep else ""
        if keep:
            columns[name] = vectors[:, j] * np.sqrt(values[j])
        info.append({"set": prefix, "axis": int(j + 1), "eigenvalue": round(float(values[j]), 4),
                     "morans_I": round(float(moran), 4), "expected_I": round(expected, 4),
                     "kept_as": name, "mst_threshold_km": round(threshold, 2), "n_sites": n})
    return pd.DataFrame(columns, index=labels), info


# --------------------------------------------------------------------------- collinearity

def vif(frame: pd.DataFrame) -> pd.Series:
    """Variance inflation factors: the diagonal of the inverse correlation matrix."""
    corr = np.corrcoef(frame.to_numpy(dtype=float), rowvar=False)
    return pd.Series(np.diag(np.linalg.pinv(corr)), index=frame.columns)


def backward_vif(frame: pd.DataFrame, threshold: float) -> tuple[list[str], pd.DataFrame]:
    """Drop the highest-VIF predictor until every VIF is at or below ``threshold``."""
    keep, history = list(frame.columns), []
    while len(keep) > 1:
        values = vif(frame[keep]).sort_values(ascending=False)
        worst, worst_value = values.index[0], float(values.iloc[0])
        stop = worst_value <= threshold
        history.append({"step": len(history) + 1, "n_predictors": len(keep),
                        "max_vif_predictor": worst, "max_vif": round(worst_value, 3),
                        "action": "stop" if stop else "drop", "retained": ";".join(keep)})
        if stop:
            break
        keep.remove(worst)
    return keep, pd.DataFrame(history)


def plot_correlation(corr: pd.DataFrame, selected: list[str], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 7))
    image = ax.imshow(corr.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(corr)))
    ax.set_xticklabels(corr.columns, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(corr)))
    ax.set_yticklabels(corr.index, fontsize=8)
    for i in range(len(corr)):
        for j in range(len(corr)):
            ax.text(j, i, f"{corr.iat[i, j]:.2f}", ha="center", va="center", fontsize=6)
    for k, name in enumerate(corr.columns):
        if name in selected:
            ax.get_xticklabels()[k].set_fontweight("bold")
            ax.get_yticklabels()[k].set_fontweight("bold")
    fig.colorbar(image, ax=ax, shrink=0.8, label="Pearson r")
    ax.set_title("Candidate predictors and geography (sites in the environmental model)\n"
                 "selected predictors in bold")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------- main

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--predictors", default=DEFAULT_PREDICTORS,
                        help=f"comma-separated predictors for the RDA (default: {DEFAULT_PREDICTORS})")
    parser.add_argument("--candidates", default=DEFAULT_CANDIDATES,
                        help="comma-separated candidate predictors for the collinearity screen")
    parser.add_argument("--coops", default=DEFAULT_COOPS,
                        help="location=CO-OPS station pairs for sites without an Ecology station "
                             f"(default: {DEFAULT_COOPS}; empty string to skip)")
    parser.add_argument("--vif-threshold", type=float, default=5.0,
                        help="VIF above which backward elimination drops a predictor (default: 5)")
    parser.add_argument("--force", action="store_true", help="re-download cached CO-OPS files")
    args = parser.parse_args()

    predictor_table = STEP07_DIR / "tables/site_predictors_ecology.tsv"
    for path in (SITE_COORDINATES, predictor_table, STEP07_DIR / "metadata.json", STEP07_SCRIPT):
        if not path.exists():
            print(f"ERROR: {path} not found; run code/07_orca_ecology_data.py first.", file=sys.stderr)
            return 1

    started = time.time()
    for directory in (RAW_DIR, TABLES_DIR, FIGURES_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    setup_logging()
    step07 = load_step07()
    step07_meta = json.loads((STEP07_DIR / "metadata.json").read_text())
    years = list(range(step07_meta["parameters"]["start_year"], step07_meta["parameters"]["end_year"] + 1))
    min_months = step07_meta["parameters"]["min_months"]
    selected = [p.strip() for p in args.predictors.split(",") if p.strip()]
    candidates = [p.strip() for p in args.candidates.split(",") if p.strip()]
    candidates += [p for p in selected if p not in candidates]

    sites = pd.read_csv(SITE_COORDINATES, sep="\t")
    env = pd.read_csv(predictor_table, sep="\t")
    env["env_source"] = np.where(env["ecology_station"].fillna("") != "", "WA Ecology", "")
    env["env_station"] = env["ecology_station"].fillna("")
    env["env_distance_km"] = env["ecology_distance_km"]
    log.info("%d sites from step 07; window %d-%d", len(env), years[0], years[-1])

    # CO-OPS temperature for sites without an Ecology station ---------------
    coops_sources = []
    for pair in filter(None, args.coops.split(",")):
        location, station = pair.split("=")
        row = env.index[env["location"] == location]
        if row.empty:
            log.warning("--coops: %s is not a site; skipped", location)
            continue
        if env.loc[row[0], "env_station"]:
            log.info("%s already has Ecology station %s; CO-OPS %s not used", location,
                     env.loc[row[0], "env_station"], station)
            continue
        daily = coops_daily(station, years, args.force, step07)
        if daily.empty:
            log.warning("%s: CO-OPS %s returned no temperature for %d-%d", location, station, years[0], years[-1])
            continue
        clim = step07.monthly_climatology(daily.assign(date=pd.to_datetime(daily["date"])))
        summary = {k: v for k, v in step07.summarize_climatology(clim, min_months).items() if k.startswith("temp_")}
        for key, value in summary.items():
            env.loc[row[0], key] = value
        env.loc[row[0], ["env_source", "env_station"]] = ["NOAA CO-OPS (temperature only)", f"CO-OPS {station}"]
        coops_sources.append({"location": location, "station": station, "n_days": int(len(daily)),
                              "url_template": COOPS_URL})
        log.info("%s: CO-OPS %s, %d daily means, summer mean %.2f C", location, station, len(daily),
                 summary.get("temp_summer_mean", np.nan))

    # Model membership and shared stations ------------------------------------
    missing = [p for p in candidates if p not in env.columns]
    if missing:
        print(f"ERROR: predictors not in step 07 output: {', '.join(missing)}", file=sys.stderr)
        return 1
    env["in_env_model"] = env[selected].notna().all(axis=1)
    groups = env[env["env_station"] != ""].groupby("env_station")["location"].apply(list).to_dict()
    env["shares_station_with"] = [";".join(x for x in groups.get(s, []) if x != loc)
                                  for loc, s in zip(env["location"], env["env_station"])]
    model = env[env["in_env_model"]].reset_index(drop=True)
    n_profiles = model.groupby("env_station").ngroups
    log.info("environmental model: %d sites, %d distinct stations; excluded: %s", len(model), n_profiles,
             ", ".join(env.loc[~env["in_env_model"], "location"]) or "none")

    # Geography -----------------------------------------------------------------
    mem_model, mem_info = dbmem(model["latitude"].to_numpy(), model["longitude"].to_numpy(),
                                model["location"].tolist(), "MEM")
    mem_all, info_all = dbmem(env["latitude"].to_numpy(), env["longitude"].to_numpy(),
                              env["location"].tolist(), "MEMall")
    env = (env.merge(mem_model, left_on="location", right_index=True, how="left")
              .merge(mem_all, left_on="location", right_index=True, how="left"))
    log.info("dbMEM: %d axes kept for the model sites, %d for all sites", mem_model.shape[1], mem_all.shape[1])

    # Collinearity --------------------------------------------------------------
    X = model[candidates].astype(float)
    if X.isna().any().any():
        dropped = X.columns[X.isna().any()].tolist()
        log.warning("candidates missing at some model sites, left out of the screen: %s", ", ".join(dropped))
        X = X.dropna(axis=1)
    if len(X) <= X.shape[1]:
        log.warning("%d sites for %d candidates: the VIF screen is rank-deficient", len(X), X.shape[1])

    family_rows = []
    for family in FAMILIES:
        cols = [f"{family}_{s}" for s in FAMILY_STATS if f"{family}_{s}" in model.columns]
        block = model[cols].astype(float).dropna(axis=1)
        corr = block.corr()
        for i, a in enumerate(corr.columns):
            for b in corr.columns[i + 1:]:
                family_rows.append({"family": family, "predictor_a": a, "predictor_b": b,
                                    "pearson_r": round(float(corr.loc[a, b]), 4)})
    family_corr = pd.DataFrame(family_rows)

    geography = model[["latitude", "longitude"]].join(mem_model.reset_index(drop=True))
    corr = pd.concat([X, geography], axis=1).corr()
    pairs = [{"predictor_a": a, "predictor_b": b, "pearson_r": round(float(corr.loc[a, b]), 4)}
             for i, a in enumerate(corr.columns) for b in corr.columns[i + 1:]]
    pairs = pd.DataFrame(pairs).assign(abs_r=lambda d: d["pearson_r"].abs()).sort_values("abs_r", ascending=False)

    survivors, history = backward_vif(X, args.vif_threshold)
    selected_in_x = [p for p in selected if p in X.columns]
    vif_rows = ([{"set": "all_candidates", "predictor": k, "vif": round(float(v), 3)} for k, v in vif(X).items()]
                + [{"set": "selected", "predictor": k, "vif": round(float(v), 3)}
                   for k, v in vif(X[selected_in_x]).items()])
    sets = {**CANDIDATE_SETS, "selected": selected_in_x, "vif_survivors": survivors}
    set_rows = []
    for name, cols in sets.items():
        cols = [c for c in cols if c in X.columns]
        if len(cols) < 2:
            continue
        values = vif(X[cols])
        set_rows.append({"candidate_set": name, "n_predictors": len(cols), "predictors": ";".join(cols),
                         "max_vif": round(float(values.max()), 3), "mean_vif": round(float(values.mean()), 3)})
    set_table = pd.DataFrame(set_rows).sort_values("max_vif")

    Z = (X - X.mean()) / X.std(ddof=1)
    eig = np.linalg.svd(Z.to_numpy(), compute_uv=False) ** 2 / (len(Z) - 1)
    pca = pd.DataFrame({"axis": [f"PC{i + 1}" for i in range(len(eig))], "eigenvalue": eig.round(4),
                        "proportion": (eig / eig.sum()).round(4),
                        "cumulative": np.cumsum(eig / eig.sum()).round(4)})
    max_selected = float(vif(X[selected_in_x]).max()) if len(selected_in_x) > 1 else 1.0
    log.info("selected %s: max VIF %.2f; VIF survivors: %s", ", ".join(selected_in_x), max_selected,
             ", ".join(survivors))
    if max_selected > args.vif_threshold:
        log.warning("selected set exceeds the VIF threshold of %g; see tables/candidate-sets.tsv",
                    args.vif_threshold)
    log.info("environmental PCA: PC1 %.1f%%, PC2 %.1f%%", 100 * pca["proportion"].iloc[0],
             100 * pca["proportion"].iloc[1])

    # Outputs ------------------------------------------------------------------
    front = ["location", "region", "latitude", "longitude", "n_samples", "in_env_model", "env_source",
             "env_station", "env_distance_km", "shares_station_with"]
    env = env[front + [c for c in env.columns if c not in front and c not in
                       ("ecology_station", "ecology_distance_km")]]
    paths = {
        "site_env_matrix": TABLES_DIR / "site-env-matrix.tsv",
        "predictor_correlation": TABLES_DIR / "predictor-correlation.tsv",
        "within_family_correlation": TABLES_DIR / "within-family-correlation.tsv",
        "vif": TABLES_DIR / "vif.tsv",
        "vif_elimination": TABLES_DIR / "vif-elimination.tsv",
        "candidate_sets": TABLES_DIR / "candidate-sets.tsv",
        "environment_pca": TABLES_DIR / "environment-pca.tsv",
        "dbmem": TABLES_DIR / "dbmem.tsv",
        "figure": FIGURES_DIR / "predictor-correlation.png",
    }
    env.to_csv(paths["site_env_matrix"], sep="\t", index=False, float_format="%.6g")
    pairs.drop(columns="abs_r").to_csv(paths["predictor_correlation"], sep="\t", index=False)
    family_corr.to_csv(paths["within_family_correlation"], sep="\t", index=False)
    pd.DataFrame(vif_rows).to_csv(paths["vif"], sep="\t", index=False)
    history.to_csv(paths["vif_elimination"], sep="\t", index=False)
    set_table.to_csv(paths["candidate_sets"], sep="\t", index=False)
    pca.to_csv(paths["environment_pca"], sep="\t", index=False)
    pd.DataFrame(mem_info + info_all).to_csv(paths["dbmem"], sep="\t", index=False)
    plot_correlation(corr, selected_in_x, paths["figure"])

    with open(OUTPUT_DIR / "metadata.json", "w") as handle:
        json.dump({
            "script": SCRIPT_NAME,
            "inputs": [SITE_COORDINATES.as_posix(), predictor_table.as_posix(),
                       (STEP07_DIR / "metadata.json").as_posix()],
            "parameters": {"predictors": selected, "candidates": candidates, "coops": args.coops,
                           "vif_threshold": args.vif_threshold, "years": years,
                           "min_months": min_months, "force": args.force},
            "sources": {"climatologies": "step 07 (WA Ecology monthly CTD casts, 0-5 m)",
                        "coops": coops_sources},
            "results": {
                "n_sites": int(len(env)),
                "n_sites_in_env_model": int(env["in_env_model"].sum()),
                "excluded_from_env_model": env.loc[~env["in_env_model"], "location"].tolist(),
                "distinct_stations_in_env_model": int(n_profiles),
                "sites_sharing_a_station": {k: v for k, v in groups.items() if len(v) > 1},
                "selected_max_vif": round(max_selected, 3),
                "vif_survivors": survivors,
                "mem_axes_env_model": mem_model.columns.tolist(),
                "mem_axes_all_sites": mem_all.columns.tolist(),
                "environment_pc1_pc2": pca["proportion"].iloc[:2].tolist(),
            },
            "produced_files": [p.as_posix() for p in paths.values()],
            "date": datetime.now(timezone.utc).isoformat(),
            "runtime_seconds": round(time.time() - started, 1),
            "software": {"python": platform.python_version(), "pandas": pd.__version__,
                         "numpy": np.__version__, "matplotlib": matplotlib.__version__},
        }, handle, indent=2, default=str)
    log.info("Wrote %s/", OUTPUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
