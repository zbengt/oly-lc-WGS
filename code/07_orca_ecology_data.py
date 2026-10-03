#!/usr/bin/env python3
"""Site-level marine climatologies from WA Ecology profiles and UW/NANOOS ORCA moorings.

Replaces the 30-day NOAA snapshot of step 04 as the environmental layer for
genotype-environment analyses (RDA). For every putative site in
``output/04_environmental_data/site-coordinates.tsv`` the script:

1. downloads the WA Dept. of Ecology yearly marine-water profile netCDFs for
   the requested years, keeps QC-passed temperature, salinity, dissolved
   oxygen and chlorophyll fluorescence in a shallow depth band, and assigns
   each site the nearest Ecology station with enough profiles in the window,
2. pulls the ORCA L3 gridded profiles from the NANOOS ERDDAP for moorings
   within range of a site (QARTOD PASS values only, same depth band),
3. builds monthly climatologies per site and source, and
4. derives one row of summary predictors per site (annual, summer, winter
   means; monthly extremes; seasonal range) for use in a later RDA step.

Ecology is the primary source because it covers every Washington site with a
single method (monthly CTD casts); ORCA is written separately as a
higher-frequency check for the sites it reaches. Coos Bay (Oregon) has no
source here and is left missing.

Usage:
    python code/07_orca_ecology_data.py [--start-year 2015] [--end-year 2018]
"""

from __future__ import annotations

import argparse
import io
import json
import logging
import math
import platform
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SCRIPT_NAME = "07_orca_ecology_data.py"
SITE_COORDINATES = Path("output/04_environmental_data/site-coordinates.tsv")
OUTPUT_DIR = Path("output/07_orca_ecology_data")
RAW_DIR = OUTPUT_DIR / "raw"
TABLES_DIR = OUTPUT_DIR / "tables"
FIGURES_DIR = OUTPUT_DIR / "figures"
LOGS_DIR = OUTPUT_DIR / "logs"

ECOLOGY_URL = (
    "https://fortress.wa.gov/ecy/ezshare/EAP/SalishSea/netCDF-files/"
    "MarineWaterProfilesAndNutrientsYear{year}.nc"
)
# data.nanoos.org redirects here; the griddap request must percent-encode [ ].
ERDDAP_BASE = "https://erddap.nanoos.org/erddap"

USER_AGENT = "oly-lc-WGS environmental data assembly (github.com/zbengt/oly-lc-WGS)"
FILL_VALUE = -99999.9

# Output variable name -> (Ecology netCDF variable, ORCA ERDDAP variable, units).
VARIABLES = {
    "temperature": ("Temp", "sea_water_temperature", "degC"),
    "salinity": ("Salinity", "sea_water_practical_salinity", "PSU"),
    "oxygen": ("DOAdjusted", "mass_concentration_of_oxygen_in_sea_water", "mg/L"),
    "chlorophyll": ("FluorAdjusted", "mass_concentration_of_chlorophyll_a_in_sea_water", "mg/m3"),
}
# Column prefix for each variable in the per-site summary tables.
PREFIX = {"temperature": "temp", "salinity": "sal", "oxygen": "do", "chlorophyll": "chl"}
ECOLOGY_QC_PASS = 2   # QC code: 0 None, 1 Fail, 2 Pass
ORCA_QC_PASS = 1      # QARTOD aggregate: 1 PASS, 2 NOT_EVALUATED, 3 SUSPECT, 4 FAIL, 9 MISSING

# ORCA moorings: ERDDAP dataset, description, buoy lat/lon from the dataset
# global attributes (same list as step 04).
ORCA_MOORINGS = [
    ("orca_hydro_twanoh", "Twanoh, southern Hood Canal", 47.375, -123.008),
    ("orca_hydro_hoodsport", "Hoodsport, Hood Canal", 47.425, -123.113),
    ("orca_hydro_dabobbay", "Dabob Bay, Hood Canal", 47.803, -122.803),
    ("orca_hydro_hansville", "Hansville, northern Puget Sound", 47.907, -122.627),
    ("orca_hydro_pointwells", "Point Wells, central Puget Sound", 47.761, -122.397),
    ("orca_hydro_carrinlet", "Carr Inlet, southern Puget Sound", 47.280, -122.730),
]

# Sites whose nearest station by straight-line distance lies across land in a
# different water body. Each maps to the station on the same passage system.
STATION_OVERRIDES = {
    # Liberty Bay drains to Port Orchard; HCB010 (14.6 km) is across the Kitsap
    # Peninsula in Hood Canal. OCH014 in Liberty Bay was not sampled 2015-2018.
    "Dogfish_Bay": "SIN001",
}

SUMMER_MONTHS = (7, 8, 9)
WINTER_MONTHS = (12, 1, 2)

log = logging.getLogger(SCRIPT_NAME)


def setup_logging() -> None:
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    log.setLevel(logging.INFO)
    for handler in (logging.FileHandler(LOGS_DIR / "pipeline.log", mode="w"),
                    logging.StreamHandler(sys.stdout)):
        handler.setFormatter(formatter)
        log.addHandler(handler)


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Great-circle distance in kilometres."""
    radius = 6371.0088
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = p2 - p1
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * radius * math.asin(math.sqrt(a))


def download(url: str, dest: Path, timeout: int = 300, retries: int = 3) -> None:
    """Stream a URL to ``dest`` via a ``.part`` file, retrying transient failures."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=timeout) as response, open(part, "wb") as out:
                while True:
                    chunk = response.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            part.rename(dest)
            return
        except urllib.error.HTTPError as error:
            last_error = error
            if error.code == 404:
                break
            time.sleep(5 * (attempt + 1))
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            last_error = error
            time.sleep(5 * (attempt + 1))
    if part.exists():
        part.unlink()
    raise RuntimeError(f"failed to fetch {url}: {last_error}")


# --------------------------------------------------------------------------- Ecology

def decode_chars(array: np.ndarray) -> list[str]:
    """Decode a netCDF (n, strlen) char array into stripped strings."""
    return [b"".join(row).decode("ascii", errors="replace").strip("\x00 ") for row in array]


def read_ecology_year(path: Path, max_depth: float) -> pd.DataFrame:
    """One row per QC-passed observation at or above ``max_depth`` metres."""
    with h5py.File(path, "r") as nc:
        stations = decode_chars(nc["Station"][:])
        station_pos = {int(n): i for i, n in enumerate(nc["station_number"][:])}
        lat, lon = nc["Latitude"][:], nc["Longitude"][:]

        # Profiles (ragged array): station and field date per profile.
        profile_pos = {int(n): i for i, n in enumerate(nc["profile_index"][:])}
        profile_station = nc["station_index"][:]
        field_date = nc["FieldDate"][:]

        depth = nc["Depth"][:]
        keep = depth <= max_depth
        obs_profile = nc["obs_index"][:][keep]
        rows = np.array([profile_pos[int(p)] for p in obs_profile], dtype=int)
        sidx = np.array([station_pos[int(profile_station[r])] for r in rows], dtype=int)

        frame = pd.DataFrame({
            "station": np.array(stations, dtype=object)[sidx],
            "station_lat": lat[sidx].astype(float),
            "station_lon": lon[sidx].astype(float),
            "profile": [f"{path.stem[-4:]}-{int(p)}" for p in obs_profile],
            "date": pd.to_datetime(field_date[rows], unit="D"),
            "depth_m": depth[keep].astype(float),
        })
        for name, (eco_var, _orca_var, _units) in VARIABLES.items():
            values = nc[eco_var][:][keep].astype(float)
            qc = nc[eco_var + "QC"][:][keep]
            ok = (qc == ECOLOGY_QC_PASS) & (values > FILL_VALUE / 10)
            frame[name] = np.where(ok, values, np.nan)
    return frame


def ecology_profiles(years: list[int], max_depth: float, force: bool) -> tuple[pd.DataFrame, list[str]]:
    """Depth-band means per Ecology profile across ``years``."""
    frames, sources = [], []
    for year in years:
        url = ECOLOGY_URL.format(year=year)
        path = RAW_DIR / "ecology" / Path(urllib.parse.urlparse(url).path).name
        if force or not path.exists():
            log.info("Downloading Ecology %s ...", year)
            download(url, path)
        frames.append(read_ecology_year(path, max_depth))
        sources.append(url)
        log.info("  Ecology %s: %d observations within %g m", year, len(frames[-1]), max_depth)
    obs = pd.concat(frames, ignore_index=True)
    names = list(VARIABLES)
    profiles = (obs.groupby(["station", "station_lat", "station_lon", "profile", "date"])[names]
                .mean().reset_index().sort_values(["station", "date"]))
    return profiles, sources


def assign_ecology_stations(sites: pd.DataFrame, profiles: pd.DataFrame,
                            max_km: float, min_profiles: int) -> pd.DataFrame:
    """Nearest station within ``max_km`` with at least ``min_profiles`` profiles, unless overridden."""
    stations = (profiles.groupby(["station", "station_lat", "station_lon"])
                .agg(n_profiles=("profile", "nunique"),
                     first_date=("date", "min"), last_date=("date", "max"))
                .reset_index())
    rows = []
    for _, site in sites.iterrows():
        dist = stations.apply(lambda s: haversine_km(site["latitude"], site["longitude"],
                                                     s["station_lat"], s["station_lon"]), axis=1)
        ranked = stations.assign(distance_km=dist.round(2)).sort_values("distance_km")
        nearest = ranked.iloc[0]
        usable = ranked[(ranked["distance_km"] <= max_km) & (ranked["n_profiles"] >= min_profiles)]
        chosen = usable.iloc[0] if not usable.empty else None
        method = "nearest"
        override = STATION_OVERRIDES.get(site["location"])
        if override is not None and override in set(ranked["station"]):
            chosen = ranked[ranked["station"] == override].iloc[0]
            method = "override"
        rows.append({
            "location": site["location"],
            "putative_site": site["putative_site"],
            "nearest_station": nearest["station"],
            "nearest_distance_km": nearest["distance_km"],
            "nearest_n_profiles": int(nearest["n_profiles"]),
            "assigned_station": chosen["station"] if chosen is not None else "",
            "assignment_method": method if chosen is not None else "none",
            "assigned_distance_km": chosen["distance_km"] if chosen is not None else np.nan,
            "assigned_n_profiles": int(chosen["n_profiles"]) if chosen is not None else 0,
            "assigned_first_date": chosen["first_date"].date().isoformat() if chosen is not None else "",
            "assigned_last_date": chosen["last_date"].date().isoformat() if chosen is not None else "",
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- ORCA

def griddap_url(dataset: str, start: str, end: str, max_depth: float) -> str:
    """CSV request for all variables plus QC flags over a time and depth window."""
    constraint = f"[({start}):1:({end})][(0.0):1:({max_depth})]"
    names = []
    for _name, (_eco, orca_var, _units) in VARIABLES.items():
        names += [orca_var + constraint, orca_var + "_qc_aggregate" + constraint]
    query = ",".join(names)
    query = query.replace("[", "%5B").replace("]", "%5D").replace("<", "%3C").replace(">", "%3E")
    return f"{ERDDAP_BASE}/griddap/{dataset}.csv?{query}"


def read_orca_csv(path: Path) -> pd.DataFrame:
    """ERDDAP CSV (row 2 holds units) -> QC-passed values per cast and depth."""
    frame = pd.read_csv(path, skiprows=[1])
    out = pd.DataFrame({
        "time": pd.to_datetime(frame["cast_start_time"], utc=True),
        "depth_m": frame["depth"].astype(float),
    })
    for name, (_eco, orca_var, _units) in VARIABLES.items():
        values = pd.to_numeric(frame[orca_var], errors="coerce")
        qc = pd.to_numeric(frame[orca_var + "_qc_aggregate"], errors="coerce")
        out[name] = values.where(qc == ORCA_QC_PASS)
    return out


def orca_daily(dataset: str, years: list[int], max_depth: float, force: bool) -> tuple[pd.DataFrame, list[dict]]:
    """Daily means of cast depth-band means for one mooring; failures are reported, not raised."""
    frames, status = [], []
    for year in years:
        path = RAW_DIR / "orca" / f"{dataset}_{year}_0-{max_depth:g}m.csv"
        url = griddap_url(dataset, f"{year}-01-01T00:00:00Z", f"{year}-12-31T23:59:59Z", max_depth)
        if force or not path.exists():
            try:
                download(url, path, timeout=600)
            except RuntimeError as error:
                message = str(error)
                # ERDDAP answers 404 when the time window is outside the record.
                state = "no data in window" if "404" in message else "failed"
                log.warning("  ORCA %s %s: %s", dataset, year, message.split(": ", 1)[-1])
                status.append({"dataset": dataset, "year": year, "status": state, "url": url})
                if state == "failed":
                    # Server-side failure (e.g. 503): later years will fail the same way.
                    break
                continue
        frame = read_orca_csv(path)
        frames.append(frame)
        status.append({"dataset": dataset, "year": year, "status": "ok", "rows": int(len(frame)),
                       "url": url})
    if not frames:
        return pd.DataFrame(), status
    obs = pd.concat(frames, ignore_index=True)
    names = list(VARIABLES)
    casts = obs.groupby("time")[names].mean().reset_index()
    casts["date"] = casts["time"].dt.tz_convert(None).dt.normalize()
    daily = casts.groupby("date")[names].mean().reset_index()
    return daily, status


# --------------------------------------------------------------------------- summaries

def monthly_climatology(frame: pd.DataFrame) -> pd.DataFrame:
    """Mean per calendar month of the per-(year, month) means, long format."""
    names = [n for n in VARIABLES if n in frame]
    work = frame.assign(year=frame["date"].dt.year, month=frame["date"].dt.month)
    by_year_month = work.groupby(["year", "month"])[names].mean().reset_index()
    rows = []
    for name in names:
        grouped = by_year_month.dropna(subset=[name]).groupby("month")[name]
        n_obs = work.dropna(subset=[name]).groupby("month")[name].size()
        for month, values in grouped:
            rows.append({"variable": name, "units": VARIABLES[name][2], "month": int(month),
                         "mean": round(float(values.mean()), 3),
                         "n_years": int(values.size), "n_obs": int(n_obs.get(month, 0))})
    return pd.DataFrame(rows)


def summarize_climatology(clim: pd.DataFrame, min_months: int) -> dict:
    """Annual/summer/winter means, monthly extremes, and seasonal range per variable."""
    out = {}
    for name in VARIABLES:
        months = clim[clim["variable"] == name].set_index("month")["mean"]
        prefix = PREFIX[name]
        complete = months.size >= min_months
        summer = months.reindex(SUMMER_MONTHS).dropna()
        winter = months.reindex(WINTER_MONTHS).dropna()
        out[f"{prefix}_n_months"] = int(months.size)
        out[f"{prefix}_annual_mean"] = round(float(months.mean()), 3) if complete else np.nan
        out[f"{prefix}_summer_mean"] = round(float(summer.mean()), 3) if summer.size >= 2 else np.nan
        out[f"{prefix}_winter_mean"] = round(float(winter.mean()), 3) if winter.size >= 2 else np.nan
        out[f"{prefix}_month_max"] = round(float(months.max()), 3) if complete else np.nan
        out[f"{prefix}_month_min"] = round(float(months.min()), 3) if complete else np.nan
        out[f"{prefix}_range"] = round(float(months.max() - months.min()), 3) if complete else np.nan
    return out


def daily_extremes(daily: pd.DataFrame) -> dict:
    """Tail statistics only a high-frequency record can give."""
    out = {}
    if daily.empty:
        return out
    temp = daily["temperature"].dropna()
    oxy = daily["oxygen"].dropna()
    sal = daily["salinity"].dropna()
    out["temp_daily_p95"] = round(float(temp.quantile(0.95)), 3) if temp.size else np.nan
    out["temp_days_above_18C_per_yr"] = (
        round(float((temp > 18).sum() / temp.size * 365.25), 1) if temp.size else np.nan)
    out["do_daily_p05"] = round(float(oxy.quantile(0.05)), 3) if oxy.size else np.nan
    out["do_days_below_4mgL_per_yr"] = (
        round(float((oxy < 4).sum() / oxy.size * 365.25), 1) if oxy.size else np.nan)
    out["sal_daily_p05"] = round(float(sal.quantile(0.05)), 3) if sal.size else np.nan
    out["n_days"] = int(len(daily))
    return out


def plot_climatologies(clim: pd.DataFrame, path: Path) -> None:
    names = list(VARIABLES)
    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex=True)
    eco = clim[clim["source"] == "Ecology"]
    # Several sites share a station; draw each station once, labelled with its sites.
    served = eco.groupby("station")["location"].unique()
    cmap = plt.get_cmap("tab20")
    for ax, name in zip(axes.flat, names):
        sub = eco[eco["variable"] == name]
        for i, (station, locations) in enumerate(served.items()):
            line = (sub[(sub["station"] == station) & (sub["location"] == locations[0])]
                    .sort_values("month"))
            if line.empty:
                continue
            ax.plot(line["month"], line["mean"], marker="o", ms=3, lw=1.2, color=cmap(i % 20),
                    label=f"{station}: {', '.join(sorted(locations))}")
        ax.set_title(f"{name} ({VARIABLES[name][2]})")
        ax.set_xticks(range(1, 13))
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("month")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", fontsize=8, frameon=False)
    fig.suptitle("Ecology monthly climatology at the station assigned to each site")
    fig.tight_layout(rect=(0, 0, 0.74, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------- main

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--start-year", type=int, default=2015,
                        help="first year of the climatology window (default: 2015)")
    parser.add_argument("--end-year", type=int, default=2018,
                        help="last year of the climatology window (default: 2018)")
    parser.add_argument("--max-depth", type=float, default=5.0,
                        help="deepest bin averaged, metres (default: 5)")
    parser.add_argument("--ecology-max-km", type=float, default=25.0,
                        help="maximum site-to-Ecology-station distance (default: 25)")
    parser.add_argument("--min-profiles", type=int, default=12,
                        help="minimum Ecology profiles in the window for a station to be assigned (default: 12)")
    parser.add_argument("--min-months", type=int, default=10,
                        help="calendar months required for annual summaries (default: 10)")
    parser.add_argument("--orca-max-km", type=float, default=30.0,
                        help="maximum site-to-ORCA-mooring distance (default: 30)")
    parser.add_argument("--skip-orca", action="store_true", help="do not query the NANOOS ERDDAP")
    parser.add_argument("--force", action="store_true", help="re-download cached raw files")
    args = parser.parse_args()

    if not SITE_COORDINATES.exists():
        print(f"ERROR: {SITE_COORDINATES} not found; run code/04_environmental_data.py first.",
              file=sys.stderr)
        return 1
    if args.end_year < args.start_year:
        print("ERROR: --end-year is before --start-year", file=sys.stderr)
        return 1

    started = time.time()
    for directory in (RAW_DIR, TABLES_DIR, FIGURES_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    setup_logging()
    years = list(range(args.start_year, args.end_year + 1))
    sites = pd.read_csv(SITE_COORDINATES, sep="\t")
    log.info("%d sites, window %d-%d, depth 0-%g m", len(sites), years[0], years[-1], args.max_depth)

    # Ecology ------------------------------------------------------------
    profiles, ecology_sources = ecology_profiles(years, args.max_depth, args.force)
    profiles_path = TABLES_DIR / "ecology_profiles.tsv"
    profiles.assign(date=profiles["date"].dt.date).round(3).to_csv(profiles_path, sep="\t", index=False)

    assignment = assign_ecology_stations(sites, profiles, args.ecology_max_km, args.min_profiles)
    assignment_path = TABLES_DIR / "site_station_assignment.tsv"

    clim_frames, predictor_rows = [], []
    for _, row in assignment.iterrows():
        base = {"location": row["location"]}
        if not row["assigned_station"]:
            log.info("%s: no Ecology station within %g km with >= %d profiles (nearest %s at %.1f km)",
                     row["location"], args.ecology_max_km, args.min_profiles,
                     row["nearest_station"], row["nearest_distance_km"])
            predictor_rows.append({**base, "ecology_station": "", "ecology_distance_km": np.nan})
            continue
        station_profiles = profiles[profiles["station"] == row["assigned_station"]]
        clim = monthly_climatology(station_profiles)
        clim.insert(0, "location", row["location"])
        clim.insert(1, "source", "Ecology")
        clim.insert(2, "station", row["assigned_station"])
        clim_frames.append(clim)
        predictor_rows.append({**base, "ecology_station": row["assigned_station"],
                               "ecology_distance_km": row["assigned_distance_km"],
                               "ecology_n_profiles": row["assigned_n_profiles"],
                               **summarize_climatology(clim, args.min_months)})
        log.info("%s: Ecology %s at %.1f km, %d profiles (%s)", row["location"],
                 row["assigned_station"], row["assigned_distance_km"], row["assigned_n_profiles"],
                 row["assignment_method"])

    # ORCA ---------------------------------------------------------------
    orca_status, orca_rows, orca_pairs = [], [], []
    server_down = False
    if not args.skip_orca:
        for dataset, description, mlat, mlon in ORCA_MOORINGS:
            near = []
            for _, site in sites.iterrows():
                distance = haversine_km(site["latitude"], site["longitude"], mlat, mlon)
                if distance <= args.orca_max_km:
                    near.append((site["location"], round(distance, 2)))
            if not near:
                continue
            log.info("ORCA %s (%s): %d sites within %g km", dataset, description, len(near),
                     args.orca_max_km)
            if server_down:
                log.warning("  skipped: ERDDAP returned server errors for the previous mooring")
                orca_status.append({"dataset": dataset, "status": "skipped (server unavailable)"})
                daily = pd.DataFrame()
            else:
                daily, status = orca_daily(dataset, years, args.max_depth, args.force)
                orca_status += status
                server_down = all(s["status"] == "failed" for s in status)
            for location, distance in near:
                orca_pairs.append({"location": location, "dataset": dataset,
                                   "mooring": description, "distance_km": distance,
                                   "n_days": int(len(daily))})
            if daily.empty:
                continue
            clim = monthly_climatology(daily)
            extremes = daily_extremes(daily)
            for location, distance in near:
                clim_frames.append(clim.assign(location=location, source="ORCA", station=dataset)
                                   [["location", "source", "station"] + list(clim.columns)])
                orca_rows.append({"location": location, "orca_dataset": dataset,
                                  "orca_distance_km": distance,
                                  **summarize_climatology(clim, args.min_months), **extremes})

    # Outputs ------------------------------------------------------------
    if orca_pairs:
        pairs = pd.DataFrame(orca_pairs)
        nearest_orca = pairs.sort_values("distance_km").drop_duplicates("location")
        assignment = assignment.merge(
            nearest_orca[["location", "dataset", "distance_km", "n_days"]].rename(columns={
                "dataset": "orca_dataset", "distance_km": "orca_distance_km",
                "n_days": "orca_n_days"}),
            on="location", how="left")
    assignment.to_csv(assignment_path, sep="\t", index=False)

    clim_all = pd.concat(clim_frames, ignore_index=True) if clim_frames else pd.DataFrame()
    clim_path = TABLES_DIR / "monthly_climatology.tsv"
    clim_all.to_csv(clim_path, sep="\t", index=False)

    predictors = sites[["location", "region", "latitude", "longitude", "n_samples"]].merge(pd.DataFrame(predictor_rows), on="location", how="left")
    predictors_path = TABLES_DIR / "site_predictors_ecology.tsv"
    predictors.to_csv(predictors_path, sep="\t", index=False)

    orca_path = TABLES_DIR / "site_summary_orca.tsv"
    pd.DataFrame(orca_rows).to_csv(orca_path, sep="\t", index=False)

    produced = [profiles_path, assignment_path, clim_path, predictors_path, orca_path]
    if not clim_all.empty:
        figure_path = FIGURES_DIR / "ecology_monthly_climatology.png"
        plot_climatologies(clim_all, figure_path)
        produced.append(figure_path)

    orca_ok = sum(1 for s in orca_status if s["status"] == "ok")
    orca_unfilled = sum(1 for s in orca_status if s["status"] != "ok"
                        and s["status"] != "no data in window")
    with open(OUTPUT_DIR / "metadata.json", "w") as handle:
        json.dump({
            "script": SCRIPT_NAME,
            "inputs": [SITE_COORDINATES.as_posix()],
            "parameters": {
                "start_year": args.start_year, "end_year": args.end_year,
                "max_depth_m": args.max_depth, "ecology_max_km": args.ecology_max_km,
                "min_profiles": args.min_profiles, "min_months": args.min_months,
                "orca_max_km": args.orca_max_km, "skip_orca": args.skip_orca, "force": args.force,
                "ecology_qc_kept": ECOLOGY_QC_PASS, "orca_qc_kept": ORCA_QC_PASS,
                "station_overrides": STATION_OVERRIDES,
                "summer_months": list(SUMMER_MONTHS), "winter_months": list(WINTER_MONTHS),
            },
            "sources": {
                "ecology_netcdf": ecology_sources,
                "orca_erddap": ERDDAP_BASE,
            },
            "orca_requests": orca_status,
            "counts": {
                "sites": int(len(sites)),
                "ecology_profiles": int(profiles["profile"].nunique()),
                "ecology_stations": int(profiles["station"].nunique()),
                "sites_with_ecology_station": int((assignment["assigned_station"] != "").sum()),
                "orca_requests_ok": orca_ok,
                "orca_requests_total": len(orca_status),
                "sites_with_orca_summary": int(len({r["location"] for r in orca_rows})),
            },
            "produced_files": [p.as_posix() for p in produced],
            "date": datetime.now(timezone.utc).isoformat(),
            "runtime_seconds": round(time.time() - started, 1),
            "software": {
                "python": platform.python_version(),
                "pandas": pd.__version__,
                "numpy": np.__version__,
                "h5py": h5py.__version__,
                "matplotlib": matplotlib.__version__,
            },
            "notes": (
                "Site coordinates are approximate centroids from step 04, not recorded collection "
                "points. Ecology stations are mid-channel monthly CTD casts and characterize the "
                "basin, not the tideflat. Coos Bay (Oregon) has no source in this step."
            ),
        }, handle, indent=2, default=str)

    log.info("Wrote %d sites (%d with an Ecology station; %d with ORCA) to %s/",
             len(sites), int((assignment["assigned_station"] != "").sum()),
             len({r["location"] for r in orca_rows}), OUTPUT_DIR)
    if orca_unfilled:
        log.warning("%d ORCA requests failed or were skipped; rerun later to fill them in "
                    "(cached years are not re-requested).", orca_unfilled)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
