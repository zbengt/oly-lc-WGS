#!/usr/bin/env python3
"""Assemble environmental data from buoys/stations near each putative sample site.

For every inferred sampling location in
``output/01_align_and_visualize/metrics/sample_metadata.tsv`` the script:

1. attaches the putative site name and approximate coordinates documented in
   the repository README,
2. finds the closest NOAA NDBC buoys and NOAA CO-OPS tide/met stations,
3. downloads recent observations (CO-OPS water temperature/meteorology and the
   NDBC realtime2 standard meteorological feed) for the nearest stations, and
4. writes tables, a per-site summary, and ``metadata.json`` to
   ``output/04_environmental_data/``.

Usage:
    python code/04_environmental_data.py [--days 30] [--radius-km 75]
"""

from __future__ import annotations

import argparse
import io
import json
import math
import platform
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pandas as pd

SCRIPT_NAME = "04_environmental_data.py"
SAMPLE_METADATA = Path("output/01_align_and_visualize/metrics/sample_metadata.tsv")
OUTPUT_DIR = Path("output/04_environmental_data")
STATIONS_DIR = OUTPUT_DIR / "stations"
OBS_DIR = OUTPUT_DIR / "observations"

NDBC_STATION_TABLE = "https://www.ndbc.noaa.gov/data/stations/station_table.txt"
NDBC_REALTIME = "https://www.ndbc.noaa.gov/data/realtime2/{station}.txt"
COOPS_MDAPI = "https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations.json?type={kind}"
COOPS_DATA = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"

USER_AGENT = "oly-lc-WGS environmental data assembly (github.com/zbengt/oly-lc-WGS)"

# Approximate coordinates for the putative sites named in README.md. Sites
# flagged ``uncertain`` there carry the same flag here: the coordinate is a
# best-guess centroid for the interpreted site, not a recorded collection point.
SITES = {
    "CS18_22_Wild_plate1": {
        "putative_site": "Central Sound wild collection, 2018 (Clam Bay / Manchester vicinity)",
        "lat": 47.5760, "lon": -122.5460, "region": "Central Puget Sound", "certain": False,
    },
    "Coos_Bay": {
        "putative_site": "Coos Bay, Oregon",
        "lat": 43.3450, "lon": -124.3170, "region": "Oregon coast", "certain": True,
    },
    "Dogfish_Bay": {
        "putative_site": "Dogfish Bay, Liberty Bay vicinity, Kitsap Peninsula",
        "lat": 47.6960, "lon": -122.6300, "region": "Central Puget Sound", "certain": True,
    },
    "FB18_Wild": {
        "putative_site": "Fidalgo Bay wild collection, 2018, Anacortes",
        "lat": 48.4880, "lon": -122.5760, "region": "Northern Puget Sound", "certain": True,
    },
    "Fidalgo_Bay": {
        "putative_site": "Fidalgo Bay, Anacortes",
        "lat": 48.4880, "lon": -122.5760, "region": "Northern Puget Sound", "certain": True,
    },
    "HC18_Triton_Wild": {
        "putative_site": "Triton Cove, Hood Canal",
        "lat": 47.6070, "lon": -122.9770, "region": "Hood Canal", "certain": True,
    },
    "LS": {
        "putative_site": "Little Skookum Inlet, southern Puget Sound",
        "lat": 47.1600, "lon": -123.0300, "region": "South Puget Sound", "certain": False,
    },
    "MB": {
        "putative_site": "Mud Bay, Eld Inlet, southern Puget Sound",
        "lat": 47.1080, "lon": -122.9770, "region": "South Puget Sound", "certain": False,
    },
    "NS18_Disco_Wild": {
        "putative_site": "Discovery Bay, north Olympic Peninsula",
        "lat": 48.0450, "lon": -122.8880, "region": "Strait of Juan de Fuca", "certain": True,
    },
    "NS18_Sequim_Wild": {
        "putative_site": "Sequim Bay, north Olympic Peninsula",
        "lat": 48.0640, "lon": -123.0300, "region": "Strait of Juan de Fuca", "certain": True,
    },
    "Ostrich_Bay": {
        "putative_site": "Ostrich Bay, Dyes Inlet, Bremerton",
        "lat": 47.5860, "lon": -122.6900, "region": "Central Puget Sound", "certain": True,
    },
    "PGB18_Wild": {
        "putative_site": "Port Gamble Bay, northern Hood Canal",
        "lat": 47.8500, "lon": -122.5800, "region": "Hood Canal", "certain": True,
    },
    "SS18_North_Bay_Wild": {
        "putative_site": "North Bay, Case Inlet, southern Puget Sound",
        "lat": 47.3800, "lon": -122.8300, "region": "South Puget Sound", "certain": True,
    },
    "Squaxin_Island": {
        "putative_site": "Squaxin Island, southern Puget Sound",
        "lat": 47.1800, "lon": -122.9200, "region": "South Puget Sound", "certain": True,
    },
    "WB": {
        "putative_site": "Westcott Bay, San Juan Island",
        "lat": 48.5850, "lon": -123.1600, "region": "San Juan Islands", "certain": False,
    },
}

# Curated reference list of NANOOS/UW ORCA moored buoys. These are the closest
# real-time oceanographic moorings to most Puget Sound sites, but their data are
# distributed through NANOOS (nvs.nanoos.org) rather than a public bulk endpoint,
# so they are recorded here as pointers rather than downloaded.
ORCA_BUOYS = [
    ("ORCA_Twanoh", "Twanoh, southern Hood Canal", 47.375, -123.008),
    ("ORCA_Hoodsport", "Hoodsport, Hood Canal", 47.425, -123.113),
    ("ORCA_DabobBay", "Dabob Bay, Hood Canal", 47.803, -122.803),
    ("ORCA_NorthBuoy", "North Buoy, Hood Canal", 47.900, -122.622),
    ("ORCA_PointWells", "Point Wells, central Puget Sound", 47.761, -122.397),
    ("ORCA_CarrInlet", "Carr Inlet, southern Puget Sound", 47.280, -122.730),
]

BLANK_LOCATIONS = {"Blank"}


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Great-circle distance in kilometres."""
    radius = 6371.0088
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = p2 - p1
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * radius * math.asin(math.sqrt(a))


def fetch(url: str, timeout: int = 60, retries: int = 3) -> str:
    """GET a URL as text, retrying transient failures."""
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            last_error = error
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"failed to fetch {url}: {last_error}")


def parse_ndbc_dms(location: str) -> tuple[float, float] | None:
    """Extract decimal lat/lon from the NDBC station-table LOCATION field."""
    parts = location.split()
    if len(parts) < 4:
        return None
    try:
        lat = float(parts[0])
        lon = float(parts[2])
    except ValueError:
        return None
    if parts[1].upper() == "S":
        lat = -lat
    if parts[3].upper() == "W":
        lon = -lon
    return lat, lon


def load_ndbc_stations() -> pd.DataFrame:
    text = fetch(NDBC_STATION_TABLE)
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("|")
        if len(fields) < 7:
            continue
        coords = parse_ndbc_dms(fields[6])
        if coords is None:
            continue
        rows.append({
            "network": "NDBC",
            "station_id": fields[0].strip(),
            "station_name": fields[4].strip(),
            "station_type": fields[2].strip(),
            "owner": fields[1].strip(),
            "lat": coords[0],
            "lon": coords[1],
        })
    return pd.DataFrame(rows)


def load_coops_stations() -> pd.DataFrame:
    """CO-OPS stations, tagged with the observation types each one reports."""
    frames = {}
    for kind, label in [
        ("watertemp", "water_temperature"),
        ("met", "meteorology"),
        ("waterlevels", "water_level"),
    ]:
        payload = json.loads(fetch(COOPS_MDAPI.format(kind=kind)))
        for station in payload.get("stations", []):
            sid = str(station.get("id"))
            entry = frames.setdefault(sid, {
                "network": "NOAA CO-OPS",
                "station_id": sid,
                "station_name": station.get("name", ""),
                "station_type": "",
                "owner": "NOAA CO-OPS",
                "lat": station.get("lat"),
                "lon": station.get("lng"),
                "products": set(),
            })
            entry["products"].add(label)
    rows = []
    for entry in frames.values():
        if entry["lat"] is None or entry["lon"] is None:
            continue
        entry["station_type"] = ",".join(sorted(entry.pop("products")))
        rows.append(entry)
    return pd.DataFrame(rows)


def nearest_stations(site: dict, stations: pd.DataFrame, radius_km: float, top_n: int) -> pd.DataFrame:
    distances = stations.apply(
        lambda row: haversine_km(site["lat"], site["lon"], row["lat"], row["lon"]), axis=1
    )
    nearby = stations.assign(distance_km=distances.round(2))
    nearby = nearby[nearby["distance_km"] <= radius_km]
    return nearby.sort_values("distance_km").head(top_n)


def download_coops(station_id: str, product: str, begin: str, end: str) -> pd.DataFrame:
    url = (
        f"{COOPS_DATA}?product={product}&application=oly-lc-WGS&begin_date={begin}"
        f"&end_date={end}&datum=MLLW&station={station_id}&time_zone=gmt"
        f"&units=metric&interval=h&format=csv"
    )
    text = fetch(url)
    if not text.strip() or text.lstrip().startswith("{") or "Error" in text[:200]:
        return pd.DataFrame()
    frame = pd.read_csv(io.StringIO(text))
    frame.columns = [c.strip().lower().replace(" ", "_") for c in frame.columns]
    return frame


def trim_ndbc_window(frame: pd.DataFrame, begin: datetime) -> pd.DataFrame:
    """Keep only realtime2 rows at or after ``begin`` (feed spans ~45 days)."""
    needed = {"YY", "MM", "DD", "hh", "mm"}
    if frame.empty or not needed.issubset(frame.columns):
        return frame
    stamps = pd.to_datetime(dict(
        year=frame["YY"], month=frame["MM"], day=frame["DD"],
        hour=frame["hh"], minute=frame["mm"],
    ), errors="coerce", utc=True)
    return frame[stamps >= begin]


def download_ndbc_realtime(station_id: str) -> pd.DataFrame:
    """Standard meteorological realtime2 feed (~45 days of observations)."""
    try:
        text = fetch(NDBC_REALTIME.format(station=station_id.upper()), retries=2)
    except RuntimeError:
        return pd.DataFrame()
    lines = text.splitlines()
    if len(lines) < 3:
        return pd.DataFrame()
    header = lines[0].lstrip("#").split()
    units = lines[1].lstrip("#").split()
    records = [line.split() for line in lines[2:] if line.strip()]
    frame = pd.DataFrame([r for r in records if len(r) == len(header)], columns=header)
    frame = frame.replace({"MM": pd.NA})
    for column in frame.columns:
        converted = pd.to_numeric(frame[column], errors="coerce")
        if converted.notna().any():
            frame[column] = converted
    frame.attrs["units"] = dict(zip(header, units))
    return frame


def summarize_numeric(frame: pd.DataFrame, column: str) -> dict:
    if frame.empty or column not in frame:
        return {}
    values = pd.to_numeric(frame[column], errors="coerce").dropna()
    if values.empty:
        return {}
    return {
        "n": int(values.size),
        "mean": round(float(values.mean()), 3),
        "min": round(float(values.min()), 3),
        "max": round(float(values.max()), 3),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--days", type=int, default=30,
                        help="days of recent observations to download (default: 30)")
    parser.add_argument("--radius-km", type=float, default=75.0,
                        help="search radius for nearby stations (default: 75)")
    parser.add_argument("--top-n", type=int, default=5,
                        help="stations to list per site and network (default: 5)")
    parser.add_argument("--no-download", action="store_true",
                        help="only build the site/station crosswalk, skip observations")
    args = parser.parse_args()

    if not SAMPLE_METADATA.exists():
        print(f"ERROR: {SAMPLE_METADATA} not found; run code/01_align_and_visualize.py first.",
              file=sys.stderr)
        return 1

    started = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    STATIONS_DIR.mkdir(parents=True, exist_ok=True)
    OBS_DIR.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(SAMPLE_METADATA, sep="\t")
    counts = metadata["location"].value_counts().to_dict()

    print("Fetching station catalogs ...")
    ndbc = load_ndbc_stations()
    coops = load_coops_stations()
    catalog = pd.concat([ndbc, coops], ignore_index=True)
    catalog.to_csv(STATIONS_DIR / "station-catalog.tsv", sep="\t", index=False)
    print(f"  NDBC stations: {len(ndbc)}   CO-OPS stations: {len(coops)}")

    site_rows, station_rows, summary_rows = [], [], []
    end = datetime.now(timezone.utc)
    begin = end - timedelta(days=args.days)
    begin_str, end_str = begin.strftime("%Y%m%d"), end.strftime("%Y%m%d")

    for location in sorted(counts):
        if location in BLANK_LOCATIONS:
            continue
        site = SITES.get(location)
        if site is None:
            print(f"  WARNING: no coordinate defined for location '{location}'; skipped.")
            continue

        site_rows.append({
            "location": location,
            "putative_site": site["putative_site"],
            "region": site["region"],
            "latitude": site["lat"],
            "longitude": site["lon"],
            "coordinate_confidence": "high" if site["certain"] else "uncertain",
            "n_samples": counts[location],
        })

        selected = pd.concat([
            nearest_stations(site, ndbc, args.radius_km, args.top_n),
            nearest_stations(site, coops, args.radius_km, args.top_n),
        ], ignore_index=True).sort_values("distance_km")

        for orca_id, orca_name, orca_lat, orca_lon in ORCA_BUOYS:
            distance = haversine_km(site["lat"], site["lon"], orca_lat, orca_lon)
            if distance <= args.radius_km:
                selected = pd.concat([selected, pd.DataFrame([{
                    "network": "NANOOS/UW ORCA",
                    "station_id": orca_id,
                    "station_name": orca_name,
                    "station_type": "moored buoy (CTD, oxygen, chlorophyll, met)",
                    "owner": "University of Washington",
                    "lat": orca_lat, "lon": orca_lon,
                    "distance_km": round(distance, 2),
                }])], ignore_index=True).sort_values("distance_km")

        selected.insert(0, "location", location)
        selected.insert(1, "putative_site", site["putative_site"])
        station_rows.append(selected)
        print(f"{location}: {len(selected)} stations within {args.radius_km:g} km")

        if args.no_download:
            continue

        # Nearest CO-OPS station reporting water temperature. Searched against
        # the full catalog rather than the per-site top-N list, because the very
        # closest CO-OPS stations often report water level only.
        # A station can advertise water temperature yet return nothing for the
        # requested window (sensor offline), so walk outwards until one delivers.
        wt_catalog = coops[coops["station_type"].str.contains("water_temperature", na=False)]
        for _, station in nearest_stations(site, wt_catalog, args.radius_km, args.top_n).iterrows():
            frame = download_coops(station["station_id"], "water_temperature", begin_str, end_str)
            if frame.empty:
                continue
            path = OBS_DIR / f"{location}_coops-{station['station_id']}_water-temperature.tsv"
            frame.to_csv(path, sep="\t", index=False)
            stats = summarize_numeric(frame, "water_temperature")
            summary_rows.append({
                "location": location, "putative_site": site["putative_site"],
                "network": "NOAA CO-OPS", "station_id": station["station_id"],
                "station_name": station["station_name"],
                "distance_km": station["distance_km"], "variable": "water_temperature_degC",
                "file": path.as_posix(), **stats,
            })
            print(f"  CO-OPS {station['station_id']} ({station['station_name']}) "
                  f"water temperature: {stats}")
            break

        # Nearest NDBC buoy with a live standard-met feed.
        for _, station in selected[selected["network"] == "NDBC"].iterrows():
            frame = trim_ndbc_window(download_ndbc_realtime(station["station_id"]), begin)
            if frame.empty:
                continue
            path = OBS_DIR / f"{location}_ndbc-{station['station_id']}_standard-met.tsv"
            frame.to_csv(path, sep="\t", index=False)
            for column, label in [("WTMP", "water_temperature_degC"),
                                  ("ATMP", "air_temperature_degC"),
                                  ("WSPD", "wind_speed_m_s")]:
                stats = summarize_numeric(frame, column)
                if stats:
                    summary_rows.append({
                        "location": location, "putative_site": site["putative_site"],
                        "network": "NDBC", "station_id": station["station_id"],
                        "station_name": station["station_name"],
                        "distance_km": station["distance_km"], "variable": label,
                        "file": path.as_posix(), **stats,
                    })
            print(f"  NDBC {station['station_id']} realtime2: {len(frame)} records")
            break

    sites_frame = pd.DataFrame(site_rows)
    sites_frame.to_csv(OUTPUT_DIR / "site-coordinates.tsv", sep="\t", index=False)

    stations_frame = pd.concat(station_rows, ignore_index=True) if station_rows else pd.DataFrame()
    stations_frame.to_csv(STATIONS_DIR / "nearby-stations.tsv", sep="\t", index=False)

    summary_frame = pd.DataFrame(summary_rows)
    summary_frame.to_csv(OUTPUT_DIR / "environmental-summary.tsv", sep="\t", index=False)

    with open(OUTPUT_DIR / "metadata.json", "w") as handle:
        json.dump({
            "script": SCRIPT_NAME,
            "input": SAMPLE_METADATA.as_posix(),
            "parameters": {
                "days": args.days,
                "radius_km": args.radius_km,
                "top_n": args.top_n,
                "no_download": args.no_download,
                "observation_window_utc": [begin.isoformat(), end.isoformat()],
            },
            "sources": {
                "ndbc_station_table": NDBC_STATION_TABLE,
                "ndbc_realtime2": NDBC_REALTIME,
                "coops_metadata_api": COOPS_MDAPI,
                "coops_data_api": COOPS_DATA,
                "nanoos_orca": "https://nvs.nanoos.org/",
            },
            "counts": {
                "sites": int(len(sites_frame)),
                "station_matches": int(len(stations_frame)),
                "observation_series": int(len(summary_frame)),
            },
            "date": datetime.now(timezone.utc).isoformat(),
            "runtime_seconds": round(time.time() - started, 1),
            "software": {
                "python": platform.python_version(),
                "pandas": pd.__version__,
            },
            "notes": (
                "Site coordinates are approximate centroids for the putative sites named in "
                "README.md, not recorded collection points. Sites marked 'uncertain' have "
                "ambiguous name prefixes and their station matches should be re-checked once "
                "the original collection records are available."
            ),
        }, handle, indent=2)

    print(f"\nWrote {len(sites_frame)} sites, {len(stations_frame)} station matches, "
          f"{len(summary_frame)} observation series to {OUTPUT_DIR}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
