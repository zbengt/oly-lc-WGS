#!/usr/bin/env python3
"""
Map of all putative Olympia oyster sampling sites.

Steps:
1. Loads the approximate site coordinates defined in
   ``code/04_environmental_data.py`` (``SITES``) and the per-location sample
   counts from ``output/01_align_and_visualize/metrics/sample_metadata.tsv``.
2. Downloads Natural Earth 10 m land polygons once and caches them in
   ``output/05_site_map/cache/``.
3. Draws a regional overview (Puget Sound to Coos Bay) and a Puget Sound
   detail panel. Sites flagged uncertain are drawn as open markers.

Outputs: ``output/05_site_map/figures/site_map.png`` (+ ``.pdf``),
``output/05_site_map/tables/site_coordinates.tsv``.

All paths are relative to the repository root; do not modify files in `data/`.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Polygon, Rectangle

LAND_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
    "master/geojson/ne_10m_land.geojson"
)

REGION_COLORS = {
    "San Juan Islands": "#1f77b4",
    "Northern Puget Sound": "#17becf",
    "Strait of Juan de Fuca": "#9467bd",
    "Hood Canal": "#2ca02c",
    "Central Puget Sound": "#ff7f0e",
    "South Puget Sound": "#d62728",
    "Oregon coast": "#8c564b",
}

OVERVIEW_EXTENT = (-125.2, -121.8, 42.8, 49.1)
DETAIL_EXTENT = (-123.45, -122.15, 46.95, 48.75)

# Label offsets (points) for the detail panel, tuned to avoid overlaps.
LABEL_OFFSETS = {
    "WB": (-8, 6, "right"),
    "Fidalgo_Bay + FB18_Wild": (8, 4, "left"),
    "NS18_Sequim_Wild": (-8, 6, "right"),
    "NS18_Disco_Wild": (8, 6, "left"),
    "PGB18_Wild": (8, 4, "left"),
    "Dogfish_Bay": (8, 2, "left"),
    "HC18_Triton_Wild": (-8, 4, "right"),
    "Ostrich_Bay": (-6, -11, "right"),
    "CS18_22_Wild_plate1": (8, -6, "left"),
    "SS18_North_Bay_Wild": (8, 2, "left"),
    "Squaxin_Island": (8, 3, "left"),
    "LS": (-8, 2, "right"),
    "MB": (8, -6, "left"),
}


def load_sites(repo_root: Path) -> dict:
    spec = importlib.util.spec_from_file_location(
        "env_data", repo_root / "code" / "04_environmental_data.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.SITES


def load_land(cache_dir: Path) -> list:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / "ne_10m_land.geojson"
    if not path.exists():
        print(f"Downloading {LAND_URL}")
        urllib.request.urlretrieve(LAND_URL, path)
    features = json.loads(path.read_text())["features"]
    rings = []
    for feature in features:
        geom = feature["geometry"]
        polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
        for poly in polys:
            rings.append(poly[0])  # exterior ring; lakes are negligible here
    return rings


def clip_rings(rings: list, extent: tuple) -> list:
    x0, x1, y0, y1 = extent
    keep = []
    for ring in rings:
        xs = [p[0] for p in ring]
        ys = [p[1] for p in ring]
        if max(xs) < x0 or min(xs) > x1 or max(ys) < y0 or min(ys) > y1:
            continue
        keep.append(ring)
    return keep


def draw_base(ax, rings: list, extent: tuple) -> None:
    ax.set_facecolor("#dbe9f4")
    for ring in clip_rings(rings, extent):
        ax.add_patch(Polygon(ring, closed=True, facecolor="#f2efe6",
                             edgecolor="#8a8a8a", linewidth=0.4))
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    mid_lat = (extent[2] + extent[3]) / 2
    ax.set_aspect(1 / math.cos(math.radians(mid_lat)))
    ax.tick_params(labelsize=7)
    ax.set_xlabel("Longitude (°E)", fontsize=8)
    ax.set_ylabel("Latitude (°N)", fontsize=8)


def plot_points(ax, table: pd.DataFrame, size_scale: float) -> None:
    for _, row in table.iterrows():
        color = REGION_COLORS.get(row["region"], "black")
        ax.scatter(
            row["lon"], row["lat"], s=size_scale * row["n_samples"],
            facecolor=color if row["certain"] else "white",
            edgecolor=color, linewidth=1.6, zorder=5,
        )


def build_table(sites: dict, metadata: Path) -> pd.DataFrame:
    counts = pd.read_csv(metadata, sep="\t")["location"].value_counts()
    rows = []
    for location, info in sites.items():
        rows.append({
            "location": location,
            "putative_site": info["putative_site"],
            "region": info["region"],
            "lat": info["lat"],
            "lon": info["lon"],
            "certain": info["certain"],
            "n_samples": int(counts.get(location, 0)),
        })
    return pd.DataFrame(rows)


def merge_colocated(table: pd.DataFrame) -> pd.DataFrame:
    """Collapse sites sharing identical coordinates into one labelled point."""
    grouped = []
    for (lat, lon), grp in table.groupby(["lat", "lon"], sort=False):
        grouped.append({
            "label": " + ".join(grp["location"]),
            "region": grp["region"].iloc[0],
            "lat": lat,
            "lon": lon,
            "certain": bool(grp["certain"].all()),
            "n_samples": int(grp["n_samples"].sum()),
        })
    return pd.DataFrame(grouped)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    out_dir = repo_root / "output" / "05_site_map"
    fig_dir = out_dir / "figures"
    tab_dir = out_dir / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    sites = load_sites(repo_root)
    table = build_table(
        sites, repo_root / "output" / "01_align_and_visualize" / "metrics" / "sample_metadata.tsv"
    )
    table.to_csv(tab_dir / "site_coordinates.tsv", sep="\t", index=False)
    points = merge_colocated(table)
    rings = load_land(out_dir / "cache")

    fig, (ax_over, ax_det) = plt.subplots(
        1, 2, figsize=(10.5, 7.5), gridspec_kw={"width_ratios": [0.375, 0.48]}
    )

    # Overview panel
    draw_base(ax_over, rings, OVERVIEW_EXTENT)
    plot_points(ax_over, points, size_scale=6)
    x0, x1, y0, y1 = DETAIL_EXTENT
    ax_over.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False,
                                edgecolor="black", linewidth=1, linestyle="--", zorder=6))
    coos = points[points["label"] == "Coos_Bay"].iloc[0]
    ax_over.annotate(f"Coos_Bay (n={coos['n_samples']})", (coos["lon"], coos["lat"]),
                     xytext=(8, 0), textcoords="offset points", fontsize=8, va="center")
    ax_over.text(-122.0, 44.2, "OREGON", fontsize=9, color="#777", ha="right")
    ax_over.text(-122.0, 46.6, "WASHINGTON", fontsize=9, color="#777", ha="right")
    ax_over.text(-124.9, 45.8, "Pacific\nOcean", fontsize=9, color="#5a7fa0",
                 style="italic", ha="left")
    ax_over.set_xticks([-125, -124, -123, -122])
    ax_over.set_title("A  Region", loc="left", fontsize=10, fontweight="bold")

    # Detail panel
    draw_base(ax_det, rings, DETAIL_EXTENT)
    plot_points(ax_det, points, size_scale=14)
    for _, row in points.iterrows():
        if row["label"] == "Coos_Bay":
            continue
        dx, dy, ha = LABEL_OFFSETS.get(row["label"], (8, 0, "left"))
        suffix = "" if row["certain"] else " ?"
        ax_det.annotate(f"{row['label']} (n={row['n_samples']}){suffix}",
                        (row["lon"], row["lat"]), xytext=(dx, dy),
                        textcoords="offset points", fontsize=7.5, ha=ha, va="center",
                        zorder=7)
    ax_det.set_title("B  Puget Sound and Salish Sea", loc="left",
                     fontsize=10, fontweight="bold")

    # Legend: regions + certainty
    handles = [
        plt.Line2D([], [], marker="o", linestyle="", markersize=8,
                   markerfacecolor=c, markeredgecolor=c, label=r)
        for r, c in REGION_COLORS.items() if r in set(points["region"])
    ]
    handles += [
        plt.Line2D([], [], marker="o", linestyle="", markersize=8,
                   markerfacecolor="white", markeredgecolor="#444",
                   label="Uncertain site (?)"),
    ]
    fig.legend(handles=handles, loc="center right", fontsize=8, frameon=False,
               title="Region", title_fontsize=9)

    fig.suptitle(
        f"Olympia oyster lc-WGS sampling sites "
        f"({len(table)} locations, {table['n_samples'].sum()} samples)",
        fontsize=12,
    )
    fig.text(0.5, 0.01,
             "Coordinates are approximate centroids of putative sites (code/04_environmental_data.py); "
             "marker area scales with sample count. Land: Natural Earth 10 m.",
             ha="center", fontsize=7, color="#555")
    fig.tight_layout(rect=(0, 0.03, 0.82, 0.96), w_pad=3)

    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"site_map.{ext}", dpi=args.dpi)
    print(f"Wrote {fig_dir / 'site_map.png'}")


if __name__ == "__main__":
    main()
