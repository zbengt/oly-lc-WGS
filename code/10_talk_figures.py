#!/usr/bin/env python3
"""Presentation figures for talks, redrawn from the committed step 05-09 tables.

The step figures are built for notebooks and reports: small text, sample-code
site names and many panels. This step redraws the ones used in talks at slide
size, with labels that stay readable when a figure fills about half of a
1920 x 1080 slide:

- plain site names (``HC18_Triton_Wild`` -> "Triton Cove"),
- one colour per region, the same in every figure (Okabe-Ito, colour-blind safe),
- sites ordered geographically, outer coast -> South -> Central -> Hood Canal
  -> Strait -> North Sound,
- figure sizes in hundreds of slide pixels (``figsize=(10.4, 6.48)`` fills a
  1040 x 648 px box), saved at 200 dpi, base font 18 pt (about 25 px on the slide).

Figures (``output/10_talk_figures/figures/``):
    site_map.png     coast overview + Salish Sea panel, sites coloured by region
    depth.png        per-oyster mean depth, xbOstLuri2 vs the old assembly
    pca.png          PCAngsd PC1 vs PC2 with the three genetic groups labelled
    admixture.png    PCAngsd admixture at K = 2 and K = 3
    fst.png          weighted pairwise Fst heat map
    env_space.png    summer temperature vs salinity range, bubble = summer chlorophyll
    rda_biplot.png   environment-only RDA biplot
    null.png         individual-level permutation nulls (free vs among sites)

Every input is a committed table, so this step runs on a laptop with numpy,
pandas and matplotlib; no genotype files are needed.

Usage:
    python code/10_talk_figures.py [--force]
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

SCRIPT_NAME = "10_talk_figures"
OUT_DIR = Path("output") / SCRIPT_NAME
FIG_DIR = OUT_DIR / "figures"
LOGS_DIR = OUT_DIR / "logs"

INPUTS = {
    "land": Path("output/05_site_map/cache/ne_10m_land.geojson"),
    "sites": Path("output/05_site_map/tables/site_coordinates.tsv"),
    "depth": Path("output/05_realign_xbOstLuri2/metrics/comparison_v081.tsv"),
    "pca_scores": Path("output/06_angsd_structure/tables/pca_scores.tsv"),
    "pca_variance": Path("output/06_angsd_structure/tables/pca_variance.tsv"),
    "admix_k2": Path("output/06_angsd_structure/tables/admixture_K2.tsv"),
    "admix_k3": Path("output/06_angsd_structure/tables/admixture_K3.tsv"),
    "fst": Path("output/06_angsd_structure/tables/fst_matrix.tsv"),
    "env": Path("output/08_environmental_predictors/tables/site-env-matrix.tsv"),
    "rda_units": Path("output/09_rda/tables/rda-unit-scores.tsv"),
    "rda_biplot": Path("output/09_rda/tables/rda-biplot-scores.tsv"),
    "rda_eigen": Path("output/09_rda/tables/rda-eigenvalues.tsv"),
    "ind_null": Path("output/09_rda/tables/individual-permutation-null.tsv"),
    "ind_tests": Path("output/09_rda/tables/individual-level-tests.tsv"),
}
STEP_FOR = {"land": "05_site_map", "sites": "05_site_map", "depth": "05_realign_xbOstLuri2",
            "env": "08_environmental_predictors"}

INK = "#14303A"
MUTED = "#4A5A60"
GRID = "#E4E0D6"
ARROW = "#9C3D10"

NAME = {
    "CS18_22_Wild_plate1": "Clam Bay", "Coos_Bay": "Coos Bay", "Dogfish_Bay": "Dogfish Bay",
    "FB18_Wild": "Fidalgo Bay '18", "Fidalgo_Bay": "Fidalgo Bay", "HC18_Triton_Wild": "Triton Cove",
    "LS": "Little Skookum", "MB": "Mud Bay", "NS18_Disco_Wild": "Discovery Bay",
    "NS18_Sequim_Wild": "Sequim Bay", "Ostrich_Bay": "Ostrich Bay", "PGB18_Wild": "Port Gamble",
    "SS18_North_Bay_Wild": "North Bay", "Squaxin_Island": "Squaxin Island", "WB": "Willapa Bay",
}
REGION = {
    "Coos_Bay": "Outer coast", "WB": "Outer coast",
    "LS": "South Sound", "MB": "South Sound", "Squaxin_Island": "South Sound", "SS18_North_Bay_Wild": "South Sound",
    "CS18_22_Wild_plate1": "Central Sound", "Dogfish_Bay": "Central Sound", "Ostrich_Bay": "Central Sound",
    "HC18_Triton_Wild": "Hood Canal", "PGB18_Wild": "Hood Canal",
    "NS18_Sequim_Wild": "Strait", "NS18_Disco_Wild": "Strait",
    "FB18_Wild": "North Sound", "Fidalgo_Bay": "North Sound",
}
REGION_COLOUR = {
    "Outer coast": "#6B4C9A", "South Sound": "#D55E00", "Central Sound": "#E69F00",
    "Hood Canal": "#009E73", "Strait": "#0072B2", "North Sound": "#56B4E9",
}
# Text in the pale North Sound blue is hard to read; labels use a darker shade.
LABEL_COLOUR = {**REGION_COLOUR, "North Sound": "#2A86B8"}
GEO_ORDER = ["WB", "Coos_Bay", "LS", "MB", "Squaxin_Island", "SS18_North_Bay_Wild",
             "CS18_22_Wild_plate1", "Dogfish_Bay", "Ostrich_Bay", "HC18_Triton_Wild", "PGB18_Wild",
             "NS18_Sequim_Wild", "NS18_Disco_Wild", "Fidalgo_Bay", "FB18_Wild"]

log = logging.getLogger(SCRIPT_NAME)


def setup_logging() -> None:
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    log.setLevel(logging.INFO)
    for handler in (logging.FileHandler(LOGS_DIR / "pipeline.log", mode="w"),
                    logging.StreamHandler(sys.stdout)):
        handler.setFormatter(formatter)
        log.addHandler(handler)


def check_inputs() -> None:
    for key, path in INPUTS.items():
        if not path.exists():
            step = STEP_FOR.get(key, path.parts[1])
            raise FileNotFoundError(f"{path} is missing; run step {step} first")


def style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 18, "axes.labelsize": 20, "axes.titlesize": 22,
        "xtick.labelsize": 17, "ytick.labelsize": 17, "legend.fontsize": 17,
        "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
        "axes.edgecolor": "#8A9499", "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": "white", "axes.facecolor": "white",
    })


def region_handles(regions: list[str], size: float = 13) -> list[Line2D]:
    return [Line2D([0], [0], marker="o", ls="", ms=size, mfc=REGION_COLOUR[r], mec="white", label=r)
            for r in regions]


def save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)
    log.info("wrote %s", path)


def pc_label(variance: pd.DataFrame, pc: int) -> str:
    pct = variance.loc[variance["PC"] == pc, "variance_pct"].iloc[0]
    return f"PC{pc} ({pct:.1f}%)"


def p_label(p: float) -> str:
    return f"p = {p:.3f}" if p >= 0.001 else "p < 0.001"


# ---------------------------------------------------------------- figures

def fig_site_map(path: Path) -> None:
    geo = json.loads(INPUTS["land"].read_text())
    rings = []
    for feature in geo["features"]:
        geom = feature["geometry"]
        polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
        for poly in polys:
            ring = np.asarray(poly[0])
            if (ring[:, 0].max() < -126 or ring[:, 0].min() > -121
                    or ring[:, 1].max() < 42.5 or ring[:, 1].min() > 49.5):
                continue
            rings.append(ring)
    sites = pd.read_csv(INPUTS["sites"], sep="\t").set_index("location")

    def basemap(ax, xlim, ylim, lat0):
        ax.set_facecolor("#DCE9F2")
        for ring in rings:
            ax.add_patch(Polygon(ring, closed=True, fc="#EFEBE1", ec="#9AA3A8", lw=0.6))
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect(1 / np.cos(np.radians(lat0)))
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("#8A9499")

    fig = plt.figure(figsize=(10.4, 6.48))
    ax_coast = fig.add_axes([0.01, 0.17, 0.24, 0.82])
    ax_sound = fig.add_axes([0.27, 0.17, 0.72, 0.82])

    basemap(ax_coast, (-125.3, -121.9), (43.0, 49.0), 46)
    for loc in ["WB", "Coos_Bay"]:
        ax_coast.scatter(sites.loc[loc, "lon"], sites.loc[loc, "lat"], s=180,
                         c=REGION_COLOUR["Outer coast"], ec="white", lw=1.5, zorder=3)
        ax_coast.text(sites.loc[loc, "lon"] + 0.2, sites.loc[loc, "lat"], NAME[loc], fontsize=17, va="center")
    ax_coast.add_patch(plt.Rectangle((-123.35, 46.95), 1.15, 1.75, fill=False, ls="--", lw=1.5, ec=INK))
    ax_coast.text(-124.9, 45.4, "Pacific\nOcean", fontsize=16, style="italic", color="#4C6E8A")

    basemap(ax_sound, (-123.45, -122.15), (47.0, 48.62), 47.8)
    # (dx, dy, ha) label offsets in degrees; the two Fidalgo collections share one point
    offsets = {
        "FB18_Wild": (0.04, 0.0, "left"), "NS18_Sequim_Wild": (-0.03, 0.07, "center"),
        "NS18_Disco_Wild": (0.04, -0.06, "left"), "PGB18_Wild": (0.04, 0, "left"),
        "HC18_Triton_Wild": (-0.04, 0.05, "right"), "Dogfish_Bay": (0.04, 0.02, "left"),
        "Ostrich_Bay": (-0.03, -0.06, "right"), "CS18_22_Wild_plate1": (0.04, -0.03, "left"),
        "SS18_North_Bay_Wild": (0.04, 0, "left"), "Squaxin_Island": (0.04, 0.01, "left"),
        "LS": (-0.04, 0, "right"), "MB": (0.04, -0.04, "left"),
    }
    for loc, (dx, dy, ha) in offsets.items():
        lon, lat = sites.loc[loc, "lon"], sites.loc[loc, "lat"]
        ax_sound.scatter(lon, lat, s=260, c=REGION_COLOUR[REGION[loc]], ec="white", lw=1.8, zorder=3)
        label = "Fidalgo Bay (×2)" if loc == "FB18_Wild" else NAME[loc]
        ax_sound.text(lon + dx, lat + dy, label, fontsize=18, ha=ha, va="center", fontweight="bold")
    fig.legend(handles=region_handles(list(REGION_COLOUR), 14), loc="lower center", ncol=3,
               frameon=False, fontsize=16, columnspacing=1.2, handletextpad=0.3, bbox_to_anchor=(0.5, 0.0))
    save(fig, path)


def fig_depth(path: Path) -> None:
    depth = pd.read_csv(INPUTS["depth"], sep="\t")
    depth = depth[~depth["location"].str.startswith("Blank")]
    fig, ax = plt.subplots(figsize=(7.6, 6.0))
    fig.subplots_adjust(left=0.15, bottom=0.15, right=0.97, top=0.95)
    ax.plot([0, 9], [0, 9], ls="--", color="#8A9499", lw=1.5)
    ax.text(6.2, 5.6, "1 : 1", color=MUTED, fontsize=18)
    for region, colour in REGION_COLOUR.items():
        sub = depth[depth["location"].map(REGION) == region]
        ax.scatter(sub["meandepth_v081"], sub["meandepth_xbOstLuri2"], s=110, c=colour, ec="white", lw=1)
    ax.set_xlim(0, 9)
    ax.set_ylim(0, 9)
    ax.set_xlabel("Mean depth, old assembly (×)")
    ax.set_ylabel("Mean depth, xbOstLuri2 (×)")
    ax.grid(color=GRID)
    ax.legend(handles=region_handles(list(REGION_COLOUR)), loc="upper left", frameon=False, fontsize=15)
    save(fig, path)


def fig_pca(path: Path) -> None:
    scores = pd.read_csv(INPUTS["pca_scores"], sep="\t")
    variance = pd.read_csv(INPUTS["pca_variance"], sep="\t")
    fig, ax = plt.subplots(figsize=(11.2, 6.0))
    fig.subplots_adjust(left=0.13, bottom=0.15, right=0.70, top=0.97)
    for region, colour in REGION_COLOUR.items():
        sub = scores[scores["location"].map(REGION) == region]
        ax.scatter(sub["PC1"], sub["PC2"], s=150, c=colour, ec="white", lw=1, alpha=0.95)
    ax.set_xlabel(pc_label(variance, 1))
    ax.set_ylabel(pc_label(variance, 2))
    ax.grid(color=GRID)
    # Group labels sit in empty space next to each cluster (positions in PC units).
    ax.text(0.23, -0.015, "Outer coast", ha="center", fontsize=19, fontweight="bold",
            color=REGION_COLOUR["Outer coast"])
    ax.text(0.03, 0.185, "Hood Canal–Strait\n+ Mud Bay", ha="left", va="center", fontsize=19,
            fontweight="bold", color="#007A59")
    ax.text(-0.04, -0.115, "Central–South Sound", ha="left", fontsize=19, fontweight="bold", color="#B04A00")
    ax.set_ylim(-0.13, 0.22)
    ax.legend(handles=region_handles(list(REGION_COLOUR)), loc="center left", bbox_to_anchor=(1.02, 0.5),
              frameon=False, fontsize=17)
    save(fig, path)


def fig_admixture(path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11.0, 6.4))
    fig.subplots_adjust(left=0.1, right=0.99, top=0.97, bottom=0.3, hspace=0.12)
    rank = {loc: i for i, loc in enumerate(GEO_ORDER)}
    for ax, k in zip(axes, [2, 3]):
        q = pd.read_csv(INPUTS[f"admix_k{k}"], sep="\t")
        q = q.assign(order=q["location"].map(rank)).sort_values("order", kind="stable").reset_index(drop=True)
        comps = [c for c in q.columns if c.startswith("V")]
        # Colour each ancestry component by the group it dominates, so colours are stable across reruns.
        if k == 2:
            coast = max(comps, key=lambda c: q.loc[q["location"] == "WB", c].mean())
            stack = [(next(c for c in comps if c != coast), "#E69F00"), (coast, "#0072B2")]
        else:
            anchors = [("LS", "#E69F00"), ("HC18_Triton_Wild", "#009E73"), ("WB", "#6B4C9A")]
            stack = [(max(comps, key=lambda c, a=a: q.loc[q["location"] == a, c].mean()), col)
                     for a, col in anchors]
        bottom = np.zeros(len(q))
        for comp, colour in stack:
            ax.bar(range(len(q)), q[comp], bottom=bottom, width=1.0, color=colour, lw=0)
            bottom += q[comp].to_numpy()
        ends = q.groupby("location", sort=False).size().cumsum().to_numpy()
        for end in ends[:-1]:
            ax.axvline(end - 0.5, color="white", lw=2.5)
        ax.set_xlim(-0.5, len(q) - 0.5)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_ylabel(f"K = {k}", rotation=0, ha="right", va="center", fontsize=20)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        if k == 3:
            locs = q["location"].unique()
            mids = (np.r_[0, ends[:-1]] + ends - 1) / 2
            ax.set_xticks(mids)
            ax.set_xticklabels([NAME[loc] for loc in locs], rotation=50, ha="right", fontsize=17)
            for tick, loc in zip(ax.get_xticklabels(), locs):
                tick.set_color(LABEL_COLOUR[REGION[loc]])
            ax.tick_params(length=0)
    save(fig, path)


def fig_fst(path: Path) -> None:
    fst = pd.read_csv(INPUTS["fst"], sep="\t", index_col=0).loc[GEO_ORDER, GEO_ORDER]
    values = fst.to_numpy(dtype=float)
    np.fill_diagonal(values, np.nan)
    fig, ax = plt.subplots(figsize=(8.2, 6.48))
    fig.subplots_adjust(left=0.24, bottom=0.25, right=0.93, top=0.98)
    image = ax.imshow(values, cmap="YlOrBr", vmin=0.03, vmax=0.145)
    names = [NAME[loc] for loc in GEO_ORDER]
    ax.set_xticks(range(len(names)))
    ax.set_yticks(range(len(names)))
    ax.set_xticklabels(names, rotation=55, ha="right", fontsize=15)
    ax.set_yticklabels(names, fontsize=15)
    for tick, loc in zip(ax.get_xticklabels() + ax.get_yticklabels(), GEO_ORDER * 2):
        tick.set_color(LABEL_COLOUR[REGION[loc]])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    bar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.02)
    bar.set_label("Weighted Fst", fontsize=18)
    bar.ax.tick_params(labelsize=15)
    # Rule off the two outer-coast sites
    ax.axhline(1.5, color=INK, lw=2)
    ax.axvline(1.5, color=INK, lw=2)
    save(fig, path)


def fig_env_space(path: Path) -> None:
    env = pd.read_csv(INPUTS["env"], sep="\t")
    env = env[env["in_env_model"].astype(str) == "True"]
    stations = (env.groupby("env_station")
                .agg(temp=("temp_summer_mean", "first"), sal=("sal_range", "first"),
                     chl=("chl_summer_mean", "first"), locs=("location", list))
                .reset_index())
    shared = {("CS18_22_Wild_plate1", "Dogfish_Bay", "Ostrich_Bay"): "Clam, Dogfish & Ostrich Bays",
              ("FB18_Wild", "Fidalgo_Bay"): "Fidalgo Bay (×2)"}
    offsets = {  # (dx °C, dy PSU, ha)
        "Willapa Bay": (-0.3, 0, "right"), "Little Skookum": (-0.3, 0, "right"),
        "Triton Cove": (0.3, 0.15, "left"), "North Bay": (-0.2, 0.45, "center"), "Mud Bay": (0.3, 0, "left"),
        "Clam, Dogfish & Ostrich Bays": (0.3, -0.45, "left"), "Squaxin Island": (-0.4, -0.65, "center"),
        "Port Gamble": (0.4, 0.65, "center"), "Discovery Bay": (0, -0.65, "center"),
        "Sequim Bay": (-0.2, 0, "right"), "Fidalgo Bay (×2)": (-0.1, 0.6, "center"),
    }
    fig, ax = plt.subplots(figsize=(10.4, 5.94))
    fig.subplots_adjust(left=0.1, bottom=0.16, right=0.98, top=0.95)
    for _, row in stations.iterrows():
        label = shared.get(tuple(sorted(row["locs"])), NAME[row["locs"][0]])
        ax.scatter(row["temp"], row["sal"], s=40 + row["chl"] * 14, c=REGION_COLOUR[REGION[row["locs"][0]]],
                   ec="white", lw=1.5, alpha=0.95, zorder=3)
        dx, dy, ha = offsets.get(label, (0.3, 0, "left"))
        ax.text(row["temp"] + dx, row["sal"] + dy, label, fontsize=17, ha=ha, va="center")
    ax.set_xlabel("Summer mean temperature (°C)")
    ax.set_ylabel("Salinity range (PSU)")
    ax.grid(color=GRID)
    ax.set_xlim(9.3, 21)
    ax.set_ylim(0.5, 11.5)
    for chl in [5, 20, 50]:
        ax.scatter([], [], s=40 + chl * 14, c="#B8B2A4", ec="white", label=f"{chl} mg m⁻³")
    ax.legend(title="Summer chlorophyll", loc="upper left", frameon=False, labelspacing=1.3,
              fontsize=15, title_fontsize=16, borderpad=1)
    save(fig, path)


def fig_rda_biplot(path: Path) -> None:
    units = pd.read_csv(INPUTS["rda_units"], sep="\t")
    arrows = pd.read_csv(INPUTS["rda_biplot"], sep="\t")
    eigen = pd.read_csv(INPUTS["rda_eigen"], sep="\t").set_index("axis")
    fig, ax = plt.subplots(figsize=(10.8, 6.48))
    fig.subplots_adjust(left=0.1, bottom=0.14, right=0.98, top=0.97)
    ax.axhline(0, color="#B8B2A4", lw=1, ls=":")
    ax.axvline(0, color="#B8B2A4", lw=1, ls=":")
    for _, row in units.iterrows():
        ax.scatter(row["RDA1"], row["RDA2"], s=200, c=REGION_COLOUR[REGION[row["unit"]]], ec="white",
                   lw=1.5, zorder=3)
    # The three Central Sound sites share one Ecology station and overlap; they get one label.
    offsets = {"Squaxin_Island": (-1, -6, "right"), "HC18_Triton_Wild": (3, 3.5, "left"),
               "MB": (3, -3, "left"), "PGB18_Wild": (3, 0, "left"), "NS18_Disco_Wild": (-3, 0, "right"),
               "NS18_Sequim_Wild": (-3, 0, "right"), "Fidalgo_Bay": (3, 1.5, "left"),
               "FB18_Wild": (3, -2, "left"), "LS": (3, 1, "left"), "SS18_North_Bay_Wild": (3, 0, "left"),
               "WB": (0, 5, "center")}
    for _, row in units.iterrows():
        if row["unit"] in offsets:
            dx, dy, ha = offsets[row["unit"]]
            ax.text(row["RDA1"] + dx, row["RDA2"] + dy, NAME[row["unit"]], fontsize=16, ha=ha, va="center")
    central = units[units["unit"].isin(["CS18_22_Wild_plate1", "Dogfish_Bay", "Ostrich_Bay"])]
    ax.text(central["RDA1"].mean(), central["RDA2"].mean() + 8.5, "Clam, Dogfish\n& Ostrich Bays",
            fontsize=16, ha="center", va="center")
    nice = {"temp_summer_mean": "Summer temp.", "sal_range": "Salinity range", "chl_summer_mean": "Summer chl."}
    scale = 60
    for _, row in arrows.iterrows():
        x, y = row["RDA1"] * scale, row["RDA2"] * scale
        ax.annotate("", xy=(x, y), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", color=ARROW, lw=2.5, mutation_scale=22))
        if row["term"] == "chl_summer_mean":
            tx, ty = -17, -6
        else:
            tx, ty = x * 1.12, y * 1.12 + (2 if row["term"] == "temp_summer_mean" else 0)
        ax.text(tx, ty, nice.get(row["term"], row["term"]), color=ARROW, fontsize=18, fontweight="bold",
                ha="center", va="center")
    ax.set_xlabel(f"RDA1 ({eigen.loc['RDA1', 'prop_total'] * 100:.1f}% of total variance)")
    ax.set_ylabel(f"RDA2 ({eigen.loc['RDA2', 'prop_total'] * 100:.1f}%)")
    ax.set_xlim(-45, 100)
    ax.set_ylim(-55, 72)
    ax.legend(handles=region_handles(list(REGION_COLOUR)), loc="lower right", frameon=False, fontsize=15)
    save(fig, path)


def fig_null(path: Path) -> None:
    null = pd.read_csv(INPUTS["ind_null"], sep="\t")
    tests = pd.read_csv(INPUTS["ind_tests"], sep="\t")
    env = tests[tests["model"] == "env"].set_index("permutation")
    observed = env["F"].iloc[0]
    p_free = env.loc["free (pseudoreplicated)", "p"]
    p_sites = env.loc["among sites", "p"]
    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    fig.subplots_adjust(left=0.1, bottom=0.22, right=0.98, top=0.95)
    bins = np.linspace(0.9, 2.0, 56)
    ax.hist(null.loc[null["permutation"] == "free", "F"], bins=bins, color="#E69F00", alpha=0.85)
    ax.hist(null.loc[null["permutation"] == "among sites", "F"], bins=bins, color="#2E6E7E", alpha=0.8)
    ax.axvline(observed, color=INK, lw=3)
    ax.text(observed + 0.02, 230, f"observed\nF = {observed:.2f}", fontsize=16, va="top")
    ax.text(1.06, 395, f"Oysters shuffled freely\n{p_label(p_free)}", va="top", fontsize=16,
            color="#9A6400", fontweight="bold")
    ax.text(1.17, 175, f"Whole sites shuffled\n{p_label(p_sites)}", va="top", fontsize=16,
            color="#1F5560", fontweight="bold")
    ax.set_xlabel("Pseudo-F under the null")
    ax.set_ylabel("Permutations")
    ax.set_yticks([0, 200, 400])
    save(fig, path)


FIGURES = {
    "site_map.png": fig_site_map, "depth.png": fig_depth, "pca.png": fig_pca,
    "admixture.png": fig_admixture, "fst.png": fig_fst, "env_space.png": fig_env_space,
    "rda_biplot.png": fig_rda_biplot, "null.png": fig_null,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--force", action="store_true", help="redraw figures that already exist")
    args = parser.parse_args()

    start = time.time()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    setup_logging()
    check_inputs()
    style()

    written = []
    for name, draw in FIGURES.items():
        path = FIG_DIR / name
        if path.exists() and not args.force:
            log.info("skip %s (exists; use --force to redraw)", path)
        else:
            draw(path)
        written.append(str(path))

    metadata = {
        "script": f"code/{SCRIPT_NAME}.py",
        "date": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "runtime_seconds": round(time.time() - start, 1),
        "parameters": {"force": args.force, "dpi": 200, "base_font_pt": 18},
        "inputs": {key: str(path) for key, path in INPUTS.items()},
        "software": {"python": sys.version.split()[0], "numpy": np.__version__,
                     "pandas": pd.__version__, "matplotlib": matplotlib.__version__},
        "outputs": {"figures": written, "log": str(LOGS_DIR / "pipeline.log")},
    }
    (OUT_DIR / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    log.info("wrote %s", OUT_DIR / "metadata.json")


if __name__ == "__main__":
    main()
