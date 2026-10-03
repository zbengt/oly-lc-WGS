# 04 — Environmental data near putative sample sites

Buoy and shore-station observations assembled for every inferred sampling
location in the lc-WGS dataset, produced by
[`code/04_environmental_data.py`](../../code/04_environmental_data.py).

Site names and coordinates follow the sites documented in the top-level
[`README.md`](../../README.md). **Coordinates are approximate centroids for each
site, not recorded collection points.**

## Contents

| File | Contents |
| --- | --- |
| `site-coordinates.tsv` | One row per location: putative site, region, lat/lon, sample count |
| `stations/station-catalog.tsv` | Full NDBC + NOAA CO-OPS station catalog with coordinates (as fetched) |
| `stations/nearby-stations.tsv` | Every station within the search radius of each site, with distance |
| `observations/*_coops-*_water-temperature.tsv` | Hourly CO-OPS water temperature for the nearest station returning data |
| `observations/*_ndbc-*_standard-met.tsv` | NDBC realtime2 standard meteorological records, trimmed to the requested window |
| `environmental-summary.tsv` | Per-site, per-variable n/mean/min/max with the source station and file |
| `metadata.json` | Run parameters, sources, counts, software versions |

## Data sources

- **NOAA NDBC** — station table and `realtime2` standard meteorological feed
  (<https://www.ndbc.noaa.gov/>). Includes moored buoys and coastal/shore met
  stations; not every station reports water temperature.
- **NOAA CO-OPS Tides & Currents** — station metadata API and data getter
  (<https://api.tidesandcurrents.noaa.gov/>), used for hourly water temperature.
- **NANOOS / UW ORCA moored buoys** — Twanoh, Hoodsport, Dabob Bay, Hansville,
  Point Wells, Carr Inlet. These are the closest oceanographic moorings to most
  Puget Sound sites and carry CTD, oxygen, chlorophyll, nitrate, and met
  sensors. They are listed in `nearby-stations.tsv` as pointers and are **not**
  downloaded by this step, but their data *are* openly available from the NANOOS
  ERDDAP server (<https://data.nanoos.org/erddap/>, profiles back to 2005) — see
  [`docs/environmental-data-access-plan.md`](../../docs/environmental-data-access-plan.md)
  for the access details and the plan to ingest them.

## Current snapshot

Run parameters: 30-day window ending
2026-10-03 (UTC),
75 km search radius,
15 sites, 200 station
matches, 44 observation series.

Closest station that actually returned water temperature during the window:

| Location | Putative site | Closest station reporting water temperature | Network | Distance (km) | Mean (°C) | Min–Max (°C) |
| --- | --- | --- | --- | ---: | ---: | ---: |
| `CS18_22_Wild_plate1` | Central Sound wild collection, 2018 (Clam Bay / Manchester vicinity) | Tacoma (9446484) | NOAA CO-OPS | 35.8 | 13.74 | 13.4–14.6 |
| `Coos_Bay` | Coos Bay, Oregon | Port Orford (9431647) | NOAA CO-OPS | 69.0 | 12.05 | 10.2–14.8 |
| `Dogfish_Bay` | Dogfish Bay, Liberty Bay vicinity, Kitsap Peninsula | Port Townsend (9444900) | NOAA CO-OPS | 47.2 | 12.55 | 11.7–13.7 |
| `FB18_Wild` | Fidalgo Bay wild collection, 2018, Anacortes | Port Townsend (9444900) | NOAA CO-OPS | 44.0 | 12.54 | 11.7–13.7 |
| `Fidalgo_Bay` | Fidalgo Bay, Anacortes | Port Townsend (9444900) | NOAA CO-OPS | 44.0 | 12.54 | 11.7–13.7 |
| `HC18_Triton_Wild` | Triton Cove, Hood Canal | Tacoma (9446484) | NOAA CO-OPS | 56.8 | 13.74 | 13.4–14.6 |
| `LS` | Little Skookum Inlet, southern Puget Sound | Tacoma (9446484) | NOAA CO-OPS | 48.1 | 13.74 | 13.4–14.6 |
| `MB` | Mud Bay, Eld Inlet, southern Puget Sound | Tacoma (9446484) | NOAA CO-OPS | 46.1 | 13.74 | 13.4–14.6 |
| `NS18_Disco_Wild` | Discovery Bay, north Olympic Peninsula | Port Townsend, WA (ptww1) | NDBC | 12.0 | 12.54 | 11.7–13.6 |
| `NS18_Sequim_Wild` | Sequim Bay, north Olympic Peninsula | Port Townsend, WA (ptww1) | NDBC | 20.7 | 12.54 | 11.7–13.6 |
| `Ostrich_Bay` | Ostrich Bay, Dyes Inlet, Bremerton | Tacoma (9446484) | NOAA CO-OPS | 41.2 | 13.74 | 13.4–14.6 |
| `PGB18_Wild` | Port Gamble Bay, northern Hood Canal | Port Townsend (9444900) | NOAA CO-OPS | 32.0 | 12.54 | 11.7–13.7 |
| `SS18_North_Bay_Wild` | North Bay, Case Inlet, southern Puget Sound | Tacoma (9446484) | NOAA CO-OPS | 33.8 | 13.74 | 13.4–14.6 |
| `Squaxin_Island` | Squaxin Island, southern Puget Sound | Tacoma (9446484) | NOAA CO-OPS | 39.5 | 13.74 | 13.4–14.6 |
| `WB` | Stony Point, Willapa Bay, Washington coast | Toke Point, WA (tokw1) | NDBC | 20.8 | 16.16 | 13.8–19.0 |

Several sites share a station: the small inlets these oysters come from have no
in-water NOAA sensor of their own, so the nearest reporting station is a basin
away (Tacoma serves the South Sound sites, Port Townsend the north Sound and
Hood Canal sites; Toke Point, inside Willapa Bay, serves `WB`). Treat these
as basin-scale context, not site conditions — the
ORCA moorings above are far closer for Hood Canal and South Sound, and
Washington Department of Ecology's monthly CTD stations are within 4-21 km of
every Puget Sound site here. Both are openly accessible; see
[`docs/environmental-data-access-plan.md`](../../docs/environmental-data-access-plan.md).

## Reproducing

```bash
python code/04_environmental_data.py --days 30 --radius-km 75
```

Add `--no-download` to rebuild only the site/station crosswalk. The script
queries live NOAA endpoints, so re-running produces a new observation window.
