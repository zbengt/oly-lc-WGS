# Environmental data access plan — NANOOS and WA Dept. of Ecology

Follow-up to [`output/04_environmental_data/`](../output/04_environmental_data/README.md),
which found that no NOAA in-water sensor sits in any of the inlets these oysters
were collected from. This document records what closer data actually exist, a
specific inquiry to send to NANOOS, and the plan for pulling WA Dept. of Ecology
marine monitoring data.

Everything marked *verified* below was checked against the live endpoint on
2026-09-04.

---

## 0. Correction to the step-04 write-up

The step-04 README stated that NANOOS/UW ORCA mooring data are "not available
from a public bulk endpoint." **That is wrong.** NANOOS runs a public ERDDAP
server at <https://data.nanoos.org/erddap/> that serves the full ORCA record
with no authentication:

| Dataset | Type | Coverage (verified) | Casts | Max depth |
| --- | --- | --- | ---: | ---: |
| `orca_hydro_twanoh` | griddap | 2005-01-19 → 2026-09-04 | 37,585 | 35 m |
| `orca_hydro_hoodsport` | griddap | 2005-10-20 → 2026-09-04 | 30,277 | 120 m |
| `orca_hydro_hansville` | griddap | 2005-11-09 → 2026-06-24 | 16,087 | 100 m |
| `orca_hydro_dabobbay` | griddap | 2010-06-01 → **2025-05-02** | 12,419 | 105 m |
| `orca_hydro_carrinlet` | griddap | 2010-10-13 → 2026-09-04 | 19,377 | 47 m |
| `orca_hydro_pointwells` | griddap | 2010-04-01 → 2026-09-04 | 15,002 | 100 m |

Profile variables: temperature, practical salinity, conductivity, sigma-theta,
dissolved oxygen (mass concentration and % saturation), chlorophyll-*a*,
nitrate, PAR, turbidity — each with a QC aggregate flag, gridded to 0.25 m bins.
Matching surface meteorology is served as tabledap (`orca_met_*`), and
`ORCA_*` tabledap datasets add air/sea pCO2 at some moorings.

License (from the dataset metadata): free use and redistribution, not intended
for legal use.

Two mechanics worth writing down, both of which cost time to find:

- The profile data are **griddap**, not tabledap; `tabledap/orca_hydro_*`
  returns 404. Request form is
  `griddap/<id>.csv?<var>[(<start>):(<end>)][(<mindepth>):(<maxdepth>)]`.
- The server sits behind Tomcat, which rejects raw `<`, `>`, `[`, `]` in the
  request target. Percent-encode them (`%3C`, `%3E`, `%5B`, `%5D`) or the request
  fails with `HTTP 400 – Invalid character found in the request target`.

Verified sample pull (Twanoh, 2025-08-28, 4.5 m): 16.06 °C, 28.20 PSU,
6.53 mg/L O₂, 7.34 mg/m³ chlorophyll-*a*.

**Consequence:** the ORCA data no longer need an access request. The NANOOS
inquiry below is narrower — it is about interpretation, gaps, and the one
mooring whose record stops.

---

## 1. Specific request to NANOOS / NWEM

**Who.** The ORCA moorings are operated by the Northwest Environmental Moorings
(NWEM) group at UW-APL (<https://nwem.apl.washington.edu/>), not by NANOOS
itself — NANOOS distributes the data. The ERDDAP metadata lists Seth Travis
(`setht1@uw.edu`) as creator/publisher and John Mickett (`mickett@uw.edu`) as a
contributor. Send to NWEM directly and copy the NANOOS help address from
<https://nanoos.org/>; NANOOS is the right cc for anything about the ERDDAP
service itself rather than the moorings.

**What to ask.** Five questions, all of them things the public metadata does not
answer and all of them decision-relevant for pairing environmental data to the
lc-WGS population structure:

1. **Dabob Bay stopped.** `orca_hydro_dabobbay` ends 2025-05-02 while the other
   five moorings run to the present. Is the mooring recovered, is the gap a
   processing backlog, and is a return planned? Dabob Bay is our closest proxy
   for the Hood Canal sites.
2. **Which QC level to use.** The L3 gridded product carries `*_qc_aggregate`
   flags. Which flag values are safe to keep for multi-year climatology, and is
   there a documented flag-value → meaning table beyond the CF attributes?
3. **Pre-2021 record.** Twanoh and Hoodsport go back to 2005, but the tabledap
   `ORCA_*` datasets start in 2021. Is the 2005–2020 portion of the gridded
   product processed to the same standard, or were there sensor/calibration
   changes we should treat as a break in the series?
4. **Depth to use for an intertidal/shallow-subtidal organism.** Olympia oysters
   sit at roughly 0–2 m MLLW. Is the shallowest reliable bin from a winched CTD
   the 0.25–1 m range, or is there a near-surface bin the group treats as
   unreliable (wave contamination, wiper interference)?
5. **Anything closer.** Do NWEM or NANOOS partners hold non-ORCA data — pilot
   moorings, shore stations, tribal or shellfish-grower sensors served through
   NVS — inside Case Inlet, Totten/Little Skookum, Eld Inlet, Dyes Inlet,
   Port Gamble Bay, Fidalgo Bay, or Westcott Bay? These are the sites where our
   nearest reporting station is currently 30–57 km away.

**Draft message** — review and send yourself; I have not sent anything.

> Subject: ORCA mooring data for an Olympia oyster population-genomics study — QC and coverage questions
>
> Hi Seth (cc: John, NANOOS),
>
> I'm working on a low-coverage whole-genome sequencing study of Olympia oyster
> (*Ostrea lurida*) populations across Puget Sound — roughly 110 samples from 15
> sites including Case Inlet, Totten/Little Skookum, Eld Inlet, Squaxin Island,
> Dyes Inlet, Dogfish Bay, Port Gamble Bay, Triton Cove, Discovery and Sequim
> bays, and Fidalgo Bay. We're pairing the genomic structure with environmental
> context, and the ORCA moorings are by far the closest oceanographic records to
> most of our sites.
>
> I've pulled the L3 gridded profiles from the NANOOS ERDDAP and they're exactly
> what we need. Five questions before we build the analysis on them:
>
> 1. `orca_hydro_dabobbay` ends 2025-05-02 while the other moorings run current.
>    Is that a recovery, a processing backlog, or a longer outage — and is a
>    redeployment planned? It's our closest proxy for the Hood Canal sites.
> 2. For multi-year climatology, which `*_qc_aggregate` values do you recommend
>    keeping, and is there documentation of the flag scheme beyond the CF
>    attributes?
> 3. The gridded record starts in 2005 at Twanoh and Hoodsport. Is the pre-2021
>    portion processed to the same standard, or are there sensor or calibration
>    changes we should treat as breaks in the series?
> 4. Our animals sit at about 0–2 m MLLW. What's the shallowest bin you'd trust
>    from the winched CTD casts?
> 5. Are there any other data you serve or know of — pilot moorings, shore
>    stations, tribal or grower-operated sensors — inside Case Inlet, Totten or
>    Little Skookum Inlet, Eld Inlet, Dyes Inlet, Port Gamble Bay, Fidalgo Bay,
>    or Westcott Bay? Those are the sites where our nearest reporting station is
>    still 30+ km away.
>
> Happy to share what we build back, and we'll cite the NWEM/NANOOS data as
> requested in the dataset metadata.
>
> Thanks,
> Steven Roberts
> University of Washington, School of Aquatic and Fishery Sciences

**Ask before sending.** Confirm the current NWEM contact — the addresses above
come from dataset metadata, which can lag staffing.

---

## 2. WA Dept. of Ecology marine monitoring — investigated

Ecology's Long-term Marine Water Quality Monitoring Program flies monthly CTD
casts to fixed stations across the Salish Sea, Willapa Bay, and Grays Harbor.
**This is the closest routine in-water sampling to our sites, by a wide margin.**

### Access paths (verified)

| Path | URL | Notes |
| --- | --- | --- |
| Yearly netCDF, 1999–present | `https://fortress.wa.gov/ecy/ezshare/EAP/SalishSea/netCDF-files/MarineWaterProfilesAndNutrientsYear<YYYY>.nc` | Direct download, no auth. 2025 file is 39 MB, netCDF-4/HDF5. **Verified: downloaded and parsed.** |
| EIM database (study `MarineWater`) | <https://apps.ecology.wa.gov/eim/search/default.aspx> | Web search/download UI; needed for discrete samples not in the netCDF |
| Program page | <https://ecology.wa.gov/research-data/monitoring-assessment/puget-sound-and-marine-monitoring/water-column-data> | Documents the netCDF and EIM routes |
| Ferry transects | Salish and Victoria Clipper THREDDS catalogs linked from the program page | Underway surface data, central Sound tracks |

Note: `apps.ecology.wa.gov` returns 403 to a default curl user agent — set a
browser user agent. The fortress.wa.gov netCDF path has no such restriction.

### Structure of the netCDF (verified against the 2025 file)

Ragged-array CF profile format: 39 stations, 266 profiles, 40,222 observation
rows. Per-observation variables with QA/QC/QF flag triplets:

- **Sensor:** `Pres`, `Depth`, `Temp` (°C), `Cond`, `Salinity` (PSU),
  `Density`, `DOAdjusted` (mg/L), `DOSatAdjusted`, `Turb`, `FluorAdjusted`
  (mg/m³), `Xmiss_25cm`, `BatC`, `PAR`
- **Discrete nutrients:** `NO3`, `NO2`, `NH4`, `PO4`, `SiOH4`
- **Indexing:** `Station` (char), `Latitude`, `Longitude`, `station_index`,
  `row_size`, `profile_index`, `FieldDate` (days since 1970-01-01, UTC)

Missing values are `-99999.9`, not NaN — filter explicitly. Reading requires
`h5py`, `netCDF4`, or `xarray`; `scipy.io.netcdf` will not open these
(HDF5-backed netCDF-4).

### Nearest Ecology station to each site (verified, computed from the 2025 file)

| Location | Putative site | Nearest Ecology station | km | vs. nearest NOAA station |
| --- | --- | --- | ---: | ---: |
| `Squaxin_Island` | Squaxin Island | DNA001 | 4.3 | 39.5 |
| `MB` | Mud Bay, Eld Inlet | BUD005 | 4.9 | 46.1 |
| `Ostrich_Bay` | Ostrich Bay, Dyes Inlet | SIN001 | 5.4 | 41.2 |
| `LS` | Little Skookum Inlet | OAK004 | 6.9 | 48.1 |
| `SS18_North_Bay_Wild` | North Bay, Case Inlet | HCB007 | 7.7 | 33.8 |
| `CS18_22_Wild_plate1` | Central Sound | SIN001 | 7.7 | 35.8 |
| `HC18_Triton_Wild` | Triton Cove, Hood Canal | HCB003 | 8.0 | 56.8 |
| `PGB18_Wild` | Port Gamble Bay | ADM003 | 8.0 | 32.0 |
| `NS18_Disco_Wild` | Discovery Bay | PTH005 | 10.2 | 12.0 |
| `Dogfish_Bay` | Dogfish Bay | OCH014 | 3.7 | 47.2 |
| `NS18_Sequim_Wild` | Sequim Bay | ADM002 | 19.6 | 20.7 |
| `FB18_Wild` / `Fidalgo_Bay` | Fidalgo Bay | RSR837 | 19.8 | 44.0 |
| `WB` | Westcott Bay | SJF000 | 21.2 | 55.2 |
| `Coos_Bay` | Coos Bay, OR | — | 346 | 69.0 |

Every Puget Sound site gains a closer in-water record, most of them
dramatically — and unlike the NOAA stations, these are full-depth profiles with
oxygen, chlorophyll, and nutrients rather than a single surface temperature.

### Caveats

- **Monthly, not continuous.** One cast per station per month, weather
  permitting. Good for seasonal climatology and among-site contrasts; useless
  for extremes, tidal-band variability, or short thermal events. Pair with the
  ORCA moorings where a site is near one.
- **Coverage is uneven.** DNA001 has 7 profiles in the 2025 file; some
  parameters are fill-valued in a given cast (the DNA001 cast inspected had
  temperature and salinity but no DO or fluorescence).
- **Deep-water stations.** The stations are mid-channel, not intertidal. They
  characterize the basin the oysters sit in, not the tideflat.
- **Coos Bay is out of scope** — Ecology is Washington only. That site needs an
  Oregon source (see below).

---

## 3. Recommended next steps

1. **Add a step 05 that pulls both sources.** `code/05_orca_ecology_data.py`:
   ORCA gridded profiles from NANOOS ERDDAP for the moorings within range, plus
   Ecology yearly netCDFs subset to the nearest station per site. Restrict to a
   shallow depth band (0–5 m) matched to oyster habitat, filter on QC flags and
   the `-99999.9` fill, and emit per-site monthly climatologies of temperature,
   salinity, dissolved oxygen, and chlorophyll.
2. **Decide the window.** The samples are 2018 collections (plus undated
   `LS`/`MB`/`WB`/Squaxin/Ostrich/Dogfish/Coos/Fidalgo sets). Both sources cover
   2018, so a climatology for the years preceding collection is available — this
   is more defensible than the current 30-day snapshot, which describes 2026
   conditions for animals sampled in 2018.
3. **Resolve the four uncertain sites first.** `CS18_22_Wild_plate1`, `LS`,
   `MB`, and `WB` drive station assignment. Confirming them against the original
   collection records changes which Ecology station each one gets.
4. **Send the NANOOS/NWEM message** once the contact is confirmed.
5. **Find an Oregon source for Coos Bay.** Candidates to check, none verified
   yet: OSU's Charleston/OIMB shore station, the South Slough NERR system-wide
   monitoring program (NERR SWMP data are public and Coos Bay's South Slough is
   a NERR site — likely the best match), and Oregon DEQ ambient monitoring.

## 4. Note on the step-04 outputs

The NOAA-based tables in `output/04_environmental_data/` remain valid as a
basin-scale reference and as a station crosswalk. They should not be the
environmental layer of the final analysis — Ecology plus ORCA should be, with
NOAA retained only where neither has a nearby station.
