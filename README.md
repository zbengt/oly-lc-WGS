# oly-lc-WGS

Low-coverage whole-genome sequencing (lc-WGS) workflow for Olympia oyster samples.
The repository contains Python pipelines for alignment, BAM-based diversity
summaries, and variant-level quality summaries, with outputs written to
step-matched subdirectories under `output/`.

## Repository layout

| Path | Contents |
| --- | --- |
| `code/` | Analysis scripts run in numeric order (`01_` → `03_`) |
| `data/` | Raw reads and reference assets used as read-only inputs |
| `output/` | Generated alignments, tables, figures, logs, and metadata |
| `INSTRUCTIONS.md` | Project execution conventions for agents and contributors |

Additional details for each analysis step are documented in
[`code/README.md`](code/README.md).

## Workflow overview

| Step | Script | Purpose | Key outputs |
| --- | --- | --- | --- |
| 1 | `code/01_align_and_visualize.py` | Align paired-end FASTQs to the Olympia oyster reference, build a joint VCF/PLINK dataset, and visualize genetic connectedness | `output/01_align_and_visualize/alignments/`, `variants/`, `metrics/`, `figures/genetic_connectedness.png` |
| 2 | `code/02_bam_summary.py` | Summarize mismatch, heterozygosity, IBS, and PCA directly from BAMs without rerunning variant calling | `output/02_bam_summary/tables/`, `figures/bam_connectedness.png` |
| 3 | `code/03_variant_summary.py` | Generate VCF/PLINK-based variant quality and diversity summaries | `output/03_variant_summary/` |
| 4 | `code/04_environmental_data.py` | Match each putative sampling site to nearby NOAA buoys/stations and download recent observations | `output/04_environmental_data/` |

## Requirements

Run commands from the repository root and keep `data/` read-only.

- Python 3 with `numpy` and `pandas`
- `pysam` for `code/02_bam_summary.py`
- External tools used by the workflows:
  - step 01: `bwa`, `samtools`, `bcftools`, `plink`
  - step 02: `samtools`
  - step 03: `bcftools`, `plink`
  - step 04: network access to NOAA NDBC and NOAA CO-OPS endpoints

## Typical usage

```bash
python code/01_align_and_visualize.py --threads 32 --threads-per-sample 4
python code/02_bam_summary.py --num-sites 500
python code/03_variant_summary.py --threads 32
python code/04_environmental_data.py --days 30 --radius-km 75
```

All outputs are written with relative paths so results remain reproducible across
systems.

## Results to date

The checked-in outputs currently document completed work for the alignment and
BAM-summary stages.

### Current dataset snapshot

| Metric | Value |
| --- | --- |
| Samples discovered in `sample_metadata.tsv` | 112 |
| Inferred locations/categories | 16 total |
| Blank/control samples | 2 |
| Biological locations | 15 |
| Sorted BAMs produced | 112 |
| Joint variant files produced by step 01 | `all_samples.raw.vcf.gz`, `all_samples.filtered.vcf.gz` |

### Sample counts by inferred location

| Location | Putative site (complete words) | Samples |
| --- | --- | ---: |
| Blank | Negative control (no oyster tissue) | 2 |
| CS18_22_Wild_plate1 | Central Sound wild collection, 2018 (Clam Bay / Manchester vicinity) — uncertain | 7 |
| Coos_Bay | Coos Bay, Oregon (outside Puget Sound) | 7 |
| Dogfish_Bay | Dogfish Bay, Liberty Bay vicinity, Kitsap Peninsula, Puget Sound | 8 |
| FB18_Wild | Fidalgo Bay wild collection, 2018, Anacortes, northern Puget Sound | 5 |
| Fidalgo_Bay | Fidalgo Bay, Anacortes, northern Puget Sound | 5 |
| HC18_Triton_Wild | Triton Cove, Hood Canal, 2018 wild collection | 8 |
| LS | Little Skookum Inlet, southern Puget Sound — uncertain | 7 |
| MB | Mud Bay, Eld Inlet, southern Puget Sound — uncertain | 8 |
| NS18_Disco_Wild | Discovery Bay, north Olympic Peninsula, 2018 wild collection | 8 |
| NS18_Sequim_Wild | Sequim Bay, north Olympic Peninsula, 2018 wild collection | 8 |
| Ostrich_Bay | Ostrich Bay, Dyes Inlet, Bremerton, central Puget Sound | 8 |
| PGB18_Wild | Port Gamble Bay, 2018 wild collection, northern Hood Canal | 8 |
| SS18_North_Bay_Wild | North Bay, Case Inlet, southern Puget Sound, 2018 wild collection | 8 |
| Squaxin_Island | Squaxin Island, southern Puget Sound | 7 |
| WB | Westcott Bay, San Juan Island — uncertain | 8 |

Site names in the second column are inferred from the sample-name prefixes, not
from a curated collection sheet. Rows flagged *uncertain* have prefixes that do
not map unambiguously to a single site and should be confirmed against the
original collection records before use.

### Summary metrics from committed outputs

These values come from `output/01_align_and_visualize/metrics/coverage_summary.tsv`
and `output/02_bam_summary/tables/`.

| Summary | Current value |
| --- | --- |
| Mean genome coverage across non-blank samples | 0.480 covered fraction (48.0%) |
| Mean genome-wide depth across non-blank samples | 5.44× |
| Mean location-level BAM heterozygosity | 0.00824 |
| Mean location-level BAM mismatch rate | 0.00530 |
| Highest location-level heterozygosity | `FB18_Wild` (0.01457) |
| Highest location-level mean depth at sampled BAM sites | `Ostrich_Bay` (14.91) |
| Highest location-level mismatch rate | `HC18_Triton_Wild` (0.00979) |
| Highest per-sample heterozygosity in BAM summary | `FB18_Wild_11` (0.02283) |

### Environmental context

`code/04_environmental_data.py` pairs each putative site with the buoys and
shore stations around it and pulls recent observations from NOAA NDBC and NOAA
CO-OPS. Results, including per-site water-temperature summaries and the full
station crosswalk, are in
[`output/04_environmental_data/`](output/04_environmental_data/README.md).

No NOAA in-water sensor sits inside the small inlets these oysters come from, so
the closest reporting station is typically a basin away (Tacoma for the South
Sound sites, Port Townsend for the north Sound and Hood Canal sites). The
NANOOS/UW ORCA moorings are much closer for Hood Canal and South Sound and are
listed in the station crosswalk, but their data are not downloadable from a
public bulk endpoint and so are recorded as pointers only.

### Generated figures

Step 01 connectedness figure:

![Genetic connectedness across samples](output/01_align_and_visualize/figures/genetic_connectedness.png)

Step 02 BAM-based connectedness figure:

![BAM-based connectedness across samples](output/02_bam_summary/figures/bam_connectedness.png)

### Current status of later-stage summaries

The variant-summary script is present in `code/03_variant_summary.py`, but the
checked-in repository outputs currently center on steps 01 and 02; no committed
`output/03_variant_summary/` results are present yet.
