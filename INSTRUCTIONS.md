# Project conventions

How the scripts in this repository are organized, run, and extended. Written
for contributors and for automated agents (Copilot, Claude, and similar) working
in the repo. For the scientific context see `README.md`; for per-step inputs,
outputs, and flags see `code/README.md`.

## Layout

```
oly-lc-WGS/
├── code/      analysis scripts, numbered in execution order
├── data/      raw FASTQs and the reference genome (read-only, not tracked)
├── docs/      plans and notes that are not tied to a single step
└── output/    one subdirectory per script, named to match
```

| Script | Writes to | Depends on |
| --- | --- | --- |
| `code/01_align_and_visualize.py` | `output/01_align_and_visualize/` | `data/raw/`, `data/genome/` |
| `code/02_bam_summary.py` | `output/02_bam_summary/` | step 01 BAMs, sample sheet, reference `.fai` |
| `code/03_variant_summary.py` | `output/03_variant_summary/` | step 01 filtered VCF, PLINK dataset, sample sheet |
| `code/04_environmental_data.py` | `output/04_environmental_data/` | step 01 sample sheet; live NOAA endpoints |

Steps 02, 03, and 04 are independent of one another and only need step 01's
outputs. Step 04 needs network access and no bioinformatics tools.

## Running

- Run from the repository root: `python code/NN_name.py [flags]`. Steps 01 to
  03 change directory to the repo root themselves; step 04 resolves paths
  relative to the current directory, so it must be launched from the root.
- Python 3.8 or newer with `numpy`, `pandas`, and `matplotlib`. Step 02 also
  needs `pysam`.
- External tools on `PATH`: `bwa`, `samtools`, `bcftools`, `plink` (step 01);
  `samtools` (step 02); `bcftools`, `plink` (step 03). Each script checks for
  its tools at startup and exits with the missing names.
- Thread counts are capped at 50 via `--threads`. The production run used
  `--threads 32 --threads-per-sample 4` on the UW Hyak cluster.

## Inputs and outputs

- Never write to `data/`. Step 01 symlinks the genome into
  `output/01_align_and_visualize/reference/` and builds indices there.
- Read inputs only from `data/` or from an earlier step's `output/` directory.
  Write only to your own step's directory.
- Use paths relative to the repository root in code and in committed tables.
  Absolute paths belong only in logs.
- Sample identity flows from file names. Step 01 derives `sample_id` and
  `location` from FASTQ prefixes (see `data/README.md`) and writes
  `output/01_align_and_visualize/metrics/sample_metadata.tsv`. Every later
  step reads locations from that file rather than re-deriving them.
- Locations beginning with `Blank` are negative controls. They are kept unless
  `--skip-blanks` is passed to step 01; downstream summaries should exclude or
  flag them.

## Metadata and logging

Every step writes `metadata.json` in its output directory with at least:
script name, ISO-8601 UTC timestamp, runtime in seconds, all parameters, and
the paths of the files it produced (relative to the repo root). Step 04 also
records source URLs and library versions; new steps should do the same.

Steps 01 to 03 log to `logs/pipeline.log` and to stdout through the `logging`
module. Both `metadata.json` and the log are overwritten on each run.

## Reruns and idempotence

- Scripts skip any expensive product that already exists (BAMs, VCFs, PLINK
  files, `samtools stats` reports) and pass `--force` to recompute.
- Random sampling must be seeded. Step 02 exposes `--seed` (default 42) and
  records it in `metadata.json`.
- A missing input is an error: raise `FileNotFoundError` naming the path, do
  not create placeholder outputs, and say which step must run first. Step 04
  prints the missing path and exits with status 1.
- Do not write versioned or timestamped copies. Rerunning a step replaces its
  outputs in place; the git history is the record of earlier results.

## What to commit

- Commit scripts, documentation, small tables, figures, logs, and
  `metadata.json`.
- Do not commit BAMs, VCFs, PLINK binaries, reference indices, `tmp/`
  directories, or `__pycache__/`. `output/01_align_and_visualize/.gitignore`
  and the root `.gitignore` already exclude these.
- Update `code/README.md` when a step's inputs, outputs, or flags change, and
  `README.md` when committed results change.

## Adding a step

1. Number it after the last existing script (`code/05_<name>.py`) and create
   `output/05_<name>/`.
2. Take inputs from `data/` or earlier `output/` directories via `argparse`
   flags with the repo-relative defaults shown in the existing scripts.
3. Reuse the patterns already in place: dependency check, `logging` to file and
   stdout, skip-if-exists with `--force`, `metadata.json` on completion.
4. Add `from __future__ import annotations` so `X | Y` type hints run on
   Python 3.8 and 3.9.
5. Document the step in `code/README.md` and add a row to the workflow table in
   `README.md`.

## Naming

- Scripts: `NN_snake_case.py`.
- Output files: lowercase, hyphen- or underscore-separated, no spaces. The
  existing steps mix `snake_case` (01 to 03) and `kebab-case` (04); either is
  acceptable, but be consistent within a step.
