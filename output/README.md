# output/

Generated results, one subdirectory per script in `code/`, named to match:
`code/NN_name.py` writes to `output/NN_name/`. Each step directory carries a
`metadata.json` (parameters, inputs, runtime, produced files) and a
`logs/pipeline.log`. See `code/README.md` for what each step produces and how
to run it.

## What is tracked

Small, human-readable results are committed; large binary intermediates are
not. `output/01_align_and_visualize/.gitignore` excludes `reference/`,
`alignments/`, and `variants/`. `tmp/` directories are scratch space and should
not be committed.

| Directory | Committed | Not committed |
| --- | --- | --- |
| `01_align_and_visualize/` | `metrics/` (sample sheet, coverage summary), `figures/`, `logs/`, `metadata.json` | `reference/` (symlink + indices), `alignments/` (112 sorted BAMs + `.bai`), `variants/` (joint VCFs, PLINK dataset, PCA, IBS) |
| `02_bam_summary/` | `stats/samtools_stats/` (per-sample `samtools stats`), `tables/`, `positions/`, `figures/`, `logs/`, `metadata.json` | `tmp/` |
| `03_variant_summary/` | Nothing yet. The script exists but has not completed a run against the joint VCF | |
| `04_environmental_data/` | `site-coordinates.tsv`, `stations/`, `observations/`, `environmental-summary.tsv`, `metadata.json`, `README.md` | |

## Reading the current results

- Coverage per sample: `01_align_and_visualize/metrics/coverage_summary.tsv`.
  Non-blank samples average about 5.4x depth. The `bam` column holds absolute
  cluster paths from the run that produced it.
- BAM-based diversity: `02_bam_summary/tables/`. These rest on 500 randomly
  sampled reference positions, of which only a handful are polymorphic, so the
  heterozygosity, mismatch, IBS, and PCA values are provisional. Treat them as a
  QC pass, not a population-genetic result.
- Environmental context: `04_environmental_data/README.md` explains the station
  matching and its limits. `docs/environmental-data-access-plan.md` records the
  closer NANOOS ORCA and WA Ecology sources planned for a later step.

## Rerunning

Scripts skip work whose outputs already exist and accept `--force` to redo it.
`metadata.json` and `logs/pipeline.log` are overwritten on every run, so commit
or copy them before rerunning if the previous record matters.
