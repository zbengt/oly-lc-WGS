# oly-lc-WGS Code Directory

This directory contains scripted workflows executed in ascending numeric order
according to the project conventions outlined in `INSTRUCTIONS.md`.

## 01_align_and_visualize.py

Low-coverage WGS alignment and genetic connectedness summary pipeline.

- **Inputs**
  - Paired-end FASTQ files in `data/raw/`
  - Olympia oyster reference genome `data/genome/Olurida_v081.fa`
- **Execution**
  - Run from the repository root:  
    `python code/01_align_and_visualize.py --threads 32 --threads-per-sample 4`
  - Uses up to 50 CPUs and respects a 200 MB memory ceiling by default.
- **Outputs**
  - `output/01_align_and_visualize/alignments/` sorted BAM + index files per sample
  - `output/01_align_and_visualize/variants/` joint VCFs and PLINK matrices
  - `output/01_align_and_visualize/metrics/` sample sheet and coverage summaries
  - `output/01_align_and_visualize/figures/genetic_connectedness.png`
  - `output/01_align_and_visualize/metadata.json` documenting parameters, runtime, and produced files
- **Notes**
  - Sample locations are inferred from FASTQ prefixes and propagated into all summaries.
  - Blank/control samples can be skipped with `--skip-blanks`.
  - All outputs use relative paths to maintain reproducibility.

## 02_bam_summary.py

BAM-based variation and connectedness overview (avoids joint VCF generation).

- **Inputs**
  - Sorted BAMs from step 01 (`output/01_align_and_visualize/alignments/*.sorted.bam`)
  - Sample metadata table (`output/01_align_and_visualize/metrics/sample_metadata.tsv`)
  - Reference index (`output/01_align_and_visualize/reference/Olurida_v081.fa.fai`)
- **Execution**
  - Run from the repository root:  
    `python code/02_bam_summary.py --num-sites 500`
- **Outputs**
  - `output/02_bam_summary/stats/` per-sample `samtools stats` reports
  - `output/02_bam_summary/tables/` mismatch, heterozygosity, IBS, and PCA summaries
  - `output/02_bam_summary/figures/bam_connectedness.png`
  - `output/02_bam_summary/metadata.json` documenting parameters, runtime, and outputs
- **Notes**
  - Randomly samples reference sites (default 500) and uses BAM base counts to estimate variation.
  - Keeps memory usage low by streaming each BAM sequentially; adjust `--num-sites` as needed.

## 03_variant_summary.py

VCF/PLINK-based variant quality assessment (requires joint VCF from step 01).

- **Inputs**
  - Filtered multi-sample VCF (`output/01_align_and_visualize/variants/all_samples.filtered.vcf.gz`)
  - PLINK dataset from step 01 (`output/01_align_and_visualize/variants/plink_dataset.*`)
  - Sample metadata table (`output/01_align_and_visualize/metrics/sample_metadata.tsv`)
  - Reference genome (`data/genome/Olurida_v081.fa`) for Ts/Tv in `bcftools stats`
- **Execution**
  - Run from the repository root:  
    `python code/03_variant_summary.py --threads 32`
- **Outputs**
  - `output/03_variant_summary/stats/` containing `bcftools` and PLINK statistics
  - `output/03_variant_summary/tables/` per-sample and per-location summaries
  - `output/03_variant_summary/figures/` bar charts for heterozygosity and call rate
  - `output/03_variant_summary/metadata.json` logging parameters, runtime, and outputs
- **Notes**
  - Generates a PLINK `--within` cluster file based on sample locations to stratify allele frequencies.
  - Respects thread limits (≤50) and reuses existing outputs unless `--force` is provided.

## 04_environmental_data.py

Environmental context for each putative sampling site from nearby buoys and shore stations.

- **Inputs**
  - Sample metadata table (`output/01_align_and_visualize/metrics/sample_metadata.tsv`)
  - Site coordinates for the putative locations documented in the top-level `README.md` (defined in the script)
  - Live NOAA endpoints: NDBC station table and `realtime2` feed; CO-OPS metadata and data APIs
- **Execution**
  - Run from the repository root:  
    `python code/04_environmental_data.py --days 30 --radius-km 75`
  - `--no-download` builds only the site/station crosswalk
- **Outputs**
  - `output/04_environmental_data/site-coordinates.tsv` putative site, region, coordinates, confidence, sample count
  - `output/04_environmental_data/stations/` full station catalog and per-site nearby-station matches
  - `output/04_environmental_data/observations/` CO-OPS hourly water temperature and NDBC standard-met records
  - `output/04_environmental_data/environmental-summary.tsv` per-site, per-variable n/mean/min/max
  - `output/04_environmental_data/metadata.json` parameters, sources, runtime, software versions
- **Notes**
  - Site coordinates are approximate centroids for interpreted sites, not recorded collection points; ambiguous locations (`CS18_22_Wild_plate1`, `LS`, `MB`, `WB`) are flagged `uncertain`.
  - A station may advertise water temperature but return nothing for the window, so the script walks outwards until one delivers data.
  - NANOOS/UW ORCA moorings are listed as pointers only; their data are not available from a public bulk endpoint.
