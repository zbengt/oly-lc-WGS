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
  - Each BAM carries an `@RG` line with `SM:<sample_id>`, so bcftools and PLINK name
    samples by ID. BAMs aligned before this was added lack the tag and are named by
    file path in the VCF; the script warns about them and maps the names back for
    plotting, but re-align with `--force` (or `samtools addreplacerg`) to fix the VCF.
  - The IBS heatmap plots PLINK's `.mibs` values directly (proportion of alleles
    shared); an earlier version inverted them.
  - `--fix-read-groups` adds the missing `@RG` line in place to existing BAMs with
    `samtools addreplacerg`, in parallel, before variant calling.
  - If the filtered VCF carries file paths as sample names (from BAMs aligned
    without read groups), the script rewrites its header to sample IDs with
    `bcftools reheader` and regenerates the PLINK outputs and figure.
  - `--force-variants` re-runs mpileup/call/filter and everything downstream;
    `--force-plink` re-runs only PLINK and the figure from the existing filtered
    VCF. Plain `--force` also re-aligns and re-indexes the reference.
  - With `--skip-blanks`, blanks already present in an existing VCF are dropped
    from the PLINK PCA/IBS via `--keep`.
  - `--min-mean-depth` (default 1.0x) drops samples below that genome-wide depth in
    `metrics/coverage_summary.tsv` from the PCA/IBS; the excluded IDs are recorded in
    `metadata.json`. `HC18_Triton_Wild_10` (0.26x) is the only current sample affected.
  - **PLINK contig subset.** PLINK 1.9 accepts at most 65,280 distinct
    chromosome codes even with `--allow-extra-chr`, and `Olurida_v081` has
    159,429 contigs (158,839 carry filtered variants), so the VCF is subset to
    long contigs before conversion. Contig lengths are read from the reference
    index `output/01_align_and_visualize/reference/Olurida_v081.fa.fai` (not the
    VCF header); contigs with length `>= --plink-min-contig-length` are kept,
    sorted longest-first, and truncated to `--plink-max-contigs`. The selection
    is written to `variants/plink_contigs.txt` as a three-column
    `bcftools -R` regions file (`contig<TAB>1<TAB>length`), the count is logged
    with the genome fraction retained, and `metadata.json` records the
    thresholds plus `plink_contigs_selected`. If the list on disk differs from
    the current selection the subset VCF is rebuilt, so changing the flags or
    the reference never silently reuses a stale `subset_for_plink.vcf.gz`.
  - **Choosing the threshold.** The assembly is fragmented (N50 12.9 kb).
    Retained by `--plink-min-contig-length` for `Olurida_v081`:

    | min length | contigs | bp | genome | filtered variants |
    |---:|---:|---:|---:|---:|
    | 100 kb (old default) | 32 | 3.8 Mb | 0.33 % | 160 k |
    | 50 kb | 894 | 57.7 Mb | 5.1 % | 2.6 M |
    | 20 kb (default) | 11,501 | 360 Mb | 31.6 % | 16.2 M |
    | 10 kb | 35,320 | 692 Mb | 60.6 % | 31.4 M |
    | none | 159,429 | 1.14 Gb | 100 % | 56.1 M (exceeds PLINK limit) |

    The default is 20 kb: it samples roughly a third of the genome across
    ~11.5k independent contigs, stays well under the PLINK cap, and keeps the
    subset VCF and PCA tractable. 10 kb is also valid (35k contigs) if more
    sites are wanted; anything below ~7 kb exceeds the cap. Contigs are in any
    case only a proxy for independence; ANGSD/PCAngsd on the BAMs would avoid
    both the hard-call genotypes and the PLINK contig limit.
  - Note: the committed `figures/genetic_connectedness.png` was produced before
    this selection logic was fixed and rests on a single contig (`Contig19646`,
    7,865 SNPs); see the note under *Generated figures* in the root `README.md`.
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
  - Site coordinates are approximate centroids for each site, not recorded collection points.
  - A station may advertise water temperature but return nothing for the window, so the script walks outwards until one delivers data.
  - NANOOS/UW ORCA moorings are listed as pointers only; their data are openly served from the NANOOS ERDDAP and ingesting them is planned separately (see `docs/environmental-data-access-plan.md`).

## 05_realign_xbOstLuri2.Rmd

Re-alignment of every sample to the chromosome-level NCBI reference
GCA_061535525.1 (xbOstLuri2, USDA-ARS / PSRF, released 2026-09-21), replacing the
Olurida_v081 alignments from step 01 for all downstream genetics.

- **Inputs**
  - Paired-end FASTQs in `data/raw/` (restored from the gannet backup by the `restore-fastq` chunk if missing)
  - `data/genome/GCA_061535525.1_xbOstLuri2_USDA-ARS_PSRF_primary_genomic.fna` (from NCBI Datasets)
  - Sample sheet from step 01 (`output/01_align_and_visualize/metrics/sample_metadata.tsv`)
- **Execution**
  - Open in RStudio on Hyak and run the chunks in order, or `Rscript -e 'rmarkdown::render("code/05_realign_xbOstLuri2.Rmd")'`.
  - `restore-fastq`, `index-genome` and `align` each submit a SLURM job (`coenv` / `cpu-g2`) and return; `align`
    is a 112-task array, 24 at a time, chained to the other two with `--dependency=afterok`. The `status` chunk
    reports progress; run the summary and metadata chunks after every sample has a BAM.
  - Tools come from `/mmfs1/gscratch/srlab/containers/srlab-R4.4-bioinformatics-container-c3d3116.sif`
    (bwa 0.7.19, samtools 1.20); nothing needs to be installed.
- **Outputs**
  - `output/05_realign_xbOstLuri2/alignments/` sorted, duplicate-marked BAM + index per sample (gitignored)
  - `output/05_realign_xbOstLuri2/reference/` genome symlink plus BWA, faidx and dict indices (gitignored)
  - `output/05_realign_xbOstLuri2/metrics/` per-sample `flagstat`, `coverage` and `markdup` reports,
    `coverage_summary.tsv`, `coverage_by_location.tsv`, `comparison_v081.tsv`
  - `output/05_realign_xbOstLuri2/figures/depth_v081_vs_xbOstLuri2.png`
  - `output/05_realign_xbOstLuri2/metadata.json` and `logs/`
- **Notes**
  - Each BAM carries `@RG` with `ID`, `SM`, `LB` and `PL` set to the sample ID.
  - Duplicates are marked (flag 1024) with `samtools markdup`, not removed; step 01 never marked them.
  - A sample is skipped when its `metrics/per_sample/<sample>.coverage.tsv` is newer than its BAM; BAMs are
    written as `.part` and renamed on success, so interrupted tasks are redone on the next submission.

## 06_angsd_structure.Rmd

Population structure, admixture, diversity and pairwise Fst from genotype
likelihoods on the step 05 xbOstLuri2 BAMs, replacing the hard-call PCA of
step 01 (7,865 SNPs on one contig) with a genome-wide, low-coverage-appropriate
analysis.

- **Inputs**
  - Sorted, duplicate-marked BAMs from step 05 (`output/05_realign_xbOstLuri2/alignments/`)
  - Step 05 reference FASTA and `.fai`; chromosome list from `data/genome/GCA_061535525.1_sequence_report.jsonl`
  - `output/05_realign_xbOstLuri2/metrics/coverage_summary.tsv` for sample selection (blanks and samples
    under 2x mean depth are excluded: currently the two blanks and `HC18_Triton_Wild_10`, leaving 109)
- **Execution**
  - Run the chunks in order from RStudio, or `Rscript -e 'rmarkdown::render("code/06_angsd_structure.Rmd")'`.
  - `angsd-gl` (44 windows of 25 Mb over the 10 chromosomes), `pcangsd`, `saf` (15 populations) and `fst`
    (105 pairs) each submit SLURM jobs on `coenv` / `cpu-g2` and return; `pcangsd` and `fst` chain to their
    prerequisites with `--dependency=afterok`. `status` reports progress; run the R chunks when done.
  - Tools come from the conda env `/mmfs1/gscratch/srlab/sr320/miniforge3/envs/angsd` (ANGSD 0.940,
    PCAngsd 1.36.4); rebuild with `mamba create -n angsd -c conda-forge -c bioconda angsd pcangsd`.
- **Outputs**
  - `output/06_angsd_structure/samples/` sample table, BAM lists per population, chromosome list, windows, pairs, GL parameters
  - `output/06_angsd_structure/gl/` per-window and merged Beagle and `mafs` files (gitignored)
  - `output/06_angsd_structure/pca/` PCAngsd covariance matrix, site mask, admixture Q and P for K = 2 to 6
  - `output/06_angsd_structure/saf/` per-population SAF (gitignored); `sfs/` folded SFS and `thetaStat` summaries
  - `output/06_angsd_structure/fst/` per-pair global Fst (`*.fst.txt`)
  - `output/06_angsd_structure/tables/` PCA scores and variance, admixture per K, diversity by population,
    Fst pairs and matrix, SNPs per chromosome
  - `output/06_angsd_structure/figures/` `pca.png`, `admixture.png`, `fst_heatmap.png`
  - `output/06_angsd_structure/metadata.json` and `logs/`
- **Notes**
  - Filters for every ANGSD run: MAPQ >= 20, base quality >= 20, unique proper pairs, `-remove_bads 1`
    (drops the duplicates marked in step 05), BAQ with `-C 50`, site covered in >= 80% of samples, total
    depth between a third and twice the summed mean depth of the samples in the run.
  - Genotype likelihoods use `-GL 1` with `-SNP_pval 1e-6` and minor allele frequency >= 0.05; SFS use all
    sites and are folded (the reference is not an ancestral sequence).
  - Only the 10 chromosomes (995 of 1,029 Mb) are analysed; unplaced scaffolds are skipped.
