# data/

Read-only inputs for the Olympia oyster (*Ostrea lurida*) low-coverage WGS
workflow. Nothing in this directory is modified by the scripts in `code/`, and
the large files are not tracked in git (see the root `.gitignore`). Only this
README is committed.

## Layout

| Path | Contents | Tracked |
| --- | --- | --- |
| `raw/` | Paired-end Illumina FASTQs, one R1/R2 pair per sample | No |
| `genome/Olurida_v081.fa` | Olympia oyster reference assembly (v081, ~1.14 Gb, highly fragmented: contig IDs run past 680,000) | No |

The full working tree, with `raw/`, `genome/`, and all gitignored outputs
(BAMs, VCFs, reference indices), is archived on the Roberts Lab web server:
<https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/oly-lc-WGS/>
(`bu.sh` is the rsync script that pushes it there). The November 2025 runs
were executed on UW Hyak under `/mmfs1/gscratch/scrubbed/sr320/github/oly-lc-WGS/`,
a scrubbed filesystem, so that copy may no longer exist.

The current working copy is on the lab server raven at
`/home/shared/8TB_HDD_03/sr320/oly-lc-WGS/`, restored from the gannet archive in
September 2026, with the tools in the conda environment `oly-lc-wgs`
(samtools, bcftools, PLINK, ANGSD, PCAngsd, pysam). Activate it with
`conda activate oly-lc-wgs` before running any step there.

## FASTQ naming

`code/01_align_and_visualize.py` discovers samples from file names, so the
naming convention matters:

```
<sample_id>_S<n>_L004_R1_001.fastq.gz
<sample_id>_S<n>_L004_R2_001.fastq.gz
```

- `sample_id` is everything before the `_S<n>_` sequencer index.
- `location` is `sample_id` with its final `_`-separated token removed.
- Any location beginning with `Blank` is treated as a negative control and can be
  excluded with `--skip-blanks`.

| File | `sample_id` | `location` |
| --- | --- | --- |
| `CS18_22_Wild_plate1_A7_S47_L004_R1_001.fastq.gz` | `CS18_22_Wild_plate1_A7` | `CS18_22_Wild_plate1` |
| `LS_10a_S12_L004_R1_001.fastq.gz` | `LS_10a` | `LS` |
| `Blank_plate1_S92_L004_R1_001.fastq.gz` | `Blank_plate1` | `Blank` |

The current dataset holds 112 samples across 16 inferred locations (15 sites
plus blanks). The derived sample sheet is written to
`output/01_align_and_visualize/metrics/sample_metadata.tsv` and reused by every
later step. Putative site names and their confidence are documented in the
root `README.md`; `LS`, `MB`, `WB`, and `CS18_22_Wild_plate1` are not yet
confirmed against collection records.

## Reference indices

Step 01 does not write indices next to the FASTA here. It symlinks the genome
into `output/01_align_and_visualize/reference/` and builds the BWA, `faidx`, and
`dict` indices there, keeping this directory untouched.
