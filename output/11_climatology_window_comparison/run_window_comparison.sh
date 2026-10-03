#!/bin/bash
#SBATCH --job-name=oly-env-window
#SBATCH --account=coenv
#SBATCH --partition=compute
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/mmfs1/gscratch/scrubbed/sr320/github/oly-lc-WGS/output/11_climatology_window_comparison/logs/job_%j.out
#SBATCH --chdir=/mmfs1/gscratch/scrubbed/sr320/github/oly-lc-WGS
# Steps 07 and 08 for 2021-2024 with redirected output, then the comparison.
set -euo pipefail
py=/mmfs1/gscratch/srlab/sr320/miniforge3/envs/angsd/bin/python
echo "start $(date -u +%FT%TZ) host $(hostname)"
${py} code/07_orca_ecology_data.py --start-year 2021 --end-year 2024 --output-dir output/07_orca_ecology_data_2021-2024
${py} code/08_environmental_predictors.py --step07-dir output/07_orca_ecology_data_2021-2024 --output-dir output/08_environmental_predictors_2021-2024
${py} code/11_climatology_window_comparison.py
echo "end $(date -u +%FT%TZ)"
