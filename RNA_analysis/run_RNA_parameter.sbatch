#!/usr/bin/env bash
#SBATCH --job-name=run_RNA
#SBATCH --partition=day
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=08:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

set -Eeuo pipefail
export LC_ALL=C

# ===== user settings =====
JSON="/path/to/parameter.json"

CONDA_SH="/path/to/miniforge3/etc/profile.d/conda.sh"
CONDA_ENV="astro"

# ===== run =====
source "$CONDA_SH"
conda activate "$CONDA_ENV"

CPU="${SLURM_CPUS_PER_TASK:-32}"
export OMP_NUM_THREADS="$CPU"
export PYTHONNOUSERSITE=1

ASTRO "$JSON"