#!/usr/bin/env bash
set -Eeuo pipefail
export LC_ALL=C
export PYTHONUNBUFFERED=1


PIPELINE_PY="scripts/run_sgvfinder2_pipeline.py"

DB_PREFIX="/path/to/sgvfinder2_database_prefix"
FASTQ_DIR="/path/to/fastq_dir"

OUT_ROOT="./sgvfinder2_output"
ICRA_ROOT="${OUT_ROOT}/icra"
SMP_DIR="${OUT_ROOT}/smp"
FINAL_DIR="${OUT_ROOT}/final"

MIN_SAMP_CUTOFF=2

mkdir -p "$ICRA_ROOT" "$SMP_DIR" "$FINAL_DIR"

python "$PIPELINE_PY" \
  --db_prefix "$DB_PREFIX" \
  --fastq_dir "$FASTQ_DIR" \
  --icra_root "$ICRA_ROOT" \
  --smp_dir "$SMP_DIR" \
  --final_dir "$FINAL_DIR" \
  --min_samp_cutoff "$MIN_SAMP_CUTOFF" \
  --write_csv