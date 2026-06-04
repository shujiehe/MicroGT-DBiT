#!/usr/bin/env bash
set -Eeuo pipefail
export LC_ALL=C

# Edit the paths below before running.
R1="/path/to/non_host_R1.fq.gz"
R2="/path/to/non_host_R2.fq.gz"
WHITELIST="/path/to/spatial_barcodes.txt"
REF="/path/to/bacteria_all.fa"

LINKER2="ATCCACGTGCTTGAGAGGCCAGAGCATTCG"

DEMUX_PY="scripts/demultiplexer.py"
MAP_PY="scripts/mapping.py"
COUNT_PY="scripts/count_contig.py"

OUTDIR="./spatial_dna_output"
DEMUX_DIR="${OUTDIR}/01_demultiplex"
MAP_DIR="${OUTDIR}/02_mapping"
COUNT_DIR="${OUTDIR}/03_count_contig"

mkdir -p "$DEMUX_DIR" "$MAP_DIR" "$COUNT_DIR"

python "$DEMUX_PY" \
  --R1 "$R1" \
  --R2 "$R2" \
  --whitelist "$WHITELIST" \
  --linker2 "$LINKER2" \
  --outdir "$DEMUX_DIR"

python "$MAP_PY" \
  --ref "$REF" \
  --fq "${DEMUX_DIR}/combine.fq" \
  --outdir "$MAP_DIR" \
  --threads 16 \
  --mapq 10

python "$COUNT_PY" \
  --bam "${MAP_DIR}/align.sorted.bam" \
  --barcodes "$WHITELIST" \
  --outdir "$COUNT_DIR" \
  --threads 16