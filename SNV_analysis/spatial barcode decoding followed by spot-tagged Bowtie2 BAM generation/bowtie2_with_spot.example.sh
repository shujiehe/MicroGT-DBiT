#!/usr/bin/env bash
set -Eeuo pipefail
export LC_ALL=C

# Input files
R1="/path/to/non_host_R1.fq.gz"
R2="/path/to/non_host_R2.fq.gz"
IDX_PREFIX="/path/to/bowtie2_index/bacteria"

BARCODE_MAP="/path/to/spatial_barcodes.csv"
POSITIONS="/path/to/position.txt"

DECODE_PY="scripts/decode_r2_barcode_to_spot.py"
TAG_PY="scripts/tag_sam_with_spot.py"

# Output directory
OUTDIR="./bowtie2_with_spot"
mkdir -p "$OUTDIR"

# Parameters
T=24
S=4
MAPQ_MIN=0
EXCL_FLAGS=$((0x4 + 0x100 + 0x800))

MAP_TSV="${OUTDIR}/qname_to_spot.tsv.gz"
BAM="${OUTDIR}/bacteria.R1.bowtie2.sorted.bam"

python "$DECODE_PY" \
  --r2 "$R2" \
  --barcode-map "$BARCODE_MAP" \
  --positions "$POSITIONS" \
  --out "$MAP_TSV" \
  --x-slice 32:40 \
  --y-slice 70:78 \
  --max-mismatch 1

bowtie2 \
  -x "$IDX_PREFIX" \
  -U "$R1" \
  -p "$T" \
  --very-sensitive \
  --reorder \
| python "$TAG_PY" \
    --map "$MAP_TSV" \
| samtools view -@ "$S" -b -F "$EXCL_FLAGS" -q "$MAPQ_MIN" - \
| samtools sort -@ "$S" -o "$BAM" -

samtools index "$BAM"
samtools flagstat "$BAM" > "${OUTDIR}/flagstat.txt"
samtools idxstats "$BAM" > "${OUTDIR}/idxstats.txt"