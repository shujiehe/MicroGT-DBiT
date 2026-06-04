#!/usr/bin/env bash
set -Eeuo pipefail
export LC_ALL=C

# Input files
R1="/path/to/bacteria_R1.fq.gz"
IDX_PREFIX="/path/to/bowtie2_index/bacteria"

# Output directory
OUTDIR="./bowtie2_no_spot"
mkdir -p "$OUTDIR"

# Parameters
T=16
S=4
MAPQ_MIN=0
EXCL_FLAGS=$((0x4 + 0x100 + 0x800))

BAM="${OUTDIR}/bacteria.R1.bowtie2.sorted.bam"
LOG="${OUTDIR}/bowtie2.log"

bowtie2 \
  -x "$IDX_PREFIX" \
  -U "$R1" \
  -p "$T" \
  --very-sensitive \
  2> "$LOG" \
| samtools view -@ "$S" -b -F "$EXCL_FLAGS" -q "$MAPQ_MIN" - \
| samtools sort -@ "$S" -o "$BAM" -

samtools index "$BAM"
samtools flagstat "$BAM" > "${OUTDIR}/flagstat.txt"
samtools idxstats "$BAM" > "${OUTDIR}/idxstats.txt"