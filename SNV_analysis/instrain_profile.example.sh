#!/usr/bin/env bash
set -Eeuo pipefail
export LC_ALL=C

# Input files
BAM="/path/to/bacteria.R1.bowtie2.sorted.bam"
REF_FA="/path/to/bacterial_reference.fa"
STB="/path/to/bacterial_reference.stb"

# Output directory
OUTDIR="./instrain_profile"
mkdir -p "$OUTDIR"

# Parameters
T=16

MIN_MAPQ=40
MIN_READ_ANI=0.92
MIN_COV=5
MIN_FREQ=0.05
FDR=1e-6

inStrain profile \
  "$BAM" "$REF_FA" \
  -o "$OUTDIR" \
  -p "$T" \
  --pairing_filter all_reads \
  --stb "$STB" \
  --min_mapq "$MIN_MAPQ" \
  --min_read_ani "$MIN_READ_ANI" \
  --min_cov "$MIN_COV" \
  --min_freq "$MIN_FREQ" \
  --fdr "$FDR" \
  --skip_mm_profiling