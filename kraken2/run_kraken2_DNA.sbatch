#!/usr/bin/env bash
#SBATCH --partition=day
#SBATCH --job-name=kraken2_DNA
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 24:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err

set -Eeuo pipefail
export LC_ALL=C

# Paths
R1="/path/to/sample_R1.fastq.gz"
REF="/path/to/host_reference.fa"
KRAKEN_DB_STD="/path/to/kraken2/Standard"

BWA="/path/to/bwa-mem2"
SAMTOOLS="/path/to/samtools"
KRAKEN2="/path/to/kraken2"

# Parameters
T=${SLURM_CPUS_PER_TASK:-32}
CONF=0.05
HITGROUPS=2
OUTDIR="./kraken2_out"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Temporary files
TMPPFX="${SLURM_TMPDIR:-/tmp}/sort.$(hostname).${SLURM_JOB_ID:-$$}.$RANDOM"
trap 'rm -f "$TMPPFX"* non_host_R1.fq.gz 2>/dev/null || true' EXIT

# Host depletion
"$BWA" mem -K 1000000 -t "$T" "$REF" "$R1" \
| "$SAMTOOLS" view -@ $((T/2)) -b -f 4 -F 0x900 - \
| "$SAMTOOLS" sort -n -@ $((T/2)) -T "$TMPPFX" - \
| "$SAMTOOLS" fastq -@ $((T/2)) - \
| gzip -c > non_host_R1.fq.gz

# Kraken2 classification
"$KRAKEN2" \
  --db "$KRAKEN_DB_STD" \
  --threads "$T" \
  --use-names \
  --gzip-compressed \
  --preload \
  --confidence "$CONF" \
  --minimum-hit-groups "$HITGROUPS" \
  --report std.report \
  --output /dev/null \
  non_host_R1.fq.gz

rm -f non_host_R1.fq.gz