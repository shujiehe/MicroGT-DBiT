#!/usr/bin/env bash
#SBATCH --partition=day
#SBATCH --job-name=kraken2_RNA
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 24:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err

set -Eeuo pipefail
export LC_ALL=C

# Paths
R1="/path/to/sample_R1.trim.fastq.gz"
STAR_INDEX="/path/to/STAR_index"
KRAKEN_DB_STD="/path/to/kraken2/Standard"

STAR_BIN="/path/to/STAR"
KRAKEN2="/path/to/kraken2"

# Parameters
T=${SLURM_CPUS_PER_TASK:-32}
CONF=0.05
HITGROUPS=2
OUTDIR="./kraken2_out"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Host depletion
TMPBASE="${SLURM_TMPDIR:-/tmp}"
OUTTMP="$(mktemp -u -p "$TMPBASE" star_tmp.${SLURM_JOB_ID:-$$}.XXXXXX)"
trap '[[ -d "$OUTTMP" ]] && rm -rf "$OUTTMP"' EXIT

"$STAR_BIN" \
  --runThreadN "$T" \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$R1" \
  --readFilesCommand zcat \
  --outSAMtype None \
  --outReadsUnmapped Fastx \
  --outTmpDir "$OUTTMP" \
  --genomeLoad NoSharedMemory

gzip -f Unmapped.out.mate1
mv Unmapped.out.mate1.gz non_host_R1.fq.gz

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