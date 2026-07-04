# kraken2

This module performs taxonomic classification of sequencing reads using [Kraken2](https://ccb.jhu.edu/software/kraken2/). It includes host depletion prior to classification to remove contaminating host reads.

## DNA Pipeline Overview

```
Raw R1 FASTQ
        │
        ▼
BWA-MEM2 alignment to host genome
        │
        ▼
samtools (keep unmapped reads only)
        │
        ▼
Kraken2 classification
        │
        ▼
std.report  (taxonomic classification report)
```

## Scripts

### `run_kraken2_DNA.sbatch`
SLURM batch script for taxonomic classification of **DNA** reads.

Steps:
1. Aligns reads to the host genome using BWA-MEM2
2. Extracts only unmapped reads (non-host reads) with samtools
3. Runs Kraken2 on the host-depleted reads

Kraken2 parameters:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--confidence` | 0.05 | Minimum confidence score for a classification |
| `--minimum-hit-groups` | 2 | Minimum number of distinct k-mer groups required |
| Database | Standard | Kraken2 standard database |


**Submit jobs:**
```bash
sbatch run_kraken2_DNA.sbatch
```

Edit the scripts to set paths for your host genome, Kraken2 database, and input FASTQ files before submitting.

## Input Files

| File | Description |
|------|-------------|
| `*_R1.fastq.gz` | Read 1 containing genomic sequences from raw paired-end sequencing reads |
| Host genome FASTA | Reference genome for host depletion (e.g., mouse or human) |
| BWA-MEM2 index | Pre-built index of the host genome |
| Kraken2 database | Standard or custom Kraken2 database directory |

## Output Files

| File | Description |
|------|-------------|
| `std.report` | Kraken2 report: per-taxon read counts and percentages |
| Host-depleted FASTQ | Temporary host-depleted R1 FASTQ generated from BWA-MEM2-unmapped reads; removed after Kraken2 classification |

## Dependencies

- BWA-MEM2
- samtools
- Kraken2 (with a built database, e.g., standard)
- SLURM workload manager (for `.sbatch` submission)


## RNA Pipeline Overview

```
Trimmed R1 FASTQ
        │
        ▼
STAR alignment to host genome
        │
        ▼
STAR outputs unmapped reads
        │
        ▼
Kraken2 classification
        │
        ▼
std.report  (taxonomic classification report)
```

## Scripts

### `run_kraken2_RNA.sbatch`
SLURM batch script for taxonomic classification of **RNA** reads.

Steps:
1. Aligns reads to the host genome using STAR
2. Keeps STAR-unmapped reads as host-depleted reads
3. Runs Kraken2 on the host-depleted reads

Kraken2 parameters:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--confidence` | 0.05 | Minimum confidence score for a classification |
| `--minimum-hit-groups` | 2 | Minimum number of distinct k-mer groups required |
| Database | Standard | Kraken2 standard database |


**Submit jobs:**
```bash
sbatch run_kraken2_RNA.sbatch
```

Edit the scripts to set paths for your host genome, Kraken2 database, and input FASTQ files before submitting.

## Input Files

| File | Description |
|------|-------------|
| `*_R1.trim.fastq.gz` | Read 1 containing genomic sequences from raw paired-end sequencing reads |
| Host genome FASTA | Reference genome for host depletion (e.g., mouse or human) |
| STAR index | Pre-built index of the host genome |
| Kraken2 database | Standard or custom Kraken2 database directory |

## Output Files

| File | Description |
|------|-------------|
| `std.report` | Kraken2 report: per-taxon read counts and percentages |
| Host-depleted FASTQ | Temporary host-depleted R1 FASTQ generated from STAR-unmapped reads; removed after Kraken2 classification |

## Dependencies

- STAR
- Kraken2 (with a built database, e.g., standard)
- SLURM workload manager (for `.sbatch` submission)

## Notes

- Adjust `--confidence` and `--minimum-hit-groups` thresholds depending on your reference database and desired specificity.
