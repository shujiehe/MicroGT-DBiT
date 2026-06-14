# kraken2

This module performs taxonomic classification of sequencing reads using [Kraken2](https://ccb.jhu.edu/software/kraken2/). It includes host depletion prior to classification to remove contaminating host reads.

## Pipeline Overview

```
Raw R1 + R2 FASTQ
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

### `run_kraken2_RNA.sbatch`
SLURM batch script for taxonomic classification of **RNA** reads. Follows the same host-depletion and classification steps as the DNA script.

**Submit jobs:**
```bash
sbatch run_kraken2_DNA.sbatch
sbatch run_kraken2_RNA.sbatch
```

Edit the scripts to set paths for your host genome, Kraken2 database, and input FASTQ files before submitting.

## Input Files

| File | Description |
|------|-------------|
| `*_R1.fastq.gz` / `*_R2.fastq.gz` | Raw paired-end sequencing reads |
| Host genome FASTA | Reference genome for host depletion (e.g., mouse or human) |
| BWA-MEM2 index | Pre-built index of the host genome |
| Kraken2 database | Standard or custom Kraken2 database directory |

## Output Files

| File | Description |
|------|-------------|
| `std.report` | Kraken2 report: per-taxon read counts and percentages |
| `classified.fastq.gz` | Reads that were classified by Kraken2 |
| `unclassified.fastq.gz` | Reads that could not be classified |
| Host-depleted FASTQ | Intermediate non-host reads passed to Kraken2 |

## Dependencies

- BWA-MEM2
- samtools
- Kraken2 (with a built database, e.g., standard)
- SLURM workload manager (for `.sbatch` submission)

## Notes

- The non-host reads produced here are also the input to `DNA_analysis/` for spatial mapping.
- Adjust `--confidence` and `--minimum-hit-groups` thresholds depending on your reference database and desired specificity.
