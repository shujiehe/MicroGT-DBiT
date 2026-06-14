# RNA_analysis

This module processes spatially-barcoded RNA sequencing data from MicroGT-DBiT experiments. It decodes dual spatial barcodes, extracts UMIs, aligns transcripts to a combined host-plus-bacterial genome, and generates a UMI-deduplicated gene-by-spot expression matrix.

## Pipeline Overview

```
R1 (transcript reads) + R2 (spatial barcodes + UMI)
        │
        ▼
ASTRO (configured via RNA_parameter.json)
   ├── Barcode decoding  (8bp + linker + linker + 8bp structure)
   ├── UMI extraction    (24 bp fixed anchor + 10 bp random UMI)
   ├── STAR alignment    (host_plus_bac GTF)
   └── UMI deduplication
        │
        ▼
Gene × Spot expression matrix
```

RNA processing is handled by the external **ASTRO** tool. This folder contains the configuration and SLURM submission script needed to run it.

## Files

### `RNA_parameter.json`
JSON configuration file that defines all parameters for the ASTRO RNA pipeline.

Key fields:

| Parameter | Description |
|-----------|-------------|
| `R1` | Path to transcript read FASTQ (R1) |
| `R2` | Path to barcode/UMI read FASTQ (R2) |
| `spatial_barcodes` | Path to `spatial_barcodes.txt` whitelist |
| `UMI_structure` | `CAAGCGTTGGCTTCTCGCATCT_10` — 24 bp fixed anchor + 10 bp random UMI |
| `barcode_structure` | `8bp_LINKER_LINKER_8bp_LINKER` — dual 8 bp barcodes with flanking linkers |
| `GTF` | Combined host + bacterial annotation file (`host_plus_bac.gtf`) |
| `STAR_index` | Path to STAR genome index directory |
| `output` | Output directory for expression matrices and intermediate files |
| `steps` | Processing steps to run (1–7) |

Edit this file to point to your data before running.

### `run_RNA_parameter.sbatch`
SLURM batch script that submits the ASTRO job to an HPC cluster.

Resource requirements:
- **CPUs:** 32
- **Memory:** 180 GB
- **Time limit:** 8 hours
- **Conda environment:** `astro`

**Submit the job:**
```bash
sbatch run_RNA_parameter.sbatch
```

**Or run interactively:**
```bash
conda activate astro
ASTRO RNA_parameter.json
```

## Input Files

| File | Description |
|------|-------------|
| `*_R1.fastq.gz` | Read 1 — RNA transcript reads |
| `*_R2.fastq.gz` | Read 2 — spatial barcodes and UMI |
| `spatial_barcodes.txt` | Barcode whitelist: `BARCODE_B+BARCODE_A  X  Y` |
| `host_plus_bac.gtf` | Combined host + bacterial GTF annotation (in `test_data/built_gtf/`) |
| STAR index | Pre-built STAR genome index |

## Output Files

| File | Description |
|------|-------------|
| `matrix.mtx` / `barcodes.tsv` / `features.tsv` | Sparse gene × spot count matrix (10x-compatible format) |
| STAR alignment BAM | Aligned reads with spatial barcode tags |
| Log files | Per-step processing statistics |

## Dependencies

- **ASTRO** (installed in conda environment `astro`)
- STAR aligner
- SLURM workload manager (for `.sbatch` submission)

## Notes

- The `host_plus_bac.gtf` annotation in `test_data/built_gtf/` combines host genome gene annotations with bacterial gene models, allowing simultaneous quantification of both host and microbial transcripts in a single alignment step.
- UMI deduplication is applied to correct for PCR amplification bias.
