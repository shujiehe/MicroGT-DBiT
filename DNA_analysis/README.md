# DNA_analysis

This module processes spatially-barcoded DNA sequencing data from MicroGT-DBiT experiments. It decodes spatial barcodes from R2 reads, maps the resulting reads to a bacterial reference genome, and produces a contig-by-spot abundance matrix.

## Pipeline Overview

```
R1 + R2 FASTQ (non-host reads)
        │
        ▼
demultiplexer.py   ──►  combine.fq   (reads tagged with spatial coordinates)
        │
        ▼
mapping.py         ──►  align.sorted.bam
        │
        ▼
count_contig.py    ──►  expmat_contig.tsv  (contig × spot matrix)
```

The full pipeline can be run end-to-end using `run_spatial_dna_pipeline.sh`.

## Scripts

### `demultiplexer.py`
Extracts spatial barcodes from R2 reads and assigns each read to a tissue spot.

- Locates the "L2" linker sequence (30 bp anchor) in R2 to find barcode positions
- Extracts barcode B (8 bp) and barcode A (8 bp) flanking the linker
- Matches the combined B+A barcode against `spatial_barcodes.txt` with up to 1 mismatch allowed
- Encodes the matched spot coordinates into the read name as `@read_id|:_:|X_Y`
- Outputs `combine.fq` containing only reads that matched a valid spatial barcode

**Usage:**
```bash
python demultiplexer.py \
  --r1 <R1.fastq.gz> \
  --r2 <R2.fastq.gz> \
  --barcode <spatial_barcodes.txt> \
  --output <combine.fq>
```

### `mapping.py`
Maps demultiplexed reads to a bacterial reference genome using BWA-MEM2.

- Aligns `combine.fq` to the reference genome
- Filters alignments by mapping quality (MAPQ ≥ 10)
- Sorts and indexes the output BAM file
- Preserves the spatial spot tag in read names throughout

**Usage:**
```bash
python mapping.py \
  --input <combine.fq> \
  --reference <reference.fasta> \
  --output <align.sorted.bam> \
  --threads <N>
```

### `count_contig.py`
Counts aligned reads per contig per spatial spot to produce the abundance matrix.

- Iterates over the sorted BAM, parsing the `X_Y` spot coordinate from each read name
- Accumulates read counts per (contig, spot) pair
- Outputs a tab-separated matrix with contigs as rows and spot coordinates as columns

**Usage:**
```bash
python count_contig.py \
  --bam <align.sorted.bam> \
  --output <expmat_contig.tsv>
```

### `run_spatial_dna_pipeline.sh`
Bash wrapper that runs all three steps sequentially with configurable parameters. Edit the variables at the top of the script before running.

**Usage:**
```bash
bash run_spatial_dna_pipeline.sh
```

## Input Files

| File | Description |
|------|-------------|
| `*_R1.fastq.gz` | Read 1 — DNA fragment reads (host depleted) |
| `*_R2.fastq.gz` | Read 2 — contains spatial barcodes (B + linker + A) |
| `spatial_barcodes.txt` | Whitelist: `BARCODE_B+BARCODE_A  X  Y` |
| Reference FASTA | Bacterial reference genome/contigs for mapping |

## Output Files

| File | Description |
|------|-------------|
| `combine.fq` | Demultiplexed reads with spot coordinates in read names |
| `align.sorted.bam` | Sorted, filtered BAM retaining spatial coordinate tags |
| `expmat_contig.tsv` | Contig × spot abundance matrix (tab-separated) |
| `demux.log` | Barcode matching statistics |

## Dependencies

- Python 3 with `regex`, `pandas`
- BWA-MEM2
- samtools
