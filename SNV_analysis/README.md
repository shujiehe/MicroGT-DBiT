# SNV_analysis

This module detects single nucleotide variants (SNVs) in bacterial genomes with spatial resolution. It provides two workflows depending on whether spatial barcode decoding is needed.

## Workflows

### Workflow 1 — Spatial barcode decoding + spot-tagged BAM
**Folder:** `spatial barcode decoding followed by spot-tagged Bowtie2 BAM generation/`

Use this when your BAM does not yet have spatial spot assignments. Decodes spatial barcodes from R2 and tags each alignment with its spot coordinate before variant calling.

```
R1 + R2 FASTQ
        │
        ▼
decode_r2_barcode_to_spot.py  ──►  qname_to_spot.tsv.gz
        │
        ▼
Bowtie2 mapping (with --reorder flag)
+ tag_sam_with_spot.py (pipe)
        │
        ▼
bacteria.R1.bowtie2.sorted.bam  (ZS:Z:<spot> tag on each alignment)
        │
        ▼
inStrain profile  ──►  SNV profiles per spot
```

### Workflow 2 — Standard Bowtie2 mapping (no spatial barcodes)
**Folder:** `standard Bowtie2 mapping without spatial barcode decoding/`

Use this for bulk (non-spatial) SNV analysis, or when spatial information is not required.

```
R1 FASTQ
    │
    ▼
Bowtie2 mapping
    │
    ▼
sorted BAM  ──►  inStrain profile
```

---

## Scripts

### `decode_r2_barcode_to_spot.py`
Extracts spatial barcodes from R2 reads and maps each read to a tissue spot.

- Reads barcode B (positions 33–40: 8bp) and barcode A (positions 71–78: 8bp) from R2
- Combines barcode B + barcode A and matches it to the spatial barcode whitelist
- Keeps only spots listed in the position file
- Outputs a gzipped read-to-spot table: `qname_to_spot.tsv.gz`

**Usage:**
```bash
python decode_r2_barcode_to_spot.py \
  --r2 <R2.fastq.gz> \
  --barcode-map <barcode_map.csv> \
  --position <position.txt> \
  --out qname_to_spot.tsv.gz
```

### `tag_sam_with_spot.py`
Tags SAM records with their spatial spot using the mapping produced by `decode_r2_barcode_to_spot.py`.

- Reads SAM records from stdin
- Reads `qname_to_spot.tsv.gz` in the same read order
- Appends a `ZS:Z:<spot>` tag to records with `status=OK`
- Discards reads with failed or ambiguous barcode assignments
- Requires the SAM stream and qname-to-spot table to have the same read order
- Bowtie2 must be run with `--reorder` to preserve read order for streaming

**Usage (as part of Bowtie2 pipe):**
```bash
bowtie2 --reorder -x <index> -U <R1.fastq.gz> \
  | python tag_sam_with_spot.py --map qname_to_spot.tsv.gz \
  | samtools sort -o bacteria.R1.bowtie2.sorted.bam
```

See `bowtie2_with_spot.example.sh` for a complete runnable example.

### `bowtie2_with_spot.example.sh`
End-to-end example script for Workflow 1 (spatial barcode decoding + spot-tagged BAM).

### `bowtie2_no_spot.example.sh`
Example script for Workflow 2 (standard Bowtie2 mapping, no spatial barcodes).

### `instrain_profile.example.sh`
Example inStrain command for SNV profiling from a sorted BAM.

Key parameters:
| Parameter | Value | Description |
|-----------|-------|-------------|
| `--min_mapq` | 40 | Minimum mapping quality |
| `--min_read_ani` | 0.92 | Minimum average nucleotide identity per read |
| `--min_coverage` | 5 | Minimum coverage depth for SNV profiling |

**Usage:**
```bash
inStrain profile \
  <sorted.bam> <reference.fasta> \
  -o <output_dir> \
  -p 16 \
  --pairing_filter all_reads \
  --stb <contig_to_species.stb> \
  --min_mapq 40 \
  --min_read_ani 0.92 \
  --min_cov 5 \
  --min_freq 0.05 \
  --fdr 1e-6 \
  --skip_mm_profiling
```

---

## Input Files

| File | Description |
|------|-------------|
| `*_R1.fastq.gz` | Read 1 — sequencing reads for alignment |
| `*_R2.fastq.gz` | Read 2 — contains spatial barcodes (Workflow 1 only) |
| `barcode_map.csv` | Barcode-to-coordinate mapping CSV (Workflow 1 only) |
| `position.txt` | Valid tissue spot coordinates (Workflow 1 only) |
| Bowtie2 index | Pre-built Bowtie2 index for the bacterial reference |
| Reference FASTA | Bacterial reference genome used by inStrain profiling |
| `*.stb` | Contig-to-species mapping file used by inStrain |

## Output Files

| File | Description |
|------|-------------|
| `qname_to_spot.tsv.gz` | Read-to-spot mapping table generated from R2 barcode decoding (Workflow 1 only) |
| `*.bowtie2.sorted.bam` | Sorted BAM file from Bowtie2 alignment; contains `ZS:Z:<spot>` spatial tags in Workflow 1 |
| inStrain output directory | inStrain profiling results, including SNV tables, coverage, and species-level microdiversity metrics |

## Dependencies

- Python 3 with `pandas`
- Bowtie2
- samtools
- inStrain
