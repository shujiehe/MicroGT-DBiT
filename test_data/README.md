# test_data

Example data for testing and validating the MicroGT-DBiT pipeline. All files are small subsets intended for quick end-to-end testing, not biological analysis.

## Directory Structure

```
test_data/
├── barcode/               # Spatial barcode whitelist
├── position/              # Tissue spot coordinate file
├── built_gtf/             # Gene annotation file
├── spatial_DNA_data/      # Test DNA FASTQ files
└── spatial_RNA_data/      # Test RNA FASTQ files
```

## Contents

### `barcode/spatial_barcodes.txt`
Whitelist of valid spatial barcodes mapping each barcode combination to a tissue spot coordinate.

Format (tab-separated):
```
BARCODE_B+BARCODE_A    X    Y
AACGTGATAACGTGAT       50   1
...
```

- ~768 barcodes arranged in a 48 × 16 spatial grid
- Used by `DNA_analysis/demultiplexer.py` and `SNV_analysis/decode_r2_barcode_to_spot.py`

### `position/position.txt`
List of tissue spot coordinates to retain after demultiplexing. Spots not in this file are discarded, allowing analysis to be restricted to tissue-covered regions.

Format: comma- or space-separated `XxY` coordinate pairs (1000+ spots).

### `built_gtf/host_plus_bac.gtf`
Combined GTF annotation file containing gene models for both the host genome and bacterial genomes. Used by STAR during RNA alignment to quantify both host and microbial transcripts simultaneously.

### `spatial_DNA_data/`
| File | Description |
|------|-------------|
| `DNA_test_R1.fastq.gz` | Test DNA read 1 (fragment reads) |
| `DNA_test_R2.fastq.gz` | Test DNA read 2 (spatial barcodes) |

These reads are used to test the full `DNA_analysis/` pipeline.

### `spatial_RNA_data/`
| File | Description |
|------|-------------|
| `RNA_test_R1.fastq.gz` | Test RNA read 1 (transcript reads) |
| `RNA_test_R2.fastq.gz` | Test RNA read 2 (spatial barcodes + UMI) |

These reads are used to test the `RNA_analysis/` pipeline.

## Running the Test

**DNA pipeline test:**
```bash
cd DNA_analysis/
bash run_spatial_dna_pipeline.sh \
  --r1 ../test_data/spatial_DNA_data/DNA_test_R1.fastq.gz \
  --r2 ../test_data/spatial_DNA_data/DNA_test_R2.fastq.gz \
  --barcode ../test_data/barcode/spatial_barcodes.txt \
  --reference <your_reference.fasta>
```

**RNA pipeline test:**
Update `RNA_analysis/RNA_parameter.json` to point to the test FASTQ files and GTF, then:
```bash
cd RNA_analysis/
conda activate astro
ASTRO RNA_parameter.json
```

## Notes

- Sample names in the test data (e.g., `P.Colon`, `D.Colon`, `DSS_Colon`) reflect colon tissue sections from a murine model — these are example labels and not required by the pipeline.
- A pre-built reference genome and Bowtie2/BWA-MEM2 index are required separately and are not included here due to file size.
