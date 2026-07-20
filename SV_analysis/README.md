# SV_analysis

This module detects structural variants (SVs) — primarily variable-coverage regions and deletions — in bacterial genomes across multiple samples using [SGVFinder2](https://github.com/segallab/SGVFinder2).

## Pipeline Overview

```
R1 FASTQ files (one per sample)
        │
        ▼
run_sgvfinder2_pipeline.py
   ├── Per-sample:  Bowtie2 mapping → delta-depth calling → sample map
   └── Collection:  Cross-sample SV calling
        │
        ▼
variable_sgv.pkl   (variable-coverage SV regions)
deletion_sgv.pkl   (high-confidence deletion SV regions)
```

## Scripts

### `run_sgvfinder2_pipeline.py`
Python wrapper that orchestrates per-sample processing and cross-sample SV calling using the SGVFinder2 library.

**Per-sample steps (run for each input FASTQ):**
1. Bowtie2 alignment to the bacterial reference
2. Delta-depth calculation per genomic position
3. Generation of a per-sample map file

**Collection step (run once across all samples):**
- Combines sample maps and calls SVs using SGVFinder2's collection mode

Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sensitivity` | very-sensitive | Bowtie2 sensitivity setting used by SGVFinder2 |
| `max_spacing` | 10 | Maximum gap between positions in an SV region |
| `dels_detect_thresh` | 0.25 | Depth threshold fraction for deletion detection |

**Usage:**
```bash
python run_sgvfinder2_pipeline.py \
  --db_prefix <sgvfinder2_database_prefix> \
  --fastq_dir <fastq_dir> \
  --icra_root <output_dir>/icra \
  --smp_dir <output_dir>/smp \
  --final_dir <output_dir>/final \
  --min_samp_cutoff 2 \
  --write_csv
```

### `run_sgvfinder2.sh`
Shell script example for running `run_sgvfinder2_pipeline.py` with typical parameters.

**Usage:**
```bash
bash run_sgvfinder2.sh
```

## Input Files

| File | Description |
|------|-------------|
| `*_R1.fastq.gz` | Read 1 FASTQ files (one per sample) |
| SGVFinder2 database prefix | Pre-built SGVFinder2 database prefix for bacterial genome(s) |

## Output Files

| File | Description |
|------|-------------|
| `variable_sgv.pkl` | Python pickle: table of variable coverage SV regions |
| `deletion_sgv.pkl` | Python pickle: table of high-confidence deletion SV regions |
| Per-sample map files | Intermediate depth/delta-depth files (one per sample) |

Loading output in Python:
```python
import pandas as pd
variable_sgv = pd.read_pickle("variable_sgv.pkl")
deletion_sgv = pd.read_pickle("deletion_sgv.pkl")
```

## Dependencies

- Python 3 with `sgvfinder2`, `pandas`
- Bowtie2
