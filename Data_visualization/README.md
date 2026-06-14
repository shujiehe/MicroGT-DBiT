# Data_visualization

R scripts for generating publication-quality spatial plots and statistical analyses of microbial community data from MicroGT-DBiT experiments.

## Scripts

### `species_spatial_plot.R`
Generates spatial heatmaps showing the abundance of each bacterial species across tissue coordinates.

**What it does:**
- Reads a feature-by-spot count matrix (TSV format)
- Parses spot positions from a positions file (`XxY` coordinate format)
- For each species/feature, generates a color-coded spatial heatmap
- Saves one PNG per species at 600 DPI

**Usage:**
```r
Rscript species_spatial_plot.R \
  --matrix <expmat_contig.tsv> \
  --positions <position.txt> \
  --output-dir species_spatial_plots/
```

**Output:** `species_spatial_plots/<species_name>.png`

---

### `DNA_species_colocalization_sim9.R`
Tests whether pairs of bacterial species co-occur more (or less) than expected by chance, using a null model randomization approach.

**What it does:**
- Loads the species-by-spot presence/absence or abundance matrix
- Computes pairwise co-occurrence statistics between all species pairs
- Runs EcoSimR null model randomization (N_SWAPS = 25,000 swap iterations, N_NULL = 50 null communities)
- Applies FDR correction to p-values
- Outputs heatmaps of pairwise colocalization scores and significance

**Key parameters (editable at top of script):**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N_SWAPS` | 25,000 | Number of matrix swaps per null model |
| `N_NULL` | 50 | Number of null communities to generate |
| FDR threshold | 0.05 | Significance cutoff after correction |

**Output:**
- Pairwise colocalization heatmap (PNG)
- Table of species pairs with colocalization scores and FDR-corrected p-values

---

### `RNA_spatial_clustersing.R`
Performs spatial clustering and visualization of RNA expression data across tissue spots.

**What it does:**
- Reads a gene-by-spot RNA expression matrix
- Clusters spots based on gene expression patterns
- Generates spatial plots colored by cluster assignment

---

## Input Files

| File | Description |
|------|-------------|
| `expmat_contig.tsv` | Contig/species × spot count matrix from `DNA_analysis/` |
| `position.txt` | Tissue spot coordinates (XxY format) |
| RNA expression matrix | Gene × spot matrix from `RNA_analysis/` |

## Output Files

| File | Description |
|------|-------------|
| `species_spatial_plots/*.png` | Per-species spatial abundance heatmaps (600 DPI) |
| Colocalization heatmap PNG | Pairwise species colocalization significance map |
| Colocalization table | Species pairs with scores and FDR p-values |
| Cluster spatial plot | RNA-based spot cluster assignments on tissue |

## Dependencies

Install required R packages:
```r
install.packages(c("data.table", "ggplot2", "ragg", "EcoSimR"))
```

| Package | Purpose |
|---------|---------|
| `data.table` | Fast reading and manipulation of large matrices |
| `ggplot2` | Spatial heatmap plotting |
| `ragg` | High-quality PNG output at 600 DPI |
| `EcoSimR` | Null model randomization for colocalization testing |
