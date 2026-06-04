library(Seurat)
library(Matrix)
library(ggplot2)
library(qs)

set.seed(123)

# ============================================================
# User settings
# ============================================================
sample_tbl <- data.frame(
  sampleNames = c(
    "Sample1",
    "Sample2",
    "Sample3"
  ),
  matrix_file = c(
    "/path/to/Sample1/Filtered_matrix.sparse.rds",
    "/path/to/Sample2/Filtered_matrix.sparse.rds",
    "/path/to/Sample3/Filtered_matrix.sparse.rds"
  ),
  stringsAsFactors = FALSE
)

OUTDIR <- "./host_harmony_umap"
SPATIAL_OUTDIR <- file.path(OUTDIR, "spatial_by_sample")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SPATIAL_OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Read one sample
# ============================================================
read_one_sample <- function(rds_path, sampleName) {
  mtx <- readRDS(rds_path)
  
  obj <- CreateSeuratObject(counts = mtx, assay = "RNA", min.cells = 1, min.features = 0)
  obj <- RenameCells(obj, new.names = paste0(sampleName, "#", colnames(obj)))
  obj$sampleNames <- sampleName
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^[Mm][Tt]-")
  
  spot_raw <- sub("^.*?#", "", colnames(obj))
  sp <- strsplit(spot_raw, "x", fixed = TRUE)
  
  obj$X <- suppressWarnings(as.integer(vapply(sp, `[`, "", 1)))
  obj$Y <- suppressWarnings(as.integer(vapply(sp, `[`, "", 2)))
  
  obj
}

# ============================================================
# Read samples
# ============================================================
objs <- lapply(seq_len(nrow(sample_tbl)), function(i) {
  read_one_sample(
    rds_path = sample_tbl$matrix_file[i],
    sampleName = sample_tbl$sampleNames[i]
  )
})

names(objs) <- sample_tbl$sampleNames

# ============================================================
# Merge samples
# ============================================================
obj_merged <- merge(x = objs[[1]], y = objs[-1])

rm(objs)
gc()

# ============================================================
# SCT, PCA, Harmony, UMAP, clustering
# ============================================================
DefaultAssay(obj_merged) <- "RNA"
obj_merged <- SCTransform(obj_merged, vars.to.regress = "percent.mt", verbose = FALSE)

DefaultAssay(obj_merged) <- "SCT"
obj_merged <- RunPCA(obj_merged, assay = "SCT", npcs = 50, verbose = FALSE)
pElbow <- ElbowPlot(obj_merged, ndims = 50)

ggsave(file.path(OUTDIR, "ElbowPlot.png"), pElbow, width = 6, height = 4, dpi = 600)

obj_merged <- IntegrateLayers(
  object = obj_merged,
  method = HarmonyIntegration,
  group.by.vars = "sampleNames",
  normalization.method = "SCT",
  orig.reduction = "pca",
  new.reduction = "harmony"
)

obj_merged <- RunUMAP(obj_merged, reduction = "harmony", dims = 1: 10)
obj_merged <- FindNeighbors(obj_merged, reduction = "harmony", dims = 1: 10, graph.name = "SCT_snn_harmony")
obj_merged <- FindClusters(obj_merged, graph.name = "SCT_snn_harmony", resolution = 0.4)

# ============================================================
# UMAP plots
# ============================================================
p1 <- DimPlot(obj_merged, reduction = "umap", group.by = "sampleNames", label = TRUE, repel = TRUE, pt.size = 0.1)
p2 <- DimPlot(obj_merged, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 0.1)
ggsave(file.path(OUTDIR, "umap_by_sample.png"),  p1, width = 7, height = 6, dpi = 600)
ggsave(file.path(OUTDIR, "umap_by_cluster.png"), p2, width = 7, height = 6, dpi = 600)

qsave(obj_merged, file.path(OUTDIR, "host_harmony_umap.qs"))

# ============================================================
# Spatial plots by sample
# ============================================================
df_all <- obj_merged@meta.data
df_all$spot <- rownames(df_all)

df_all <- data.frame(
  spot = df_all$spot,
  sample = df_all$sampleNames,
  X = as.numeric(df_all$X),
  Y = as.numeric(df_all$Y),
  ctype = as.character(df_all$seurat_clusters),
  stringsAsFactors = FALSE
)

df_all <- df_all[
  is.finite(df_all$X) &
    is.finite(df_all$Y) &
    !is.na(df_all$sample) &
    !is.na(df_all$ctype),
]

gut_palette <- c(
  "0" = "#93CC82",
  "1" = "#4D97CD",
  "2" = "#B379B4",
  "3" = "#DB6968",
  "4" = "#E8C15C"
)

ctype_order <- c("0", "1", "2", "3", "4")

df_all$ctype <- factor(as.character(df_all$ctype), levels = ctype_order)

plot_one_sample <- function(dat, sample_name) {
  ggplot(dat, aes(x = X, y = Y, color = ctype)) +
    geom_point(shape = 16, size = 3.2, alpha = 0.95) +
    scale_color_manual(
      values = gut_palette,
      breaks = ctype_order,
      drop = FALSE,
      name = "Cell type"
    ) +
    scale_y_reverse() +
    coord_equal() +
    labs(title = sample_name) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      legend.position = "right",
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 18),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

samples <- sort(unique(df_all$sample))

for (s in samples) {
  dat_s <- df_all[df_all$sample == s, ]
  
  p <- plot_one_sample(dat_s, s)
  
  ggsave(
    file.path(SPATIAL_OUTDIR, sprintf("spatial_%s.png", s)),
    plot = p,
    width = 8.6,
    height = 8.6,
    dpi = 600
  )
}