library(Matrix)
library(data.table)
library(parallel)
library(EcoSimR)
library(ggplot2)
library(scales)
library(ragg)

# ============================================================
# User settings
# ============================================================
sample_names <- c(
  "P.Colon",
  "D.Colon",
  "DSS_P.Colon",
  "DSS_D.Colon"
)

matrix_files <- c(
  "/path/to/P.Colon/speciesxspot.tsv",
  "/path/to/D.Colon/speciesxspot.tsv",
  "/path/to/DSS_P.Colon/speciesxspot.tsv",
  "/path/to/DSS_D.Colon/speciesxspot.tsv"
)

out_dir <- "./species_colocalization_sim9"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Parameters
# ============================================================
SPOT_TOTAL_MIN <- 25

PRESENT_COUNT_MIN <- 2
PRESENT_RA_MIN <- 0.01

SPOT_RICHNESS_MIN <- 2
SPECIES_PRESENT_SPOTS_MIN <- 10

N_SWAPS <- 25000
N_NULL <- 50
SEED <- 123

FDR_CUTOFF <- 0.05

TOP_N_HEATMAP <- 100
HEATMAP_Z_CAP <- 10

n_workers <- as.integer(Sys.getenv("N_THREADS", "1"))
if (is.na(n_workers) || n_workers < 1) {
  n_workers <- 1
}

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)

# ============================================================
# Helper functions
# ============================================================
read_species_spot_matrix <- function(matrix_file) {
  dt <- fread(
    matrix_file,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    data.table = TRUE,
    nThread = n_workers,
    showProgress = FALSE
  )
  
  setnames(dt, 1, "species")
  dt <- dt[!is.na(species) & trimws(as.character(species)) != ""]
  dt[, species := as.character(species)]
  
  spot_cols <- setdiff(names(dt), "species")
  
  dt[, (spot_cols) := lapply(.SD, function(x) {
    x <- suppressWarnings(as.numeric(x))
    x[is.na(x)] <- 0
    x
  }), .SDcols = spot_cols]
  
  dt <- dt[, lapply(.SD, sum, na.rm = TRUE), by = species, .SDcols = spot_cols]
  
  mat <- as.matrix(dt[, ..spot_cols])
  rownames(mat) <- dt$species
  storage.mode(mat) <- "double"
  
  mat
}

make_upper_pairs <- function(p) {
  n_pairs <- p * (p - 1) / 2
  ii <- integer(n_pairs)
  jj <- integer(n_pairs)
  
  pos <- 1L
  for (i in seq_len(p - 1L)) {
    n <- p - i
    idx <- pos:(pos + n - 1L)
    ii[idx] <- i
    jj[idx] <- (i + 1L):p
    pos <- pos + n
  }
  
  list(i = ii, j = jj)
}

jaccard_stats_vec <- function(X, ii, jj) {
  X <- Matrix::drop0(Matrix::Matrix(X > 0, sparse = TRUE))
  
  species_present_spots <- Matrix::colSums(X)
  cooccur_mat <- Matrix::crossprod(X)
  
  observed_cooccur <- as.numeric(cooccur_mat[cbind(ii, jj)])
  union_spots <- species_present_spots[ii] + species_present_spots[jj] - observed_cooccur
  
  jaccard_similarity <- observed_cooccur / pmax(union_spots, 1)
  jaccard_distance <- 1 - jaccard_similarity
  
  list(
    observed_cooccur = observed_cooccur,
    jaccard_similarity = jaccard_similarity,
    jaccard_distance = jaccard_distance
  )
}

save_plot_png <- function(plot_obj, filename, width, height, res = 600) {
  ragg::agg_png(
    filename = filename,
    width = width,
    height = height,
    units = "in",
    res = res
  )
  print(plot_obj)
  dev.off()
}

# ============================================================
# Per-sample analysis
# ============================================================
run_one_sample <- function(sample_name, matrix_file) {
  sample_out_dir <- file.path(out_dir, sample_name)
  dir.create(sample_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  count_mat <- read_species_spot_matrix(matrix_file)
  
  raw_species <- nrow(count_mat)
  raw_spots <- ncol(count_mat)
  raw_spot_total <- colSums(count_mat)
  
  keep_spot_total <- raw_spot_total >= SPOT_TOTAL_MIN
  
  count_spot_filtered <- count_mat[, keep_spot_total, drop = FALSE]
  spot_total_filtered <- colSums(count_spot_filtered)
  
  ra_mat <- sweep(count_spot_filtered, 2, spot_total_filtered, "/")
  present_all <- (count_spot_filtered >= PRESENT_COUNT_MIN) & (ra_mat >= PRESENT_RA_MIN)
  present_all[is.na(present_all)] <- FALSE
  
  keep_species <- rowSums(present_all) >= SPECIES_PRESENT_SPOTS_MIN
  keep_spots <- colSums(present_all[keep_species, , drop = FALSE]) > SPOT_RICHNESS_MIN
  
  for (iter_i in seq_len(20)) {
    keep_species_new <- rowSums(present_all[, keep_spots, drop = FALSE]) >= SPECIES_PRESENT_SPOTS_MIN
    keep_spots_new <- colSums(present_all[keep_species_new, , drop = FALSE]) > SPOT_RICHNESS_MIN
    
    if (identical(keep_species_new, keep_species) && identical(keep_spots_new, keep_spots)) {
      break
    }
    
    keep_species <- keep_species_new
    keep_spots <- keep_spots_new
  }
  
  present_final <- present_all[keep_species, keep_spots, drop = FALSE]
  
  Xobs <- Matrix::Matrix(t(present_final * 1L), sparse = TRUE)
  rownames(Xobs) <- colnames(present_final)
  colnames(Xobs) <- rownames(present_final)
  
  n_spots_final <- nrow(Xobs)
  n_species_final <- ncol(Xobs)
  
  final_species_present_spots <- Matrix::colSums(Xobs)
  
  # ============================================================
  # Added output: species prevalence after QC
  # ============================================================
  count_final_spots <- count_mat[colnames(Xobs), rownames(Xobs), drop = FALSE]
  
  species_prevalence_dt <- data.table(
    sample = sample_name,
    species = colnames(Xobs),
    total_count_in_final_spots = as.numeric(rowSums(count_final_spots)),
    present_spots = as.numeric(final_species_present_spots),
    prevalence = as.numeric(final_species_present_spots / n_spots_final)
  )
  
  species_prevalence_dt[, mean_count_per_present_spot := total_count_in_final_spots / present_spots]
  setorder(species_prevalence_dt, -present_spots, -total_count_in_final_spots)
  
  fwrite(
    species_prevalence_dt,
    file.path(sample_out_dir, "species_prevalence_after_qc.tsv"),
    sep = "\t",
    quote = FALSE
  )
  
  taxa <- colnames(Xobs)
  p <- length(taxa)
  
  pair_idx <- make_upper_pairs(p)
  ii <- pair_idx$i
  jj <- pair_idx$j
  n_pairs <- length(ii)
  
  obs_stats <- jaccard_stats_vec(Xobs, ii, jj)
  
  d_obs <- obs_stats$jaccard_distance
  observed_cooccur <- obs_stats$observed_cooccur
  observed_jaccard_similarity <- obs_stats$jaccard_similarity
  
  speciesData0 <- (t(as.matrix(Xobs)) > 0) * 1L
  
  rep_ids <- seq_len(N_NULL)
  rep_seeds <- SEED + rep_ids
  
  worker_one_null <- function(rep_id) {
    set.seed(rep_seeds[rep_id])
    
    M <- speciesData0
    
    for (swap_i in seq_len(N_SWAPS)) {
      M <- EcoSimR::sim9_single(M)
    }
    
    X_random <- Matrix::Matrix(t(M) > 0, sparse = TRUE)
    colnames(X_random) <- taxa
    
    jaccard_stats_vec(X_random, ii, jj)$jaccard_distance
  }
  
  d_null_list <- parallel::mclapply(
    rep_ids,
    worker_one_null,
    mc.cores = min(n_workers, N_NULL),
    mc.preschedule = FALSE
  )
  
  d_null_mat <- do.call(cbind, d_null_list)
  
  null_mean_distance <- rowMeans(d_null_mat)
  null_sd_distance <- apply(d_null_mat, 1, sd)
  
  z_score <- (null_mean_distance - d_obs) / null_sd_distance
  z_score[!is.finite(z_score)] <- NA_real_
  
  p_value <- 2 * pnorm(-abs(z_score))
  p_value[!is.finite(p_value)] <- NA_real_
  
  FDR <- p.adjust(p_value, method = "BH")
  
  direction <- ifelse(
    !is.na(FDR) & FDR < FDR_CUTOFF & !is.na(z_score) & z_score > 0,
    "positive_colocalization",
    ifelse(
      !is.na(FDR) & FDR < FDR_CUTOFF & !is.na(z_score) & z_score < 0,
      "negative_colocalization",
      "not_significant"
    )
  )
  
  all_pairs <- data.table(
    sample = sample_name,
    species1 = taxa[ii],
    species2 = taxa[jj],
    present_spots1 = as.numeric(final_species_present_spots[ii]),
    present_spots2 = as.numeric(final_species_present_spots[jj]),
    prevalence1 = as.numeric(final_species_present_spots[ii] / n_spots_final),
    prevalence2 = as.numeric(final_species_present_spots[jj] / n_spots_final),
    observed_cooccur_spots = observed_cooccur,
    observed_jaccard_similarity = observed_jaccard_similarity,
    observed_jaccard_distance = d_obs,
    null_mean_jaccard_distance = null_mean_distance,
    null_sd_jaccard_distance = null_sd_distance,
    z_score = z_score,
    p_value = p_value,
    FDR = FDR,
    direction = direction
  )
  
  all_pairs[, abs_z_score := abs(z_score)]
  setorder(all_pairs, FDR, -abs_z_score)
  
  fwrite(
    all_pairs,
    file.path(sample_out_dir, "all_pairwise_species_colocalization.tsv"),
    sep = "\t",
    quote = FALSE
  )
  
  significant_positive_n <- sum(all_pairs$direction == "positive_colocalization", na.rm = TRUE)
  significant_negative_n <- sum(all_pairs$direction == "negative_colocalization", na.rm = TRUE)
  significant_total_n <- significant_positive_n + significant_negative_n
  
  z_mat <- matrix(NA_real_, nrow = n_species_final, ncol = n_species_final)
  rownames(z_mat) <- taxa
  colnames(z_mat) <- taxa
  
  z_mat[cbind(ii, jj)] <- z_score
  z_mat[cbind(jj, ii)] <- z_score
  diag(z_mat) <- NA_real_
  
  z_mat_dt <- as.data.table(z_mat, keep.rownames = "species")
  
  fwrite(
    z_mat_dt,
    file.path(sample_out_dir, "pairwise_zscore_matrix.tsv"),
    sep = "\t",
    quote = FALSE
  )
  
  top_species <- species_prevalence_dt$species[
    seq_len(min(TOP_N_HEATMAP, nrow(species_prevalence_dt)))
  ]
  
  if (length(top_species) >= 3) {
    heat_pairs <- all_pairs[species1 %in% top_species & species2 %in% top_species]
    
    z_top <- matrix(NA_real_, nrow = length(top_species), ncol = length(top_species))
    rownames(z_top) <- top_species
    colnames(z_top) <- top_species
    
    sig_top <- matrix("", nrow = length(top_species), ncol = length(top_species))
    rownames(sig_top) <- top_species
    colnames(sig_top) <- top_species
    
    i1 <- match(heat_pairs$species1, top_species)
    i2 <- match(heat_pairs$species2, top_species)
    
    z_top[cbind(i1, i2)] <- heat_pairs$z_score
    z_top[cbind(i2, i1)] <- heat_pairs$z_score
    
    sig_idx <- !is.na(heat_pairs$FDR) &
      heat_pairs$FDR < FDR_CUTOFF &
      is.finite(heat_pairs$z_score)
    
    sig_top[cbind(i1[sig_idx], i2[sig_idx])] <- "*"
    sig_top[cbind(i2[sig_idx], i1[sig_idx])] <- "*"
    
    diag(z_top) <- NA_real_
    diag(sig_top) <- ""
    
    z_cluster <- z_top
    z_cluster[is.na(z_cluster) | !is.finite(z_cluster)] <- 0
    diag(z_cluster) <- 0
    
    hc <- hclust(dist(z_cluster), method = "ward.D2")
    species_order <- hc$labels[hc$order]
    
    z_top[z_top > HEATMAP_Z_CAP] <- HEATMAP_Z_CAP
    z_top[z_top < -HEATMAP_Z_CAP] <- -HEATMAP_Z_CAP
    
    heatmap_dt <- as.data.table(as.table(z_top))
    setnames(heatmap_dt, c("species1", "species2", "z_score"))
    
    sig_dt <- as.data.table(as.table(sig_top))
    setnames(sig_dt, c("species1", "species2", "sig_label"))
    
    heatmap_dt <- merge(
      heatmap_dt,
      sig_dt,
      by = c("species1", "species2"),
      all.x = TRUE,
      sort = FALSE
    )
    
    heatmap_dt[is.na(sig_label), sig_label := ""]
    
    heatmap_dt[, species1 := factor(species1, levels = species_order)]
    heatmap_dt[, species2 := factor(species2, levels = rev(species_order))]
    
    p_heat <- ggplot(heatmap_dt, aes(x = species1, y = species2, fill = z_score)) +
      geom_tile(color = "grey90", linewidth = 0.1) +
      geom_text(
        data = heatmap_dt[sig_label != ""],
        aes(label = sig_label),
        color = "black",
        size = 2.0,
        fontface = "bold"
      ) +
      scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-HEATMAP_Z_CAP, HEATMAP_Z_CAP),
        breaks = c(-HEATMAP_Z_CAP, 0, HEATMAP_Z_CAP),
        labels = c(
          paste0("-", HEATMAP_Z_CAP),
          "0",
          as.character(HEATMAP_Z_CAP)
        ),
        oob = scales::squish,
        na.value = "grey95",
        name = "Z-score"
      ) +
      coord_fixed() +
      theme_bw(base_family = "Arial") +
      labs(
        title = paste0(sample_name, " whole-region species co-localization"),
        subtitle = paste0("Top ", length(top_species), " species by prevalence; * FDR < 0.05"),
        x = NULL,
        y = NULL
      ) +
      theme(
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 5),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8)
      )
    
    save_plot_png(
      p_heat,
      file.path(sample_out_dir, "top_species_pairwise_zscore_heatmap.png"),
      width = 10,
      height = 9
    )
  }
  
  data.table(
    sample = sample_name,
    raw_species = raw_species,
    raw_spots = raw_spots,
    spots_after_spot_total_filter = ncol(count_spot_filtered),
    final_spots = n_spots_final,
    final_species = n_species_final,
    total_pairs_tested = n_pairs,
    significant_positive_pairs = significant_positive_n,
    significant_negative_pairs = significant_negative_n,
    total_significant_pairs = significant_total_n
  )
}

# ============================================================
# Run all samples
# ============================================================
all_summary <- rbindlist(
  lapply(seq_along(sample_names), function(i) {
    run_one_sample(sample_names[i], matrix_files[i])
  }),
  fill = TRUE
)

fwrite(
  all_summary,
  file.path(out_dir, "all_samples_summary.tsv"),
  sep = "\t",
  quote = FALSE
)