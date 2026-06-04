library(data.table)
library(ggplot2)
library(scales)
library(grid)
library(ragg)

# ============================================================
# User settings
# ============================================================
input_matrix_tsv <- "/path/to/feature_by_spot_matrix.tsv"
species_list_txt <- "/path/to/species_list.txt"
position_file <- "/path/to/position.txt"

out_dir <- "./species_spatial_plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Plot settings
# ============================================================
POINT_SIZE <- 1.72
BG_COLOR <- "white"
BACKGROUND_SPOT_COLOR <- "#FFE1B0"

COUNT_COLORS <- c("#FFC788", "#FFAA68", "#F37953", "#D9483A", "#B2262C", "#971822")
COUNT_VALUES <- scales::rescale(c(0, 0.22, 0.44, 0.66, 0.80, 0.90))

PLOT_W <- 4.2
PLOT_H <- 3.6
PNG_DPI <- 600

COLOR_MAX_MODE <- "median10"

# ============================================================
# Helper functions
# ============================================================
safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

map_features_to_species <- function(feature_vec, species_vec) {
  species_sorted <- species_vec[order(nchar(species_vec), decreasing = TRUE)]
  species_for_feature <- rep(NA_character_, length(feature_vec))
  names(species_for_feature) <- feature_vec
  
  for (sp in species_sorted) {
    idx <- is.na(species_for_feature) & (feature_vec == sp | startsWith(feature_vec, paste0(sp, "_")))
    species_for_feature[idx] <- sp
  }
  
  species_for_feature
}

read_position_spots <- function(position_file) {
  pos_txt <- paste(readLines(position_file, warn = FALSE), collapse = "\n")
  pos_all <- unique(unlist(regmatches(pos_txt, gregexpr("\\d+x\\d+", pos_txt))))
  sp_all <- strsplit(pos_all, "x", fixed = TRUE)
  
  data.table(
    spot = pos_all,
    A = as.numeric(vapply(sp_all, `[[`, "", 1)),
    B = as.numeric(vapply(sp_all, `[[`, "", 2))
  )
}

get_color_max <- function(x, mode = COLOR_MAX_MODE) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  xpos <- x[x > 0]
  
  if (length(xpos) == 0) return(1)
  
  if (mode == "median10") {
    vmax <- 10 * median(xpos, na.rm = TRUE)
  } else if (mode == "p99") {
    vmax <- as.numeric(quantile(x, 0.99, na.rm = TRUE))
  } else if (mode == "max") {
    vmax <- max(x, na.rm = TRUE)
  } else {
    vmax <- 10 * median(xpos, na.rm = TRUE)
  }
  
  if (!is.finite(vmax) || vmax <= 0) vmax <- max(xpos, na.rm = TRUE)
  vmax
}

nice_ceiling <- function(x) {
  x <- as.numeric(x)
  if (!is.finite(x) || x <= 0) return(1)
  
  pow10 <- 10^floor(log10(x))
  scaled <- x / pow10
  
  nice_scaled <- if (scaled <= 1) 1 else if (scaled <= 2) 2 else if (scaled <= 5) 5 else 10
  nice_scaled * pow10
}

format_nice_label <- function(x) {
  if (!is.finite(x)) return("NA")
  
  if (x >= 1000 && x %% 1000 == 0) {
    paste0(x / 1000, "k")
  } else if (x >= 1000) {
    format(x, scientific = FALSE, trim = TRUE, big.mark = "")
  } else if (x >= 1) {
    format(round(x), scientific = FALSE, trim = TRUE)
  } else {
    sprintf("%.2f", x)
  }
}

plot_one_species <- function(full_spot_dt, matrix_value_dt, species_name) {
  plot_dt <- copy(full_spot_dt)
  plot_dt[, count := NA_real_]
  
  m <- match(matrix_value_dt$spot, plot_dt$spot)
  keep <- !is.na(m)
  plot_dt$count[m[keep]] <- matrix_value_dt$count[keep]
  
  vmax <- nice_ceiling(get_color_max(plot_dt$count))
  
  p <- ggplot(plot_dt, aes(x = A, y = B, color = count)) +
    geom_point(shape = 16, size = POINT_SIZE) +
    scale_color_gradientn(
      colours = COUNT_COLORS,
      values = COUNT_VALUES,
      limits = c(0, vmax),
      breaks = c(0, vmax),
      labels = c("0", format_nice_label(vmax)),
      oob = scales::squish,
      na.value = BACKGROUND_SPOT_COLOR,
      name = NULL
    ) +
    scale_y_reverse() +
    coord_equal(expand = TRUE, clip = "off") +
    labs(x = NULL, y = NULL, title = species_name) +
    theme_void(base_family = "Arial") +
    theme(
      text = element_text(family = "Arial", face = "plain"),
      plot.title = element_text(hjust = 0.5, size = 13, face = "plain"),
      plot.background = element_rect(fill = BG_COLOR, colour = NA),
      panel.background = element_rect(fill = BG_COLOR, colour = NA),
      legend.title = element_blank(),
      legend.text = element_text(size = 10, family = "Arial", face = "plain"),
      legend.position = "right",
      plot.margin = margin(t = 8, r = 12, b = 8, l = 8)
    ) +
    guides(color = guide_colorbar(barheight = unit(30, "mm"), barwidth = unit(3, "mm")))
  
  ggsave(
    file.path(out_dir, paste0(safe_name(species_name), "_spatial.png")),
    p,
    width = PLOT_W,
    height = PLOT_H,
    units = "in",
    dpi = PNG_DPI,
    device = ragg::agg_png
  )
}

# ============================================================
# Read input files
# ============================================================
full_spot_dt <- read_position_spots(position_file)

dt <- fread(input_matrix_tsv, sep = "\t", header = TRUE, check.names = FALSE)
feature_col <- names(dt)[1]
setnames(dt, feature_col, "feature")

spot_cols <- setdiff(names(dt), "feature")

dt[, (spot_cols) := lapply(.SD, function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
}), .SDcols = spot_cols]

species_list <- fread(species_list_txt, header = FALSE, col.names = "species")$species
species_list <- trimws(species_list)
species_list <- unique(species_list[species_list != "" & !is.na(species_list)])

dt[, species := map_features_to_species(feature, species_list)]

# ============================================================
# Plot all species in species_list
# ============================================================
summary_list <- vector("list", length(species_list))

for (i in seq_along(species_list)) {
  sp <- species_list[i]
  dt_sp <- dt[species == sp]
  
  if (nrow(dt_sp) > 0) {
    sp_sum_dt <- dt_sp[, lapply(.SD, sum, na.rm = TRUE), .SDcols = spot_cols]
    count_vec <- as.numeric(unlist(sp_sum_dt, use.names = FALSE))
    matched_features <- nrow(dt_sp)
    total_count <- sum(count_vec, na.rm = TRUE)
  } else {
    count_vec <- rep(0, length(spot_cols))
    matched_features <- 0L
    total_count <- 0
  }
  
  matrix_value_dt <- data.table(spot = spot_cols, count = count_vec)
  plot_one_species(full_spot_dt, matrix_value_dt, sp)
  
  summary_list[[i]] <- data.table(
    species = sp,
    matched_features = matched_features,
    total_count = total_count
  )
}

plot_summary_dt <- rbindlist(summary_list)

fwrite(
  plot_summary_dt,
  file.path(out_dir, "species_spatial_plot_summary.tsv"),
  sep = "\t",
  quote = FALSE
)