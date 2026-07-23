#!/usr/bin/env Rscript

# ==============================================================================
# Script: 04_analysis/07_spatial_agreement_rho_max.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   Compute pairwise spatial Spearman correlations between ρ_max maps
#   (one per indicator × month). The C1 vegetation mask (Et > 0 any month
#   AND mean kNDVI > 0.025) is applied before sampling so only vegetated
#   pixels contribute. Output feeds the heatmap in fig4.
#
# Input:
#   NetCDF rasters: outputs/intermediate/max_spearman/max_correlation_kndvi_<v>.nc
#   Vegetation mask C1: outputs/intermediate/vegetation_mask_c1.tif
#
# Output:
#   outputs/tables/spatial_correlations_rho_max.csv
#   Columns: Pair, Month, rho, pval
#
# Dependencies:
#   terra, dplyr, tidyr
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

dir_vi <- if (exists("paths") && !is.null(paths$intermediate)) {
  file.path(paths$intermediate, "max_spearman")
} else {
  file.path("outputs", "intermediate", "max_spearman")
}

file_veg_mask <- if (exists("paths") && !is.null(paths$veg_mask_c1)) {
  paths$veg_mask_c1
} else {
  file.path("outputs", "intermediate", "vegetation_mask_c1.tif")
}

out_file <- if (exists("paths") && !is.null(paths$tables)) {
  file.path(paths$tables, "spatial_correlations_rho_max.csv")
} else {
  file.path("outputs", "tables", "spatial_correlations_rho_max.csv")
}

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

vars <- c("Ep", "Et", "ED", "SMrz", "SMs")

labels <- c(
  Ep   = "AED",
  Et   = "Et",
  ED   = "ED",
  SMrz = "SMrz",
  SMs  = "SMs"
)

meses <- c(
  "January", "February", "March", "April", "May", "June",
  "July", "August", "September", "October", "November", "December"
)

n_samp         <- 50000
extent_longlat <- ext(-170, 170, -65, 90)

# ------------------------------------------------------------------------------
# Load and resample vegetation mask
# ------------------------------------------------------------------------------

if (!file.exists(file_veg_mask))
  stop("Vegetation mask not found: ", file_veg_mask, call. = FALSE)

raster_files <- setNames(
  file.path(dir_vi, paste0("max_correlation_kndvi_", vars, ".nc")),
  vars
)
missing_files <- raster_files[!file.exists(raster_files)]
if (length(missing_files) > 0)
  stop("Missing input raster(s):\n",
       paste(" ", missing_files, collapse = "\n"), call. = FALSE)

message("Resampling vegetation mask to raster grid...")
reference_raster <- rast(raster_files[[vars[1]]])[[1]]
veg_mask_raw     <- rast(file_veg_mask)
veg_mask_res     <- crop(resample(veg_mask_raw, reference_raster, method = "near"),
                         extent_longlat)

# ------------------------------------------------------------------------------
# Helper: Spearman rho and p-value for one pair
# ------------------------------------------------------------------------------

pair_corr <- function(df, xvar, yvar) {
  ct <- suppressWarnings(
    cor.test(df[[xvar]], df[[yvar]], method = "spearman")
  )
  data.frame(
    Pair = paste0(labels[[xvar]], " vs ", labels[[yvar]]),
    rho  = round(unname(ct$estimate), 2),
    pval = ct$p.value
  )
}

# ------------------------------------------------------------------------------
# Main loop
# ------------------------------------------------------------------------------

message("Computing spatial correlations (rho_max, vegmask) for all months.")

df_all_corr <- data.frame()

for (mes_input in meses) {
  message("  Processing: ", mes_input)
  m <- match(mes_input, meses)

  layers <- lapply(vars, function(v) {
    r     <- rast(raster_files[[v]])
    if (is.na(crs(r))) crs(r) <- "EPSG:4326"
    r_cor <- r[[grep("max_correlation", names(r))]]
    if (nlyr(r_cor) < m)
      stop("Raster for '", v, "' has fewer than ", m, " layers.", call. = FALSE)
    lyr <- crop(r_cor[[m]], extent_longlat)
    # Apply vegetation mask: non-vegetated pixels → NA
    lyr[veg_mask_res == 0] <- NA
    names(lyr) <- v
    lyr
  })

  R <- rast(layers)

  samp <- terra::spatSample(
    R,
    size   = n_samp,
    method = "random",
    na.rm  = TRUE,
    as.df  = TRUE,
    xy     = FALSE
  )
  colnames(samp) <- vars
  samp <- drop_na(samp)

  pairs  <- combn(vars, 2, simplify = FALSE)
  df_mes <- bind_rows(lapply(pairs, function(p) pair_corr(samp, p[1], p[2])))
  df_mes$Month <- mes_input

  df_all_corr <- bind_rows(df_all_corr, df_mes[, c("Pair", "Month", "rho", "pval")])
}

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

write.csv(df_all_corr, out_file, row.names = FALSE)
message("Saved: ", out_file)
