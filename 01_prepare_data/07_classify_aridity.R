#!/usr/bin/env Rscript

# ==============================================================================
# Script: 01_prepare_data/07_classify_aridity.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   Classify the mean annual Aridity Index (AI) raster into five aridity zones
#   following the thresholds of Zomer et al. (2022). The classified raster is
#   cropped to latitudes >= -65° and saved as a GeoTIFF. An Equal Earth
#   reprojection is also produced for visual inspection.
#
# Input:
#   NetCDF file containing the mean annual AI for the study period:
#
#     data/processed/aridity/ai_1982_2022_period.nc
#
#   The AI is computed from CRU TS precipitation and PET; see
#   data/raw/AI_CRU/ for the source data and generation scripts.
#
# Output:
#   data/processed/aridity/AI_1982_2022_Zomer_classes.tif
#
#   A categorical raster with five aridity classes:
#     1 = Hyper-arid    (AI < 0.03)
#     2 = Arid          (0.03 <= AI < 0.20)
#     3 = Semi-arid     (0.20 <= AI < 0.50)
#     4 = Dry sub-humid (0.50 <= AI < 0.65)
#     5 = Humid         (AI >= 0.65)
#
# Notes:
#   - Thresholds follow Zomer et al. (2022), where hyper-arid is defined as
#     AI < 0.03. This differs from the IPCC classification (AI < 0.05).
#   - Intervals are left-closed: [lower, upper).
#   - The raster is cropped at 65° S to exclude Antarctica.
#   - The Equal Earth reprojection (EPSG:8857) is used for cartographic
#     display only and is not saved to disk.
#   - AI values outside [0, 1] are not masked here; filtering is applied
#     downstream in analysis scripts.
#
# References:
#   Zomer, R.J., et al. (2022). Version 3 of the Global Aridity Index and
#   Potential Evapotranspiration Database. Scientific Data, 9, 409.
#   https://doi.org/10.1038/s41597-022-01493-1
#
# Dependencies:
#   terra
#
# Repository location:
#   01_prepare_data/07_classify_aridity.R
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
})

# ------------------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------------------

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) {
  source(config_file)
}

input_file <- if (exists("paths") && !is.null(paths$aridity)) {
  file.path(paths$aridity, "ai_1982_2022_period.nc")
} else {
  file.path("data", "processed", "aridity", "ai_1982_2022_period.nc")
}

output_file <- if (exists("paths") && !is.null(paths$aridity)) {
  file.path(paths$aridity, "AI_1982_2022_Zomer_classes.tif")
} else {
  file.path("data", "processed", "aridity", "AI_1982_2022_Zomer_classes.tif")
}

# Zomer et al. (2022) thresholds.
rcl_zomer <- matrix(c(
  -Inf, 0.03, 1,
  0.03, 0.20, 2,
  0.20, 0.50, 3,
  0.50, 0.65, 4,
  0.65,  Inf, 5
), ncol = 3, byrow = TRUE)

class_labels <- c("Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid", "Humid")

cols_ai <- c("#8C510A", "#D8B365", "#F6E8C3", "#C7EAE5", "#01665E")

# ------------------------------------------------------------------------------
# Load and classify
# ------------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

message("Classifying Aridity Index.")
message("Input:  ", input_file)
message("Output: ", output_file)

ai <- rast(input_file)[[1]]

ai_class <- classify(ai, rcl_zomer, right = FALSE)

levels(ai_class) <- data.frame(
  value = 1:5,
  class = class_labels
)

# Crop to exclude Antarctica (lat >= -65°).
e      <- ext(ai_class)
e[3]   <- -65
ai_class_65S <- crop(ai_class, e)

# ------------------------------------------------------------------------------
# Save output
# ------------------------------------------------------------------------------

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

writeRaster(ai_class_65S, output_file, overwrite = TRUE)

message("Saved: ", output_file)

# ------------------------------------------------------------------------------
# Visual check (Equal Earth reprojection — not saved)
# ------------------------------------------------------------------------------

r_8857 <- project(ai_class_65S, "EPSG:8857", method = "near")

plot(r_8857, col = cols_ai, plg = list(title = "Aridity Index (Equal Earth)"))

legend(
  "bottomleft",
  legend = class_labels,
  fill   = cols_ai,
  bty    = "n",
  title  = "Aridity class"
)
