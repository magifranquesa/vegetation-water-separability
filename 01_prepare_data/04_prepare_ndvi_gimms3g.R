#!/usr/bin/env Rscript

# ==============================================================================
# Script: 01_prepare_data/04_prepare_ndvi_gimms3g.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   Prepare the monthly maximum NDVI file from the GIMMS NDVI3g+ dataset
#   (Pinzon et al. 2023). The pipeline has two steps:
#
#   Step 1 (this script): Merge annual NC4 files into a single NetCDF covering
#     the full time series. The raw files contain fortnightly (15-day) composites.
#     NDVI values are rescaled from integer (x10000) to float.
#
#   Step 2 (CDO, run after this script): Aggregate fortnightly composites to
#     monthly maxima using:
#       cdo monmax <merged_file> <output_file>
#
# Data source:
#   Pinzon, J. E., Pak, E. W., Tucker, C. J., Bhatt, U. S., Frost, G. V., &
#   Macander, M. J. (2023). Global Vegetation Greenness (NDVI) from AVHRR
#   GIMMS-3G+, 1981-2022 (Version 1). ORNL DAAC.
#   https://doi.org/10.3334/ORNLDAAC/2187
#
# Input:
#   data/raw/NDVI/GIMMS3G/*.nc4   — annual GIMMS NDVI3g+ files (fortnightly)
#
# Output (Step 1):
#   data/raw/NDVI/ndvi3g_global_1982_2022.nc   — merged fortnightly time series
#
# Output (Step 2, CDO):
#   data/raw/NDVI/NDVI3g_monthlymax_1982_2022.nc   — monthly maximum composites
#
# Dependencies:
#   ncdf4
#   R/helpers_ncdf.R   (included in this repository)
#
# Usage:
#   Run from the repository root:
#     Rscript 01_prepare_data/04_prepare_ndvi_gimms3g.R
#   Then:
#     cdo monmax data/raw/NDVI/ndvi3g_global_1982_2022.nc \
#                data/raw/NDVI/NDVI3g_monthlymax_1982_2022.nc
# ==============================================================================

suppressPackageStartupMessages(library(ncdf4))
source(file.path("R", "function_write_ncdf.R"))

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

input_dir  <- file.path("data", "raw", "NDVI", "GIMMS3G")
merged_out <- file.path("data", "raw", "NDVI", "ndvi3g_global_1982_2022.nc")

if (!dir.exists(input_dir)) {
  stop("Input directory not found: ", input_dir,
       "\nDownload GIMMS NDVI3g+ NC4 files from https://doi.org/10.3334/ORNLDAAC/2187",
       " and place them in: ", input_dir, call. = FALSE)
}

nc_files <- list.files(input_dir, pattern = "\\.nc4$", full.names = TRUE)

if (length(nc_files) == 0) {
  stop("No .nc4 files found in: ", input_dir, call. = FALSE)
}

nc_files <- sort(nc_files)
message("Found ", length(nc_files), " NC4 files.")

# ------------------------------------------------------------------------------
# Step 1: Read first file to define structure
# ------------------------------------------------------------------------------

message("Reading structure from: ", basename(nc_files[1]))
nc_example <- nc_open(nc_files[1])

lats      <- ncvar_get(nc_example, "lat")
lons      <- ncvar_get(nc_example, "lon")
nc_data   <- ncvar_get(nc_example, "ndvi") / 10000
time      <- ncvar_get(nc_example, "time")
time_bnds <- ncvar_get(nc_example, "time_bnds")

nc_close(nc_example)

var_list   <- list(
  name      = "ndvi",
  name.long = "maximum normalized difference vegetation index (NDVI) over composite period",
  unit      = "1",
  range     = c(-0.3, 1)
)
lon_def    <- list(name = "lon", unit = "degrees_east",  values = lons, name.long = "longitude coordinate")
lat_def    <- list(name = "lat", unit = "degrees_north", values = lats, name.long = "latitude coordinate")
climatology <- list(cell.methods = "maximum value for the composite period", time.bounds = time_bnds)

# ------------------------------------------------------------------------------
# Step 1: Write first file and open for incremental writing
# ------------------------------------------------------------------------------

message("Creating output file: ", merged_out)
write_ncdf(nc_data,
           file.out    = merged_out,
           var         = var_list,
           lon         = lon_def,
           lat         = lat_def,
           time        = list(values = time),
           climatology = climatology,
           crs         = list(epsg = "4326"),
           overwrite   = TRUE)

nc_out       <- nc_open(merged_out, write = TRUE)
time_counter <- length(time)

# ------------------------------------------------------------------------------
# Step 1: Append remaining files
# ------------------------------------------------------------------------------

for (f in nc_files[-1]) {
  message("  Appending: ", basename(f))
  nc        <- nc_open(f)
  nc_data   <- ncvar_get(nc, "ndvi") / 10000
  time      <- ncvar_get(nc, "time")
  time_bnds <- ncvar_get(nc, "time_bnds")
  nc_close(nc)

  ncvar_put(nc_out, "ndvi", nc_data,
            start = c(1, 1, time_counter + 1),
            count = c(-1, -1, length(time)))
  ncvar_put(nc_out, "time", time,
            start = time_counter + 1,
            count = length(time))
  ncvar_put(nc_out, "time_bounds", time_bnds,
            start = c(1, time_counter + 1),
            count = c(2, length(time)))

  time_counter <- time_counter + length(time)
}

nc_close(nc_out)
message("Step 1 complete. Merged file: ", merged_out)

# ------------------------------------------------------------------------------
# Step 2: Monthly maximum via CDO
# ------------------------------------------------------------------------------

monthly_out <- file.path("data", "raw", "NDVI", "NDVI3g_monthlymax_1982_2022.nc")

if (Sys.which("cdo") == "") {
  message("\nCDO not found. Run manually:")
  message("  cdo monmax ", merged_out, " ", monthly_out)
} else {
  message("\nStep 2: Computing monthly maximum (cdo monmax)...")
  cmd <- paste("cdo monmax", shQuote(merged_out), shQuote(monthly_out))
  ret <- system(cmd)
  if (ret == 0) {
    message("Step 2 complete. Monthly max file: ", monthly_out)
    message("Removing intermediate merged file...")
    file.remove(merged_out)
  } else {
    warning("CDO returned non-zero exit code. Check output.")
  }
}
