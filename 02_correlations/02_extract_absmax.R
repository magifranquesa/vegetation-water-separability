#!/usr/bin/env Rscript

# ==============================================================================
# Script: 02_correlations/02_extract_absmax.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   Extract, for each grid cell, calendar month, and hydroclimatic indicator, the
#   Spearman correlation coefficient with the largest absolute magnitude across
#   timescales. The selected coefficient is retained with its original
#   sign and is referred to in the manuscript as rho* (ρ*). The corresponding
#   timescale is stored as scale*.
#
# Input:
#   NetCDF files containing monthly Spearman correlations between kNDVI and each
#   hydroclimatic indicator at individual timescales:
#
#     spearman_correlation_kndvi_<indicator>_scale_<scale>.nc
#
#   Expected indicators:
#     Ep   : atmospheric evaporative demand (AED)
#     Et   : plant transpiration
#     ED   : evaporation deficit
#     SMrz : root-zone soil moisture
#     SMs  : surface soil moisture
#
#   Expected accumulation timescales:
#     1, 3, 6, 9, and 12 months
#
# Output:
#   One NetCDF file per hydroclimatic indicator:
#
#     abs_max_correlation_kndvi_<indicator>.nc
#
#   Each output file contains:
#     abs_max_correlation   : signed Spearman correlation with largest |rho|
#     abs_max_significance  : significance code associated with selected scale
#     abs_max_p_value       : p-value associated with selected scale
#     abs_max_scale         : accumulation timescale at which abs_max_correlation occurs
#
# Notes:
#   - rho* is selected as the correlation with the largest absolute magnitude
#     across timescales, while retaining its sign.
#   - rho* is not the maximum positive correlation. The maximum positive
#     correlation (rho_max), used for the spatial agreement heatmap, is computed
#     separately.
#   - scale* stores the actual timescale in months, not the index of
#     the selected scale.
#   - In case of exact ties in |rho| across timescales, which.max()
#     returns the first occurrence; because the scales are ordered increasingly,
#     the shortest tied timescale is selected.
#   - The output NetCDF files preserve the spatial dimension order used in the
#     correlation files: lon, lat, month.
#   - Ep is used as the GLEAM variable representing AED.
#   - Total evaporation (E) is not processed here because it is used only to
#     compute ED and is not analysed as an independent hydroclimatic indicator.
#
# Dependencies:
#   ncdf4
#   abind
#
# Repository location:
#   02_correlations/02_extract_absmax.R
# ==============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(abind)
})

# ------------------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------------------

# If an R/config.R file exists, paths can be defined there. Otherwise, the
# relative defaults below are used.
config_file <- file.path("R", "config.R")
if (file.exists(config_file)) {
  source(config_file)
}

# Default paths. Replace these with project-specific paths if not using config.
input_dir <- if (exists("paths") && !is.null(paths$correlations_spearman)) {
  paths$correlations_spearman
} else {
  file.path("outputs", "intermediate", "correlations_spearman")
}

output_dir <- if (exists("paths") && !is.null(paths$absmax_spearman)) {
  paths$absmax_spearman
} else {
  file.path("outputs", "intermediate", "absmax_spearman")
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

indicator_list <- c("Ep", "Et", "ED", "SMrz", "SMs")
scales <- c(1, 3, 6, 9, 12)
n_scales <- length(scales)

input_file_template <- "spearman_correlation_kndvi_%s_scale_%d.nc"
output_file_template <- "abs_max_correlation_kndvi_%s.nc"

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

read_month_dimension <- function(nc) {
  dim_names <- names(nc$dim)

  if ("month" %in% dim_names) {
    vals <- ncvar_get(nc, "month")
  } else if ("time" %in% dim_names) {
    vals <- ncvar_get(nc, "time")
  } else {
    vals <- seq_len(12)
  }

  vals
}

check_required_variables <- function(nc, file_path) {
  required_vars <- c("correlation", "significance", "p_value")
  missing_vars <- setdiff(required_vars, names(nc$var))

  if (length(missing_vars) > 0) {
    stop(
      "Missing variable(s) in ", file_path, ": ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
}

read_correlation_file <- function(file_path, expected_dims = NULL) {
  if (!file.exists(file_path)) {
    stop("Input file not found: ", file_path, call. = FALSE)
  }

  nc <- nc_open(file_path)
  on.exit(nc_close(nc), add = TRUE)

  check_required_variables(nc, file_path)

  correlation <- ncvar_get(nc, "correlation")
  significance <- ncvar_get(nc, "significance")
  p_value <- ncvar_get(nc, "p_value")

  if (!is.null(expected_dims) && !all(dim(correlation) == expected_dims)) {
    stop(
      "Unexpected dimensions in ", file_path, ". Expected ",
      paste(expected_dims, collapse = " x "), " but found ",
      paste(dim(correlation), collapse = " x "), ".",
      call. = FALSE
    )
  }

  list(
    correlation = correlation,
    significance = significance,
    p_value = p_value
  )
}

write_absmax_netcdf <- function(
  output_file,
  lon_vals,
  lat_vals,
  month_vals,
  indicator,
  scales,
  abs_max_correlation,
  abs_max_significance,
  abs_max_p_value,
  abs_max_scale
) {
  fill_value_double <- -9999
  fill_value_integer <- -9999L

  # Replace NA values explicitly before writing.
  corr_out <- abs_max_correlation
  sig_out <- abs_max_significance
  p_out <- abs_max_p_value
  scale_out <- abs_max_scale

  corr_out[is.na(corr_out)] <- fill_value_double
  sig_out[is.na(sig_out)] <- fill_value_double
  p_out[is.na(p_out)] <- fill_value_double
  scale_out[is.na(scale_out)] <- fill_value_integer

  dim_lon <- ncdim_def("lon", "degrees_east", lon_vals)
  dim_lat <- ncdim_def("lat", "degrees_north", lat_vals)
  dim_month <- ncdim_def("month", "calendar_month", month_vals)

  var_abs_max_corr <- ncvar_def(
    "abs_max_correlation",
    "1",
    list(dim_lon, dim_lat, dim_month),
    fill_value_double,
    longname = "Signed Spearman correlation coefficient with largest absolute magnitude across accumulation windows",
    prec = "double"
  )

  var_abs_max_sig <- ncvar_def(
    "abs_max_significance",
    "1",
    list(dim_lon, dim_lat, dim_month),
    fill_value_double,
    longname = "Significance code associated with selected accumulation window",
    prec = "double"
  )

  var_abs_max_p <- ncvar_def(
    "abs_max_p_value",
    "1",
    list(dim_lon, dim_lat, dim_month),
    fill_value_double,
    longname = "P-value associated with selected accumulation window",
    prec = "double"
  )

  var_abs_max_sc <- ncvar_def(
    "abs_max_scale",
    "months",
    list(dim_lon, dim_lat, dim_month),
    fill_value_integer,
    longname = "Accumulation window at which abs_max_correlation occurs",
    prec = "integer"
  )

  nc_out <- nc_create(
    output_file,
    list(var_abs_max_corr, var_abs_max_sig, var_abs_max_p, var_abs_max_sc),
    force_v4 = TRUE
  )
  on.exit(nc_close(nc_out), add = TRUE)

  ncvar_put(nc_out, "abs_max_correlation", corr_out)
  ncvar_put(nc_out, "abs_max_significance", sig_out)
  ncvar_put(nc_out, "abs_max_p_value", p_out)
  ncvar_put(nc_out, "abs_max_scale", scale_out)

  # Global attributes.
  ncatt_put(
    nc_out, 0, "title",
    "Absolute-maximum monthly Spearman correlations between kNDVI and hydroclimatic indicators"
  )
  ncatt_put(
    nc_out, 0, "project",
    "Aridity and timescale bound the separability of vegetation water stress"
  )
  ncatt_put(nc_out, 0, "period", "1982-2022")
  ncatt_put(nc_out, 0, "indicator", indicator)
  ncatt_put(
    nc_out, 0, "method",
    "For each grid cell and calendar month, the Spearman correlation with the largest absolute magnitude across accumulation windows was selected while retaining its sign."
  )
  ncatt_put(nc_out, 0, "accumulation_windows_months", paste(scales, collapse = ", "))
  ncatt_put(nc_out, 0, "dimension_order", "lon, lat, month")
  ncatt_put(
    nc_out, 0, "significance_codes",
    "-2 = negative significant; -1 = negative non-significant; 1 = positive non-significant; 2 = positive significant"
  )
  ncatt_put(nc_out, 0, "created_by", "02_extract_absmax_spearman_correlations.R")

  # Variable attributes.
  ncatt_put(nc_out, "abs_max_correlation", "description", "rho*: signed Spearman correlation coefficient with largest absolute magnitude across accumulation windows")
  ncatt_put(nc_out, "abs_max_significance", "description", "Significance code corresponding to the selected rho* accumulation window")
  ncatt_put(nc_out, "abs_max_p_value", "description", "P-value corresponding to the selected rho* accumulation window")
  ncatt_put(nc_out, "abs_max_scale", "description", "scale*: actual accumulation window in months corresponding to rho*")
}

# ------------------------------------------------------------------------------
# Main processing loop
# ------------------------------------------------------------------------------

message("Starting absolute-maximum Spearman correlation extraction.")
message("Input directory:  ", input_dir)
message("Output directory: ", output_dir)

for (indicator in indicator_list) {
  message("\nProcessing indicator: ", indicator)

  file_template <- function(scale) {
    file.path(input_dir, sprintf(input_file_template, indicator, scale))
  }

  # Read coordinates and dimensions from the first scale file.
  reference_file <- file_template(scales[1])
  if (!file.exists(reference_file)) {
    stop("Reference input file not found: ", reference_file, call. = FALSE)
  }

  nc_ref <- nc_open(reference_file)
  lon_vals <- ncvar_get(nc_ref, "lon")
  lat_vals <- ncvar_get(nc_ref, "lat")
  month_vals <- read_month_dimension(nc_ref)
  nc_close(nc_ref)

  dims <- c(length(lon_vals), length(lat_vals), length(month_vals))

  message("  Expected dimensions [lon, lat, month]: ", paste(dims, collapse = " x "))

  abs_max_correlation <- array(NA_real_, dims)
  abs_max_significance <- array(NA_real_, dims)
  abs_max_p_value <- array(NA_real_, dims)
  abs_max_scale <- array(NA_integer_, dims)

  correlation_list <- vector("list", n_scales)
  significance_list <- vector("list", n_scales)
  p_value_list <- vector("list", n_scales)

  # Read all accumulation timescales.
  for (i in seq_along(scales)) {
    scale <- scales[i]
    input_file <- file_template(scale)

    message("  Reading accumulation timescale: ", scale, " months")

    data <- read_correlation_file(input_file, expected_dims = dims)

    correlation_list[[i]] <- data$correlation
    significance_list[[i]] <- data$significance
    p_value_list[[i]] <- data$p_value
  }

  # Build 4D arrays: [scale_index, lon, lat, month].
  correlation_array <- abind::abind(correlation_list, along = 0)
  significance_array <- abind::abind(significance_list, along = 0)
  p_value_array <- abind::abind(p_value_list, along = 0)

  expected_4d_dims <- c(n_scales, dims)
  if (!all(dim(correlation_array) == expected_4d_dims)) {
    stop(
      "Unexpected dimensions in correlation_array for ", indicator,
      ". Expected ", paste(expected_4d_dims, collapse = " x "),
      " but found ", paste(dim(correlation_array), collapse = " x "), ".",
      call. = FALSE
    )
  }

  # Select the correlation with the largest absolute magnitude across
  # accumulation timescales for each grid cell and calendar month.
  #
  # In case of ties, which.max() returns the first occurrence, corresponding to
  # the shortest accumulation timescale among tied |rho| values.
  for (lon in seq_len(dims[1])) {
    if (lon %% 100 == 0) {
      message("  Processed longitude index ", lon, " / ", dims[1])
    }

    for (lat in seq_len(dims[2])) {
      for (month in seq_len(dims[3])) {
        scale_values <- correlation_array[, lon, lat, month]

        if (all(is.na(scale_values))) {
          abs_max_correlation[lon, lat, month] <- NA_real_
          abs_max_significance[lon, lat, month] <- NA_real_
          abs_max_p_value[lon, lat, month] <- NA_real_
          abs_max_scale[lon, lat, month] <- NA_integer_
        } else {
          idx <- which.max(abs(scale_values))

          abs_max_correlation[lon, lat, month] <- correlation_array[idx, lon, lat, month]
          abs_max_significance[lon, lat, month] <- significance_array[idx, lon, lat, month]
          abs_max_p_value[lon, lat, month] <- p_value_array[idx, lon, lat, month]
          abs_max_scale[lon, lat, month] <- scales[idx]
        }
      }
    }
  }

  output_file <- file.path(output_dir, sprintf(output_file_template, indicator))

  write_absmax_netcdf(
    output_file = output_file,
    lon_vals = lon_vals,
    lat_vals = lat_vals,
    month_vals = month_vals,
    indicator = indicator,
    scales = scales,
    abs_max_correlation = abs_max_correlation,
    abs_max_significance = abs_max_significance,
    abs_max_p_value = abs_max_p_value,
    abs_max_scale = abs_max_scale
  )

  message("Completed absolute-maximum extraction for indicator: ", indicator)
  message("Output written to: ", output_file)
}

message("\nAll indicators processed successfully.")
