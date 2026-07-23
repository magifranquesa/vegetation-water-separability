#!/usr/bin/env Rscript

# ==============================================================================
# Script: 02_correlations/01_monthly_spearman.R
# Project: Aridity and timescale bound the separability of vegetation water stress
# Purpose: Compute monthly pixel-wise Spearman correlations between kNDVI and
#          hydroclimatic indicators at multiple accumulation timescales.
#
# Description:
#   For each hydroclimatic indicator, timescale, grid cell, and
#   calendar month, this script computes Spearman's rank correlation between
#   monthly kNDVI and the corresponding hydroclimatic series over 1982–2022.
#   Both series are linearly detrended before correlation. The script writes
#   one NetCDF file per indicator and timescale, containing:
#     - correlation: signed Spearman correlation coefficient
#     - significance: sign/significance code
#     - p_value: Spearman correlation p-value
#     - n_obs: number of valid paired observations
#
# Dimension convention:
#   Calculations are performed row-wise as [lat, lon, month] for efficiency.
#   Before writing, arrays are transposed to [lon, lat, month] so that the
#   NetCDF output follows the same dimension order as the input datasets.
#   This incorporates the dimension correction previously applied downstream.
#
# Significance codes:
#    2  positive significant correlation (p <= alpha)
#    1  positive non-significant correlation
#   -1  negative non-significant correlation
#   -2  negative significant correlation (p <= alpha)
#
# Inputs:
#   - kNDVI NetCDF file with variable `kndvi`
#   - Accumulated hydroclimatic NetCDF files named as:
#       {indicator}_GLEAM_v4.2a_MO_1982-2022_scale_{scale}.nc
#
# Outputs:
#   - spearman_correlation_kndvi_{indicator}_scale_{scale}.nc
#
# Notes:
#   - `Ep` is used here as the GLEAM potential evaporation variable representing
#     atmospheric evaporative demand (AED).
#   - `Et` is plant transpiration.
#   - `ED` is evaporation deficit, computed upstream as E - AED for each
#     timescale.
#   - Total evaporation `E` is not analysed as an independent hydroclimatic
#     indicator in the final manuscript and is therefore not included by default.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(parallel)
  library(pracma)
  library(stats)
})

cat("Starting parallel computation on Linux/WSL...\n\n")

# ==============================================================================
# 1. User configuration
# ==============================================================================

# Edit these paths as needed, or replace this block by sourcing R/config.R.
config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

dir_accumulated <- if (exists("paths") && !is.null(paths$accumulated)) {
  paths$accumulated
} else {
  file.path("outputs", "intermediate", "accumulated")
}

file_kndvi <- if (exists("paths") && !is.null(paths$kndvi)) {
  paths$kndvi
} else {
  file.path("outputs", "intermediate", "kndvi", "kndvi.nc")
}

dir_output <- if (exists("paths") && !is.null(paths$correlations_spearman)) {
  paths$correlations_spearman
} else {
  file.path("outputs", "intermediate", "correlations_spearman")
}

# Resolve vegetation and hydroclimate sources from the `sources` registry
# (defined in R/config.R or an override config). Falls back to the
# paths-based defaults above so that the script remains self-contained.
veg_file <- if (exists("sources") && !is.null(sources$veg)) sources$veg$file else file_kndvi
veg_var  <- if (exists("sources") && !is.null(sources$veg)) sources$veg$var  else "kndvi"

hydro_dir      <- if (exists("sources") && !is.null(sources$hydro)) sources$hydro$dir      else dir_accumulated
hydro_file_tpl <- if (exists("sources") && !is.null(sources$hydro)) sources$hydro$file_tpl else "{var}_GLEAM_v4.2a_MO_1982-2022_scale_{scale}.nc"
hydro_var_tpl  <- if (exists("sources") && !is.null(sources$hydro)) sources$hydro$var_tpl  else "{var}"

# Hydroclimatic indicators to analyse.
# Use Ep for AED; Et for transpiration; ED for evaporation deficit.
# indicators <- c("Ep", "Et", "ED", "SMrz", "SMs")
indicators <- c("Ep")
# Accumulation timescales in months.
scales <- c(3, 6, 9, 12)

# Statistical settings.
alpha <- 0.05
minimum_valid_pairs <- 5

# Parallel settings. Forking is memory-efficient on Linux/WSL.
requested_cores <- 24
available_cores <- max(1, parallel::detectCores(logical = FALSE))
ncores <- min(requested_cores, max(1, available_cores - 1))

if (!dir.exists(dir_output)) dir.create(dir_output, recursive = TRUE)

message("Starting monthly Spearman correlation analysis")
message("  Vegetation source: ", if (exists("sources") && !is.null(sources$veg)) sources$veg$name else "kNDVI-GIMMS3g (default)")
message("  Hydro source:      ", if (exists("sources") && !is.null(sources$hydro)) sources$hydro$name else "GLEAM-v4.2a (default)")
message("  Input directory: ", hydro_dir)
message("  Output directory: ", dir_output)
message("  Indicators: ", paste(indicators, collapse = ", "))
message("  Scales: ", paste(scales, collapse = ", "))
message("  Cores: ", ncores)

# ==============================================================================
# 2. Helper functions
# ===============================================================================

read_nc_variable <- function(nc_file, variable_name) {
  if (!file.exists(nc_file)) {
    stop("Input file not found: ", nc_file, call. = FALSE)
  }
  nc <- nc_open(nc_file)
  on.exit(nc_close(nc), add = TRUE)
  ncvar_get(nc, variable_name)
}

linear_detrend <- function(x) {
  # Linear detrend via pracma::detrend.
  # Detrending is applied before removing missing pairs.
  if (all(is.na(x))) return(x)
  pracma::detrend(x, tt = "linear")
}

compute_spearman_pixel <- function(x, y, alpha = 0.05, minimum_valid_pairs = 5) {
  x <- linear_detrend(x)
  y <- linear_detrend(y)
  
  valid <- !is.na(x) & !is.na(y)
  x <- x[valid]
  y <- y[valid]
  n_valid <- length(x)
  
  out <- list(correlation = NA_real_, significance = NA_real_, p_value = NA_real_, n_obs = n_valid)
  
  if (n_valid < minimum_valid_pairs || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(out)
  }
  
  ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate)
  pval <- ct$p.value
  
  sig_code <- if (rho >= 0 && pval <= alpha) {
    2
  } else if (rho >= 0) {
    1
  } else if (rho < 0 && pval <= alpha) {
    -2
  } else {
    -1
  }
  
  list(correlation = rho, significance = sig_code, p_value = pval, n_obs = n_valid)
}

compute_month_by_latitude <- function(i, kndvi_month, indicator_month, lon_dim, alpha, minimum_valid_pairs) {
  corr_row <- rep(NA_real_, lon_dim)
  sig_row  <- rep(NA_real_, lon_dim)
  p_row    <- rep(NA_real_, lon_dim)
  n_row    <- rep(NA_real_, lon_dim)
  
  for (j in seq_len(lon_dim)) {
    res <- compute_spearman_pixel(
      x = kndvi_month[j, i, ],
      y = indicator_month[j, i, ],
      alpha = alpha,
      minimum_valid_pairs = minimum_valid_pairs
    )
    
    corr_row[j] <- res$correlation
    sig_row[j]  <- res$significance
    p_row[j]    <- res$p_value
    n_row[j]    <- res$n_obs
  }
  
  list(corr = corr_row, sig = sig_row, p = p_row, n = n_row)
}

write_correlation_netcdf <- function(output_file, lon, lat, months,
                                     correlation_latlon, significance_latlon,
                                     p_value_latlon, n_obs_latlon,
                                     indicator, scale, alpha,
                                     minimum_valid_pairs) {
  # Convert from [lat, lon, month] to [lon, lat, month].
  correlation_lonlat  <- aperm(correlation_latlon,  c(2, 1, 3))
  significance_lonlat <- aperm(significance_latlon, c(2, 1, 3))
  p_value_lonlat      <- aperm(p_value_latlon,      c(2, 1, 3))
  n_obs_lonlat        <- aperm(n_obs_latlon,        c(2, 1, 3))
  
  dim_lon <- ncdim_def("lon", "degrees_east", lon)
  dim_lat <- ncdim_def("lat", "degrees_north", lat)
  dim_month <- ncdim_def("month", "month", months)
  
  var_cor <- ncvar_def(
    "correlation", "1", list(dim_lon, dim_lat, dim_month), -9999,
    longname = "Spearman rank correlation coefficient between monthly kNDVI and hydroclimatic indicator"
  )
  var_sig <- ncvar_def(
    "significance", "1", list(dim_lon, dim_lat, dim_month), -9999,
    longname = "Correlation sign and significance code"
  )
  var_p <- ncvar_def(
    "p_value", "1", list(dim_lon, dim_lat, dim_month), -9999,
    longname = "P-value of Spearman rank correlation"
  )
  var_n <- ncvar_def(
    "n_obs", "1", list(dim_lon, dim_lat, dim_month), -9999,
    longname = "Number of valid paired observations"
  )
  
  nc <- nc_create(output_file, list(var_cor, var_sig, var_p, var_n), force_v4 = TRUE)
  on.exit(nc_close(nc), add = TRUE)
  
  ncvar_put(nc, var_cor, correlation_lonlat)
  ncvar_put(nc, var_sig, significance_lonlat)
  ncvar_put(nc, var_p, p_value_lonlat)
  ncvar_put(nc, var_n, n_obs_lonlat)
  
  ncatt_put(nc, 0, "title", "Monthly Spearman correlations between kNDVI and hydroclimatic indicators")
  ncatt_put(nc, 0, "project", "Aridity and timescale bound the separability of vegetation water stress")
  ncatt_put(nc, 0, "indicator", indicator)
  ncatt_put(nc, 0, "accumulation_window_months", scale)
  ncatt_put(nc, 0, "analysis_period", "1982-2022")
  ncatt_put(nc, 0, "calendar_months", "1=January, ..., 12=December")
  ncatt_put(nc, 0, "method", "Spearman rank correlation after linear detrending of both variables")
  ncatt_put(nc, 0, "alpha", alpha)
  ncatt_put(nc, 0, "minimum_valid_pairs", minimum_valid_pairs)
  ncatt_put(nc, 0, "dimension_order", "lon, lat, month")
  ncatt_put(nc, 0, "significance_codes", "2 positive significant; 1 positive non-significant; -1 negative non-significant; -2 negative significant")
  ncatt_put(nc, 0, "created_by", basename(sys.frame(1)$ofile %||% "01_compute_monthly_spearman_correlations.R"))
  ncatt_put(nc, 0, "date_created", as.character(Sys.time()))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# ==============================================================================
# 3. Load kNDVI
# ===============================================================================

message("Loading vegetation index: ", veg_var, " from ", veg_file, "...")
kndvi_nc <- nc_open(veg_file)
kndvi_data <- ncvar_get(kndvi_nc, veg_var)
lon <- ncvar_get(kndvi_nc, "lon")
lat <- ncvar_get(kndvi_nc, "lat")
nc_close(kndvi_nc)

lon_dim <- length(lon)
lat_dim <- length(lat)
time_dim <- dim(kndvi_data)[3]

message("Loaded dimensions:")
message("  lon_dim = ", lon_dim)
message("  lat_dim = ", lat_dim)
message("  time_dim = ", time_dim)

if (time_dim %% 12 != 0) {
  warning("The kNDVI time dimension is not a multiple of 12. Check monthly indexing.")
}

# ==============================================================================
# 4. Main processing loop
# ===============================================================================


for (indicator in indicators) {
  message("\nProcessing indicator: ", indicator)
  
  for (scale in scales) {
    message("  Accumulation timescale: ", scale, " month(s)")
    
    fname          <- gsub("\\{scale\\}", scale,
                           gsub("\\{var\\}", indicator, hydro_file_tpl))
    indicator_file <- file.path(hydro_dir, fname)
    var_name       <- gsub("\\{var\\}", indicator, hydro_var_tpl)

    indicator_data <- read_nc_variable(indicator_file, var_name)
    
    if (!all(dim(indicator_data)[1:2] == c(lon_dim, lat_dim))) {
      stop(
        "Spatial dimensions do not match kNDVI for ", indicator_file,
        ". Expected [", lon_dim, ", ", lat_dim, ", time], got [",
        paste(dim(indicator_data), collapse = ", "), "].",
        call. = FALSE
      )
    }
    
    if (dim(indicator_data)[3] != time_dim) {
      stop(
        "Time dimension does not match kNDVI for ", indicator_file,
        ". Expected ", time_dim, ", got ", dim(indicator_data)[3], ".",
        call. = FALSE
      )
    }
    
    # Arrays are stored during calculation as [lat, lon, month].
    correlation  <- array(NA_real_, c(lat_dim, lon_dim, 12))
    significance <- array(NA_real_, c(lat_dim, lon_dim, 12))
    p_value      <- array(NA_real_, c(lat_dim, lon_dim, 12))
    n_obs        <- array(NA_real_, c(lat_dim, lon_dim, 12))
    
    for (month in 1:12) {
      message("    Calendar month: ", month)
      
      month_indices <- seq(month, time_dim, by = 12)
      kndvi_month <- kndvi_data[, , month_indices]
      indicator_month <- indicator_data[, , month_indices]
      
      if (.Platform$OS.type == "unix" && ncores > 1) {
        results <- parallel::mclapply(
          X = seq_len(lat_dim),
          FUN = compute_month_by_latitude,
          kndvi_month = kndvi_month,
          indicator_month = indicator_month,
          lon_dim = lon_dim,
          alpha = alpha,
          minimum_valid_pairs = minimum_valid_pairs,
          mc.cores = ncores
        )
      } else {
        results <- lapply(
          X = seq_len(lat_dim),
          FUN = compute_month_by_latitude,
          kndvi_month = kndvi_month,
          indicator_month = indicator_month,
          lon_dim = lon_dim,
          alpha = alpha,
          minimum_valid_pairs = minimum_valid_pairs
        )
      }
      
      for (i in seq_len(lat_dim)) {
        correlation[i, , month]  <- results[[i]]$corr
        significance[i, , month] <- results[[i]]$sig
        p_value[i, , month]      <- results[[i]]$p
        n_obs[i, , month]        <- results[[i]]$n
      }
      
      rm(kndvi_month, indicator_month, results)
      gc(verbose = FALSE)
    }
    
    output_file <- file.path(
      dir_output,
      paste0("spearman_correlation_kndvi_", indicator, "_scale_", scale, ".nc")
    )
    
    write_correlation_netcdf(
      output_file = output_file,
      lon = lon,
      lat = lat,
      months = 1:12,
      correlation_latlon = correlation,
      significance_latlon = significance,
      p_value_latlon = p_value,
      n_obs_latlon = n_obs,
      indicator = indicator,
      scale = scale,
      alpha = alpha,
      minimum_valid_pairs = minimum_valid_pairs
    )
    
    message("  Saved: ", output_file)
    rm(indicator_data, correlation, significance, p_value, n_obs)
    gc(verbose = FALSE)
  }
}

message("\nMonthly Spearman correlation analysis completed.")
