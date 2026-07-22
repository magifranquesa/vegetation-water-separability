# ==============================================================================
# R/config_era5.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Canonical paths and metadata for the consolidated ERA5-Land raw data.
# Source this file instead of hard-coding paths, so that a future move needs a
# single edit here.
#
#   source(file.path("R", "config_era5.R"))
#
# See data/raw/ERA5_land/MANIFEST.md for why each file uses the product it does.
# ==============================================================================

era5_dir <- file.path("data", "raw", "ERA5_land")

era5 <- list(
  # --- accumulated fluxes: product 'monthly_averaged_reanalysis_by_hour_of_day'
  #     at 00:00 (step 24 h). The plain monthly product is corrupted for these
  #     variables from September 2022 (ECMWF known issue 2024-01).
  transpiration = file.path(era5_dir, "era5land_transpiration_mnth_1981-2022.nc"),
  totalevap     = file.path(era5_dir, "era5land_totalevap_mnth_1981-2022.nc"),
  ssrd          = file.path(era5_dir, "era5land_ssrd_mnth_1981-2022.nc"),

  # --- instantaneous states: product 'monthly_averaged_reanalysis'. The by-hour
  #     product at 00:00 would give only the 00 UTC snapshot, not a monthly mean.
  soillayers    = file.path(era5_dir, "era5land_soillayers_moda_1981-2022.nc"),
  dewpoint      = file.path(era5_dir, "era5land_dewpoint_moda_1981-2022.nc"),
  wind          = file.path(era5_dir, "era5land_wind_moda_1981-2022.nc"),

  # --- derived monthly means of daily extremes. Originally 1980-2022 (516 months);
  #     the leading year was trimmed with CDO so that they now start in January 1981
  #     like every other input. Note the time axis still carries the 1980 reference
  #     origin ("days since 1980-01-01"), which is only the epoch, not the first month.
  tmin          = file.path(era5_dir, "era5land_tmin_1981-2022.nc"),
  tmax          = file.path(era5_dir, "era5land_tmax_1981-2022.nc"),

  # --- time-invariant
  geopotential  = file.path(era5_dir, "era5land_geopotential.nc")
)

# Variable name INSIDE each file. Renaming the files did not rename the variables.
# 'transpiration' still holds the variable 'evabs': ERA5-Land has a documented
# known issue by which evabs / evaow / evavt hold each other's values, so the
# variable labelled "evaporation from bare soil" actually contains transpiration.
# Verified physically in 06_robustness_era5/checks/verify_evabs.R.
era5_var <- list(
  transpiration = "evabs",
  totalevap     = "e",
  ssrd          = "ssrd",
  soillayers    = c("swvl1", "swvl2", "swvl3", "swvl4"),
  dewpoint      = "d2m",
  wind          = c("u10", "v10"),
  tmin          = "t2m",          # both temperature files carry the name 't2m'
  tmax          = "t2m",
  geopotential  = "z"
)

# All files now start in January 1981, so month t is at index t everywhere: no offset.
# Kept as an explicit table (rather than dropped) because the temperature files DID
# start in 1980 until they were trimmed, and reading them one year out is an error that
# raises no exception: the seasonal cycle would look perfect and only the interannual
# variability -- the very thing the analysis uses -- would be wrong. The consumers assert
# the alignment against these values, so a future file with a different period is caught.
era5_t_offset <- list(
  transpiration = 0L, totalevap = 0L, ssrd = 0L,
  soillayers = 0L, dewpoint = 0L, wind = 0L,
  tmin = 0L, tmax = 0L, geopotential = 0L
)

# Common grid: 3600 lon (0 .. 359.9) x 1801 lat (90 .. -90), node-registered.
# The GLEAM/kNDVI target grid is 3600 x 1800, cell-centred, -180/180.
ERA5_NLON <- 3600L
ERA5_NLAT <- 1801L
ERA5_NTIME <- 504L      # 1981-01 .. 2022-12, the series every script works on
