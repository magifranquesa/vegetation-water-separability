# ==============================================================================
# R/config.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Central path configuration for all R scripts in the project.
# All scripts source this file if it exists; otherwise they fall back to the
# same relative paths defined here.
#
# To use on a remote server with a different directory layout, set absolute
# paths below and source this file before running any script.
#
# Run all scripts from the repository root.
# ==============================================================================

paths <- list(

  # --------------------------------------------------------------------------
  # Raw and processed input data
  # --------------------------------------------------------------------------

  # Mean aridity index NetCDF (distributed with the repository)
  aridity_index  = file.path("data", "processed", "aridity", "ai_1982_2022_period.nc"),

  # Directory for aridity outputs (ai_1982_2022_period.nc, Zomer classes TIF)
  aridity        = file.path("data", "processed", "aridity"),

  # World administrative boundaries in Equal Earth projection
  world_equal_earth   = file.path("data", "external", "admin_equal_earth_clean.shp"),

  # --------------------------------------------------------------------------
  # Intermediate outputs
  # --------------------------------------------------------------------------

  # Root of the intermediate outputs directory
  intermediate          = file.path("outputs", "intermediate"),

  # kNDVI resampled to the GLEAM 0.1 degree grid
  kndvi                 = file.path("outputs", "intermediate", "kndvi", "kndvi.nc"),

  # Concatenated GLEAM Et time series (monthly 1981-2022)
  et_gleam              = file.path("outputs", "intermediate", "concatenated",
                                    "Et_GLEAM_v4.2a_MO_1981-2022_NCO.nc"),

  # Vegetated domain mask C1 (Et > 0 any month AND mean kNDVI > 0.025)
  veg_mask_c1           = file.path("outputs", "intermediate",
                                    "vegetation_mask_c1.tif"),

  # Accumulated hydroclimatic variables (5 indicators x 5 timescales)
  accumulated            = file.path("outputs", "intermediate", "accumulated"),

  # Monthly Spearman correlations (25 NetCDF files)
  correlations_spearman = file.path("outputs", "intermediate", "correlations_spearman"),

  # Absolute-maximum Spearman correlation across timescales (rho*)
  absmax_spearman       = file.path("outputs", "intermediate", "absmax_spearman"),

  # Maximum positive Spearman correlation across timescales
  max_spearman          = file.path("outputs", "intermediate", "max_spearman"),

  # --------------------------------------------------------------------------
  # Final outputs
  # --------------------------------------------------------------------------

  figures = file.path("outputs", "figures"),
  tables  = file.path("outputs", "tables")
)
