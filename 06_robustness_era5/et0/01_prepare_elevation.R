# ==============================================================================
# 06_robustness_era5/et0/01_prepare_elevation.R
#
# Build the elevation grid (z, in metres) required by the FAO-56 pressure and
# clear-sky-radiation terms, from the ERA5-Land invariant geopotential field.
#
# Rationale: the elevation must be consistent with the model orography on which
# the meteorological variables were produced, and must sit on exactly the same
# 0.1 deg ERA5-Land grid. The invariant geopotential satisfies both, with no
# reprojection or resampling.
#
# Input : ERA5-Land invariant geopotential NetCDF (path in R/config_era5.R).
# Output: elevation_era5land.rds  -- a [longitude x latitude] matrix of metres,
#         aligned cell-by-cell with the meteorological files.
#
# Run from the repository root:  Rscript 06_robustness_era5/et0/01_prepare_elevation.R
# ==============================================================================

library(ncdf4)

here <- file.path("06_robustness_era5", "et0")
source(file.path(here, "fao56_functions.R"))

source(file.path("R", "config_era5.R"))   # canonical paths

geopotential_nc <- era5$geopotential
out_rds         <- file.path(here, "elevation_era5land.rds")

g <- 9.80665                                       # standard gravity [m s-2]

# ------------------------------------------------------------------------------
nc  <- nc_open(geopotential_nc)
var <- data_vars(nc)[1]                            # usually "z" (units m2 s-2)
message("Reading geopotential variable: ", var)

# Geopotential is time-invariant; read the first time slice if a time dim exists.
d <- nc$var[[var]]$size
if (length(d) == 3) {
  geop <- ncvar_get(nc, var, start = c(1, 1, 1), count = c(-1, -1, 1))
} else {
  geop <- ncvar_get(nc, var)
}
nc_close(nc)

# This file carries _FillValue = -32767; ncdf4 usually maps it to NA, but guard
# explicitly so ocean cells never leak in as spurious elevations.
geop[!is.finite(geop) | geop <= -32760] <- NA

elevation <- geop / g                              # geopotential -> metres (orography)
# Genuine below-sea-level land (Dead Sea, Caspian, etc.) is kept: it is
# physically correct for the atmospheric-pressure term. Only ocean stays NA.

message(sprintf("Elevation grid: %d x %d cells | range [%.0f, %.0f] m | %.1f%% land",
                nrow(elevation), ncol(elevation),
                min(elevation, na.rm = TRUE), max(elevation, na.rm = TRUE),
                100 * mean(!is.na(elevation))))

saveRDS(elevation, out_rds)
message("Saved: ", out_rds)
