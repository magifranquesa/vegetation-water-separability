# ==============================================================================
# 06_robustness_era5/et0/fao56_functions.R
#
# FAO-56 Penman-Monteith reference evapotranspiration (ET0), monthly step.
# Reference: Allen, Pereira, Raes & Smith (1998), FAO Irrigation and Drainage
# Paper 56. Equation numbers below refer to that document.
#
# All functions are vectorised: arguments may be scalars, vectors or matrices,
# as long as their shapes are mutually compatible. Units follow FAO-56:
#   temperature [degC], vapour pressure [kPa], radiation [MJ m-2 day-1],
#   wind [m s-1], elevation [m], resulting ET0 [mm day-1].
# ==============================================================================

# Saturation vapour pressure [kPa] at air temperature T [degC]        (eq. 11)
svp <- function(T) 0.6108 * exp(17.27 * T / (T + 237.3))

# Slope of the saturation vapour pressure curve [kPa degC-1]          (eq. 13)
svp_slope <- function(Tmean) 4098 * svp(Tmean) / (Tmean + 237.3)^2

# Atmospheric pressure [kPa] from elevation z [m]                     (eq. 7)
atm_pressure <- function(z) 101.3 * ((293 - 0.0065 * z) / 293)^5.26

# Psychrometric constant [kPa degC-1] from pressure P [kPa]           (eq. 8)
psy_constant <- function(P) 0.000665 * P

# Convert wind speed measured at 10 m to the 2 m reference height     (eq. 47)
wind_10m_to_2m <- function(u10, v10) {
  sqrt(u10^2 + v10^2) * (4.87 / log(67.8 * 10 - 5.42))   # factor = 0.748
}

# Extraterrestrial radiation Ra [MJ m-2 day-1]                        (eq. 21)
#   lat_rad : latitude in radians (negative in the southern hemisphere)
#   J       : day of the year (1-365), use the mid-month day for monthly ET0
extraterrestrial_radiation <- function(lat_rad, J) {
  Gsc <- 0.0820                                    # solar constant [MJ m-2 min-1]
  dr  <- 1 + 0.033 * cos(2 * pi * J / 365)         # inverse Earth-Sun distance (eq. 23)
  dec <- 0.409 * sin(2 * pi * J / 365 - 1.39)      # solar declination          (eq. 24)
  ws  <- acos(pmin(pmax(-tan(lat_rad) * tan(dec), -1), 1))  # sunset hour angle  (eq. 25)
  (24 * 60 / pi) * Gsc * dr *
    (ws * sin(lat_rad) * sin(dec) + cos(lat_rad) * cos(dec) * sin(ws))
}

# Clear-sky solar radiation Rso [MJ m-2 day-1]                        (eq. 37)
clear_sky_radiation <- function(Ra, z) (0.75 + 2e-5 * z) * Ra

# Net radiation Rn [MJ m-2 day-1]                             (eqs. 38, 39, 40)
net_radiation <- function(Rs, Rso, Tmax, Tmin, ea, albedo = 0.23) {
  Rns   <- (1 - albedo) * Rs                       # net shortwave              (eq. 38)
  sigma <- 4.903e-9                                # Stefan-Boltzmann [MJ K-4 m-2 day-1]
  TmaxK <- Tmax + 273.16
  TminK <- Tmin + 273.16
  ratio <- pmin(Rs / Rso, 1)                       # relative shortwave, capped at 1
  Rnl   <- sigma * ((TmaxK^4 + TminK^4) / 2) *
    (0.34 - 0.14 * sqrt(pmax(ea, 0))) *
    (1.35 * ratio - 0.35)                          # net longwave               (eq. 39)
  Rns - Rnl
}

# FAO-56 Penman-Monteith reference ET0 [mm day-1]                     (eq. 6)
#   Tmax, Tmin, Tdew : monthly-mean daily max / min / dewpoint temperature [degC]
#   u2               : wind speed at 2 m [m s-1]
#   Rs               : incoming solar radiation [MJ m-2 day-1]
#   Ra               : extraterrestrial radiation [MJ m-2 day-1]
#   z                : elevation [m]
#   G                : soil heat flux [MJ m-2 day-1]; 0 is standard for monthly ET0
et0_penman_monteith <- function(Tmax, Tmin, Tdew, u2, Rs, Ra, z,
                                G = 0, albedo = 0.23) {
  Tmean <- (Tmax + Tmin) / 2
  delta <- svp_slope(Tmean)
  P     <- atm_pressure(z)
  gamma <- psy_constant(P)
  es    <- (svp(Tmax) + svp(Tmin)) / 2             # saturation vp             (eq. 12)
  ea    <- svp(Tdew)                               # actual vp from dewpoint   (eq. 14)
  Rso   <- clear_sky_radiation(Ra, z)
  Rn    <- net_radiation(Rs, Rso, Tmax, Tmin, ea, albedo)

  num <- 0.408 * delta * (Rn - G) +
    gamma * (900 / (Tmean + 273)) * u2 * (es - ea)
  den <- delta + gamma * (1 + 0.34 * u2)
  et0 <- num / den
  et0[et0 < 0] <- 0                                # ET0 is non-negative by definition
  et0
}

# ------------------------------------------------------------------------------
# NetCDF helpers
# ------------------------------------------------------------------------------

# Names of the 3-D data variables (valid_time, latitude, longitude) in an
# open ncdf4 object, ignoring bookkeeping variables such as `number`/`expver`.
data_vars <- function(nc) {
  vars <- names(nc$var)
  vars[vapply(vars, function(v) nc$var[[v]]$ndims == 3, logical(1))]
}

# Number of days in a given calendar month.
days_in_month <- function(year, month) {
  first <- as.Date(sprintf("%04d-%02d-01", year, month))
  nextm <- seq(first, by = "month", length.out = 2)[2]
  as.integer(nextm - first)
}
