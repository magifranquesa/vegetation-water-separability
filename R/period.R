# ==============================================================================
# R/period.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Period-window helper shared by the variance-partitioning scripts.
#
# The analysis normally runs over the full record (1982-2022). Setting a period
# restricts every per-calendar-month series to the years of that window. This
# also makes the linear detrending WINDOW-LOCAL, which is required when
# comparing windows: if both forcings dry over time, a trend in the means would
# inflate the correlation between levels without any change in the coupling of
# anomalies.
#
# Two ways to set it (the environment variable wins, so the committed default
# always reproduces the published full-record results):
#
#   VWS_PERIOD=1982-2001 Rscript 03_variance_partitioning/01_varpart.R
#   VWS_PERIOD=2002-2022 Rscript 03_variance_partitioning/01_varpart.R
#
#   ...or edit PERIOD_DEFAULT below.
#
# Every output path built with PERIOD_SFX is kept separate from the full-record
# results. This matters beyond tidiness: 01_varpart.R caches latitude stripes as
# RDS and SKIPS any stripe whose file already exists. Without a distinct suffix
# a window run would silently reuse the full-record stripes and report the wrong
# numbers with no error.
# ==============================================================================

PERIOD_DEFAULT <- NULL       # NULL = full record. Do not commit a window here.

YEAR_START <- 1982L
YEAR_END   <- 2022L

.vws_env_period <- Sys.getenv("VWS_PERIOD")

PERIOD <- if (nzchar(.vws_env_period)) {
  .p <- suppressWarnings(as.integer(strsplit(.vws_env_period, "[-_:to]+")[[1]]))
  .p <- .p[!is.na(.p)]
  if (length(.p) != 2L)
    stop("VWS_PERIOD must look like '1982-2001'; got: ", .vws_env_period, call. = FALSE)
  .p
} else {
  PERIOD_DEFAULT
}

if (!is.null(PERIOD)) {
  if (PERIOD[1] < YEAR_START || PERIOD[2] > YEAR_END)
    stop(sprintf("PERIOD %d-%d falls outside the record %d-%d",
                 PERIOD[1], PERIOD[2], YEAR_START, YEAR_END), call. = FALSE)
  if (PERIOD[1] >= PERIOD[2])
    stop(sprintf("PERIOD must be increasing; got %d-%d", PERIOD[1], PERIOD[2]),
         call. = FALSE)
}

# Suffix appended to every output path. Empty for the full record, so existing
# filenames and downstream figure scripts are untouched.
PERIOD_SFX <- if (is.null(PERIOD)) "" else sprintf("_%d_%d", PERIOD[1], PERIOD[2])

# Indices (1-based, into the year axis) of the years kept by the current window.
period_year_index <- function() {
  if (is.null(PERIOD)) {
    seq_len(YEAR_END - YEAR_START + 1L)
  } else {
    (PERIOD[1] - YEAR_START + 1L):(PERIOD[2] - YEAR_START + 1L)
  }
}

period_n_years <- function() length(period_year_index())

period_label <- function() {
  if (is.null(PERIOD)) sprintf("%d-%d (full record)", YEAR_START, YEAR_END)
  else sprintf("%d-%d", PERIOD[1], PERIOD[2])
}

message(sprintf("[period] window: %s  (%d years)  suffix: '%s'",
                period_label(), period_n_years(), PERIOD_SFX))
