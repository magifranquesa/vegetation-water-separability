#!/usr/bin/env Rscript

# ==============================================================================
# Script: 04_analysis/05_pc1_variance_explained.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   Extract the EXACT within-block PC1 variance-explained fraction used by the
#   global variance partitioning, for the Methods sentence:
#     "PC1 captured X% (supply) and Y% (demand) of within-block variance".
#
#   It replicates the per-cell, per-calendar-month block PCA of
#   03_variance_partitioning/01_varpart.R (reduce_block): same linear detrend
#   (pracma::detrend, tt="linear"), same prcomp(center=TRUE, scale.=TRUE), same
#   block columns:
#       supply (SUMINISTRO) = SMs + SMrz, 5 accumulation scales -> 10 series
#       demand (DEMANDA)    = Ep,         5 accumulation scales ->  5 series
#   The ONLY difference: instead of taking PC1 scores for the partition, it
#   records the proportion of variance carried by PC1 (and PC1+PC2).
#
#   *** It does NOT recompute the partition. Inputs are READ-ONLY. ***
#
#   Sampling: a latitude-band sample over the vegetated domain (all vegetated
#   cells within several representative bands, capped per band), × 12 calendar
#   months. Area weighting uses cos(lat), which is proportional to the EPSG:4326
#   geodesic cell area on a regular lon/lat grid (so the area-weighted MEAN is
#   identical to a terra::cellSize-weighted mean).
#
# Reports per block (supply, demand):
#   - median and p5–p95 of the PC1 variance fraction across cell×month
#   - area-weighted mean (cos-lat weights)
#   - PC1+PC2 fraction (context)
#
# Inputs (READ-ONLY):
#   outputs/intermediate/kndvi/kndvi.nc            (grid axes only)
#   outputs/intermediate/vegetation_mask_c1.tif    (vegetated cells)
#   outputs/intermediate/accumulated/{var}_GLEAM_v4.2a_MO_1982-2022_scale_{sc}.nc
#
# Output:
#   outputs/varpart_global/pc1_variance_explained.csv
#
# Run from repo root:
#   Rscript 04_analysis/05_pc1_variance_explained.R
#
# Dependencies: ncdf4, terra, pracma
# ==============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(terra)
  library(pracma)
})

t_start <- proc.time()

# ------------------------------------------------------------------------------
# Configuration (block/scale definitions)
# ------------------------------------------------------------------------------

SCALES      <- c(1, 3, 6, 9, 12)
BLOCKS      <- list(SUMINISTRO = c("SMs", "SMrz"), DEMANDA = "Ep")
BLOCK_NAMES <- names(BLOCKS)
N_TIME      <- 492L          # 1982–2022: 41 yr × 12 mo
file_tpl    <- "{var}_GLEAM_v4.2a_MO_1982-2022_scale_{scale}.nc"

# Sampling controls
BAND_DEG       <- 3                                   # latitude band thickness
BAND_CENTERS   <- c(-50, -35, -22, -10, 2, 15, 28, 42, 55, 68)
CELLS_PER_BAND <- 5000L
set.seed(42L)

# Absolute repo root so the script works regardless of the session's working
# directory (e.g. when run via source() from the console).
project_root <- "."

cfg <- file.path(project_root, "R", "config.R")
if (file.exists(cfg)) source(cfg)

# Use a config path only if it actually exists locally; otherwise the default.
pick <- function(p, default) if (!is.null(p) && file.exists(p)) p else default

file_kndvi    <- pick(if (exists("paths")) paths$kndvi else NULL,
                      file.path(project_root, "outputs", "intermediate", "kndvi", "kndvi.nc"))
dir_accum     <- pick(if (exists("paths")) paths$accumulated else NULL,
                      file.path(project_root, "outputs", "intermediate", "accumulated"))
file_veg_mask <- pick(if (exists("paths")) paths$veg_mask_c1 else NULL,
                      file.path(project_root, "outputs", "intermediate", "vegetation_mask_c1.tif"))

out_csv <- file.path(project_root, "outputs", "varpart_global", "pc1_variance_explained.csv")
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(file_kndvi), file.exists(file_veg_mask), dir.exists(dir_accum))

# mtime record (read-only protection)
series_keys <- unlist(lapply(BLOCK_NAMES, function(bn)
  unlist(lapply(BLOCKS[[bn]], function(ind) paste0(ind, "_s", SCALES)))))
series_files <- setNames(lapply(BLOCK_NAMES, function(bn)
  do.call(c, lapply(BLOCKS[[bn]], function(ind)
    setNames(file.path(dir_accum,
      vapply(SCALES, function(sc)
        gsub("\\{scale\\}", sc, gsub("\\{var\\}", ind, file_tpl)), character(1))),
      paste0(ind, "_s", SCALES))))), BLOCK_NAMES)
all_series_files <- do.call(c, series_files)
mt_before <- file.info(c(file_kndvi, file_veg_mask, as.character(all_series_files)))$mtime

# ------------------------------------------------------------------------------
# Helpers (detrend)
# ------------------------------------------------------------------------------

detrend_series <- function(x) {
  if (sum(!is.na(x)) < 5L) return(rep(NA_real_, length(x)))
  pracma::detrend(x, tt = "linear")
}

# Same complete.cases / sd>0 / prcomp(scale.) logic as reduce_block, but returns
# the variance fractions instead of the PC1 scores.
block_pca_var <- function(mat) {
  ok <- complete.cases(mat)
  if (sum(ok) < 10L) return(NULL)
  mat_ok <- mat[ok, , drop = FALSE]
  col_sd <- apply(mat_ok, 2L, sd)
  keep   <- !is.na(col_sd) & col_sd > 0
  if (!any(keep)) return(NULL)
  mat_ok <- mat_ok[, keep, drop = FALSE]
  if (ncol(mat_ok) == 0L) return(NULL)
  pca <- prcomp(mat_ok, center = TRUE, scale. = TRUE)
  v   <- pca$sdev^2
  tot <- sum(v)
  c(pc1  = v[1L] / tot,
    pc12 = sum(v[seq_len(min(2L, length(v)))]) / tot,
    ncol = ncol(mat_ok))
}

# ------------------------------------------------------------------------------
# Reference grid + vegetation mask
# ------------------------------------------------------------------------------

message("Loading grid axes and vegetation mask...")
nc0     <- nc_open(file_kndvi)
lon_all <- ncvar_get(nc0, "lon")
lat_all <- ncvar_get(nc0, "lat")
nc_close(nc0)
lon_n <- length(lon_all)

r_mask <- terra::rast(file_veg_mask)
mask_lat_range   <- c(terra::ext(r_mask)$ymin, terra::ext(r_mask)$ymax)
domain_rows      <- which(lat_all >= mask_lat_range[1] & lat_all <= mask_lat_range[2])
lat_domain_start <- min(domain_rows)
lat_domain_end   <- max(domain_rows)
lat_domain_count <- lat_domain_end - lat_domain_start + 1L
lat_domain       <- lat_all[lat_domain_start:lat_domain_end]

arr_terra <- as.array(r_mask)[, , 1]
if (lat_domain[1] < lat_domain[2]) {
  veg_mask <- t(arr_terra[nrow(arr_terra):1, ])
} else {
  veg_mask <- t(arr_terra)
}
veg_mask[is.na(veg_mask)] <- 0L
stopifnot(nrow(veg_mask) == lon_n, ncol(veg_mask) == lat_domain_count)
message(sprintf("  Domain lat n=%d [%.2f, %.2f]; vegetated cells=%d",
                lat_domain_count, min(lat_domain), max(lat_domain), sum(veg_mask > 0)))

# ------------------------------------------------------------------------------
# Process latitude bands
# ------------------------------------------------------------------------------

results <- list()

for (bc in BAND_CENTERS) {
  lo <- bc - BAND_DEG / 2
  hi <- bc + BAND_DEG / 2
  jdx <- which(lat_domain >= lo & lat_domain <= hi)        # domain lat indices
  if (length(jdx) == 0L) { message(sprintf("  band %+0.1f: no domain rows", bc)); next }

  veg_sub <- veg_mask[, jdx, drop = FALSE]
  vidx    <- which(veg_sub > 0, arr.ind = TRUE)            # col1=li, col2=local j
  if (nrow(vidx) == 0L) { message(sprintf("  band %+0.1f: no vegetated cells", bc)); next }

  if (nrow(vidx) > CELLS_PER_BAND)
    vidx <- vidx[sample.int(nrow(vidx), CELLS_PER_BAND), , drop = FALSE]

  li_vec   <- vidx[, 1L]
  j_glob   <- jdx[vidx[, 2L]]                              # domain index
  lat_vec  <- lat_domain[j_glob]
  row0     <- lat_domain_start + min(j_glob) - 1L          # nc global start row
  loc_vec  <- (lat_domain_start + j_glob - 1L) - row0 + 1L # local row in slab
  n_rows   <- (max(j_glob) - min(j_glob)) + 1L
  n_cell   <- length(li_vec)
  message(sprintf("  band %+5.1f° (lat %.1f..%.1f): %d cells, slab rows %d..%d",
                  bc, min(lat_vec), max(lat_vec), n_cell, row0, row0 + n_rows - 1L))

  # Read each of the 15 series for this band, extract sampled-cell time series.
  S <- vector("list", length(series_keys)); names(S) <- series_keys
  for (bn in BLOCK_NAMES) for (ind in BLOCKS[[bn]]) for (sc in SCALES) {
    key  <- paste0(ind, "_s", sc)
    ncg  <- nc_open(as.character(series_files[[bn]][key]))
    slab <- ncvar_get(ncg, ind, start = c(1L, row0, 1L),
                      count = c(lon_n, n_rows, N_TIME))
    nc_close(ncg)
    S[[key]] <- t(vapply(seq_len(n_cell),
                         function(k) slab[li_vec[k], loc_vec[k], ],
                         numeric(N_TIME)))
    rm(slab)
  }
  gc(verbose = FALSE)

  # Per cell × calendar month: block PCA variance fractions.
  for (k in seq_len(n_cell)) {
    for (mo in 1:12) {
      midx <- seq(mo, N_TIME, by = 12L)
      sup_mat <- do.call(cbind, lapply(BLOCKS$SUMINISTRO, function(ind)
        do.call(cbind, lapply(SCALES, function(sc)
          detrend_series(S[[paste0(ind, "_s", sc)]][k, midx])))))
      dem_mat <- do.call(cbind, lapply(BLOCKS$DEMANDA, function(ind)
        do.call(cbind, lapply(SCALES, function(sc)
          detrend_series(S[[paste0(ind, "_s", sc)]][k, midx])))))

      vs <- block_pca_var(sup_mat)
      vd <- block_pca_var(dem_mat)
      if (!is.null(vs))
        results[[length(results) + 1L]] <- data.frame(
          block = "supply", lat = lat_vec[k],
          pc1 = vs["pc1"], pc12 = vs["pc12"], ncol = vs["ncol"])
      if (!is.null(vd))
        results[[length(results) + 1L]] <- data.frame(
          block = "demand", lat = lat_vec[k],
          pc1 = vd["pc1"], pc12 = vd["pc12"], ncol = vd["ncol"])
    }
  }
  rm(S); gc(verbose = FALSE)
}

df <- do.call(rbind, results)
rownames(df) <- NULL
message(sprintf("\nCollected %d block×cell×month PCA fits.", nrow(df)))

# ------------------------------------------------------------------------------
# Summaries per block
# ------------------------------------------------------------------------------

NCOL_FULL <- c(supply = 10L, demand = 5L)

summ <- function(b) {
  d  <- df[df$block == b, ]
  # Restrict to non-degenerate fits (>=2 retained columns) so the within-block
  # PC1 fraction is well defined; report how many were dropped.
  d2 <- d[d$ncol >= 2L, ]
  w  <- cos(d2$lat * pi / 180)
  q  <- quantile(d2$pc1, c(0.05, 0.50, 0.95), na.rm = TRUE)
  q12<- quantile(d2$pc12, c(0.05, 0.50, 0.95), na.rm = TRUE)
  data.frame(
    block          = b,
    n              = nrow(d2),
    n_degenerate   = sum(d$ncol < 2L),
    pct_full_ncol  = round(mean(d$ncol == NCOL_FULL[b]) * 100, 1),
    pc1_median     = q[2], pc1_p5 = q[1], pc1_p95 = q[3],
    pc1_area_mean  = sum(w * d2$pc1) / sum(w),
    pc12_median    = q12[2], pc12_area_mean = sum(w * d2$pc12) / sum(w)
  )
}

out <- rbind(summ("supply"), summ("demand"))
rownames(out) <- NULL
write.csv(out, out_csv, row.names = FALSE)

cat("\n================================================================\n")
cat("=== Within-block PC1 variance explained (EXACT partition PCA)\n")
cat("================================================================\n\n")
for (i in seq_len(nrow(out))) {
  r <- out[i, ]
  cat(sprintf("  %-7s | n=%d (full-ncol %.1f%%, degenerate dropped=%d)\n",
              r$block, r$n, r$pct_full_ncol, r$n_degenerate))
  cat(sprintf("          PC1 : median %.1f%%  (p5–p95 %.1f–%.1f%%)  | area-weighted mean %.1f%%\n",
              r$pc1_median*100, r$pc1_p5*100, r$pc1_p95*100, r$pc1_area_mean*100))
  cat(sprintf("          PC1+PC2: median %.1f%%  | area-weighted mean %.1f%%\n\n",
              r$pc12_median*100, r$pc12_area_mean*100))
}

cat("  >>> Methods sentence (area-weighted mean):\n")
cat(sprintf("      \"PC1 captured %.0f%% (supply) and %.0f%% (demand) of within-block variance\"\n",
            out$pc1_area_mean[out$block=="supply"]*100,
            out$pc1_area_mean[out$block=="demand"]*100))
cat(sprintf("  CSV: %s\n", out_csv))

# ------------------------------------------------------------------------------
# File-protection confirmation
# ------------------------------------------------------------------------------

mt_after <- file.info(c(file_kndvi, file_veg_mask, as.character(all_series_files)))$mtime
cat("\n=== File protection ===\n")
cat(sprintf("  All %d input files UNCHANGED: %s\n",
            length(mt_before),
            if (identical(mt_before, mt_after)) "YES (read-only respected)" else "*** NO — CHECK ***"))

cat(sprintf("\nElapsed: %.0f s\n", (proc.time() - t_start)["elapsed"]))
cat("Done.\n")
