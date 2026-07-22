#!/usr/bin/env Rscript

# ==============================================================================
# Script: 03_variance_partitioning/checks/verify_varpart_identities.R
#
# Two-part verification of the global 2-block varpart NetCDF before
# any biological interpretation.
#
#  PART 1 — Arithmetic identities (decisive):
#    By construction of vegan::varpart the following must hold exactly:
#      indiv_SUMINISTRO = unique_SUMINISTRO + shared
#      indiv_DEMANDA    = unique_DEMANDA    + shared
#      total_r2         = unique_SUMINISTRO + unique_DEMANDA + shared
#    Checked on 10,000 random valid cell×months (tolerance 1e-5).
#    If these hold, layers are correctly labelled relative to each other.
#
#  PART 2 — Recalculation from source data (3 cells):
#    Re-runs varpart from raw kNDVI + GLEAM series for 3 specific cells
#    and compares each fraction against the NC value. Requires local data
#    files; skipped with instructions if not found.
#
# Run from repo root:
#   Rscript 03_variance_partitioning/checks/verify_varpart_identities.R
# ==============================================================================

suppressPackageStartupMessages(library(ncdf4))

file_nc      <- file.path("outputs", "varpart_global", "varpart_global_2blocks.nc")
file_nc_test <- file.path("outputs", "varpart_global", "varpart_global_2blocks_test.nc")
file_dom_csv <- file.path("outputs", "varpart_global", "dominance_summary_2blocks.csv")
file_dom_test<- file.path("outputs", "varpart_global", "dominance_summary_2blocks_test.csv")

FV  <- -9999.0
TOL <- 1e-5
N_SAMPLE <- 10000L

set.seed(42L)

# ==============================================================================
# Load global NC
# ==============================================================================

if (!file.exists(file_nc)) stop("Global NC not found: ", file_nc)
nc  <- nc_open(file_nc)
VARS <- c("unique_SUMINISTRO", "unique_DEMANDA", "indiv_SUMINISTRO", "indiv_DEMANDA",
          "shared", "total_r2")
arr <- list()
for (v in VARS) {
  a <- ncvar_get(nc, v); a[abs(a - FV) < 1] <- NA_real_; arr[[v]] <- a
}
nc_close(nc)

valid_idx <- which(!is.na(arr[["total_r2"]]))
cat(sprintf("Global NC: %d valid cell×months\n", length(valid_idx)))

# ==============================================================================
# PART 1 — Arithmetic identities
# ==============================================================================

cat("\n\n================================================================\n")
cat("=== PART 1: Arithmetic identities (10,000 random cell×months)\n")
cat("================================================================\n\n")

samp <- sample(valid_idx, min(N_SAMPLE, length(valid_idx)))

su  <- arr[["unique_SUMINISTRO"]][samp]
de  <- arr[["unique_DEMANDA"]][samp]
isu <- arr[["indiv_SUMINISTRO"]][samp]
ide <- arr[["indiv_DEMANDA"]][samp]
sh  <- arr[["shared"]][samp]
tot <- arr[["total_r2"]][samp]

check1 <- abs(isu - (su + sh))          # indiv_SUM = unique_SUM + shared
check2 <- abs(ide - (de + sh))          # indiv_DEM = unique_DEM + shared
check3 <- abs(tot - (su + de + sh))     # total_r2  = unique_SUM + unique_DEM + shared

pct1 <- mean(check1 < TOL) * 100
pct2 <- mean(check2 < TOL) * 100
pct3 <- mean(check3 < TOL) * 100

cat(sprintf("  Tolerance: %.0e\n\n", TOL))
cat(sprintf("  Identity 1: indiv_SUMINISTRO = unique_SUMINISTRO + shared\n"))
cat(sprintf("    Max residual : %.2e\n", max(check1)))
cat(sprintf("    %% cells pass : %.4f%%  %s\n",
            pct1, if (pct1 > 99.9) "PASS ✓" else "FAIL ✗"))

cat(sprintf("\n  Identity 2: indiv_DEMANDA = unique_DEMANDA + shared\n"))
cat(sprintf("    Max residual : %.2e\n", max(check2)))
cat(sprintf("    %% cells pass : %.4f%%  %s\n",
            pct2, if (pct2 > 99.9) "PASS ✓" else "FAIL ✗"))

cat(sprintf("\n  Identity 3: total_r2 = unique_SUMINISTRO + unique_DEMANDA + shared\n"))
cat(sprintf("    Max residual : %.2e\n", max(check3)))
cat(sprintf("    %% cells pass : %.4f%%  %s\n",
            pct3, if (pct3 > 99.9) "PASS ✓" else "FAIL ✗"))

cat("\n  Summary of residuals (all three combined):\n")
all_resid <- c(check1, check2, check3)
qr <- quantile(all_resid, probs=c(0.50,0.90,0.95,0.99,1.00))
for (i in seq_along(qr))
  cat(sprintf("    p%-3.0f : %.2e\n", as.numeric(names(qr)[i])*100, qr[i]))

if (pct1 > 99.9 && pct2 > 99.9 && pct3 > 99.9) {
  cat("\n  VERDICT: All three identities hold. Layers are consistently\n")
  cat("  labelled relative to each other — no cross-assignment detected.\n")
} else {
  cat("\n  VERDICT: *** ONE OR MORE IDENTITIES FAIL — LABEL CROSS SUSPECTED ***\n")
  cat("  Check which identity fails to pinpoint which layers are swapped.\n")
}

# ==============================================================================
# PART 2 — Recalculation from source data (3 cells)
# ==============================================================================

cat("\n\n================================================================\n")
cat("=== PART 2: Recalculation from source data (3 cells)\n")
cat("================================================================\n\n")

if (file.exists("R/config.R"))           source("R/config.R")
if (file.exists("R/config.R")) source("R/config.R")

file_kndvi <- if (exists("paths") && !is.null(paths$kndvi)) paths$kndvi else
  file.path("outputs", "intermediate", "kndvi", "kndvi.nc")
dir_accum  <- if (exists("paths") && !is.null(paths$accumulated)) paths$accumulated else
  file.path("outputs", "intermediate", "accumulated")

SCALES     <- c(1, 3, 6, 9, 12)
BLOCKS     <- list(SUMINISTRO = c("SMs", "SMrz"), DEMANDA = "Ep")
BLOCK_NAMES<- names(BLOCKS)
N_TIME     <- 492L
file_tpl   <- "{var}_GLEAM_v4.2a_MO_1982-2022_scale_{scale}.nc"

has_kndvi <- file.exists(file_kndvi)
has_gleam <- file.exists(file.path(dir_accum,
  gsub("\\{scale\\}",1,gsub("\\{var\\}","SMs",file_tpl))))

if (!has_kndvi || !has_gleam) {
  cat("  Source data not found locally. To run Part 2, either:\n")
  cat("  (a) Run the upstream pipeline (01_prepare_data, 02_correlations) first, or\n")
  cat("  (b) rsync the kNDVI and GLEAM files to local outputs/intermediate/\n\n")
  cat(sprintf("  kNDVI  : %s  [%s]\n", file_kndvi,
              if (has_kndvi) "FOUND" else "NOT FOUND"))
  cat(sprintf("  GLEAM  : %s/SMs_...  [%s]\n", dir_accum,
              if (has_gleam) "FOUND" else "NOT FOUND"))
  cat("\n  Skipping Part 2.\n")
} else {
  suppressPackageStartupMessages({
    library(pracma); library(vegan)
  })

  # ── Helper functions ───────────────────────────────────────────────────────
  detrend_series <- function(x) {
    if (sum(!is.na(x)) < 5L) return(rep(NA_real_, length(x)))
    pracma::detrend(x, tt = "linear")
  }
  reduce_block <- function(mat) {
    ok <- complete.cases(mat)
    if (sum(ok) < 10L) return(NULL)
    mat_ok <- mat[ok, , drop = FALSE]
    col_sd <- apply(mat_ok, 2L, sd)
    keep   <- !is.na(col_sd) & col_sd > 0
    if (!any(keep)) return(NULL)
    mat_ok <- mat_ok[, keep, drop = FALSE]
    pca    <- prcomp(mat_ok, center=TRUE, scale.=TRUE)
    sc     <- matrix(NA_real_, nrow(mat), 1L)
    sc[ok, ] <- pca$x[, 1L, drop = FALSE]
    sc
  }
  do_varpart_full <- function(y_dt, bpcs) {
    if (any(sapply(bpcs, is.null))) return(NULL)
    ok  <- complete.cases(cbind(y_dt, do.call(cbind, bpcs)))
    if (sum(ok) < 10L) return(NULL)
    y_s <- y_dt[ok]
    pfx <- c("SUM","DEM")
    npcs <- mapply(function(b,p) {
      m <- b[ok,,drop=FALSE]; colnames(m) <- paste0(p,"_PC",seq_len(ncol(m))); m
    }, bpcs, pfx, SIMPLIFY=FALSE)
    df   <- as.data.frame(do.call(cbind, npcs))
    fmls <- lapply(npcs, function(m)
      as.formula(paste("~", paste(colnames(m), collapse=" + "))))
    vp   <- vegan::varpart(y_s, fmls[[1]], fmls[[2]], data=df)
    adj  <- grep("Adj\\.R", colnames(vp$part$fract), value=TRUE)[1L]
    ind  <- vp$part$indfract
    fr   <- vp$part$fract
    list(
      unique_SUMINISTRO = ind[1L, adj],
      unique_DEMANDA    = ind[2L, adj],
      shared            = ind[3L, adj],
      indiv_SUMINISTRO  = fr[1L, adj],
      indiv_DEMANDA     = fr[2L, adj],
      total_r2          = fr[nrow(fr), adj],
      n_obs             = sum(ok)
    )
  }

  # ── Load NC axes ────────────────────────────────────────────────────────────
  nc0     <- nc_open(file_kndvi)
  lon_all <- ncvar_get(nc0, "lon")
  lat_all <- ncvar_get(nc0, "lat")
  nc_close(nc0)

  # Load mask for domain derivation
  file_veg_mask <- if (exists("paths") && !is.null(paths$veg_mask_c1)) paths$veg_mask_c1 else
    file.path("outputs", "intermediate", "vegetation_mask_c1.tif")
  suppressPackageStartupMessages(library(terra))
  r_mask <- terra::rast(file_veg_mask)
  mlr    <- c(terra::ext(r_mask)$ymin, terra::ext(r_mask)$ymax)
  dr     <- which(lat_all >= mlr[1L] & lat_all <= mlr[2L])
  lat_domain_start <- min(dr)

  nc_lat <- nc_open(file_nc); lat_nc <- ncvar_get(nc_lat, "lat"); nc_close(nc_lat)

  # ── Pick 3 cells: find valid cells in specific lon/lat regions ──────────────
  # South America (lon≈-60°, lat≈-10°), Africa (lon≈25°, lat≈-5°), Australia (lon≈135°, lat≈-25°)
  targets <- list(
    list(region="S.America  ", lon_tgt=-60, lat_tgt=-10, month=7L),
    list(region="Africa     ", lon_tgt= 25, lat_tgt= -5, month=4L),
    list(region="Australia  ", lon_tgt=135, lat_tgt=-25, month=1L)
  )

  cat(sprintf("  %-12s  %8s  %8s  %8s  %8s  %8s  %8s  %5s  %s\n",
              "Region", "uniq_SUM", "uniq_DEM", "shared", "indivSUM", "indivDEM",
              "total_r2", "n_obs", "match"))
  cat(sprintf("  %s\n", strrep("-", 100)))

  for (tgt in targets) {
    li   <- which.min(abs(lon_all - tgt$lon_tgt))
    la_nc<- which.min(abs(lat_nc  - tgt$lat_tgt))   # index into lat_nc (=lat_domain)
    la_d <- la_nc                                     # same: NC lat = lat_domain
    mo   <- tgt$month

    # NC values
    nc_vals <- list()
    nc_tmp <- nc_open(file_nc)
    for (v in c("unique_SUMINISTRO","unique_DEMANDA","shared",
                "indiv_SUMINISTRO","indiv_DEMANDA","total_r2")) {
      val <- ncvar_get(nc_tmp, v, start=c(li, la_nc, mo), count=c(1,1,1))
      nc_vals[[v]] <- if (abs(val - FV) < 1) NA_real_ else val
    }
    nc_close(nc_tmp)

    if (is.na(nc_vals[["total_r2"]])) {
      cat(sprintf("  %-12s  lon=%.1f lat=%.1f  mo=%d — no data in NC\n",
                  tgt$region, lon_all[li], lat_nc[la_nc], mo)); next
    }

    # Source data: load 1-cell time series
    nc_lat_row <- lat_domain_start + la_d - 1L
    midx       <- seq(mo, N_TIME, by=12L)

    nc_k <- nc_open(file_kndvi)
    y_raw <- ncvar_get(nc_k, "kndvi", start=c(li,nc_lat_row,1L), count=c(1,1,N_TIME))[1,1,]
    nc_close(nc_k)
    y_dt <- detrend_series(y_raw[midx])

    bm <- list()
    for (bn in BLOCK_NAMES)
      bm[[bn]] <- do.call(cbind, lapply(BLOCKS[[bn]], function(ind)
        do.call(cbind, lapply(SCALES, function(sc) {
          fn <- gsub("\\{scale\\}",sc,gsub("\\{var\\}",ind,file_tpl))
          nc_g <- nc_open(file.path(dir_accum, fn))
          x <- ncvar_get(nc_g, ind, start=c(li,nc_lat_row,1L), count=c(1,1,N_TIME))[1,1,]
          nc_close(nc_g)
          detrend_series(x[midx])
        }))))
    bpcs <- lapply(BLOCK_NAMES, function(bn) reduce_block(bm[[bn]]))
    names(bpcs) <- BLOCK_NAMES

    recalc <- tryCatch(do_varpart_full(y_dt, bpcs), error=function(e) NULL)
    if (is.null(recalc)) {
      cat(sprintf("  %-12s  lon=%.1f lat=%.1f  mo=%d — recalc returned NULL\n",
                  tgt$region, lon_all[li], lat_nc[la_nc], mo)); next
    }

    # Compare
    chk_vars <- c("unique_SUMINISTRO","unique_DEMANDA","shared",
                  "indiv_SUMINISTRO","indiv_DEMANDA","total_r2")
    max_diff <- max(abs(sapply(chk_vars, function(v) recalc[[v]] - nc_vals[[v]])))
    match_str <- if (max_diff < 1e-4) "OK ✓" else sprintf("DIFF %.1e", max_diff)

    cat(sprintf("  %-12s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %5d  %s\n",
                tgt$region,
                nc_vals$unique_SUMINISTRO, nc_vals$unique_DEMANDA, nc_vals$shared,
                nc_vals$indiv_SUMINISTRO,  nc_vals$indiv_DEMANDA,  nc_vals$total_r2,
                recalc$n_obs, match_str))
    cat(sprintf("  %-12s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %5s  (recalc)\n",
                "",
                recalc$unique_SUMINISTRO, recalc$unique_DEMANDA, recalc$shared,
                recalc$indiv_SUMINISTRO,  recalc$indiv_DEMANDA,  recalc$total_r2, ""))
    cat(sprintf("  %-12s  %8s  %8s  %8s  %8s  %8s  %8s\n",
                "  max|Δ|",
                formatC(abs(recalc$unique_SUMINISTRO - nc_vals$unique_SUMINISTRO), format="e", digits=1),
                formatC(abs(recalc$unique_DEMANDA    - nc_vals$unique_DEMANDA),    format="e", digits=1),
                formatC(abs(recalc$shared            - nc_vals$shared),            format="e", digits=1),
                formatC(abs(recalc$indiv_SUMINISTRO  - nc_vals$indiv_SUMINISTRO),  format="e", digits=1),
                formatC(abs(recalc$indiv_DEMANDA     - nc_vals$indiv_DEMANDA),     format="e", digits=1),
                formatC(abs(recalc$total_r2          - nc_vals$total_r2),          format="e", digits=1)))
    cat(sprintf("  %s\n", strrep("-", 100)))
  }
}
