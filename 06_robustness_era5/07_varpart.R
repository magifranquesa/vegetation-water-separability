#!/usr/bin/env Rscript

# ==============================================================================
# Script: 06_robustness_era5/07_varpart.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
# Purpose:
#   Two-block variance partitioning of interannual kNDVI between a soil-water
#   supply block (SMs + SMrz) and an atmospheric demand block (Ep), with
#   permutation-based significance tests, per grid cell and calendar month over
#   the global vegetated domain. ERA5-Land inputs (robustness assessment).
#
#   Per cell x month it computes the seven Adj.R2 fractions and three p-values:
#     p_full   — full model RDA permutation test (vegan::anova.cca)
#     p_soil   — partial RDA unique supply block (conditioned on demand)
#     p_demand — partial RDA unique demand block (conditioned on supply)
#   P-values are raw (no FDR correction; applied post-hoc over the full domain).
#
# Outputs:
#   NetCDF : outputs/varpart_global/varpart_signif_global_2blocks_era5.nc
#   CSV    : outputs/varpart_global/dominance_signif_2blocks_era5.csv
#   RDS    : outputs/varpart_global/stripes_signif_era5/
#
# Workflow:
#   1. Run with TEST_MODE <- TRUE on lat [-30,-20] (cost probe + diagnostics).
#   2. Once timings are acceptable, set TEST_MODE <- FALSE for global run.
#
# Run from repo root:
#   Rscript 06_robustness_era5/07_varpart.R
#
# Dependencies: ncdf4, terra, pracma, vegan, parallel
# ==============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(terra)
  library(pracma)
  library(vegan)
  library(parallel)
})

t_start <- proc.time()

# ==============================================================================
# 1. Configuration
# ==============================================================================

# ── Permutation parameter ────────────────────────────────────────────────────
NPERM <- 999L       # raw p-value resolution = 1/(NPERM+1) = 0.001

# ── Run mode ─────────────────────────────────────────────────────────────────
TEST_MODE      <- FALSE       # FALSE → global run
TEST_LAT_RANGE <- c(-30, -20) # only used when TEST_MODE = TRUE

N_CORES    <- 90L
STRIPE_DEG <- 10

SCALES      <- c(1, 3, 6, 9, 12)
BLOCKS      <- list(SUMINISTRO = c("SMs", "SMrz"), DEMANDA = "Ep")
BLOCK_NAMES <- names(BLOCKS)
N_TIME      <- 492L          # 1982–2022: 41 yr × 12 mo

if (file.exists(file.path("R", "config.R")))            source(file.path("R", "config.R"))
if (file.exists(file.path("R", "config_era5.R"))) source(file.path("R", "config_era5.R"))

file_kndvi    <- if (exists("paths") && !is.null(paths$kndvi)) paths$kndvi else
  file.path("outputs", "intermediate", "kndvi", "kndvi.nc")
dir_accum     <- if (exists("paths") && !is.null(paths$accumulated_era5)) paths$accumulated_era5 else
  file.path("outputs", "intermediate", "accumulated_era5")
file_veg_mask <- if (exists("paths") && !is.null(paths$veg_mask_c1)) paths$veg_mask_c1 else
  file.path("outputs", "intermediate", "vegetation_mask_c1.tif")
file_tpl <- "{var}_era5land_monthly_1982-2022_scale_{scale}.nc"

# ── Output paths ────────────────────
sfx         <- if (TEST_MODE) "_test_era5" else "_era5"
dir_stripes <- file.path("outputs", "varpart_global", paste0("stripes_signif", sfx))
out_nc      <- file.path("outputs", "varpart_global",
                          paste0("varpart_signif_global_2blocks", sfx, ".nc"))
out_dom_csv <- file.path("outputs", "varpart_global",
                          paste0("dominance_signif_2blocks", sfx, ".csv"))
dir.create(dir_stripes, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_nc), recursive = TRUE, showWarnings = FALSE)

# ── Variable metadata: 7 fractions + 3 p-values ──────────
VAR_META <- list(
  unique_SUMINISTRO = "Unique Adj.R2 for soil moisture block (SMs+SMrz, PC1)",
  unique_DEMANDA    = "Unique Adj.R2 for atmospheric demand block (Ep, PC1)",
  indiv_SUMINISTRO  = "Marginal Adj.R2 for soil moisture block alone",
  indiv_DEMANDA     = "Marginal Adj.R2 for atmospheric demand block alone",
  shared            = "Shared (joint) Adj.R2 between SUMINISTRO and DEMANDA",
  total_r2          = "Total explained variance, full 2-block model (Adj.R2)",
  unexplained       = "Unexplained fraction: 1 - max(total_r2, 0)",
  p_full            = paste0("P-value: full model RDA permutation test (NPERM=",NPERM,")"),
  p_soil            = paste0("P-value: partial RDA unique SUMINISTRO | Condition(DEMANDA) (NPERM=",NPERM,")"),
  p_demand          = paste0("P-value: partial RDA unique DEMANDA | Condition(SUMINISTRO) (NPERM=",NPERM,")")
)
VAR_NAMES <- names(VAR_META)

# ==============================================================================
# 2. Helper functions
# ==============================================================================

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
  n_drop <- sum(!keep)
  if (n_drop > 0L) mat_ok <- mat_ok[, keep, drop = FALSE]
  if (ncol(mat_ok) == 0L) return(NULL)
  pca          <- prcomp(mat_ok, center = TRUE, scale. = TRUE)
  scores       <- matrix(NA_real_, nrow(mat), 1L)
  scores[ok, ] <- pca$x[, 1L, drop = FALSE]
  attr(scores, "had_const_col") <- n_drop > 0L
  scores
}

# do_varpart_signif: Adj.R² fractions, plus 3 p-values.
#
# RDA matrix interface:
#   vegan::rda(Y, X)    → Y ~ X          (no covariable)
#   vegan::rda(Y, X, W) → Y ~ X | W      (partial, Condition on W)
#
# The full-model R²adj from rda(Y, X1+X2) MUST equal total_r2 from varpart
# (both use the same Ezekiel Adj.R² formula). Verified in TEST_MODE diagnostic.
#
# The partial p-values test significance of the unique fractions:
#   p_soil   → H0: unique_SUMINISTRO = 0
#   p_demand → H0: unique_DEMANDA    = 0

do_varpart_signif <- function(y_detrended, block_pcs) {
  if (any(sapply(block_pcs, is.null))) return(NULL)
  ok <- complete.cases(cbind(y_detrended, do.call(cbind, block_pcs)))
  if (sum(ok) < 10L) return(NULL)
  y_sub <- y_detrended[ok]
  pfx   <- c("SUM", "DEM")
  npcs  <- mapply(function(b, p) {
    m <- b[ok, , drop = FALSE]; colnames(m) <- paste0(p, "_PC", seq_len(ncol(m))); m
  }, block_pcs, pfx, SIMPLIFY = FALSE)
  df    <- as.data.frame(do.call(cbind, npcs))
  fmlas <- lapply(npcs, function(m)
    as.formula(paste("~", paste(colnames(m), collapse = " + "))))

  # ── Adj.R² fractions (vegan::varpart) ─────────────
  vp      <- suppressWarnings(vegan::varpart(y_sub, fmlas[[1]], fmlas[[2]], data = df))
  adj_col <- grep("Adj\\.R", colnames(vp$part$fract), value = TRUE)[1L]
  indf    <- vp$part$indfract
  tot     <- vp$part$fract[nrow(vp$part$fract), adj_col]
  uniq    <- setNames(indf[1:2, adj_col], paste0("unique_", BLOCK_NAMES))
  indiv   <- setNames(vp$part$fract[1:2, adj_col], paste0("indiv_", BLOCK_NAMES))

  # ── Permutation tests via RDA (matrix interface) ──────────────────────────
  # suppressWarnings: vegan emits "no unconstrained (residual) component" for
  # univariate Y because there are 0 unconstrained ordination AXES (rank(Y)=1
  # leaves no room after 1 constrained axis). The residual VARIANCE exists and
  # anova.cca uses it correctly in permutations. Verified cosmetic:
  # rho(R², p_full) = −0.993 on 594 arctic fits, unc_inertia > 0 in all cells.
  y_mat    <- matrix(y_sub, ncol = 1L)
  sum_mat  <- npcs[["SUMINISTRO"]]
  dem_mat  <- npcs[["DEMANDA"]]

  rda_full <- suppressWarnings(vegan::rda(y_mat, cbind(sum_mat, dem_mat)))
  rda_soil <- suppressWarnings(vegan::rda(y_mat, sum_mat, dem_mat))
  rda_dem  <- suppressWarnings(vegan::rda(y_mat, dem_mat, sum_mat))

  safe_pval <- function(rda_obj) {
    tryCatch({
      a <- suppressWarnings(vegan::anova.cca(rda_obj, permutations = NPERM))
      a[["Pr(>F)"]][1L]
    }, error = function(e) NA_real_)
  }

  p_full   <- safe_pval(rda_full)
  p_soil   <- safe_pval(rda_soil)
  p_demand <- safe_pval(rda_dem)

  c(uniq, indiv,
    shared      = indf[3L, adj_col],
    total_r2    = tot,
    unexplained = 1 - max(tot, 0),
    n_obs       = sum(ok),
    p_full      = p_full,
    p_soil      = p_soil,
    p_demand    = p_demand)
}

# ==============================================================================
# 3. Load reference grids and verify vegetation mask
# ==============================================================================

message("Loading kNDVI reference grid...")
nc0     <- nc_open(file_kndvi)
lon_all <- ncvar_get(nc0, "lon")
lat_all <- ncvar_get(nc0, "lat")
nc_close(nc0)
res_deg <- abs(lat_all[2] - lat_all[1])
message("  kNDVI full grid: lon=", length(lon_all), "  lat=", length(lat_all),
        "  res=", res_deg, "°")

message("Loading vegetation mask (C1 GeoTIFF): ", file_veg_mask)
if (!file.exists(file_veg_mask))
  stop("Vegetation mask not found: ", file_veg_mask)

r_kref <- terra::rast(file_kndvi)[[1]]
r_mask <- terra::rast(file_veg_mask)

mask_lat_range   <- c(terra::ext(r_mask)$ymin, terra::ext(r_mask)$ymax)
domain_rows      <- which(lat_all >= mask_lat_range[1] & lat_all <= mask_lat_range[2])
lat_domain_start <- min(domain_rows)
lat_domain_end   <- max(domain_rows)
lat_domain_count <- lat_domain_end - lat_domain_start + 1L
lat_domain       <- lat_all[lat_domain_start:lat_domain_end]
message(sprintf("  Analysis domain : lat n=%d  [%.3f, %.3f]  (kNDVI rows %d–%d)",
                lat_domain_count, min(lat_domain), max(lat_domain),
                lat_domain_start, lat_domain_end))
stopifnot("lat_domain_count must be positive"        = lat_domain_count > 0L)
stopifnot("lat_domain_count must equal mask nrow()"  = lat_domain_count == nrow(r_mask))

r_kref_crop <- terra::crop(r_kref, terra::ext(r_mask))
geom_ok <- isTRUE(all.equal(as.vector(terra::ext(r_kref_crop)),
                              as.vector(terra::ext(r_mask)), tolerance = 1e-4)) &&
           isTRUE(all.equal(terra::res(r_kref_crop), terra::res(r_mask), tolerance = 1e-6)) &&
           ncol(r_kref_crop) == ncol(r_mask) && nrow(r_kref_crop) == nrow(r_mask)
if (!geom_ok)
  stop(sprintf("Cropped kNDVI domain (%d×%d) does NOT match mask (%d×%d).",
               ncol(r_kref_crop), nrow(r_kref_crop), ncol(r_mask), nrow(r_mask)))
message(sprintf("  Geometry check : kNDVI domain %d×%d matches mask %d×%d — OK",
                ncol(r_kref_crop), nrow(r_kref_crop), ncol(r_mask), nrow(r_mask)))

arr_terra <- as.array(r_mask)[,,1]
if (lat_domain[1] < lat_domain[2]) {
  veg_mask <- t(arr_terra[nrow(arr_terra):1, ])
} else {
  veg_mask <- t(arr_terra)
}
veg_mask[is.na(veg_mask)] <- 0L
if (nrow(veg_mask) != length(lon_all) || ncol(veg_mask) != lat_domain_count)
  stop(sprintf("veg_mask dims %d×%d ≠ expected %d×%d",
               nrow(veg_mask), ncol(veg_mask), length(lon_all), lat_domain_count))
message(sprintf("  veg_mask        : %d × %d  (lon × lat_domain) — OK",
                nrow(veg_mask), ncol(veg_mask)))
message("  Vegetated cells : ", sum(veg_mask > 0), " / ", prod(dim(veg_mask)))

message("\nSanity check: mask orientation vs kNDVI...")
nc_kchk  <- nc_open(file_kndvi)
kndvi_t1 <- ncvar_get(nc_kchk, "kndvi",
                       start = c(1L, lat_domain_start, 1L),
                       count = c(length(lon_all), lat_domain_count, 1L))
nc_close(nc_kchk)
k_valid        <- !is.na(kndvi_t1)
m_veg          <- veg_mask > 0
n_veg          <- sum(m_veg);  n_nonveg <- sum(!m_veg)
pct_veg_kok    <- round(sum(m_veg  &  k_valid)  / n_veg    * 100, 1)
pct_nonveg_kna <- round(sum(!m_veg & !k_valid)  / n_nonveg * 100, 1)
message(sprintf("  mask=1 ∩ kNDVI valid : %7d / %7d  (%5.1f%%)",
                sum(m_veg & k_valid), n_veg, pct_veg_kok))
message(sprintf("  mask=0 ∩ kNDVI NA    : %7d / %7d  (%5.1f%%)",
                sum(!m_veg & !k_valid), n_nonveg, pct_nonveg_kna))
if (pct_veg_kok < 80 || pct_nonveg_kna < 40)
  stop("!! Mask orientation suspect. Check row-flip logic.")
message("  Concordance OK — proceeding.")
rm(kndvi_t1, k_valid, m_veg)

# ==============================================================================
# 4. Define latitude stripes
# ==============================================================================

STRIPE_CELLS  <- round(STRIPE_DEG / res_deg)
stripe_starts <- seq(1L, lat_domain_count, by = STRIPE_CELLS)

if (TEST_MODE) {
  test_la   <- which(lat_domain >= TEST_LAT_RANGE[1] & lat_domain <= TEST_LAT_RANGE[2])
  stripe_starts <- stripe_starts[stripe_starts >= min(test_la) &
                                   stripe_starts <= max(test_la)]
  message("\nTEST MODE: ", length(stripe_starts), " stripe(s) in lat_domain [",
          TEST_LAT_RANGE[1], ", ", TEST_LAT_RANGE[2], "]")
  message("  NPERM = ", NPERM, "  (cost probe)")
} else {
  message("\nGLOBAL MODE: ", length(stripe_starts), " stripes × ~", STRIPE_DEG,
          "° lat  |  ", N_CORES, " cores  |  NPERM=", NPERM)
}

# ==============================================================================
# 5. Process stripes
# ==============================================================================

lon_n <- length(lon_all)

for (s in seq_along(stripe_starts)) {
  la_s <- stripe_starts[s]
  la_e <- min(la_s + STRIPE_CELLS - 1L, lat_domain_count)
  la_n <- la_e - la_s + 1L

  rds_file <- file.path(dir_stripes,
                         sprintf("stripe_la%04d_%04d.rds", la_s, la_e))

  if (file.exists(rds_file)) {
    message(sprintf("[%d/%d] stripe la%04d–%04d (%.2f°–%.2f°) already done — skipping",
                    s, length(stripe_starts), la_s, la_e,
                    lat_domain[la_s], lat_domain[la_e]))
    next
  }

  message(sprintf("\n[%d/%d] stripe la%04d–%04d  lat %.2f°–%.2f°  (NPERM=%d)",
                  s, length(stripe_starts), la_s, la_e,
                  lat_domain[la_s], lat_domain[la_e], NPERM))
  t_s <- proc.time()

  veg_s   <- veg_mask[, la_s:la_e, drop = FALSE]
  veg_idx <- which(veg_s > 0, arr.ind = TRUE)
  n_veg   <- nrow(veg_idx)
  message("  Vegetated cells in stripe: ", n_veg)

  if (n_veg == 0L) {
    saveRDS(list(data = data.frame(), counts = setNames(rep(0L, 5),
      c("skip_kndvi","skip_vp","skip_error","const_drop","success"))), rds_file)
    next
  }

  nc_start_lat <- lat_domain_start + la_s - 1L

  message("  Loading kNDVI stripe (row ", nc_start_lat, ")...")
  nc_k         <- nc_open(file_kndvi)
  kndvi_stripe <- ncvar_get(nc_k, "kndvi",
                             start = c(1L, nc_start_lat, 1L),
                             count = c(lon_n, la_n, N_TIME))
  nc_close(nc_k)

  message("  Loading ERA5-Land stripe...")
  gleam_stripe <- list()
  for (bn in BLOCK_NAMES)
    for (ind in BLOCKS[[bn]])
      for (sc in SCALES) {
        fname <- gsub("\\{scale\\}", sc, gsub("\\{var\\}", ind, file_tpl))
        nc_g  <- nc_open(file.path(dir_accum, fname))
        arr   <- ncvar_get(nc_g, ind,
                           start = c(1L, nc_start_lat, 1L),
                           count = c(lon_n, la_n, N_TIME))
        nc_close(nc_g)
        gleam_stripe[[paste0(ind,"_s",sc)]] <- arr
      }
  rm(arr)

  veg_li       <- veg_idx[, 1L]
  veg_la_local <- veg_idx[, 2L]
  veg_la_glob  <- veg_la_local + la_s - 1L
  veg_lon_vals <- lon_all[veg_li]
  veg_lat_vals <- lat_domain[veg_la_glob]

  message("  mclapply over ", n_veg, " cells × 12 months on ", N_CORES, " cores...")

  cell_results <- parallel::mclapply(
    seq_len(n_veg),
    function(i) {
      tryCatch({
        li  <- veg_li[i];        la  <- veg_la_local[i]
        lag <- veg_la_glob[i];   lv  <- veg_lon_vals[i];  ltv <- veg_lat_vals[i]
        rows <- list()
        cnts <- c(skip_kndvi=0L, skip_vp=0L, skip_error=0L, const_drop=0L, success=0L)

        for (mo in 1:12) {
          midx <- seq(mo, N_TIME, by = 12L)
          y_dt <- detrend_series(kndvi_stripe[li, la, midx])

          if (sum(!is.na(y_dt)) < 10L) {
            cnts["skip_kndvi"] <- cnts["skip_kndvi"] + 1L; next
          }

          bm <- lapply(BLOCK_NAMES, function(bn)
            do.call(cbind, lapply(BLOCKS[[bn]], function(ind)
              do.call(cbind, lapply(SCALES, function(sc)
                detrend_series(gleam_stripe[[paste0(ind,"_s",sc)]][li, la, midx]))))))
          names(bm) <- BLOCK_NAMES

          bpcs <- lapply(BLOCK_NAMES, function(bn) reduce_block(bm[[bn]]))
          names(bpcs) <- BLOCK_NAMES

          if (any(sapply(bpcs, function(b) isTRUE(attr(b, "had_const_col")))))
            cnts["const_drop"] <- cnts["const_drop"] + 1L

          res <- tryCatch(do_varpart_signif(y_dt, bpcs), error = function(e) NULL)

          if (is.null(res)) { cnts["skip_vp"] <- cnts["skip_vp"] + 1L; next }

          cnts["success"] <- cnts["success"] + 1L
          rows[[length(rows)+1]] <-
            c(li=li, la=lag, lon=lv, lat=ltv, month=mo, res)
        }
        list(rows=rows, counts=cnts)

      }, error = function(e) {
        list(rows=list(),
             counts=c(skip_kndvi=0L, skip_vp=0L, skip_error=1L,
                      const_drop=0L, success=0L))
      })
    },
    mc.cores     = N_CORES,
    mc.preschedule = TRUE
  )

  rows_flat  <- unlist(lapply(cell_results, `[[`, "rows"), recursive = FALSE)
  cnts_mat   <- do.call(rbind, lapply(cell_results, `[[`, "counts"))
  stripe_cnt <- colSums(cnts_mat, na.rm = TRUE)

  df_stripe <- if (length(rows_flat) > 0)
    as.data.frame(do.call(rbind, rows_flat))
  else
    data.frame()

  elapsed <- round((proc.time() - t_s)["elapsed"])
  message(sprintf("  Done in %ds — ok=%d  skip_kndvi=%d  skip_vp=%d  err=%d  const_drop=%d",
                  elapsed, stripe_cnt["success"], stripe_cnt["skip_kndvi"],
                  stripe_cnt["skip_vp"], stripe_cnt["skip_error"], stripe_cnt["const_drop"]))

  # ── TEST_MODE: spatial footprint ──────────────────────────────────────────
  if (TEST_MODE && nrow(df_stripe) > 0) {
    lons <- as.numeric(df_stripe$lon)
    lats <- as.numeric(df_stripe$lat)
    uc   <- unique(df_stripe[, c("lon", "lat")])
    message(sprintf("\n  [TEST] Spatial footprint  (%d unique cells)", nrow(uc)))
    message(sprintf("    lon : %.2f° – %.2f°", min(lons), max(lons)))
    message(sprintf("    lat : %.2f° – %.2f°", min(lats), max(lats)))
  }

  # ── TEST_MODE: verify that RDA Adj.R² reproduces varpart fractions ─────────
  # Checks full-model and marginal Adj.R²: all three should be identical to
  # varpart fractions (different denominator for partial RDA — NOT compared here).
  if (TEST_MODE && nrow(df_stripe) > 0) {
    message("\n  [TEST] Verifying RDA Adj.R² vs varpart fractions on 3 cells...")
    ok_rows <- df_stripe[!is.na(df_stripe$total_r2), ]
    n_check <- min(3L, nrow(ok_rows))

    for (ci in seq_len(n_check)) {
      r        <- ok_rows[ci, ]
      li_c     <- as.integer(r$li)
      la_loc_c <- as.integer(r$la) - la_s + 1L
      mo_c     <- as.integer(r$month)
      midx_c   <- seq(mo_c, N_TIME, by = 12L)

      y_dt_c <- detrend_series(kndvi_stripe[li_c, la_loc_c, midx_c])
      bm_c   <- lapply(BLOCK_NAMES, function(bn)
        do.call(cbind, lapply(BLOCKS[[bn]], function(ind)
          do.call(cbind, lapply(SCALES, function(sc)
            detrend_series(gleam_stripe[[paste0(ind,"_s",sc)]][li_c, la_loc_c, midx_c]))))))
      names(bm_c) <- BLOCK_NAMES
      bpcs_c <- lapply(BLOCK_NAMES, function(bn) reduce_block(bm_c[[bn]]))
      names(bpcs_c) <- BLOCK_NAMES

      ok_c     <- complete.cases(cbind(y_dt_c, do.call(cbind, bpcs_c)))
      y_sub_c  <- y_dt_c[ok_c]
      sum_c    <- bpcs_c[["SUMINISTRO"]][ok_c, , drop=FALSE]
      dem_c    <- bpcs_c[["DEMANDA"]][ok_c, , drop=FALSE]
      y_mat_c  <- matrix(y_sub_c, ncol=1L)

      rda_full_c <- tryCatch(vegan::rda(y_mat_c, cbind(sum_c, dem_c)), error=function(e) NULL)
      rda_sum_c  <- tryCatch(vegan::rda(y_mat_c, sum_c),                error=function(e) NULL)
      rda_dem_c  <- tryCatch(vegan::rda(y_mat_c, dem_c),                error=function(e) NULL)

      rda_tot <- if (!is.null(rda_full_c)) RsquareAdj(rda_full_c)$adj.r.squared else NA
      rda_iSU <- if (!is.null(rda_sum_c))  RsquareAdj(rda_sum_c)$adj.r.squared  else NA
      rda_iDE <- if (!is.null(rda_dem_c))  RsquareAdj(rda_dem_c)$adj.r.squared  else NA

      message(sprintf("\n    Cell %d: lon=%+.2f lat=%+.2f month=%d  n=%d",
                      ci, r$lon, r$lat, mo_c, as.integer(r$n_obs)))
      message(sprintf("      total_r2        : varpart=%+.6f  rda=%+.6f  diff=%.2e",
                      r$total_r2,         rda_tot, abs(r$total_r2         - rda_tot)))
      message(sprintf("      indiv_SUMINISTRO: varpart=%+.6f  rda=%+.6f  diff=%.2e",
                      r$indiv_SUMINISTRO, rda_iSU, abs(r$indiv_SUMINISTRO - rda_iSU)))
      message(sprintf("      indiv_DEMANDA   : varpart=%+.6f  rda=%+.6f  diff=%.2e",
                      r$indiv_DEMANDA,    rda_iDE, abs(r$indiv_DEMANDA    - rda_iDE)))
    }
    message()
  }

  saveRDS(list(data = df_stripe, counts = stripe_cnt), rds_file)
  message("  Saved: ", basename(rds_file))

  rm(kndvi_stripe, gleam_stripe, cell_results, rows_flat)
  gc(verbose = FALSE)
}

# ==============================================================================
# 5b. TEST_MODE diagnostics (read from RDS — no recomputation)
# ==============================================================================

if (TEST_MODE) {
  message("\n", strrep("=", 60))
  message("[TEST_MODE] Post-stripe diagnostics")
  message(strrep("=", 60))

  rds_all <- sort(list.files(dir_stripes, "^stripe_la.*\\.rds$", full.names = TRUE))
  df_diag <- do.call(rbind, Filter(Negate(is.null), lapply(rds_all, function(f) {
    d <- readRDS(f)$data; if (nrow(d) == 0) NULL else d
  })))

  # ── Reference timing ──────────────────────────────────────────────────────
  t_stripe <- round((proc.time() - t_start)["elapsed"])
  message(sprintf("\n  Timing:"))
  message(sprintf("    This run (NPERM=%d) : %d s", NPERM, t_stripe))
  # Global extrapolation from this stripe's elapsed time
  if (nrow(veg_s <- veg_mask[, 1:1, drop=FALSE]) > 0) {   # dummy — just compute ratio
    n_stripes_global <- length(seq(1L, lat_domain_count, by = STRIPE_CELLS))
    est_global_s     <- t_stripe * n_stripes_global / length(stripe_starts)
    message(sprintf("    Estimated global    : %.0f s  (%.1f h)  [at %d cores]",
                    est_global_s, est_global_s / 3600, N_CORES))
  }

  if (!is.null(df_diag) && nrow(df_diag) > 0) {

    # ── P-value distributions ──────────────────────────────────────────────
    message(sprintf("\n  P-value distributions  (n rows = %d):", nrow(df_diag)))
    for (pv in c("p_full", "p_soil", "p_demand")) {
      vals <- as.numeric(df_diag[[pv]])
      n_ok <- sum(!is.na(vals))
      if (n_ok == 0) { message("    ", pv, " : all NA"); next }
      brks  <- c(0, 0.01, 0.05, 0.10, 0.50, 1.01)
      labs  <- c("[0,0.01]", "(0.01,0.05]", "(0.05,0.10]", "(0.10,0.50]", "(0.50,1]")
      hist_c <- table(cut(vals, breaks=brks, include.lowest=TRUE, labels=labs))
      pct_sig <- round(sum(vals <= 0.05, na.rm=TRUE) / n_ok * 100, 1)
      message(sprintf("    %-12s  (n valid=%d)  p<=0.05: %.1f%%", pv, n_ok, pct_sig))
      for (bn in names(hist_c))
        message(sprintf("      %-18s : %7d  (%5.1f%%)", bn,
                        as.integer(hist_c[bn]),
                        as.numeric(hist_c[bn]) / n_ok * 100))
    }

    # ── Coherence checks ──────────────────────────────────────────────────
    message("\n  Coherence checks (Spearman, expected negative rho: low p ↔ high R²):")

    r2   <- as.numeric(df_diag$total_r2)
    uSUM <- as.numeric(df_diag$unique_SUMINISTRO)
    uDEM <- as.numeric(df_diag$unique_DEMANDA)
    pfull  <- as.numeric(df_diag$p_full)
    psoil  <- as.numeric(df_diag$p_soil)
    pdem   <- as.numeric(df_diag$p_demand)

    ok_full <- !is.na(r2)   & !is.na(pfull)
    ok_soil <- !is.na(uSUM) & !is.na(psoil)
    ok_dema <- !is.na(uDEM) & !is.na(pdem)

    corr_full <- if (sum(ok_full) > 10)
      cor(pfull[ok_full], r2[ok_full],   method="spearman") else NA
    corr_soil <- if (sum(ok_soil) > 10)
      cor(psoil[ok_soil], uSUM[ok_soil], method="spearman") else NA
    corr_dema <- if (sum(ok_dema) > 10)
      cor(pdem[ok_dema],  uDEM[ok_dema], method="spearman") else NA

    message(sprintf("    rho(p_full,   total_r2)         = %+.4f  (expect < 0)", corr_full))
    message(sprintf("    rho(p_soil,   unique_SUMINISTRO) = %+.4f  (expect < 0)", corr_soil))
    message(sprintf("    rho(p_demand, unique_DEMANDA)    = %+.4f  (expect < 0)", corr_dema))

    if (!is.na(corr_full)  && corr_full  > 0)
      message("    !! WARN: p_full vs total_r2 correlation is POSITIVE — check RDA setup")
    if (!is.na(corr_soil)  && corr_soil  > 0)
      message("    !! WARN: p_soil vs unique_SUM correlation is POSITIVE — check partial RDA")
    if (!is.na(corr_dema)  && corr_dema  > 0)
      message("    !! WARN: p_demand vs unique_DEM correlation is POSITIVE — check partial RDA")

    # ── Cross-table: p_full threshold vs R² quartile ─────────────────────
    message("\n  Cross-table: p_full < 0.05 by total_r2 quartile")
    ok_ct <- !is.na(pfull) & !is.na(r2)
    if (sum(ok_ct) > 0) {
      qcuts <- quantile(r2[ok_ct], c(0, 0.25, 0.50, 0.75, 1))
      qlab  <- c("Q1 (low R²)", "Q2", "Q3", "Q4 (high R²)")
      r2_q  <- cut(r2[ok_ct], breaks=qcuts, include.lowest=TRUE, labels=qlab)
      sig   <- pfull[ok_ct] <= 0.05
      ct    <- table(r2_q, sig)
      for (q in rownames(ct)) {
        n_tot <- sum(ct[q,])
        n_sig <- if ("TRUE" %in% colnames(ct)) ct[q, "TRUE"] else 0L
        message(sprintf("    %-18s : %5d / %5d sig (%.0f%%)",
                        q, as.integer(n_sig), as.integer(n_tot),
                        as.integer(n_sig)/as.integer(n_tot)*100))
      }
    }
  }
}

# ==============================================================================
# 6. Assemble stripes → global NetCDF (10 variables)
# ==============================================================================

message("\n", strrep("=", 60))
message("Assembling NetCDF...")

rds_files <- sort(list.files(dir_stripes, "^stripe_la.*\\.rds$", full.names = TRUE))
if (length(rds_files) == 0) stop("No stripe RDS files found in: ", dir_stripes)

df_all <- do.call(rbind, Filter(Negate(is.null), lapply(rds_files, function(f) {
  d <- readRDS(f)$data; if (nrow(d) == 0) NULL else d
})))
message("  Total rows: ", nrow(df_all))

lon_dim <- length(lon_all); lat_dim <- lat_domain_count
out_arr <- lapply(VAR_NAMES, function(v) array(NA_real_, c(lon_dim, lat_dim, 12L)))
names(out_arr) <- VAR_NAMES

if (!is.null(df_all) && nrow(df_all) > 0) {
  idx3 <- cbind(as.integer(df_all$li), as.integer(df_all$la), as.integer(df_all$month))
  for (v in VAR_NAMES) out_arr[[v]][idx3] <- as.numeric(df_all[[v]])
}

fv        <- -9999.0
dim_lon   <- ncdim_def("lon",   "degrees_east",  lon_all)
dim_lat   <- ncdim_def("lat",   "degrees_north", lat_domain)
dim_month <- ncdim_def("month", "1",             1:12)

nc_vars <- lapply(VAR_NAMES, function(v)
  ncvar_def(v, "1", list(dim_lon, dim_lat, dim_month), fv,
            longname = VAR_META[[v]], prec = "float"))

message("  Writing: ", out_nc)
nc_out <- nc_create(out_nc, nc_vars, force_v4 = TRUE)

for (v in VAR_NAMES) {
  arr <- out_arr[[v]]; arr[is.na(arr)] <- fv
  ncvar_put(nc_out, v, arr)
}
ncatt_put(nc_out, 0, "title",
          "Global 2-block varpart + permutation tests: kNDVI ~ SUMINISTRO + DEMANDA")
ncatt_put(nc_out, 0, "method",
          "vegan::varpart + anova.cca, Adj.R2, 1 PC/block, linear detrend, 1982-2022")
ncatt_put(nc_out, 0, "blocks",
          "SUMINISTRO=SMs+SMrz, DEMANDA=Ep; 5 accumulation windows, ERA5-Land")
ncatt_put(nc_out, 0, "nperm",          as.integer(NPERM))
ncatt_put(nc_out, 0, "pvalue_note",    "Raw p-values (no FDR). Apply BH/BY correction post-hoc.")
ncatt_put(nc_out, 0, "test_mode",      as.integer(TEST_MODE))
ncatt_put(nc_out, 0, "date_created",   as.character(Sys.time()))
ncatt_put(nc_out, 0, "created_by",     "06_robustness_era5/07_varpart.R")
nc_close(nc_out)
message("  NetCDF written")

nc_check <- nc_open(out_nc)
message("\n  NC dimensions: ", paste(sapply(nc_check$dim, `[[`, "len"), collapse=" × "))
message("  NC variables : ", paste(names(nc_check$var), collapse=", "))
nc_close(nc_check)

# ==============================================================================
# 7. Dominance CSV (3-fraction classification)
# ==============================================================================

message("\nComputing dominance summary...")

if (!is.null(df_all) && nrow(df_all) > 0) {
  frac_mat <- data.frame(
    unique_SUMINISTRO = as.numeric(df_all$unique_SUMINISTRO),
    unique_DEMANDA    = as.numeric(df_all$unique_DEMANDA),
    shared            = as.numeric(df_all$shared)
  )
  dominant <- apply(frac_mat, 1, function(r) {
    if (all(is.na(r))) NA_character_ else names(r)[which.max(r)]
  })
  dom_tab <- table(dominant)
  dom_pct <- round(prop.table(dom_tab) * 100, 2)
  dom_df  <- data.frame(fraction  = names(dom_tab),
                         n_cells   = as.integer(dom_tab),
                         pct_total = as.numeric(dom_pct))
  write.csv(dom_df, out_dom_csv, row.names = FALSE)

  message("\n--- Dominance (cell×month) ---")
  for (i in seq_len(nrow(dom_df)))
    message(sprintf("  %-24s  n=%9d  (%.2f%%)",
                    dom_df$fraction[i], dom_df$n_cells[i], dom_df$pct_total[i]))
  message("  Total cell×months: ", nrow(df_all))
  message("  CSV: ", out_dom_csv)
}

# ==============================================================================
# 8. Global counters + timing
# ==============================================================================

message("\n=== Run summary ===")
all_cnt <- colSums(do.call(rbind, lapply(rds_files,
  function(f) { r <- readRDS(f)$counts; if (is.null(r)) integer(5) else r })),
  na.rm = TRUE)
for (nm in names(all_cnt))
  message(sprintf("  %-20s: %d", nm, as.integer(all_cnt[nm])))

if (!is.null(df_all) && nrow(df_all) > 0) {
  n_neg <- sum(as.numeric(df_all$shared) < 0, na.rm = TRUE)
  message(sprintf("  %-20s: %d  (%.2f%%)", "shared < 0",
                  n_neg, n_neg / nrow(df_all) * 100))

  for (pv in c("p_full", "p_soil", "p_demand")) {
    vals <- as.numeric(df_all[[pv]])
    n_s  <- sum(vals <= 0.05, na.rm=TRUE)
    n_ok <- sum(!is.na(vals))
    if (n_ok > 0)
      message(sprintf("  %-20s: %d / %d  (%.1f%% p<=0.05 raw)",
                      pv, n_s, n_ok, n_s/n_ok*100))
  }
}

t_tot <- (proc.time() - t_start)["elapsed"]
message(sprintf("  %-20s: %.0f s  (%.1f min)", "elapsed", t_tot, t_tot / 60))
message("\nNPERM = ", NPERM)
message("Done.")
