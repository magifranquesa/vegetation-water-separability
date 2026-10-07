#!/usr/bin/env Rscript

# ==============================================================================
# Script: 04_analysis/10_window_comparison.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   Compare the composition of hydroclimatic control along the aridity gradient
#   between two halves of the record (default 1982-2001 vs 2002-2022), to test
#   whether the non-separable ("joint") regime widens and whether atmospheric
#   demand gains weight over time.
#
#   Same machinery as the bottom row of Fig. 3: continuous AI axis in fine bins,
#   area-weighted, shares taken over the SIGNIFICANT subset rather than over the
#   vegetated domain. That normalisation is required here: with n~20 per window
#   the no-signal class grows, and it grows unevenly between windows, so shares
#   over the whole domain would not be comparable.
#
#   *** COMMON CELL x MONTH MASK ***
#   A cell x month enters the comparison only if the variance partitioning
#   succeeded in BOTH windows. Without this, a cell present in one window and
#   absent in the other (GIMMS3g has more gaps early in the record, especially
#   in high-latitude winter) would contribute a difference that is a change in
#   sampling, not a change in coupling. The mask is applied per cell x month,
#   so a cell dropped in January still counts in July.
#
#   Categories (as in 03_category_split.R, FDR-adjusted):
#     1 no_signal      full model not significant
#     2 soil           full sig, ONLY soil unique significant
#     3 demand         full sig, ONLY demand unique significant
#     4 coupled_pure   full sig, NEITHER unique significant   <- non-separable
#     5 both_unique    full sig, BOTH unique significant      <- separable
#
#   Two denominators are reported, both over significant area only:
#     share3 : soil + coupled_pure + demand      (Fig. 3B convention, for panels)
#     share4 : soil + demand + coupled_pure + both_unique  (full significant set)
#   share3 is what the published panel uses; share4 keeps category 5 visible,
#   which matters here because it is the opposite of category 4.
#
#   *** TWO COMPARISONS, AND WHICH ONE LEADS ***
#   The category comparison is power-limited: with n~20 only a few per cent of
#   the domain has either unique fraction significant after FDR, so coupled_pure
#   absorbs most of the significant area in BOTH windows and there is little room
#   left to detect change. The continuous Adj.R2 fractions carry no threshold and
#   are therefore the primary result; the categories are the complement.
#
# Inputs (READ-ONLY; one pair per window, produced by the varpart + FDR steps):
#   outputs/varpart_global/varpart_signif_global_2blocks_{window}.nc
#     (total_r2 + the continuous fractions unique_*, shared)
#   outputs/varpart_global/fdr_adjusted_pvalues_pooled_{window}.nc    (p_*_adj)
#   data/processed/aridity/ai_1982_2022_period.nc                     (continuous AI)
#
# Outputs (refuses to overwrite):
#   outputs/tables/window_comparison_magnitudes_by_ai.csv   <- primary
#   outputs/tables/window_comparison_shares_by_ai.csv
#   outputs/tables/window_comparison_crossover.csv
#   outputs/tables/window_comparison_by_aridity_class.csv
#
# Run from repo root, AFTER both windows have been through 01_varpart.R and the
# POOLED FDR step (02b, not 02: see FDR_MODE below and the header of 02b):
#   VWS_PERIOD=1982-2001 Rscript 03_variance_partitioning/01_varpart.R
#   VWS_PERIOD=2002-2022 Rscript 03_variance_partitioning/01_varpart.R
#   Rscript 03_variance_partitioning/02b_apply_fdr_pooled.R
#   Rscript 04_analysis/10_window_comparison.R
#
# This script does NOT read VWS_PERIOD: it reads both windows by construction.
#
# Dependencies: terra
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
})

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

# The two windows to compare, as the suffixes built by R/period.R.
WINDOWS <- c("1982_2001", "2002_2022")

ALPHA <- 0.05

# Which FDR files to threshold with.
#   "pooled"     -> fdr_adjusted_pvalues_pooled_{w}.nc   (02b_apply_fdr_pooled.R)
#   "per_window" -> fdr_adjusted_pvalues_{w}.nc          (02_apply_fdr.R)
# Use "pooled". Per-window FDR gives each window its own critical value, and on
# this split that alone turned a 4.4% shortfall in the minimum-p mass into
# "0.00% vs 4.14%" of significant unique soil. "per_window" is kept only so the
# artefact can be reproduced and reported.
FDR_MODE <- "pooled"

# Continuous Adj.R2 fractions, compared without any thresholding. These do not
# depend on test power, so they are the primary comparison: with n~20 the
# category-based one is power-limited to the point where coupled_pure absorbs
# most of the significant area in BOTH windows.
MAGNITUDES <- c("unique_SUMINISTRO", "unique_DEMANDA",
                "indiv_SUMINISTRO", "indiv_DEMANDA", "shared", "total_r2")

# ── Composition: the noise-invariant form of the question ────────────────────
# If the response carries more noise in one window, every Adj.R2 component is
# scaled by the SAME factor lambda = var(signal) / (var(signal) + var(noise)),
# because lambda is a property of y and not of the predictors. Ratios BETWEEN
# components therefore cancel it exactly.
#
# That matters because the question actually being asked is not "is more
# variance explained?" but "does demand carry more weight RELATIVE to supply?".
# The second form survives a change in instrument quality; the first does not.
# These ratios are the primary evidence on that question.
#
# Restricted to cells where total_r2 is large enough for a ratio to mean
# anything; below that the denominator approaches zero and the ratio explodes.
R2_MIN <- 0.05

# AI axis: same fine binning and plotting range as the bottom row of Fig. 3.
BW    <- 0.02
BRKS  <- seq(0, 1.6, by = BW)
MIDS  <- head(BRKS, -1) + BW / 2
XLIM  <- c(0, 0.9)

# Zomer classes, for the coarse cross-check table.
ai_breaks <- c(-Inf, 0.03, 0.20, 0.50, 0.65, Inf)
ai_labels <- c("Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid", "Humid")

CATS <- c("no_signal", "soil", "demand", "coupled_pure", "both_unique")

varpart_dir <- file.path("outputs", "varpart_global")
tables_dir  <- if (exists("paths") && !is.null(paths$tables)) paths$tables else
  file.path("outputs", "tables")
file_ai <- if (exists("paths") && !is.null(paths$aridity_index)) paths$aridity_index else
  file.path("data", "processed", "aridity", "ai_1982_2022_period.nc")

nc_varpart <- function(w)
  file.path(varpart_dir, sprintf("varpart_signif_global_2blocks_%s.nc", w))
nc_fdr <- function(w)
  file.path(varpart_dir,
            if (identical(FDR_MODE, "pooled"))
              sprintf("fdr_adjusted_pvalues_pooled_%s.nc", w)
            else
              sprintf("fdr_adjusted_pvalues_%s.nc", w))

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

out_shares    <- file.path(tables_dir, "window_comparison_shares_by_ai.csv")
out_crossover <- file.path(tables_dir, "window_comparison_crossover.csv")
out_byclass   <- file.path(tables_dir, "window_comparison_by_aridity_class.csv")
out_magn      <- file.path(tables_dir, "window_comparison_magnitudes_by_ai.csv")
out_comp      <- file.path(tables_dir, "window_comparison_composition_by_ai.csv")

for (f in c(out_shares, out_crossover, out_byclass, out_magn, out_comp))
  if (file.exists(f))
    stop("REFUSING TO OVERWRITE: ", f, "\n  Delete/rename it manually.", call. = FALSE)

for (w in WINDOWS)
  for (f in c(nc_varpart(w), nc_fdr(w)))
    if (!file.exists(f)) stop("Missing input for window ", w, ": ", f, call. = FALSE)
if (!file.exists(file_ai)) stop("Missing aridity index: ", file_ai, call. = FALSE)

mt_before <- file.info(c(sapply(WINDOWS, nc_varpart), sapply(WINDOWS, nc_fdr)))$mtime

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

set_na <- function(v) { v[v == -9999] <- NA; v }

# Five-category classification, FDR-adjusted. Unlike Fig. 3, no_signal is coded
# explicitly (1) instead of being left NA, so it can be reported: its growth
# under n~20 is exactly what motivates the significant-area denominator.
classify_month <- function(tt, pf, ps, pd) {
  valid <- !is.na(tt)
  fs    <- valid & !is.na(pf) & pf < ALPHA
  ss    <- !is.na(ps) & ps < ALPHA
  ds    <- !is.na(pd) & pd < ALPHA
  code  <- rep(NA_integer_, length(tt))
  code[valid & !fs]  <- 1L
  code[fs &  ss & !ds] <- 2L
  code[fs & !ss &  ds] <- 3L
  code[fs & !ss & !ds] <- 4L
  code[fs &  ss &  ds] <- 5L
  code
}

# Area of each category within each AI bin, for the cells selected by `keep`.
area_by_bin <- function(code, area, bin, keep, k) {
  sel <- keep & !is.na(code) & code == k
  s <- tapply(area[sel], bin[sel], sum)
  s[is.na(s)] <- 0
  as.numeric(s)
}

# Area-weighted quantiles, same convention as the violin panel of Fig. 2.
wquantile <- function(x, w, probs) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- (cumsum(w) - 0.5 * w) / sum(w)
  approx(cw, x, xout = probs, rule = 2, ties = "ordered")$y
}

# First AI bin at which the joint class overtakes soil control. Same rule as
# Fig. 3: the first upward crossing, not merely the first bin where joint > soil.
crossover_ai <- function(soil, joint, den) {
  sf <- 100 * soil  / ifelse(den > 0, den, NA)
  jf <- 100 * joint / ifelse(den > 0, den, NA)
  d  <- jf - sf
  cx <- which(d > 0 & c(FALSE, head(d, -1) <= 0) & den > 0)
  if (length(cx)) MIDS[cx[1]] else NA_real_
}

# ------------------------------------------------------------------------------
# Load grids
# ------------------------------------------------------------------------------

message("Window comparison: ", paste(WINDOWS, collapse = "  vs  "))
message("Loading rasters...")

message("  FDR mode: ", FDR_MODE)

r_tot <- lapply(WINDOWS, function(w) rast(nc_varpart(w), subds = "total_r2"))
r_pf  <- lapply(WINDOWS, function(w) rast(nc_fdr(w), subds = "p_full_adj"))
r_ps  <- lapply(WINDOWS, function(w) rast(nc_fdr(w), subds = "p_soil_adj"))
r_pd  <- lapply(WINDOWS, function(w) rast(nc_fdr(w), subds = "p_demand_adj"))
names(r_tot) <- names(r_pf) <- names(r_ps) <- names(r_pd) <- WINDOWS

# Continuous fractions, one SpatRaster per window per variable.
r_mag <- setNames(lapply(MAGNITUDES, function(v)
  setNames(lapply(WINDOWS, function(w) rast(nc_varpart(w), subds = v)), WINDOWS)),
  MAGNITUDES)

ref <- r_tot[[1]][[1]]
if (is.na(crs(ref)) || crs(ref) == "") crs(ref) <- "EPSG:4326"

for (w in WINDOWS)
  if (!all(dim(r_tot[[w]])[1:2] == dim(ref)[1:2]))
    stop("Grid mismatch between windows: ", w, call. = FALSE)

area_v <- as.vector(values(cellSize(ref, unit = "km", mask = FALSE)))

ai <- rast(file_ai)[[1]]
if (is.na(crs(ai))) crs(ai) <- "EPSG:4326"
if (xmax(ai) > 180 + 1e-6) ai <- rotate(ai)
ai_v  <- as.vector(values(resample(ai, ref, method = "bilinear")))
bin_v <- cut(ai_v, breaks = BRKS, right = FALSE, include.lowest = TRUE)
cls_v <- cut(ai_v, breaks = ai_breaks, labels = ai_labels)

# ------------------------------------------------------------------------------
# Per month: common cell x month mask, then bin along AI
# ------------------------------------------------------------------------------

share_rows <- list(); cross_rows <- list(); class_rows <- list()
magn_rows  <- list(); comp_rows <- list(); drop_log <- list()

for (mo in 1:12) {
  message(sprintf("  month %2d ...", mo))

  code <- list(); valid <- list()
  for (w in WINDOWS) {
    tt <- set_na(values(r_tot[[w]][[mo]])[, 1])
    pf <- set_na(values(r_pf[[w]][[mo]])[, 1])
    ps <- set_na(values(r_ps[[w]][[mo]])[, 1])
    pd <- set_na(values(r_pd[[w]][[mo]])[, 1])
    code[[w]]  <- classify_month(tt, pf, ps, pd)
    valid[[w]] <- !is.na(code[[w]])
  }

  # THE common mask: analysed in both windows, and with a defined AI.
  in_ai <- !is.na(bin_v)
  keep  <- valid[[WINDOWS[1]]] & valid[[WINDOWS[2]]] & in_ai

  # Two different losses, reported separately on purpose. Conflating them hides
  # the only one that threatens the comparison:
  #   - the AI restriction drops cells with no aridity index; it hits both
  #     windows identically and is not a comparability problem;
  #   - the window intersection is the one that matters, because a cell present
  #     in one window and absent in the other would contribute a difference in
  #     sampling rather than in coupling.
  a_of <- function(sel) sum(area_v[sel], na.rm = TRUE) / 1e6
  drop_log[[length(drop_log) + 1L]] <- data.frame(
    month = mo,
    w1_analysed = a_of(valid[[WINDOWS[1]]]),
    w2_analysed = a_of(valid[[WINDOWS[2]]]),
    w1_with_ai  = a_of(valid[[WINDOWS[1]]] & in_ai),
    w2_with_ai  = a_of(valid[[WINDOWS[2]]] & in_ai),
    common      = a_of(keep),
    only_w1     = a_of(valid[[WINDOWS[1]]] & !valid[[WINDOWS[2]]] & in_ai),
    only_w2     = a_of(valid[[WINDOWS[2]]] & !valid[[WINDOWS[1]]] & in_ai),
    stringsAsFactors = FALSE)

  for (w in WINDOWS) {
    a <- sapply(1:5, function(k) area_by_bin(code[[w]], area_v, bin_v, keep, k))
    colnames(a) <- CATS

    den3 <- a[, "soil"] + a[, "coupled_pure"] + a[, "demand"]
    den4 <- den3 + a[, "both_unique"]
    denA <- den4 + a[, "no_signal"]

    ok <- denA > 0 & MIDS >= XLIM[1] & MIDS <= XLIM[2]

    share_rows[[length(share_rows) + 1L]] <- data.frame(
      window = w, month = mo, AI = MIDS[ok],
      area_total_km2 = denA[ok],
      pct_no_signal_of_analysed = 100 * a[ok, "no_signal"] / denA[ok],
      soil_share3         = 100 * a[ok, "soil"]         / ifelse(den3[ok] > 0, den3[ok], NA),
      joint_share3        = 100 * a[ok, "coupled_pure"] / ifelse(den3[ok] > 0, den3[ok], NA),
      demand_share3       = 100 * a[ok, "demand"]       / ifelse(den3[ok] > 0, den3[ok], NA),
      soil_share4         = 100 * a[ok, "soil"]         / ifelse(den4[ok] > 0, den4[ok], NA),
      demand_share4       = 100 * a[ok, "demand"]       / ifelse(den4[ok] > 0, den4[ok], NA),
      joint_share4        = 100 * a[ok, "coupled_pure"] / ifelse(den4[ok] > 0, den4[ok], NA),
      both_share4         = 100 * a[ok, "both_unique"]  / ifelse(den4[ok] > 0, den4[ok], NA),
      stringsAsFactors = FALSE)

    cross_rows[[length(cross_rows) + 1L]] <- data.frame(
      window = w, month = mo,
      ai_crossover = crossover_ai(a[, "soil"], a[, "coupled_pure"], den3),
      stringsAsFactors = FALSE)

    # ── Continuous magnitudes, no thresholding ────────────────────────────────
    # Area-weighted median and quartiles per AI bin. Immune to the power problem
    # that limits the category comparison above.
    for (v in MAGNITUDES) {
      val <- set_na(values(r_mag[[v]][[w]][[mo]])[, 1])
      for (b in which(ok)) {
        sel <- keep & !is.na(bin_v) & as.integer(bin_v) == b & !is.na(val)
        if (!any(sel)) next
        q <- wquantile(val[sel], area_v[sel], c(0.25, 0.50, 0.75))
        magn_rows[[length(magn_rows) + 1L]] <- data.frame(
          window = w, month = mo, AI = MIDS[b], variable = v,
          area_km2 = sum(area_v[sel]),
          p25 = q[1], median = q[2], p75 = q[3],
          stringsAsFactors = FALSE)
      }
    }

    # ── Composition: ratios between components, immune to response noise ──────
    tt <- set_na(values(r_mag[["total_r2"]][[w]][[mo]])[, 1])
    uS <- set_na(values(r_mag[["unique_SUMINISTRO"]][[w]][[mo]])[, 1])
    uD <- set_na(values(r_mag[["unique_DEMANDA"]][[w]][[mo]])[, 1])
    iS <- set_na(values(r_mag[["indiv_SUMINISTRO"]][[w]][[mo]])[, 1])
    iD <- set_na(values(r_mag[["indiv_DEMANDA"]][[w]][[mo]])[, 1])
    sh <- set_na(values(r_mag[["shared"]][[w]][[mo]])[, 1])

    # A negative unique Adj.R2 means the block adds nothing once the other is in
    # the model; for a share it is floored at zero rather than propagated.
    uSp <- pmax(uS, 0); uDp <- pmax(uD, 0)

    base <- keep & !is.na(tt) & tt >= R2_MIN

    for (b in which(ok)) {
      sel <- base & as.integer(bin_v) == b
      if (!any(sel)) next
      a <- area_v[sel]
      den_u <- uSp[sel] + uDp[sel]
      ok_u  <- den_u > 0
      ok_i  <- iS[sel] > 0 & !is.na(iD[sel])
      comp_rows[[length(comp_rows) + 1L]] <- data.frame(
        window = w, month = mo, AI = MIDS[b],
        area_km2 = sum(a),
        # Each component as a share of the total explained variance.
        supply_of_total = wquantile(uS[sel] / tt[sel], a, 0.50),
        demand_of_total = wquantile(uD[sel] / tt[sel], a, 0.50),
        shared_of_total = wquantile(sh[sel] / tt[sel], a, 0.50),
        # THE statistic for "is demand gaining weight relative to supply?":
        # demand's share of the resolvable (unique) variance. 0.5 = parity.
        demand_weight = if (any(ok_u))
          wquantile((uDp[sel] / den_u)[ok_u], a[ok_u], 0.50) else NA_real_,
        # Same question via the marginal fractions, which are more often
        # positive and therefore a steadier ratio.
        indiv_ratio_DS = if (any(ok_i))
          wquantile((iD[sel] / iS[sel])[ok_i], a[ok_i], 0.50) else NA_real_,
        stringsAsFactors = FALSE)
    }

    # Coarse cross-check on the Zomer classes, on the same common mask.
    for (cl in ai_labels) {
      sel <- keep & !is.na(cls_v) & cls_v == cl
      ak  <- sapply(1:5, function(k) sum(area_v[sel & !is.na(code[[w]]) & code[[w]] == k],
                                         na.rm = TRUE))
      sig <- sum(ak[2:5])
      if (sig <= 0) next
      class_rows[[length(class_rows) + 1L]] <- data.frame(
        window = w, month = mo, ai_class = cl,
        area_analysed_Mkm2 = sum(ak) / 1e6,
        area_significant_Mkm2 = sig / 1e6,
        pct_no_signal_of_analysed = 100 * ak[1] / sum(ak),
        soil_of_sig         = 100 * ak[2] / sig,
        demand_of_sig       = 100 * ak[3] / sig,
        coupled_pure_of_sig = 100 * ak[4] / sig,
        both_unique_of_sig  = 100 * ak[5] / sig,
        stringsAsFactors = FALSE)
    }
  }
}

shares_df <- do.call(rbind, share_rows)
cross_df  <- do.call(rbind, cross_rows)
class_df  <- do.call(rbind, class_rows)
magn_df   <- do.call(rbind, magn_rows)
comp_df   <- do.call(rbind, comp_rows)
drop_df   <- do.call(rbind, drop_log)

write.csv(shares_df, out_shares,    row.names = FALSE)
write.csv(cross_df,  out_crossover, row.names = FALSE)
write.csv(class_df,  out_byclass,   row.names = FALSE)
write.csv(magn_df,   out_magn,      row.names = FALSE)
write.csv(comp_df,   out_comp,      row.names = FALSE)

# ------------------------------------------------------------------------------
# Console report
# ------------------------------------------------------------------------------

cat("\n----------------------------------------------------------------\n")
cat("A. Area accounting (Mkm2), the two losses separated\n")
cat("----------------------------------------------------------------\n")

# Cost of requiring a defined aridity index (identical for both windows).
drop_df$pct_kept_ai <- round(100 * drop_df$w1_with_ai / drop_df$w1_analysed, 1)
# Cost of requiring the cell x month to be analysed in BOTH windows. This is the
# one that could bias the comparison; measured within the AI domain only.
drop_df$pct_kept_common <- round(100 * drop_df$common / drop_df$w1_with_ai, 1)

print(drop_df[, c("month", "w1_analysed", "w1_with_ai", "common",
                  "only_w1", "only_w2", "pct_kept_ai", "pct_kept_common")],
      row.names = FALSE, digits = 4)

cat(sprintf("\n  Aridity-index domain      : keeps %.1f%% of the analysed area (median)\n",
            median(drop_df$pct_kept_ai)))
cat(sprintf("  Common cell x month mask  : keeps %.1f%% of that (median)\n",
            median(drop_df$pct_kept_common)))
cat(sprintf("  Asymmetry: %.3f Mkm2 only in %s, %.3f Mkm2 only in %s (annual mean)\n",
            mean(drop_df$only_w1), WINDOWS[1],
            mean(drop_df$only_w2), WINDOWS[2]))
cat("  Only the second line and the asymmetry bear on comparability. The first\n")
cat("  is the aridity domain and hits both windows identically.\n")

cat("\n----------------------------------------------------------------\n")
cat("B. Aridity index at which coupled_pure overtakes soil control\n")
cat("----------------------------------------------------------------\n")
cw <- reshape(cross_df, idvar = "month", timevar = "window", direction = "wide")
names(cw) <- sub("^ai_crossover\\.", "", names(cw))
cw$shift <- cw[[WINDOWS[2]]] - cw[[WINDOWS[1]]]
print(cw, row.names = FALSE, digits = 3)
cat(sprintf("\n  Median shift (%s - %s): %+.4f AI units\n",
            WINDOWS[2], WINDOWS[1], median(cw$shift, na.rm = TRUE)))
cat("  Negative = the non-separable regime starts at a WETTER aridity index,\n")
cat("  i.e. it has widened towards the humid end.\n")

cat("\n----------------------------------------------------------------\n")
cat("C. Composition per Zomer class, over SIGNIFICANT area, annual mean\n")
cat("----------------------------------------------------------------\n")
agg <- aggregate(cbind(pct_no_signal_of_analysed, soil_of_sig, demand_of_sig,
                       coupled_pure_of_sig, both_unique_of_sig) ~ window + ai_class,
                 data = class_df, FUN = mean)
agg$ai_class <- factor(agg$ai_class, levels = ai_labels)
agg <- agg[order(agg$ai_class, agg$window), ]
print(agg, row.names = FALSE, digits = 3)

cat("\n  Change (window 2 - window 1), percentage points of significant area:\n")
for (cl in ai_labels) {
  a1 <- agg[agg$ai_class == cl & agg$window == WINDOWS[1], ]
  a2 <- agg[agg$ai_class == cl & agg$window == WINDOWS[2], ]
  if (nrow(a1) == 0 || nrow(a2) == 0) next
  cat(sprintf("    %-14s  soil %+6.2f   demand %+6.2f   joint %+6.2f   both %+6.2f\n",
              cl,
              a2$soil_of_sig - a1$soil_of_sig,
              a2$demand_of_sig - a1$demand_of_sig,
              a2$coupled_pure_of_sig - a1$coupled_pure_of_sig,
              a2$both_unique_of_sig - a1$both_unique_of_sig))
}

cat("\n----------------------------------------------------------------\n")
cat("D. Continuous Adj.R2 fractions (no thresholding), annual mean of the\n")
cat("   area-weighted median, by aridity band\n")
cat("----------------------------------------------------------------\n")

magn_df$band <- cut(magn_df$AI, breaks = ai_breaks, labels = ai_labels)
mag_agg <- aggregate(median ~ window + variable + band, data = magn_df, FUN = mean)

for (v in MAGNITUDES) {
  cat(sprintf("\n  %s\n", v))
  cat(sprintf("    %-16s %12s %12s %10s\n", "band", WINDOWS[1], WINDOWS[2], "change"))
  for (b in ai_labels) {
    m1 <- mag_agg$median[mag_agg$variable == v & mag_agg$band == b & mag_agg$window == WINDOWS[1]]
    m2 <- mag_agg$median[mag_agg$variable == v & mag_agg$band == b & mag_agg$window == WINDOWS[2]]
    if (!length(m1) || !length(m2)) next
    cat(sprintf("    %-16s %12.4f %12.4f %+10.4f\n", b, m1, m2, m2 - m1))
  }
}

cat("\n  These carry no significance threshold, so they are unaffected by the\n")
cat("  power loss that limits sections B and C. They ARE affected by any change\n")
cat("  in the noise of the response between windows: see section E.\n")

cat("\n----------------------------------------------------------------\n")
cat("E. Composition (ratios between components) — immune to response noise\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("   Restricted to cells with total_r2 >= %.2f\n", R2_MIN))

comp_df$band <- cut(comp_df$AI, breaks = ai_breaks, labels = ai_labels)
cmp <- aggregate(cbind(supply_of_total, demand_of_total, shared_of_total,
                       demand_weight, indiv_ratio_DS) ~ window + band,
                 data = comp_df, FUN = function(z) mean(z, na.rm = TRUE))
area_by_band <- aggregate(area_km2 ~ band, data = comp_df,
                          FUN = function(z) sum(z) / 12 / 1e6)

for (b in ai_labels) {
  c1 <- cmp[cmp$band == b & cmp$window == WINDOWS[1], ]
  c2 <- cmp[cmp$band == b & cmp$window == WINDOWS[2], ]
  if (nrow(c1) == 0 || nrow(c2) == 0) next
  ab <- area_by_band$area_km2[area_by_band$band == b]
  cat(sprintf("\n  %s   (%.2f Mkm2 above the total_r2 floor)\n", b,
              if (length(ab)) ab else NA_real_))
  cat(sprintf("    %-22s %11s %11s %10s\n", "", WINDOWS[1], WINDOWS[2], "change"))
  for (nm in c("supply_of_total", "demand_of_total", "shared_of_total",
               "demand_weight", "indiv_ratio_DS"))
    cat(sprintf("    %-22s %11.4f %11.4f %+10.4f\n",
                nm, c1[[nm]], c2[[nm]], c2[[nm]] - c1[[nm]]))
}

cat("\n  demand_weight is the direct answer to \"is demand gaining weight\":\n")
cat("  demand's share of the resolvable variance, 0.5 being parity with supply.\n")
cat("  A flat value across windows refutes the hypothesis in a way that a change\n")
cat("  in instrument noise cannot produce, because such a change scales every\n")
cat("  component identically and cancels in the ratio.\n")

cat("\n  Reminder: these are two 20-year windows against each other. Neither is\n")
cat("  comparable with the 41-year full-record result, where n is twice as large\n")
cat("  and the significant fraction is correspondingly larger.\n")

cat("\nTables written:\n")
cat("  ", out_shares, "\n")
cat("  ", out_crossover, "\n")
cat("  ", out_byclass, "\n")
cat("  ", out_magn, "\n")

mt_after <- file.info(c(sapply(WINDOWS, nc_varpart), sapply(WINDOWS, nc_fdr)))$mtime
cat("\n=== File protection (inputs READ-ONLY) ===\n")
cat(sprintf("  source NetCDFs : %s\n",
            if (identical(mt_before, mt_after)) "UNCHANGED" else "*** CHANGED ***"))
cat("\n=== DONE ===\n")
