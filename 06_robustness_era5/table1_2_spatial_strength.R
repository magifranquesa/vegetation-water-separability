#!/usr/bin/env Rscript

# ==============================================================================
# Script: 06_robustness_era5/table1_2_spatial_strength.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# TWO JOBS IN ONE PASS
#   Both need the same inputs (abs_max rho*, FDR-adjusted p, vegetated mask,
#   aridity index), and those NetCDFs total ~11 GB, so they are read once.
#
#   (1) T1-FIX. The first version of T1 measured sign agreement over ALL
#       vegetated cells. That is diluted by construction: ~65% of cells are not
#       significant, rho* sits near zero there and its sign is close to random.
#       With ~30% of cells significant and agreeing almost always, and 70% noise
#       agreeing at chance, the expected value is ~65%; we observed 75-80%, i.e.
#       already better than noise. The fix is NOT to lower the threshold after
#       the fact but to measure sign where it carries information: restricted to
#       cells significant in BOTH products. The unconditioned spatial
#       correlation is kept, since it is the one metric that cannot be accused
#       of conditioning on the significance being compared. Both are reported.
#
#   (2) T2. Association strength along the aridity gradient: area-weighted mean
#       |rho*| per aridity class and indicator, in both products, plus the
#       significant area fraction within each class. Tests whether the arid
#       maximum in association strength survives the change of product.
#
#   Also emitted: vegetated area by aridity class, which is the quantity behind
#   the Abstract's "most water-controlled vegetation lies in semi-arid
#   ecosystems". NOTE: this is the PLAIN vegetated area per class. The
#   "water-controlled" restriction is a variance-partitioning concept and needs
#   T3; do not quote this table as if it were that.
#
# Aridity classes come from the SAME CRU-derived file used for GLEAM, so the
# stratification is identical for both products and cannot itself explain a
# difference.
#
# Inputs (READ-ONLY):
#   outputs/intermediate/absmax_spearman/abs_max_correlation_kndvi_<IND>.nc
#   outputs/intermediate/absmax_spearman_era5/abs_max_correlation_kndvi_<IND>_era5.nc
#   outputs/intermediate/absmax_fdr/fdr_padj_<IND>.tif
#   outputs/intermediate/absmax_fdr_era5/fdr_padj_<IND>_era5.tif
#   outputs/intermediate/vegetation_mask_c1.tif
#   data/processed/aridity/ai_1982_2022_period.nc
#
# Outputs (refuse to overwrite):
#   outputs/tables/robustness_T1_spatial_agreement_era5_v2.csv
#   outputs/tables/robustness_T2_strength_by_aridity_full.csv
#   outputs/tables/robustness_T2_strength_by_aridity_comparison.csv
#   outputs/tables/robustness_vegetated_area_by_aridity.csv
#
# Run from repo root (ideally NOT while table4_coupling.R is running: both are disk-bound):
#   Rscript 06_robustness_era5/table1_2_spatial_strength.R
#
# Dependencies: terra
# ==============================================================================

suppressPackageStartupMessages({ library(terra) })

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

dir_gleam     <- file.path("outputs", "intermediate", "absmax_spearman")
dir_era5      <- file.path("outputs", "intermediate", "absmax_spearman_era5")
dir_fdr_gleam <- file.path("outputs", "intermediate", "absmax_fdr")
dir_fdr_era5  <- file.path("outputs", "intermediate", "absmax_fdr_era5")
file_mask     <- file.path("outputs", "intermediate", "vegetation_mask_c1.tif")
file_ai       <- file.path("data", "processed", "aridity", "ai_1982_2022_period.nc")

tables_dir <- if (exists("paths") && !is.null(paths$tables)) paths$tables else
  file.path("outputs", "tables")

out_t1  <- file.path(tables_dir, "robustness_T1_spatial_agreement_era5_v2.csv")
out_t2f <- file.path(tables_dir, "robustness_T2_strength_by_aridity_full.csv")
out_t2c <- file.path(tables_dir, "robustness_T2_strength_by_aridity_comparison.csv")
out_veg <- file.path(tables_dir, "robustness_vegetated_area_by_aridity.csv")
for (f in c(out_t1, out_t2f, out_t2c, out_veg))
  if (file.exists(f)) stop("REFUSING TO OVERWRITE: ", f, call. = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.05

indicators  <- c("Ep", "Et", "ED", "SMrz", "SMs")
ind_labels  <- c(Ep = "AED", Et = "Et", ED = "ED", SMrz = "SMrz", SMs = "SMs")
months_idx  <- c(1, 4, 7, 10)
month_names <- c("1" = "January", "4" = "April", "7" = "July", "10" = "October")

ai_breaks <- c(-Inf, 0.03, 0.20, 0.50, 0.65, Inf)
ai_labels <- c("Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid", "Humid")

check <- function(fp, d) if (!file.exists(fp)) stop("Missing ", d, ": ", fp, call. = FALSE)
for (f in c(file_mask, file_ai)) check(f, "input")

# ------------------------------------------------------------------------------
# Shared geometry: vegetated domain, cell area, aridity class
# ------------------------------------------------------------------------------

ref <- rast(file.path(dir_gleam, "abs_max_correlation_kndvi_Ep.nc"),
            subds = "abs_max_correlation")[[1]]

veg_v  <- as.vector(values(resample(rast(file_mask), ref, method = "near")))
area_v <- as.vector(values(cellSize(ref, unit = "km", mask = FALSE)))

ai <- rast(file_ai)[[1]]
if (is.na(crs(ai))) crs(ai) <- "EPSG:4326"
if (xmax(ai) > 180 + 1e-6) ai <- rotate(ai)
ai_v  <- as.vector(values(resample(ai, ref, method = "bilinear")))
cls_v <- cut(ai_v, breaks = ai_breaks, labels = ai_labels, right = FALSE)

in_veg <- !is.na(veg_v) & veg_v > 0
veg_area_total <- sum(area_v[in_veg], na.rm = TRUE)

message(sprintf("Vegetated domain: %d cells, %.2f x 10^6 km2",
                sum(in_veg), veg_area_total / 1e6))
message(sprintf("  of which with valid AI: %d cells (%.2f x 10^6 km2)",
                sum(in_veg & !is.na(cls_v)),
                sum(area_v[in_veg & !is.na(cls_v)], na.rm = TRUE) / 1e6))

wsum <- function(sel) sum(area_v[sel], na.rm = TRUE)
wmean <- function(x, sel) {
  ok <- sel & is.finite(x)
  sum(x[ok] * area_v[ok], na.rm = TRUE) / sum(area_v[ok], na.rm = TRUE)
}

# ------------------------------------------------------------------------------
# Vegetated area by aridity class
# ------------------------------------------------------------------------------

veg_rows <- lapply(ai_labels, function(cl) {
  s <- in_veg & !is.na(cls_v) & cls_v == cl
  data.frame(ai_class = cl,
             area_Mkm2 = round(wsum(s) / 1e6, 3),
             pct_of_vegetated = round(100 * wsum(s) / veg_area_total, 2),
             stringsAsFactors = FALSE)
})
veg_tab <- do.call(rbind, veg_rows)
write.csv(veg_tab, out_veg, row.names = FALSE)
message("\nVegetated area by aridity class:")
print(veg_tab, row.names = FALSE)

# ------------------------------------------------------------------------------
# Main pass
# ------------------------------------------------------------------------------

t1_rows <- list()
t2_rows <- list()

for (ind in indicators) {

  f_g  <- file.path(dir_gleam, sprintf("abs_max_correlation_kndvi_%s.nc", ind))
  f_e  <- file.path(dir_era5,  sprintf("abs_max_correlation_kndvi_%s_era5.nc", ind))
  fp_g <- file.path(dir_fdr_gleam, sprintf("fdr_padj_%s.tif", ind))
  fp_e <- file.path(dir_fdr_era5,  sprintf("fdr_padj_%s_era5.tif", ind))
  for (x in list(c(f_g, "GLEAM rho*"), c(f_e, "ERA5 rho*"),
                 c(fp_g, "GLEAM padj"), c(fp_e, "ERA5 padj"))) check(x[1], x[2])

  rg <- rast(f_g, subds = "abs_max_correlation")
  re <- rast(f_e, subds = "abs_max_correlation")
  pg <- rast(fp_g); pe <- rast(fp_e)

  for (mi in months_idx) {
    mn <- month_names[[as.character(mi)]]
    message("  ", ind_labels[[ind]], " - ", mn)

    g  <- as.vector(values(rg[[mi]]));  e  <- as.vector(values(re[[mi]]))
    sg <- as.vector(values(pg[[mi]]));  se <- as.vector(values(pe[[mi]]))

    both  <- in_veg & !is.na(g) & !is.na(e)
    sig_g <- in_veg & !is.na(sg) & sg < ALPHA
    sig_e <- in_veg & !is.na(se) & se < ALPHA
    sig_b <- sig_g & sig_e & both

    # --- T1 (v2) --------------------------------------------------------------
    t1_rows[[length(t1_rows) + 1]] <- data.frame(
      indicator = ind_labels[[ind]], month = mn,
      n_cells_compared = sum(both),
      spatial_rho = round(suppressWarnings(cor(g[both], e[both], method = "spearman")), 3),
      pct_same_sign_all = round(100 * wsum(both & sign(g) == sign(e)) / wsum(both), 1),
      # the meaningful sign metric: where rho* is significant in BOTH products
      pct_same_sign_sig_both = round(100 * wsum(sig_b & sign(g) == sign(e)) / wsum(sig_b), 1),
      pct_sig_gleam = round(100 * wsum(sig_g) / veg_area_total, 1),
      pct_sig_era5  = round(100 * wsum(sig_e) / veg_area_total, 1),
      pct_sig_both  = round(100 * wsum(sig_b) / veg_area_total, 1),
      mean_abs_rho_gleam = round(mean(abs(g[both]), na.rm = TRUE), 3),
      mean_abs_rho_era5  = round(mean(abs(e[both]), na.rm = TRUE), 3),
      stringsAsFactors = FALSE)

    # --- T2 -------------------------------------------------------------------
    for (cl in ai_labels) {
      in_cl <- in_veg & !is.na(cls_v) & cls_v == cl
      if (!any(in_cl)) next
      cl_area <- wsum(in_cl)
      t2_rows[[length(t2_rows) + 1]] <- data.frame(
        indicator = ind_labels[[ind]], month = mn, ai_class = cl,
        mean_abs_rho_gleam = round(wmean(abs(g), in_cl & !is.na(g)), 3),
        mean_abs_rho_era5  = round(wmean(abs(e), in_cl & !is.na(e)), 3),
        pct_sig_gleam = round(100 * wsum(in_cl & sig_g) / cl_area, 1),
        pct_sig_era5  = round(100 * wsum(in_cl & sig_e) / cl_area, 1),
        class_area_Mkm2 = round(cl_area / 1e6, 3),
        stringsAsFactors = FALSE)
    }

    rm(g, e, sg, se, both, sig_g, sig_e, sig_b); gc(verbose = FALSE)
  }
  rm(rg, re, pg, pe); gc(verbose = FALSE)
}

t1 <- do.call(rbind, t1_rows)
t2 <- do.call(rbind, t2_rows)
write.csv(t1, out_t1,  row.names = FALSE)
write.csv(t2, out_t2f, row.names = FALSE)

# T2 comparison: averaged over the four representative months
cmp <- aggregate(cbind(mean_abs_rho_gleam, mean_abs_rho_era5,
                       pct_sig_gleam, pct_sig_era5) ~ indicator + ai_class,
                 data = t2, FUN = mean)
cmp$diff_abs_rho <- round(cmp$mean_abs_rho_era5 - cmp$mean_abs_rho_gleam, 3)
for (v in c("mean_abs_rho_gleam", "mean_abs_rho_era5"))  cmp[[v]] <- round(cmp[[v]], 3)
for (v in c("pct_sig_gleam", "pct_sig_era5"))            cmp[[v]] <- round(cmp[[v]], 1)
cmp$ai_class  <- factor(cmp$ai_class, levels = ai_labels)
cmp$indicator <- factor(cmp$indicator, levels = unname(ind_labels))
cmp <- cmp[order(cmp$indicator, cmp$ai_class), ]
write.csv(cmp, out_t2c, row.names = FALSE)

cat("\n================================================================\n")
cat("=== T1 (v2): spatial agreement, with sign measured where it means\n")
cat("================================================================\n\n")
print(t1[, c("indicator", "month", "spatial_rho",
             "pct_same_sign_all", "pct_same_sign_sig_both", "pct_sig_both")],
      row.names = FALSE)

cat("\n================================================================\n")
cat("=== T2: |rho*| by aridity class, averaged over the four months\n")
cat("================================================================\n\n")
print(cmp, row.names = FALSE)

cat("\n  >> T2 claim under test: the ARID class retains the highest |rho*|\n")
cat("     for the supply indicators, in BOTH products.\n")
cat("\nTables written:\n  ", out_t1, "\n  ", out_t2f, "\n  ", out_t2c, "\n  ", out_veg, "\n")
cat("\n=== DONE ===\n")
