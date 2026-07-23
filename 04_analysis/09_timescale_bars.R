#!/usr/bin/env Rscript

# ==============================================================================
# Script: 04_analysis/09_timescale_bars.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Builds the aridity x timescale bar TABLES (stress and non-stress regimes) that
# drive the significance bars and mean-timescale curves of Fig. 5 (a,b).
# Significance uses BH-FDR (p_adj < 0.05): a cell counts as significant when the
# sign of rho* matches the regime and fdr_padj < 0.05. scale* (abs_max_scale)
# gives the timescale at which rho* occurs; the tables aggregate perc,
# sig_perc_total and the mean timescale per aridity class.
#
# Writes the CSV tables.
#
# Inputs (READ-ONLY):
#   outputs/intermediate/absmax_spearman/abs_max_correlation_kndvi_{Ep,Et,ED,SMrz,SMs}.nc
#     (variables abs_max_correlation, abs_max_scale)
#   outputs/intermediate/absmax_fdr/fdr_padj_{...}.tif
#   data/processed/aridity/ai_1982_2022_period.nc
#   outputs/intermediate/vegetation_mask_c1.tif
#
# Outputs:
#   outputs/tables/timescale_bars_stress_fdr.csv
#   outputs/tables/timescale_bars_nonstress_fdr.csv
#
# Run from repo root (AFTER 02_correlations/04_fdr_significance_mask.R):
#   Rscript 04_analysis/09_timescale_bars.R
#
# Dependencies: terra, dplyr
# ==============================================================================

suppressPackageStartupMessages({ library(terra); library(dplyr) })

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

dir_absmax <- if (exists("paths") && !is.null(paths$absmax_spearman)) paths$absmax_spearman else
  file.path("outputs", "intermediate", "absmax_spearman")
dir_fdr <- file.path("outputs", "intermediate", "absmax_fdr")
file_ai <- if (exists("paths") && !is.null(paths$aridity_index)) paths$aridity_index else
  file.path("data", "processed", "aridity", "ai_1982_2022_period.nc")
file_veg_mask <- if (exists("paths") && !is.null(paths$veg_mask_c1)) paths$veg_mask_c1 else
  file.path("outputs", "intermediate", "vegetation_mask_c1.tif")
tables_dir <- if (exists("paths") && !is.null(paths$tables)) paths$tables else file.path("outputs", "tables")

out_stress <- file.path(tables_dir, "timescale_bars_stress_fdr.csv")
out_nonstr <- file.path(tables_dir, "timescale_bars_nonstress_fdr.csv")
for (f in c(out_stress, out_nonstr))
  if (file.exists(f)) stop("REFUSING TO OVERWRITE: ", f, call. = FALSE)

ALPHA <- 0.05
scales_vec      <- c(1, 3, 6, 9, 12)
months_selected <- c(1, 4, 7, 10)
extent_longlat  <- ext(-170, 170, -65, 90)
aridity_labels  <- c("Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid", "Humid")

fdrfile <- function(v) file.path(dir_fdr, sprintf("fdr_padj_%s.tif", v))
ncfile  <- function(v) file.path(dir_absmax, sprintf("abs_max_correlation_kndvi_%s.nc", v))

stress_config <- data.frame(
  ind       = c("Et", "SMrz", "SMs", "ED", "Ep"),
  sign_val  = c( 2L,   2L,     2L,    2L,   -2L),
  face      = c("Et+", "SMrz+", "SMs+", "ED+", "AED-"),
  ind_label = c("Et",  "SMrz",  "SMs",  "ED",  "AED"),
  stringsAsFactors = FALSE)
nonstress_config <- data.frame(
  ind       = c("Ep",  "ED"),
  sign_val  = c( 2L,   -2L),
  face      = c("AED+", "ED-"),
  ind_label = c("AED",  "ED"),
  stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# Aridity classes, veg mask, cell area
# ------------------------------------------------------------------------------

message("Setting up aridity classes, veg mask, cell area...")
ai <- rast(file_ai)[[1]]; crs(ai) <- "EPSG:4326"
ref <- rast(ncfile("Ep"), subds = "abs_max_correlation")[[1]]
ai_res <- crop(resample(ai, ref, method = "bilinear"), extent_longlat)
aridity_class <- classify(ai_res,
  rcl = matrix(c(0, 0.03, 1, 0.03, 0.20, 2, 0.20, 0.50, 3, 0.50, 0.65, 4, 0.65, Inf, 5),
               ncol = 3, byrow = TRUE), include.lowest = TRUE)
cell_area <- cellSize(aridity_class, unit = "km")
veg_res <- crop(resample(rast(file_veg_mask), ref, method = "near"), extent_longlat)

get_area <- function(mask_raster) terra::global(cell_area * mask_raster, "sum", na.rm = TRUE)$sum

veg_area_df <- bind_rows(lapply(1:5, function(cls)
  data.frame(aridity_class = cls, veg_area_km2 = get_area((aridity_class == cls) * veg_res))))
message("  veg class areas (Mkm2): ", paste(round(veg_area_df$veg_area_km2 / 1e6, 2), collapse = ", "))

# ------------------------------------------------------------------------------
# Core: area summaries with FDR significance
# ------------------------------------------------------------------------------

compute_config <- function(cfg) {
  res_scale <- list(); res_sig <- list()
  for (i in seq_len(nrow(cfg))) {
    ind <- cfg$ind[i]; sign_val <- cfg$sign_val[i]; face <- cfg$face[i]; ind_label <- cfg$ind_label[i]
    sign_str <- ifelse(sign_val > 0, "pos", "neg")
    message(sprintf("  %s (%s, %s)", face, ind, sign_str))
    r_cor_all  <- rast(ncfile(ind), subds = "abs_max_correlation")
    r_scale_all<- rast(ncfile(ind), subds = "abs_max_scale")
    r_padj_all <- rast(fdrfile(ind))
    for (m in months_selected) {
      r_cor   <- crop(r_cor_all[[m]],   extent_longlat)
      r_scale <- crop(r_scale_all[[m]], extent_longlat)
      r_padj  <- crop(r_padj_all[[m]],  extent_longlat)
      sign_ok <- if (sign_val > 0) (r_cor > 0) else (r_cor < 0)
      fdr_sig <- !is.na(r_padj) & (r_padj < ALPHA)
      for (cls in 1:5) {
        sig_mask <- sign_ok & fdr_sig & (aridity_class == cls) & (veg_res == 1)
        for (sv in scales_vec)
          res_scale[[length(res_scale) + 1L]] <- data.frame(
            indicator = ind_label, sign = sign_str, month = m, aridity_class = cls,
            scale = sv, area_km2 = get_area(sig_mask & (r_scale == sv)))
        res_sig[[length(res_sig) + 1L]] <- data.frame(
          indicator = ind_label, sign = sign_str, month = m, aridity_class = cls,
          area_sig = get_area(sig_mask))
      }
    }
  }
  df_scale <- do.call(rbind, res_scale) %>%
    group_by(indicator, sign, month, aridity_class) %>%
    mutate(total_sig = sum(area_km2), perc = ifelse(total_sig > 0, 100 * area_km2 / total_sig, 0)) %>%
    ungroup() %>% select(-total_sig)
  df_sig <- do.call(rbind, res_sig) %>%
    left_join(veg_area_df, by = "aridity_class") %>%
    mutate(sig_perc_total = 100 * area_sig / veg_area_km2)
  out <- df_scale %>%
    left_join(df_sig, by = c("indicator", "sign", "month", "aridity_class")) %>%
    mutate(month = month.name[month],
           aridity_class = aridity_labels[aridity_class]) %>%
    select(indicator, sign, month, aridity_class, scale, area_km2, perc,
           area_sig, veg_area_km2, sig_perc_total)
  out
}

message("\n=== STRESS regime (figS3 equivalent) ===")
df_s <- compute_config(stress_config)
write.csv(df_s, out_stress, row.names = FALSE)
message("  wrote: ", out_stress)

message("\n=== NON-STRESS regime (figS4 equivalent) ===")
df_n <- compute_config(nonstress_config)
write.csv(df_n, out_nonstr, row.names = FALSE)
message("  wrote: ", out_nonstr)

cat("\n=== DONE ===\n")
