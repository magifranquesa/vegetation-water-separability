#!/usr/bin/env Rscript

# ==============================================================================
# Script: 06_robustness_era5/table3_varpart_categories.R
#
# ROBUSTNESS TABLE T3 - composition of hydroclimatic control by aridity class,
# GLEAM v4.2a vs ERA5-Land.
#
# The classification is IDENTICAL to 03_variance_partitioning/03_category_split.R
# (five FDR-adjusted, area-weighted, mutually exclusive categories). The only
# thing that changes between the two runs is which pair of NetCDFs is read:
#
#   GLEAM : varpart_signif_global_2blocks.nc      + fdr_adjusted_pvalues.nc
#   ERA5  : varpart_signif_global_2blocks_era5.nc + fdr_adjusted_pvalues_era5.nc
#
# Both NetCDFs share the analysis grid (3600 x 1550 x 12), and the aridity index
# is the SAME CRU-derived file for both, so the stratification cannot itself
# explain a difference.
#
# Categories (as in 03_variance_partitioning/03_category_split.R):
#   no_signal    full model not significant
#   soil         full sig, ONLY soil unique significant
#   demand       full sig, ONLY demand unique significant
#   coupled_pure full sig, NEITHER unique significant   <- true non-separability
#   both_unique  full sig, BOTH unique significant      <- separable
#
# Output (refuses to overwrite):
#   outputs/tables/robustness_T3_varpart_categories_comparison.csv
#
# Run from repo root (AFTER both FDR files are present locally):
#   Rscript 06_robustness_era5/table3_varpart_categories.R
#
# Dependencies: terra
# ==============================================================================

suppressPackageStartupMessages({ library(terra) })

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

tables_dir <- if (exists("paths") && !is.null(paths$tables)) paths$tables else
  file.path("outputs", "tables")
file_ai <- if (exists("paths") && !is.null(paths$aridity_index)) paths$aridity_index else
  file.path("data", "processed", "aridity", "ai_1982_2022_period.nc")

vp_dir <- file.path("outputs", "varpart_global")
products <- list(
  GLEAM       = list(nc = file.path(vp_dir, "varpart_signif_global_2blocks.nc"),
                     adj = file.path(vp_dir, "fdr_adjusted_pvalues.nc")),
  `ERA5-Land` = list(nc = file.path(vp_dir, "varpart_signif_global_2blocks_era5.nc"),
                     adj = file.path(vp_dir, "fdr_adjusted_pvalues_era5.nc"))
)

out_csv <- file.path(tables_dir, "robustness_T3_varpart_categories_comparison.csv")
if (file.exists(out_csv)) stop("REFUSING TO OVERWRITE: ", out_csv, call. = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

for (p in products) for (f in unlist(p))
  if (!file.exists(f)) stop("Missing input: ", f, call. = FALSE)
if (!file.exists(file_ai)) stop("Missing aridity file: ", file_ai, call. = FALSE)

FV    <- -9999
ALPHA <- 0.05
CATS  <- c("no_signal", "soil", "demand", "coupled_pure", "both_unique")
ai_breaks <- c(-Inf, 0.03, 0.20, 0.50, 0.65, Inf)
ai_labels <- c("Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid", "Humid")

set_na <- function(v) { v[abs(v - FV) < 1] <- NA; v }

# ------------------------------------------------------------------------------
# Shared geometry (grid identical for both products)
# ------------------------------------------------------------------------------

ref <- rast(products$GLEAM$nc, subds = "shared")[[1]]
if (is.na(crs(ref)) || crs(ref) == "") crs(ref) <- "EPSG:4326"
area_v <- as.vector(values(cellSize(ref, unit = "km", mask = FALSE)))

ai <- rast(file_ai)[[1]]
if (is.na(crs(ai))) crs(ai) <- "EPSG:4326"
if (xmax(ai) > 180 + 1e-6) ai <- rotate(ai)
ai_v   <- as.vector(values(resample(ai, ref, method = "bilinear")))
ai_cls <- cut(ai_v, breaks = ai_breaks, labels = ai_labels, right = FALSE)

# ------------------------------------------------------------------------------
# Classification as in 03_variance_partitioning/03_category_split.R
# ------------------------------------------------------------------------------

classify_all <- function(nc, adj) {
  r_tot <- rast(nc,  subds = "total_r2")
  r_sh  <- rast(nc,  subds = "shared")
  r_pf  <- rast(adj, subds = "p_full_adj")
  r_ps  <- rast(adj, subds = "p_soil_adj")
  r_pd  <- rast(adj, subds = "p_demand_adj")

  code_all <- integer(0); area_all <- numeric(0); ai_all <- numeric(0)
  cls_all  <- character(0)
  for (m in 1:12) {
    tt <- set_na(values(r_tot[[m]])[, 1])
    pf <- set_na(values(r_pf[[m]])[, 1])
    ps <- set_na(values(r_ps[[m]])[, 1])
    pd <- set_na(values(r_pd[[m]])[, 1])

    valid    <- !is.na(tt)
    full_sig <- valid & !is.na(pf) & pf < ALPHA
    soil_sig <- !is.na(ps) & ps < ALPHA
    dem_sig  <- !is.na(pd) & pd < ALPHA

    code <- rep(NA_integer_, length(tt))
    code[valid & !full_sig]               <- 1L
    code[full_sig &  soil_sig & !dem_sig] <- 2L
    code[full_sig & !soil_sig &  dem_sig] <- 3L
    code[full_sig & !soil_sig & !dem_sig] <- 4L
    code[full_sig &  soil_sig &  dem_sig] <- 5L

    keep <- !is.na(code)
    code_all <- c(code_all, code[keep])
    area_all <- c(area_all, area_v[keep])
    ai_all   <- c(ai_all,   ai_v[keep])
    cls_all  <- c(cls_all,  as.character(ai_cls[keep]))
  }
  list(code = code_all, area = area_all, ai = ai_all, cls = cls_all)
}

# soil -> coupled_pure crossover, 0.02-wide AI bins
crossover <- function(d) {
  fb <- seq(0, 1, by = 0.02); mids <- head(fb, -1) + 0.01
  bin <- cut(d$ai, breaks = fb, right = FALSE, include.lowest = TRUE)
  soil_b <- tapply(d$area[d$code == 2L], bin[d$code == 2L], sum)
  pure_b <- tapply(d$area[d$code == 4L], bin[d$code == 4L], sum)
  soil_b[is.na(soil_b)] <- 0; pure_b[is.na(pure_b)] <- 0
  df <- as.numeric(pure_b) - as.numeric(soil_b)
  ok <- (as.numeric(soil_b) + as.numeric(pure_b)) > 0
  i  <- which(df > 0 & c(FALSE, head(df, -1) <= 0) & ok)
  if (length(i)) mids[i[1]] else NA_real_
}

# ------------------------------------------------------------------------------
# Run both products
# ------------------------------------------------------------------------------

rows <- list()
crossings <- list()
for (prod in names(products)) {
  message("Classifying ", prod, " ...")
  d <- classify_all(products[[prod]]$nc, products[[prod]]$adj)

  for (cl in ai_labels) {
    sel <- d$cls == cl & !is.na(d$cls)
    if (!any(sel)) next
    a <- sapply(1:5, function(k) sum(d$area[sel & d$code == k]))
    tot <- sum(a)
    rows[[length(rows) + 1L]] <- data.frame(
      product = prod, ai_class = cl,
      area_Mkm2    = round(tot / 12 / 1e6, 3),
      no_signal    = round(a[1] / tot * 100, 1),
      soil         = round(a[2] / tot * 100, 1),
      demand       = round(a[3] / tot * 100, 1),
      coupled_pure = round(a[4] / tot * 100, 1),
      both_unique  = round(a[5] / tot * 100, 1),
      stringsAsFactors = FALSE)
  }
  crossings[[prod]] <- crossover(d)
}

res <- do.call(rbind, rows)
res$ai_class <- factor(res$ai_class, levels = ai_labels)
res <- res[order(res$product, res$ai_class), ]
write.csv(res, out_csv, row.names = FALSE)

cat("\n================================================================\n")
cat("=== T3: control composition by aridity class (% of class area)\n")
cat("===   FDR-adjusted, area-weighted, 12 months\n")
cat("================================================================\n\n")
print(res, row.names = FALSE)

cat("\n--- soil -> coupled_pure crossover along the aridity gradient ---\n")
for (prod in names(crossings))
  cat(sprintf("  %-10s : AI ~ %s\n", prod,
              ifelse(is.na(crossings[[prod]]), "none", sprintf("%.2f", crossings[[prod]]))))

cat("\n  >> T3 claim under test: in BOTH products, coupled_pure is the dominant\n")
cat("     control across the semi-arid range, and the soil->coupled crossover\n")
cat("     falls at a similar AI. The GLEAM (main) analysis gives AI ~ 0.33.\n")
cat(sprintf("\nTable: %s\n", out_csv))
cat("\n=== DONE ===\n")
