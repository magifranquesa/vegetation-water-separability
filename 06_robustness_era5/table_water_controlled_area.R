#!/usr/bin/env Rscript

# ==============================================================================
# 06_robustness_era5/table_water_controlled_area.R
#
# One question: where is the water-controlled vegetated area?
#
# "Water-controlled" is taken as the paper defines the water-limited regime in
# Fig. 5a: a SIGNIFICANT and POSITIVE association with a supply variable.
# Reported for SMrz alone (the primary supply variable) and for the union of
# SMrz, SMs and Et, in both products.
#
# Output: outputs/tables/robustness_water_controlled_area_by_aridity.csv
# ==============================================================================

suppressPackageStartupMessages({ library(terra) })

dir_g  <- file.path("outputs", "intermediate", "absmax_spearman")
dir_e  <- file.path("outputs", "intermediate", "absmax_spearman_era5")
fdr_g  <- file.path("outputs", "intermediate", "absmax_fdr")
fdr_e  <- file.path("outputs", "intermediate", "absmax_fdr_era5")
f_mask <- file.path("outputs", "intermediate", "vegetation_mask_c1.tif")
f_ai   <- file.path("data", "processed", "aridity", "ai_1982_2022_period.nc")
out_csv <- file.path("outputs", "tables", "robustness_water_controlled_area_by_aridity.csv")
if (file.exists(out_csv)) stop("REFUSING TO OVERWRITE: ", out_csv, call. = FALSE)

ALPHA      <- 0.05
months_idx <- c(1, 4, 7, 10)
ai_breaks  <- c(-Inf, 0.03, 0.20, 0.50, 0.65, Inf)
ai_labels  <- c("Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid", "Humid")

ref    <- rast(file.path(dir_g, "abs_max_correlation_kndvi_Ep.nc"), subds = "abs_max_correlation")[[1]]
veg    <- as.vector(values(resample(rast(f_mask), ref, method = "near")))
area_v <- as.vector(values(cellSize(ref, unit = "km", mask = FALSE)))
ai <- rast(f_ai)[[1]]
if (is.na(crs(ai))) crs(ai) <- "EPSG:4326"
if (xmax(ai) > 180 + 1e-6) ai <- rotate(ai)
cls <- cut(as.vector(values(resample(ai, ref, method = "bilinear"))),
           breaks = ai_breaks, labels = ai_labels, right = FALSE)
in_veg <- !is.na(veg) & veg > 0

supply <- c("SMrz", "SMs", "Et")

pos_sig <- function(product, ind, mi) {
  d  <- if (product == "GLEAM") dir_g else dir_e
  fp <- if (product == "GLEAM") file.path(fdr_g, sprintf("fdr_padj_%s.tif", ind))
        else                    file.path(fdr_e, sprintf("fdr_padj_%s_era5.tif", ind))
  fnc <- if (product == "GLEAM") sprintf("abs_max_correlation_kndvi_%s.nc", ind)
         else                    sprintf("abs_max_correlation_kndvi_%s_era5.nc", ind)
  r <- as.vector(values(rast(file.path(d, fnc),
                             subds = "abs_max_correlation")[[mi]]))
  p <- as.vector(values(rast(fp)[[mi]]))
  in_veg & !is.na(r) & !is.na(p) & p < ALPHA & r > 0
}

rows <- list()
for (product in c("GLEAM", "ERA5-Land")) {
  for (mi in months_idx) {
    message(product, " - month ", mi)
    m_smrz  <- pos_sig(product, "SMrz", mi)
    m_union <- m_smrz | pos_sig(product, "SMs", mi) | pos_sig(product, "Et", mi)
    for (cl in ai_labels) {
      s <- in_veg & !is.na(cls) & cls == cl
      rows[[length(rows) + 1]] <- data.frame(
        product = product, month = mi, ai_class = cl,
        area_SMrz_pos_Mkm2  = round(sum(area_v[s & m_smrz],  na.rm = TRUE) / 1e6, 3),
        area_union_pos_Mkm2 = round(sum(area_v[s & m_union], na.rm = TRUE) / 1e6, 3),
        class_area_Mkm2     = round(sum(area_v[s], na.rm = TRUE) / 1e6, 3),
        stringsAsFactors = FALSE)
    }
    gc(verbose = FALSE)
  }
}

res <- do.call(rbind, rows)
write.csv(res, out_csv, row.names = FALSE)

agg <- aggregate(cbind(area_SMrz_pos_Mkm2, area_union_pos_Mkm2) ~ product + ai_class,
                 data = res, FUN = mean)
agg$ai_class <- factor(agg$ai_class, levels = ai_labels)
agg <- agg[order(agg$product, agg$ai_class), ]
agg$pct_of_SMrz_pos <- ave(agg$area_SMrz_pos_Mkm2, agg$product,
                           FUN = function(x) round(100 * x / sum(x), 1))
agg$pct_of_union_pos <- ave(agg$area_union_pos_Mkm2, agg$product,
                            FUN = function(x) round(100 * x / sum(x), 1))
agg$area_SMrz_pos_Mkm2  <- round(agg$area_SMrz_pos_Mkm2, 2)
agg$area_union_pos_Mkm2 <- round(agg$area_union_pos_Mkm2, 2)

cat("\n=== Water-controlled vegetated area (significant POSITIVE supply association)\n")
cat("=== mean over January, April, July, October\n\n")
print(agg, row.names = FALSE)
cat("\nWritten: ", out_csv, "\n")
