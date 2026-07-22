#!/usr/bin/env bash

# ==============================================================================
# Script: 01_prepare_data/05_compute_kndvi.sh
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   1. Compute the kernel Normalized Difference Vegetation Index (kNDVI) from
#      monthly maximum NDVI (NDVI3g) using the formula:
#
#        kNDVI = tanh(NDVI²) = (exp(2·NDVI²) - 1) / (exp(2·NDVI²) + 1)
#
#   2. Resample kNDVI from its native resolution to the 0.1° GLEAM v4.2a grid
#      using conservative (area-weighted) remapping (CDO remapcon), which
#      preserves grid-cell mean values during spatial aggregation and avoids
#      introducing artificial spatial detail compared with interpolation-based
#      methods.
#
# Input:
#   NDVI3g_monthlymax_1982_2022.nc   — monthly maximum NDVI3g, native grid
#   gleam01_grid.txt                 — GLEAM v4.2a 0.1° target grid definition
#
# Output:
#   kndvi.nc                         — kNDVI resampled to GLEAM 0.1° grid
#
#   The file kndvi.nc is the direct input to:
#     02_correlations/01_monthly_spearman.R
#
# Notes:
#   - kNDVI was proposed by Camps-Valls et al. (2021) as a theoretically
#     grounded, kernel-based vegetation index.
#   - NDVI3g is the third-generation Global Inventory Modelling and Mapping
#     Studies (GIMMS) NDVI dataset.
#   - Conservative remapping (remapcon) is preferred over bilinear (remapbil)
#     because it preserves area-integrated quantities when aggregating from a
#     coarser to a finer grid or between non-commensurate grids.
#   - The grid definition file gleam01_grid.txt is located in:
#       01_prepare_data/gleam_grid.txt
#
# References:
#   Camps-Valls, G., et al. (2021). A unified vegetation index for quantifying
#   the terrestrial biosphere. Science Advances, 7, eabc7447.
#
# Requirements:
#   CDO >= 2.0   (https://code.mpimet.mpg.de/projects/cdo)
#
# Usage:
#   Run from the directory containing the input NDVI file:
#   bash 01_prepare_data/05_compute_kndvi.sh
# ==============================================================================

set -euo pipefail

ndvi_file="data/raw/NDVI/NDVI3g_monthlymax_1982_2022.nc"
grid_file="01_prepare_data/gleam_grid.txt"
kndvi_native="kNDVI3g_monthlymax_1982_2022.nc"

OUT_DIR="outputs/intermediate/kndvi"
mkdir -p "$OUT_DIR"
kndvi_out="${OUT_DIR}/kndvi.nc"

for f in "$ndvi_file" "$grid_file"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: Required file not found: $f" >&2
    exit 1
  fi
done

# ------------------------------------------------------------------------------
# Step 1: Compute kNDVI from NDVI
# ------------------------------------------------------------------------------
echo "Step 1: Computing kNDVI = tanh(NDVI²)..."
cdo -O expr,'kndvi=(exp(2*ndvi*ndvi)-1)/(exp(2*ndvi*ndvi)+1)' \
  "$ndvi_file" \
  "$kndvi_native"

# ------------------------------------------------------------------------------
# Step 2: Resample to GLEAM 0.1° grid (conservative remapping)
# ------------------------------------------------------------------------------
echo "Step 2: Resampling to GLEAM 0.1° grid (remapcon)..."
cdo -O remapcon,"$grid_file" "$kndvi_native" "$kndvi_out"

rm -f "$kndvi_native"

echo "Done."
echo "  Output: $kndvi_out"
