#!/usr/bin/env bash

# ==============================================================================
# Script: 01_prepare_data/00_concatenate_gleam.sh
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   Concatenate annual GLEAM v4.2a NetCDF files into a single time series
#   covering 1981-2022 for each hydroclimatic variable. The 1981 starting
#   year is required so that accumulation timescales (up to 12 months) are
#   valid from January 1982 onward.
#
#   The concatenation uses NCO (ncks + ncrcat). Each annual file must first
#   have its time dimension converted to a record dimension (ncks --mk_rec_dmn)
#   before being concatenated with ncrcat.
#
# Input:
#   Annual NC files downloaded from the GLEAM v4.2a FTP server:
#     https://www.gleam.eu
#
#   Expected filename pattern per variable:
#     {var}_????_GLEAM_v4.2a_MO.nc   (one file per year)
#
#   Variables:
#     E    : total actual evaporation
#     Ep   : potential evaporation (atmospheric evaporative demand)
#     Et   : plant transpiration
#     SMrz : root-zone soil moisture
#     SMs  : surface soil moisture
#
# Output:
#   {var}_GLEAM_v4.2a_MO_1981-2022_NCO.nc   (one file per variable)
#
# Notes:
#   - Raw annual files are organized in per-variable subdirectories:
#       data/raw/GLEAM/{var}/{var}_????_GLEAM_v4.2a_MO.nc
#   - Concatenated output files are written to OUT_DIR (default:
#       outputs/intermediate/concatenated/).
#   - Temporary files (rec_*.nc) are created in OUT_DIR and removed after
#     concatenation.
#   - ED (evaporation deficit) is not a GLEAM output variable; it is computed
#     from E and Ep in the next step: 01_compute_evaporation_deficit.sh
#   - When running on a remote server, set RAW_DIR and OUT_DIR to the
#     appropriate server paths before executing.
#
# Requirements:
#   NCO >= 5.0   (https://nco.sourceforge.net)
#
# Usage:
#   Run from the repository root:
#   bash 01_prepare_data/00_concatenate_gleam.sh
# ==============================================================================

set -euo pipefail

# Paths — adjust if running on a remote server.
RAW_DIR="data/raw/GLEAM"
OUT_DIR="outputs/intermediate/concatenated"

mkdir -p "$OUT_DIR"

vars=("E" "Ep" "Et" "SMrz" "SMs")

echo "Concatenating annual GLEAM v4.2a files..."
echo "  Input:  $RAW_DIR/{var}/"
echo "  Output: $OUT_DIR/"
echo

for var in "${vars[@]}"; do

  in_dir="${RAW_DIR}/${var}"
  outfile="${OUT_DIR}/${var}_GLEAM_v4.2a_MO_1981-2022_NCO.nc"

  echo "------------------------------------------------------------"
  echo "Variable: $var → $outfile"
  echo "------------------------------------------------------------"

  if ! ls "${in_dir}/${var}_"????_GLEAM_v4.2a_MO.nc 1>/dev/null 2>&1; then
    echo "ERROR: No annual files found in: $in_dir" >&2
    exit 1
  fi

  # Convert time to record dimension in each annual file.
  echo "  Converting time to record dimension..."
  for f in "${in_dir}/${var}_"????_GLEAM_v4.2a_MO.nc; do
    ncks -O --mk_rec_dmn time "$f" "${OUT_DIR}/rec_$(basename "$f")"
  done

  # Concatenate all annual files in chronological order.
  echo "  Concatenating → $outfile"
  ncrcat -O $(ls -1 "${OUT_DIR}/rec_${var}_"????_GLEAM_v4.2a_MO.nc | sort -V) "$outfile"

  # Remove temporary record-dimension files.
  rm -f "${OUT_DIR}/rec_${var}_"????_GLEAM_v4.2a_MO.nc

  echo "  Done: $outfile"
  echo
done

echo "All variables concatenated."
