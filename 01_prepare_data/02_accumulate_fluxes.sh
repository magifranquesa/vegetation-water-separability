#!/usr/bin/env bash

# ==============================================================================
# Script: 01_prepare_data/02_accumulate_fluxes.sh
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   Compute running sums of GLEAM v4.2a monthly flux variables (E, Ep, Et, ED)
#   over accumulation windows of 1, 3, 6, 9, and 12 months.
#
#   The running sum is implemented by inverting the time axis before applying
#   cdo runsum and re-inverting afterward, so that each time step accumulates
#   over the preceding N months (look-back window, not look-ahead).
#
#   Scale 1 requires no accumulation: the original monthly values are used
#   directly, and the file is simply trimmed to the 1982-2022 period.
#
# Input:
#   outputs/intermediate/concatenated/{var}_GLEAM_v4.2a_MO_1981-2022_NCO.nc
#
#   One file per variable (Ep, Et, ED), starting January 1981 so that
#   accumulation windows are valid from January 1982 onward.
#   ED must be generated first by 01_compute_evaporation_deficit.sh.
#
# Output:
#   outputs/intermediate/accumulated/{var}_GLEAM_v4.2a_MO_1982-2022_scale_{s}.nc
#   (s = 1, 3, 6, 9, 12)
#
#   Each output file covers January 1982 – December 2022 (492 time steps).
#   The variable name inside each file matches the indicator name (e.g. "Et").
#
# Notes:
#   - Ep : potential evaporation / atmospheric evaporative demand (mm/month)
#   - Et : plant transpiration (mm/month)
#   - ED : evaporation deficit (E - Ep, computed in 01_compute_evaporation_deficit.sh)
#   - Soil moisture variables (SMrz, SMs) use running means and are processed
#     in: 03_accumulate_gleam_soil_moisture.sh
#   - The grid definition file gleam01_grid.txt describes the 0.1° global
#     grid used by GLEAM v4.2a.
#
# Requirements:
#   CDO >= 2.0   (https://code.mpimet.mpg.de/projects/cdo)
#   NCO >= 5.0   (https://nco.sourceforge.net)
#
# Usage:
#   Run from the repository root:
#   bash 01_prepare_data/02_accumulate_fluxes.sh
# ==============================================================================

set -euo pipefail

CONCAT_DIR="outputs/intermediate/concatenated"
OUT_DIR="outputs/intermediate/accumulated"

mkdir -p "$OUT_DIR"

vars=("Ep" "Et" "ED")
scales=(3 6 9 12)

echo "Processing GLEAM flux variables: ${vars[*]}"
echo "  Input:  $CONCAT_DIR/"
echo "  Output: $OUT_DIR/"
echo

for var in "${vars[@]}"; do

  infile="${CONCAT_DIR}/${var}_GLEAM_v4.2a_MO_1981-2022_NCO.nc"

  if [[ ! -f "$infile" ]]; then
    echo "ERROR: Input file not found: $infile" >&2
    exit 1
  fi

  echo "------------------------------------------------------------"
  echo "Variable: $var"
  echo "Input:    $infile"
  echo "------------------------------------------------------------"

  # ------------------------------------------------------------------
  # Scale 1: no accumulation — trim to 1982-2022 only
  # ------------------------------------------------------------------
  outfile_1="${OUT_DIR}/${var}_GLEAM_v4.2a_MO_1982-2022_scale_1.nc"
  echo "  Scale 1 → $outfile_1"
  cdo -O seltimestep,-492/-1 "$infile" "$outfile_1"

  # ------------------------------------------------------------------
  # Scales 3, 6, 9, 12: running sum with time-axis inversion
  # ------------------------------------------------------------------

  inv1="${OUT_DIR}/${var}_inv1.nc"
  echo "  Inverting time axis → $inv1"
  ncpdq -O -a -time "$infile" "$inv1"

  for s in "${scales[@]}"; do
    echo "  Scale $s → running sum..."

    inv2="${OUT_DIR}/${var}_inv2_scale_${s}.nc"
    tmpout="${OUT_DIR}/${var}_tmp_scale_${s}.nc"
    outfile="${OUT_DIR}/${var}_GLEAM_v4.2a_MO_1982-2022_scale_${s}.nc"

    cdo -O runsum,"$s" "$inv1" "$inv2"

    echo "    Re-inverting time axis..."
    ncpdq -O -a -time "$inv2" "$tmpout"

    echo "    Trimming to 1982-2022 → $outfile"
    cdo -O seltimestep,-492/-1 "$tmpout" "$outfile"

    rm -f "$inv2" "$tmpout"
  done

  rm -f "$inv1"

  echo "  Done: $var"
  echo
done

echo "All flux variables processed."
