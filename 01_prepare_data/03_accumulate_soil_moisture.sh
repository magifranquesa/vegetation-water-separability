#!/usr/bin/env bash

# ==============================================================================
# Script: 01_prepare_data/03_accumulate_soil_moisture.sh
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   Compute running means of GLEAM v4.2a monthly soil moisture variables
#   (SMrz, SMs) over accumulation timescales of 1, 3, 6, 9, and 12 months.
#
#   Running means (not sums) are used for soil moisture because these
#   variables represent storage states rather than fluxes.
#
#   The running mean is implemented by inverting the time axis before applying
#   cdo runmean and re-inverting afterward, so that each time step represents
#   the mean over the preceding N months.
#
#   Scale 1 requires no accumulation: the original monthly values are used
#   directly, and the file is simply trimmed to the 1982-2022 period.
#
# Input:
#   outputs/intermediate/concatenated/{var}_GLEAM_v4.2a_MO_1981-2022_NCO.nc
#
#   One file per variable (SMrz, SMs), starting January 1981 so that
#   accumulation timescales are valid from January 1982 onward.
#
# Output:
#   outputs/intermediate/accumulated/{var}_GLEAM_v4.2a_MO_1982-2022_scale_{s}.nc
#   (s = 1, 3, 6, 9, 12)
#
#   Each output file covers January 1982 – December 2022 (492 time steps).
#   The variable name inside each file matches the indicator name (e.g. "SMrz").
#
# Notes:
#   - SMrz : root-zone soil moisture (m³/m³)
#   - SMs  : surface soil moisture (m³/m³)
#   - Flux variables (E, Ep, Et, ED) use running sums and are processed in:
#       02_accumulate_gleam_fluxes.sh
#   - Evaporation deficit (ED) is derived from E and Ep in:
#       01_compute_evaporation_deficit.sh
#
# Requirements:
#   CDO >= 2.0   (https://code.mpimet.mpg.de/projects/cdo)
#   NCO >= 5.0   (https://nco.sourceforge.net)
#
# Usage:
#   Run from the repository root:
#   bash 01_prepare_data/03_accumulate_soil_moisture.sh
# ==============================================================================

set -euo pipefail

CONCAT_DIR="outputs/intermediate/concatenated"
OUT_DIR="outputs/intermediate/accumulated"

mkdir -p "$OUT_DIR"

vars=("SMrz" "SMs")
scales=(3 6 9 12)

echo "Processing GLEAM soil moisture variables: ${vars[*]}"
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
  # Scales 3, 6, 9, 12: running mean with time-axis inversion
  # ------------------------------------------------------------------

  inv1="${OUT_DIR}/${var}_inv1.nc"
  echo "  Inverting time axis → $inv1"
  ncpdq -O -a -time "$infile" "$inv1"

  for s in "${scales[@]}"; do
    echo "  Scale $s → running mean..."

    inv2="${OUT_DIR}/${var}_inv2_scale_${s}.nc"
    tmpout="${OUT_DIR}/${var}_tmp_scale_${s}.nc"
    outfile="${OUT_DIR}/${var}_GLEAM_v4.2a_MO_1982-2022_scale_${s}.nc"

    cdo -O runmean,"$s" "$inv1" "$inv2"

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

echo "All soil moisture variables processed."
