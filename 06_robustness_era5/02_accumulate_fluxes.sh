#!/usr/bin/env bash

# ==============================================================================
# Script: 06_robustness_era5/02_accumulate_fluxes.sh
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   ERA5-Land counterpart of 02_accumulate_gleam_fluxes.sh, for the robustness
#   replication. Computes running SUMS of the monthly flux indicators (Ep, Et, ED)
#   over accumulation windows of 1, 3, 6, 9 and 12 months.
#
#   The running sum is implemented by inverting the time axis before applying
#   cdo runsum and re-inverting afterwards, so that each time step accumulates
#   over the PRECEDING N months (look-back window, not look-ahead).
#
#   Scale 1 requires no accumulation: the monthly values are used directly and
#   the file is simply trimmed to 1982-2022.
#
#   The 1981 lead-in year exists precisely so that the 12-month window is already
#   valid in January 1982. Trimming the last 492 steps therefore yields exactly the
#   1982-2022 period, matching kNDVI.
#
# Input:
#   outputs/intermediate/concatenated_era5/{var}_era5land_monthly_1981-2022.nc
#
#   One file per variable (Ep, Et, ED), 504 monthly steps from January 1981, on the
#   GLEAM/kNDVI grid (3600 x 1800, cell-centred, -180/180). Produced by
#   06_robustness_era5/01_prepare_indicators.R.
#
#   Ep is the FAO-56 Penman-Monteith ET0 built in et0_era5land/, NOT ERA5-Land's own
#   'pev' (which is biased high and has been discarded).
#   ED = E - Ep, the same definition as the GLEAM pipeline.
#
# Output:
#   outputs/intermediate/accumulated_era5/{var}_era5land_monthly_1982-2022_scale_{s}.nc
#   (s = 1, 3, 6, 9, 12; 492 time steps each)
#
# Notes:
#   - Ep : atmospheric evaporative demand, FAO-56 ET0 (mm/month)
#   - Et : plant transpiration (mm/month). NOTE: in ERA5-Land this comes from the
#          variable labelled 'evabs', because of the documented variable rotation.
#   - ED : evaporation deficit, E - Ep (mm/month)
#   - E is NOT accumulated: it is not one of the five indicators, it only serves to
#     build ED (same convention as the GLEAM script).
#   - Soil moisture (SMrz, SMs) uses running MEANS, not sums, and is handled in
#     03_accumulate_era5_soil_moisture.sh
#
# Requirements:
#   CDO >= 2.0   (https://code.mpimet.mpg.de/projects/cdo)
#   NCO >= 5.0   (https://nco.sourceforge.net)
#
# Usage:
#   Run from the repository root:
#   bash 06_robustness_era5/02_accumulate_fluxes.sh
# ==============================================================================

set -euo pipefail

CONCAT_DIR="outputs/intermediate/concatenated_era5"
OUT_DIR="outputs/intermediate/accumulated_era5"

mkdir -p "$OUT_DIR"

vars=("Ep" "Et" "ED")
scales=(3 6 9 12)

echo "Processing ERA5-Land flux variables: ${vars[*]}"
echo "  Input:  $CONCAT_DIR/"
echo "  Output: $OUT_DIR/"
echo

for var in "${vars[@]}"; do

  infile="${CONCAT_DIR}/${var}_era5land_monthly_1981-2022.nc"

  if [[ ! -f "$infile" ]]; then
    echo "ERROR: Input file not found: $infile" >&2
    exit 1
  fi

  # The look-back window is only valid from 1982 onwards if the series really starts
  # in 1981. A file with the wrong number of steps would silently shift everything.
  nsteps=$(cdo -s ntime "$infile")
  if [[ "$nsteps" -ne 504 ]]; then
    echo "ERROR: $infile has $nsteps time steps, expected 504 (1981-01..2022-12)." >&2
    exit 1
  fi

  echo "------------------------------------------------------------"
  echo "Variable: $var"
  echo "Input:    $infile  ($nsteps steps)"
  echo "------------------------------------------------------------"

  # ------------------------------------------------------------------
  # Scale 1: no accumulation — trim to 1982-2022 only
  # ------------------------------------------------------------------
  outfile_1="${OUT_DIR}/${var}_era5land_monthly_1982-2022_scale_1.nc"
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
    outfile="${OUT_DIR}/${var}_era5land_monthly_1982-2022_scale_${s}.nc"

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

echo "All ERA5-Land flux variables processed."
