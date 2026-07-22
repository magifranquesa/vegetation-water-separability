#!/usr/bin/env bash

# ==============================================================================
# Script: 01_prepare_data/01_evaporation_deficit.sh
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   Compute the evaporation deficit (ED) as the difference between total
#   actual evaporation (E) and potential evaporation (Ep) for the full
#   1981-2022 monthly time series:
#
#     ED = E - Ep
#
#   This script must be run after 00_concatenate_gleam.sh.
#   The resulting ED file is then passed through the accumulation script
#   (02_accumulate_gleam_fluxes.sh) along with E, Ep, and Et.
#
# Input:
#   outputs/intermediate/concatenated/E_GLEAM_v4.2a_MO_1981-2022_NCO.nc
#   outputs/intermediate/concatenated/Ep_GLEAM_v4.2a_MO_1981-2022_NCO.nc
#
# Output:
#   outputs/intermediate/concatenated/ED_GLEAM_v4.2a_MO_1981-2022_NCO.nc
#
# Notes:
#   - ED is negative where actual evaporation is limited by water supply
#     (i.e. under water-stressed conditions), which is the ecologically
#     relevant signal in drylands.
#   - Positive ED values (E > Ep) can occur due to advection or model
#     artefacts; these are not masked here. See commented ncap2 command
#     at the end of this script if masking is required.
#   - The subtraction is performed on the full 1981-2022 series so that
#     the subsequent running-sum accumulation (02_accumulate_gleam_fluxes.sh)
#     operates on the same temporal coverage as the other variables.
#
# References:
#   Miralles et al. (2014). El Niño-La Niña cycle and recent trends in
#   continental evaporation. Nature Climate Change, 4, 122-126.
#
# Requirements:
#   NCO >= 5.0   (https://nco.sourceforge.net)
#
# Usage:
#   Run from the repository root:
#   bash 01_prepare_data/01_evaporation_deficit.sh
# ==============================================================================

set -euo pipefail

CONCAT_DIR="outputs/intermediate/concatenated"

e_file="${CONCAT_DIR}/E_GLEAM_v4.2a_MO_1981-2022_NCO.nc"
ep_file="${CONCAT_DIR}/Ep_GLEAM_v4.2a_MO_1981-2022_NCO.nc"
out_file="${CONCAT_DIR}/ED_GLEAM_v4.2a_MO_1981-2022_NCO.nc"

for f in "$e_file" "$ep_file"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: Input file not found: $f" >&2
    exit 1
  fi
done

echo "Computing evaporation deficit (ED = E - Ep)..."

# Step 1: Copy Ep and rename variable Ep -> E (so ncbo can subtract same-named vars).
ep_as_e="${CONCAT_DIR}/Ep_as_E.nc"
ncks -O "$ep_file" "$ep_as_e"
ncrename -O -v Ep,E "$ep_as_e"

# Step 2: Subtract. The output variable will be named E (inherited from E file).
ed_tmp="${CONCAT_DIR}/ED_tmp.nc"
ncbo -O --op_typ='-' "$e_file" "$ep_as_e" "$ed_tmp"

# Step 3: Rename variable E -> ED and add long_name attribute.
ncrename -O -v E,ED "$ed_tmp"
ncatted -O \
  -a long_name,ED,o,c,"Evaporative deficit (E - Ep)" \
  "$ed_tmp"

mv "$ed_tmp" "$out_file"

# Cleanup.
rm -f "$ep_as_e"

echo "Output: $out_file"

# ------------------------------------------------------------------------------
# Optional: mask out positive values (E > Ep).
# Uncomment if the analysis requires ED <= 0 only.
# ------------------------------------------------------------------------------
# ncap2 -O -s 'where(ED > 0) ED = 0.0f' "$out_file" "${out_file%.nc}_pos0.nc"
