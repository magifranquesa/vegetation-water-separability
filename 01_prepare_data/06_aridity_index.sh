#!/usr/bin/env bash

# ==============================================================================
# Script: 01_prepare_data/06_aridity_index.sh
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   Compute a single-layer mean Aridity Index (AI = ΣP / ΣPET) for the period
#   1982–2022 from CRU TS 4.09 monthly precipitation and potential
#   evapotranspiration grids. The result is a global 0.5° raster with one
#   value per grid cell representing the ratio of cumulative precipitation to
#   cumulative PET over the full period.
#
# Input:
#   data/raw/AI_CRU/cru_ts4.09.1901.2024.pre.dat.nc   (CRU TS precipitation)
#   data/raw/AI_CRU/cru_ts4.09.1901.2024.pet.dat.nc   (CRU TS PET, mm/day)
#
#   Source: https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.09/cruts.2503051245.v4.09/
#
# Output:
#   data/processed/aridity/ai_1982_2022_period.nc
#
# Intermediate files (written to data/raw/AI_CRU/):
#   pre_1982_2022_mmmon.nc   monthly precipitation 1982-2022 (mm/month)
#   pet_1982_2022_mmmon.nc   monthly PET 1982-2022 (mm/month)
#   pre_sum_1982_2022.nc     cumulative precipitation over the period
#   pet_sum_1982_2022.nc     cumulative PET over the period
#
# Notes:
#   - CRU PET is in mm/day; cdo muldpm converts it to mm/month before summing.
#   - The AI is computed as the ratio of period sums, not as the mean of
#     monthly ratios.
#   - Variable renaming and attribute editing use NCO (ncrename, ncatted).
#
# Requirements:
#   CDO  >= 2.0  (https://code.mpimet.mpg.de/projects/cdo)
#   NCO  >= 5.0  (https://nco.sourceforge.net)
#
# Usage:
#   bash 01_prepare_data/06_aridity_index.sh
# ==============================================================================

set -euo pipefail

RAW_DIR="data/raw/AI_CRU"
OUT_DIR="data/processed/aridity"

PRE_RAW="${RAW_DIR}/cru_ts4.09.1901.2024.pre.dat.nc"
PET_RAW="${RAW_DIR}/cru_ts4.09.1901.2024.pet.dat.nc"

PRE_MON="${RAW_DIR}/pre_1982_2022_mmmon.nc"
PET_MON="${RAW_DIR}/pet_1982_2022_mmmon.nc"
PRE_SUM="${RAW_DIR}/pre_sum_1982_2022.nc"
PET_SUM="${RAW_DIR}/pet_sum_1982_2022.nc"

AI_OUT="${OUT_DIR}/ai_1982_2022_period.nc"

mkdir -p "${OUT_DIR}"

# ------------------------------------------------------------------------------
# Step 1 — Subset precipitation to 1982-2022 (mm/month)
# ------------------------------------------------------------------------------
echo "Step 1: Subsetting precipitation..."
cdo -O seldate,1982-01-01,2022-12-31 -selname,pre "${PRE_RAW}" "${PRE_MON}"

# ------------------------------------------------------------------------------
# Step 2 — Subset PET to 1982-2022 and convert mm/day -> mm/month
# ------------------------------------------------------------------------------
echo "Step 2: Subsetting and converting PET..."
cdo -O seldate,1982-01-01,2022-12-31 -muldpm -selname,pet "${PET_RAW}" "${PET_MON}"

# ------------------------------------------------------------------------------
# Step 3 — Cumulative sums over the full period
# ------------------------------------------------------------------------------
echo "Step 3: Computing period sums..."
cdo -O timsum "${PRE_MON}" "${PRE_SUM}"
cdo -O timsum "${PET_MON}" "${PET_SUM}"

# ------------------------------------------------------------------------------
# Step 4 — AI = ΣP / ΣPET
# ------------------------------------------------------------------------------
echo "Step 4: Computing Aridity Index..."
cdo -O div "${PRE_SUM}" "${PET_SUM}" "${AI_OUT}"

# ------------------------------------------------------------------------------
# Step 5 — Rename variable and add metadata
# ------------------------------------------------------------------------------
echo "Step 5: Renaming variable and writing attributes..."
ncrename -O -v pre,ai "${AI_OUT}"

ncatted -O \
  -a long_name,ai,o,c,"Aridity index (P/PET) for 1982-2022 (ratio of sums)" \
  -a units,ai,o,c,"1" \
  -a title,global,o,c,"CRU TS4.09 Aridity Index (P/PET), 1982-2022" \
  "${AI_OUT}"

echo "Done. Output written to: ${AI_OUT}"
