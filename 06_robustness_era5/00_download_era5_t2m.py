#!/usr/bin/env python3
"""
Download hourly ERA5-Land 2m air temperature, global, 1981-2022.

Purpose
-------
Provide the input to compute monthly Tmin/Tmax for the atmospheric
evaporative demand (AED) / potential evapotranspiration.

Why hourly and not the monthly-means product?
    The `reanalysis-era5-land-monthly-means` product (used for the other
    variables: dewpoint, ssrd, wind, ...) only stores the monthly MEAN of
    each field. Daily minima/maxima cannot be recovered from a mean, so the
    hourly `reanalysis-era5-land` dataset is required. From these 24 hourly
    values per day one derives the daily min/max, and then the monthly mean
    of the daily minima (Tmin) and of the daily maxima (Tmax).

Output
------
One NetCDF file per year-month:
    data/raw/ERA5_land_monthly/era5-land_t2m_YYYY_MM.nc

Writing one file per month makes the download resumable: files that already
exist (and are non-empty) are skipped, so the script can be re-run after an
interruption. Downloads go to a `.part` temporary file that is renamed only
on success, so a half-written file is never mistaken for a complete one.

WARNING - data volume
    This is a very large download. Global hourly 0.1 deg 2m_temperature is
    several GB per month; the full 1981-2022 record is on the order of
    terabytes. Make sure the target disk has enough free space.

Requirements
------------
    pip install "cdsapi>=0.7"
A valid CDS API key. This script reads the project `.cdsapirc` (in the repo
root) explicitly, so it works even if cdsapi does not find one in $HOME.

Usage
-----
    python 06_robustness_era5/00_download_era5_t2m.py
"""

import calendar
import os
import sys
import time

import cdsapi

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
DATASET = "reanalysis-era5-land"

# Project root = one level up from this script (06_robustness_era5/..)
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Large download (order of terabytes): point this to a disk with enough space.
OUT_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "ERA5_land_monthly")
CDSAPIRC = os.path.join(PROJECT_ROOT, ".cdsapirc")

YEARS = range(1981, 2023)          # 1981 .. 2022 inclusive
MONTHS = range(1, 13)              # 01 .. 12
HOURS = [f"{h:02d}:00" for h in range(24)]   # all 24 hours (needed globally:
#   local time of Tmin/Tmax shifts with longitude, so no hour can be dropped)

MAX_RETRIES = 3
RETRY_WAIT = 60                    # seconds between retries


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def load_cdsapirc(path):
    """Parse a .cdsapirc file into (url, key). Returns (None, None) if absent."""
    if not os.path.isfile(path):
        return None, None
    conf = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or ":" not in line:
                continue
            k, v = line.split(":", 1)
            conf[k.strip()] = v.strip()
    return conf.get("url"), conf.get("key")


def build_request(year, month):
    """CDS request dict for one global year-month of hourly 2m temperature."""
    ndays = calendar.monthrange(year, month)[1]
    days = [f"{d:02d}" for d in range(1, ndays + 1)]
    return {
        "variable": ["2m_temperature"],
        "year": str(year),
        "month": f"{month:02d}",
        "day": days,
        "time": HOURS,
        "data_format": "netcdf",
        "download_format": "unarchived",
    }


def make_client():
    url, key = load_cdsapirc(CDSAPIRC)
    if url and key:
        print(f"[info] using credentials from {CDSAPIRC}")
        return cdsapi.Client(url=url, key=key)
    print("[info] .cdsapirc not found in project root; "
          "falling back to default cdsapi lookup ($HOME/.cdsapirc)")
    return cdsapi.Client()


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    client = make_client()

    tasks = [(y, m) for y in YEARS for m in MONTHS]
    total = len(tasks)
    failed = []

    for i, (year, month) in enumerate(tasks, start=1):
        fname = f"era5-land_t2m_{year}_{month:02d}.nc"
        fpath = os.path.join(OUT_DIR, fname)
        prefix = f"[{i:>3}/{total}] {fname}"

        if os.path.isfile(fpath) and os.path.getsize(fpath) > 0:
            print(f"{prefix} -> skip (already downloaded)")
            continue

        request = build_request(year, month)
        tmp = fpath + ".part"

        for attempt in range(1, MAX_RETRIES + 1):
            try:
                print(f"{prefix} -> retrieving (attempt {attempt}/{MAX_RETRIES})")
                client.retrieve(DATASET, request).download(tmp)
                os.replace(tmp, fpath)
                print(f"{prefix} -> done ({os.path.getsize(fpath) / 1e9:.2f} GB)")
                break
            except Exception as exc:                      # noqa: BLE001
                print(f"{prefix} -> ERROR: {exc}")
                if os.path.exists(tmp):
                    try:
                        os.remove(tmp)
                    except OSError:
                        pass
                if attempt < MAX_RETRIES:
                    time.sleep(RETRY_WAIT)
                else:
                    print(f"{prefix} -> giving up after {MAX_RETRIES} attempts")
                    failed.append(fname)

    print("\n==== summary ====")
    print(f"requested : {total} year-months")
    print(f"failed    : {len(failed)}")
    for f in failed:
        print(f"  - {f}")
    if failed:
        print("Re-run the script to retry the failed months "
              "(completed files are skipped).")
        sys.exit(1)


if __name__ == "__main__":
    main()
