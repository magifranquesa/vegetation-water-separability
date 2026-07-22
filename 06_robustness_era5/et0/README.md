# Atmospheric evaporative demand (ET0) from monthly ERA5-Land

Self-contained module that computes **FAO-56 Penman-Monteith reference
evapotranspiration** (ET0, used here as the atmospheric evaporative demand, AED)
from **ERA5-Land** monthly means (`reanalysis-era5-land-monthly-means`, 0.1°,
1981–2022). The FAO-56 equations are applied at **monthly** step on the native
ERA5-Land grid with `ncdf4`.

## Contents

| File | What it does |
|---|---|
| `fao56_functions.R` | Vectorised FAO-56 equations (Allen et al. 1998) + NetCDF helpers |
| `01_prepare_elevation.R` | ERA5-Land invariant geopotential → `z` (m) → `elevation_era5land.rds` |
| `02_compute_et0.R` | Monthly ET0 over the full grid → NetCDF, **native 0–360** |
| `03_regrid_to_gleam.R` | Rotates + `resample`s the native ET0 onto the **GLEAM/kNDVI grid** (−180/180) |

## Inputs

The input file paths and the variable name inside each file are defined in
`R/config_era5.R` (`era5` and `era5_var`); source that file instead of
hard-coding paths. The required ERA5-Land monthly variables are:

- 2 m dewpoint temperature — `d2m` [K]
- surface solar radiation downwards — `ssrd` [J m⁻²] (mean daily accumulation)
- 10 m wind components — `u10`, `v10` [m s⁻¹]
- 2 m temperature, monthly means of daily minimum / maximum — `t2m` [K]
- time-invariant geopotential — `z` [m² s⁻²]

All inputs share the ERA5-Land 0.1° grid in 0–360 longitude, so no grid
harmonisation is needed; `02_compute_et0.R` includes a guard that stops if any
file arrives on a different grid. Variable names are detected automatically
(`data_vars()`). Both outputs (native 0–360 and GLEAM grid) are kept; no script
deletes or overwrites files.

## Usage

Run **from the repository root**:

```bash
# 1) Once: build the elevation grid
Rscript 06_robustness_era5/et0/01_prepare_elevation.R

# 2) Compute native ET0
Rscript 06_robustness_era5/et0/02_compute_et0.R

# 3) Regrid onto the GLEAM/kNDVI grid used by the rest of the pipeline
Rscript 06_robustness_era5/et0/03_regrid_to_gleam.R
```

Outputs in `data/processed/et0_era5land/` (variable `et0`, **mm month⁻¹**,
compressed; ocean = `NA`):
- `et0_era5land_monthly_1981-2022.nc` — native ERA5-Land grid (0–360).
- `et0_era5land_monthly_1981-2022_gleamgrid.nc` — GLEAM/kNDVI grid (−180/180).

## Methods choices and assumptions

- **Method:** FAO-56 Penman-Monteith, reference crop (grass, α = 0.23),
  equations 6, 7, 8, 11–14, 21, 37–39, 47 of Allen et al. (1998).
- **Elevation:** ERA5-Land invariant geopotential (consistent with the model
  orography and without resampling), not an external DEM; converted to metres as
  `z = geopotential / g` (g = 9.80665 m s⁻²).
- **Wind:** magnitude of `u10`/`v10` brought to 2 m with the FAO-56 factor
  (eq. 47, 0.748).
- **Radiation:** `Rs = ssrd / 1e6` (MJ m⁻² day⁻¹), taking `ssrd` in the monthly
  means as the mean daily accumulation in J m⁻². Confirmed empirically: the July
  maximum reaches ≈31.5 MJ m⁻² day⁻¹, the clear-sky ceiling (0.75 × Ra), not ≈600.
- **Soil heat flux** `G = 0` (standard for monthly ET0; minor effect).
- **Monthly ET0:** daily ET0 × number of days in the month. Computing from
  monthly means slightly underestimates the sum of daily ET0 (non-linearity); a
  known and accepted limitation at monthly scale.
