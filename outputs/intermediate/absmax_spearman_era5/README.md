# outputs/intermediate/absmax_spearman_era5/

Download these files from the Zenodo archive
(https://doi.org/10.5281/zenodo.20252053) and place them in this directory:

- `abs_max_correlation_kndvi_ED_era5.nc`
- `abs_max_correlation_kndvi_Ep_era5.nc`
- `abs_max_correlation_kndvi_Et_era5.nc`
- `abs_max_correlation_kndvi_SMrz_era5.nc`
- `abs_max_correlation_kndvi_SMs_era5.nc`

ERA5-Land counterpart of the signed Spearman correlation of greatest absolute
magnitude across accumulation timescales (rho*), computed against ERA5-Land
indicators with the same method as the GLEAM rho*. Regenerating them requires the
full ERA5-Land branch (large raw downloads + reference-ET computation + the
correlation step), so they are archived instead.

The `_era5` suffix distinguishes them from the GLEAM rho* files, which carry the
same base name but live in `outputs/intermediate/absmax_spearman/`. Read by the
ERA5 FDR mask and the robustness tables (Supplementary Tables 1-4).
