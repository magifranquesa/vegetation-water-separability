# outputs/varpart_global/

Download these files from the Zenodo archive
(https://doi.org/10.5281/zenodo.20252053) and place them in this directory:

- `varpart_signif_global_2blocks.nc`        (GLEAM)
- `varpart_signif_global_2blocks_era5.nc`   (ERA5-Land)

Global two-block variance partitioning (supply = SMs + SMrz, demand = Ep) with a
per-cell, per-month permutation test, one file per hydroclimatic product. Each is
the result of a multi-hour permutation run; because the p-values are
permutation-based, archiving the files pins the exact values behind the figures
and tables rather than leaving them to a fresh random draw.

Read by the FDR step (`03_variance_partitioning/02_apply_fdr.R` and its ERA5
counterpart), the aridity/area analyses in `04_analysis/`, and Figs 2, 3 and S3.
The `fdr_adjusted_pvalues*.nc` files that also live here are regenerated cheaply
from these by the FDR step and are not archived.
