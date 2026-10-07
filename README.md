# Vegetation–water separability

Code to reproduce the analysis and figures of:

> **Aridity and timescale bound the separability of vegetation water stress.**
> Franquesa et al. (2026). *[Journal — DOI to be added on acceptance].*

The study asks whether soil-water supply and atmospheric demand can be
distinguished as controls on interannual vegetation activity, and shows that this
separability is bounded by aridity and by the timescale over which hydroclimatic
conditions are integrated. This repository contains the full analysis pipeline —
from raw data to every figure and table in the paper, including the ERA5-Land
robustness assessment reported in the supplementary materials.

---

## Repository structure

```
R/                        Shared configuration and helpers
  config.R                Relative paths for all inputs and outputs
  config_era5.R           ERA5-Land raw-input definitions
  helpers_ncdf.R          NetCDF writing helper

01_prepare_data/          Raw data -> analysis-ready fields
                          (concatenate GLEAM, evaporation deficit, accumulate to
                          1, 3, 6, 9 and 12 months, kNDVI, aridity index,
                          vegetated mask)
02_correlations/          Monthly Spearman correlations (kNDVI vs indicators),
                          rho* and rho_max extraction, BH-FDR significance mask
03_variance_partitioning/ Two-block variance partitioning + FDR + category split
04_analysis/              Derived metrics: aridity stratification, supply-demand
                          coupling by timescale, spatial agreement, area summaries
05_figures_science/       Figs. 1-4 and figs. S1-S5 (one script per figure)
06_robustness_era5/       ERA5-Land replication (tables S1 to S4):
  et0/                    Atmospheric evaporative demand (FAO-56 Penman-Monteith)
  checks/                 Physical verification of ERA5-Land variables

data/
  processed/aridity/      Derived aridity index and Zomer classes (distributed)
  external/               World boundaries (Equal Earth) for maps
outputs/                  Created when the pipeline runs (git-ignored)
```

Each figure maps to one script in `05_figures_science/`; outputs are written to
`outputs/figures_science/`:

| Figure | Script | Figure | Script |
|---|---|---|---|
| Fig. 1 | `fig1_science.R` | fig. S1 | `figS1_science.R` |
| Fig. 2 | `fig2_science.R` | fig. S2 | `figS2_science.R` |
| Fig. 3 | `fig3_science.R` | fig. S3 | `figS3_science.R` |
| Fig. 4 | `fig4_science.R` | fig. S4 | `figS4_science.R` |
| | | fig. S5 | `figS5_science.R` |

The figure scripts of an earlier submission, with a different figure numbering,
are preserved under the git tag `nature-version`.

---

## Data sources

Raw inputs are publicly available from their providers and are **not** included
here (they are large and better obtained from source):

| Dataset | Version | Source |
|---|---|---|
| GLEAM (soil moisture, evaporation, potential evaporation) | v4.2a | Zenodo, https://doi.org/10.5281/zenodo.14724263 |
| GIMMS NDVI3g+ | — | ORNL DAAC, https://doi.org/10.3334/ORNLDAAC/2187 |
| CRU TS (precipitation, PET) | 4.09 | https://crudata.uea.ac.uk/cru/data/hrg/ |
| ERA5-Land (robustness) | — | Muñoz-Sabater et al. (2021), *ESSD* 13, 4349-4383; Copernicus CDS |

GLEAM v4.2a has been superseded by v4.3a on the GLEAM server; the version used
here is permanently archived on Zenodo (link above). According to the GLEAM4.2
README, the monthly product there is distributed as one file per variable
covering the whole record; in that case `01_prepare_data/00_concatenate_gleam.sh`
(which joins per-year files) is not needed, and the files only have to be
restricted to 1981–2022.

Precomputed intermediate outputs are archived on Zenodo:
https://doi.org/10.5281/zenodo.20252053. Downloading them lets you skip the two
most expensive parts of the pipeline — the pixel-wise monthly Spearman
correlations and the multi-hour variance partitioning — as well as the entire
ERA5-Land branch, whose raw inputs are large. Place each file in the directory
shown below; these directories already exist in the repository, each with a short
README.

| Zenodo file(s) | Destination directory |
|---|---|
| `kndvi.nc` | `outputs/intermediate/kndvi/` |
| `abs_max_correlation_kndvi_{ED,Ep,Et,SMrz,SMs}.nc` | `outputs/intermediate/absmax_spearman/` |
| `max_correlation_kndvi_{ED,Ep,Et,SMrz,SMs}.nc` | `outputs/intermediate/max_spearman/` |
| `abs_max_correlation_kndvi_{ED,Ep,Et,SMrz,SMs}_era5.nc` | `outputs/intermediate/absmax_spearman_era5/` |
| `varpart_signif_global_2blocks.nc` | `outputs/varpart_global/` |
| `varpart_signif_global_2blocks_era5.nc` | `outputs/varpart_global/` |

The ERA5-Land rho* files carry an `_era5` suffix so they do not collide with the
identically named GLEAM files. With every file in place, the remaining steps are
fast: run the FDR masks (`02_correlations/04_fdr_significance_mask.R` and
`06_robustness_era5/06_fdr_significance_mask.R`), the FDR adjustment and category
split (`03_variance_partitioning/02_apply_fdr.R`, `03_category_split.R`, and the
ERA5 counterpart `06_robustness_era5/08_apply_fdr.R`), and then the analysis and
figure scripts. None of the raw GLEAM or ERA5-Land data, and neither of the two
heavy steps, are needed on that path.

Small derived files are shipped directly in the repository: the aridity index and
Zomer classes and the map boundary shapefile under `data/`; the vegetated mask
(`outputs/intermediate/vegetation_mask_c1.tif`); and the supply–demand coupling
table (`outputs/tables/supply_demand_coupling.csv`), which drives Fig. 4C and
would otherwise need the large accumulated GLEAM fields to rebuild. The
timescale-bar tables for Fig. 4A,B are not shipped: `04_analysis/09_timescale_bars.R`
regenerates them from the archived rho* files.

---

## Requirements

- **R** >= 4.2, with the packages listed in `install.R`
  (`terra`, `sf`, `ncdf4`, `ggplot2`, `patchwork`, `scico`, `ggridges`, `hexbin`,
  `vegan`, `pracma`, `dplyr`, `tidyr`, `readr`, `scales`, `abind`).
- **CDO** and **NCO** for the data-preparation shell scripts in `01_prepare_data/`.
- For the ERA5-Land download: Python with the `cdsapi` client.

Install the R packages:

```bash
Rscript install.R
```

For a version-pinned environment (recommended for exact reproducibility):

```r
install.packages("renv")
renv::init()       # discovers dependencies from the scripts
renv::snapshot()   # writes renv.lock with the exact installed versions
```

---

## Running the pipeline

Run everything **from the repository root** (scripts resolve paths relative to
it via `R/config.R`). Stages are ordered; run them in sequence.

```bash
# 1. Prepare data (GLEAM + GIMMS + CRU -> analysis-ready fields)
bash 01_prepare_data/00_concatenate_gleam.sh
#    ... through 08_vegetation_mask.R (see folder, run in numbered order)

# 2. Correlations and significance
Rscript 02_correlations/01_monthly_spearman.R
Rscript 02_correlations/02_extract_absmax.R
Rscript 02_correlations/03_extract_max.R
Rscript 02_correlations/04_fdr_significance_mask.R

# 3. Variance partitioning
Rscript 03_variance_partitioning/01_varpart.R      # heavy: permutation tests
Rscript 03_variance_partitioning/02_apply_fdr.R
Rscript 03_variance_partitioning/03_category_split.R

# 4. Derived metrics
Rscript 04_analysis/01_rhostar_areas.R             # ... run 01-08 in order

# 5. Figures
Rscript 05_figures_science/fig1_science.R          # ... fig1-fig4, figS1-figS5
```

The ERA5-Land robustness assessment (tables S1 to S4) is a parallel
track under `06_robustness_era5/`, run in numbered order after the ERA5-Land raw
data have been obtained.

### Shortcut: from the archived intermediates to the figures

With the Zenodo files placed as listed above (and the small files already shipped
in the repository), stages 01–03 and the entire ERA5-Land branch can be skipped.
Only these fast steps are needed to reproduce every main and supplementary figure:

```bash
# 1. Regenerate the fast FDR intermediates (seconds each)
Rscript 02_correlations/04_fdr_significance_mask.R      # GLEAM -> absmax_fdr/
Rscript 03_variance_partitioning/02_apply_fdr.R         # GLEAM -> fdr_adjusted_pvalues.nc
Rscript 03_variance_partitioning/03_category_split.R
Rscript 06_robustness_era5/06_fdr_significance_mask.R   # ERA5  -> absmax_fdr_era5/
Rscript 06_robustness_era5/08_apply_fdr.R               # ERA5  -> fdr_adjusted_pvalues_era5.nc

# 2. Table that feeds Fig. 4A,B (reads the archived rho*)
Rscript 04_analysis/09_timescale_bars.R

# 3. Figures, any order
Rscript 05_figures_science/fig1_science.R    # ... fig2-fig4, figS1-figS5 (*_science.R)
```

figs. S2 and S4 read small precomputed tables shipped in the
repository (`outputs/tables/proportions_vegetated_fdr.csv` and
`spatial_correlations_rho_max.csv`), so both render directly at step 3; to rebuild
those tables from the archived correlations run
`04_analysis/08_vegetated_area_proportions_fdr.R` and
`04_analysis/07_spatial_agreement_rho_max.R` respectively.

Tables S1 to S4 reproduce from the same archived data after step 1:
`table1_2_spatial_strength.R`, `table3_varpart_categories.R` and
`table_water_controlled_area.R` in `06_robustness_era5/`. The two exceptions are
`04_analysis/06_supply_demand_coupling.R` (GLEAM) and `06_robustness_era5/table4_coupling.R`
(ERA5, table S4): both need the raw accumulated fields, which are too
large to archive. The GLEAM coupling table they feed
(`outputs/tables/supply_demand_coupling.csv`, Fig. 4C) is shipped in the
repository; the rest of that path requires rebuilding the accumulated fields from
the raw data via stages 01 and 06.

> **Note on compute.** The variance partitioning (`03_variance_partitioning/01_varpart.R`
> and its ERA5 counterpart) runs a permutation test per grid cell and month over
> the global vegetated domain and is designed for a multi-core machine.

Outputs are written under `outputs/` (figures in `outputs/figures_science/`,
tables in `outputs/tables/`), which is git-ignored.

---

## Citation

If you use this code, please cite the paper (above) and this repository.

## License

MIT — see [LICENSE](LICENSE).
