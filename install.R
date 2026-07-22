# Install the CRAN packages required to run this repository (latest versions).
#
# For an exact, version-pinned environment, use renv instead: the renv.lock in
# the repository recreates the exact package versions with
#   renv::restore()
# and you do not need this script. See the README for details.
#
# To install with this script instead, run once:
#   Rscript install.R

pkgs <- c(
  "terra",      # raster / NetCDF handling
  "sf",         # vector data, map projections
  "ncdf4",      # low-level NetCDF I/O
  "ggplot2",    # figures
  "patchwork",  # figure composition
  "scico",      # perceptual colour scales
  "ggridges",   # ridgeline plots
  "hexbin",     # hexbin density (Fig. 4)
  "vegan",      # variance partitioning (RDA, permutation tests)
  "pracma",     # detrending
  "dplyr", "tidyr", "readr", "scales", "abind"  # data wrangling
)

installed <- rownames(installed.packages())
to_install <- setdiff(pkgs, installed)

if (length(to_install)) {
  message("Installing: ", paste(to_install, collapse = ", "))
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  message("All required packages are already installed.")
}

# Base packages used but shipped with R (no installation needed):
#   grid, parallel, stats
