#' Write a 3D array to a NetCDF file (CF-1.5 compliant)
#'
#' @param data        3D array [lon, lat, time].
#' @param file.out    Output file path.
#' @param var         List with elements: name, name.long, unit, range.
#' @param lon         List with elements: name, name.long, unit, values.
#' @param lat         List with elements: name, name.long, unit, values.
#' @param time        List with element: values (numeric, days since origin).
#' @param climatology List with elements: cell.methods, time.bounds. NULL if unused.
#' @param crs         List with element: epsg (character). Only "4326" is fully supported.
#' @param overwrite   Logical. Overwrite existing file? Default FALSE.
#' @param comment     Character. Optional global summary attribute.
#' @param vers        Character. Optional version string.
#' @param meta        List. Optional non-standard variable attributes.
#'
#' @note Time units are hardcoded to "days since 1982-01-01 00:00:00".

write_ncdf <- function(data,
                       file.out,
                       var = list(name = NULL, name.long = NULL, unit = NULL, range = NULL),
                       lon = list(name = NULL, name.long = NULL, unit = NULL, values = NULL),
                       lat = list(name = NULL, name.long = NULL, unit = NULL, values = NULL),
                       time = list(values = NULL),
                       climatology = list(cell.methods = NULL, time.bounds = NULL),
                       crs = list(epsg = NULL),
                       overwrite = FALSE,
                       comment = NULL,
                       vers = NULL,
                       meta = NULL) {

  require(ncdf4)

  if (length(dim(data)) != 3) stop("Input data must have 3 dimensions.")

  dim_lon <- ncdim_def(lon$name, lon$unit, lon$values, longname = lon$name.long, unlim = FALSE)
  dim_lat <- ncdim_def(lat$name, lat$unit, lat$values, longname = lat$name.long, unlim = FALSE)
  dim_tme <- ncdim_def("time", "days since 1982-01-01 00:00:00", time$values,
                       longname = "Time", calendar = "standard", unlim = TRUE)

  out.nc.var <- ncvar_def(var$name, var$unit, list(dim_lon, dim_lat, dim_tme),
                          missval = NA, longname = var$name.long, compression = 9)
  out.nc.crs_def <- ncvar_def("crs", "", list(),
                              longname = "Coordinate Reference System definition", prec = "char")

  vars_list <- list(out.nc.crs_def, out.nc.var)

  if (!is.null(climatology)) {
    out.nc.bounds_def <- ncdim_def("bnds", "", c(1, 2),
                                   longname = "Time bounds", create_dimvar = TRUE)
    out.nc.time_bounds_def <- ncvar_def("time_bounds", "days since 1982-01-01 00:00:00",
                                        list(out.nc.bounds_def, out.nc.var$dim[[3]]),
                                        longname = "Time bounds of the temporal sub-intervals")
    vars_list <- list(out.nc.crs_def, out.nc.time_bounds_def, out.nc.var)
  }

  if (file.exists(file.out)) {
    if (!overwrite) stop("`file.out` already exists. Set `overwrite = TRUE` to replace it.")
    file.remove(file.out)
  }

  out.nc <- nc_create(file.out, vars_list)

  ncatt_put(out.nc, var$name, "grid_mapping", "crs")
  ncatt_put(out.nc, var$name, "valid_range", var$range)
  if (!is.null(climatology)) ncatt_put(out.nc, var$name, "cell_methods", climatology$cell.methods)

  ncatt_put(out.nc, "crs", "long_name", "CRS definition")
  ncatt_put(out.nc, "crs", "semi_major_axis", 6378137.0)
  ncatt_put(out.nc, "crs", "inverse_flattening", 298.257222101)
  if (!is.null(crs) && crs$epsg == "4326") {
    ncatt_put(out.nc, "crs", "grid_mapping_name", "latitude_longitude")
    ncatt_put(out.nc, "crs", "longitude_of_prime_meridian", "0.0")
  }

  ncatt_put(out.nc, "lon", "axis", "X")
  ncatt_put(out.nc, "lat", "axis", "Y")
  ncatt_put(out.nc, "time", "axis", "T")
  if (!is.null(climatology)) {
    ncatt_put(out.nc, "time", "bounds", "time_bounds")
  }

  ncatt_put(out.nc, 0, "conventions", "CF-1.5")
  if (!is.null(comment)) ncatt_put(out.nc, 0, "summary", comment)
  if (!is.null(vers))    ncatt_put(out.nc, 0, "version", vers)
  ncatt_put(out.nc, 0, "creation_date", date())

  ntme <- length(time$values)
  if (ntme < 1000) {
    ncvar_put(out.nc, var$name, data)
  } else {
    splits <- split(seq_len(ntme), ceiling(seq_along(seq_len(ntme)) / 1000))
    for (i in seq_along(splits)) ncvar_put(out.nc, var$name, data[splits[[i]]])
  }

  if (!is.null(climatology)) ncvar_put(out.nc, "time_bounds", climatology$time.bounds)

  nc_close(out.nc)
}
