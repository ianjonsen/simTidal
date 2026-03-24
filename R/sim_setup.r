#' @title Pre-simulation setup
#'
#' @description load required rasters, receiver locations
#'
#' @author Ian Jonsen \email{jonsen@stochastic-qc.org}
#'
#' @param config - path to config.R script containing file.paths for required
#'   & optional environmental layers (see Details)
#' @param month - one or more month names matching the u/v subdirectories
#'   produced by \code{process_month()} in \code{create_envt.R}.
#'   Must be a subset of \code{c("May", "June", "July", "August")}.
#'   Multiple months are automatically sorted into calendar order.
#'   E.g. \code{"July"}, \code{c("July", "August")},
#'   or \code{c("May", "June", "July", "August")}.
#' @param year  - 2-digit year suffix, e.g. \code{"22"}
#'
#' @details Loads u and v as lazy file-backed SpatRasters built directly from
#'   the per-day files written by \code{process_month()}. When multiple months
#'   are supplied, the day-files from each month are concatenated in calendar
#'   order into a single file-backed stack — no large assembled file is needed
#'   or written. \code{terra::rast(file_vector)} reads only the specific tiles
#'   needed at each \code{extract()} call.
#'
#'   \code{fvcom.origin} and \code{fvcom_step_secs} are derived from the first
#'   two u layer names so neither needs to be set manually in \code{mpar}.
#'
#' @importFrom terra rast nlyr
#' @export
#'
sim_setup <- function(config = config,
                      month  = "July",
                      year   = "22") {

  valid_months <- c("May", "June", "July", "August")
  month <- match.arg(month, choices = valid_months, several.ok = TRUE)

  ## Enforce calendar order regardless of how the user supplied the vector
  month <- valid_months[sort(match(month, valid_months))]

  suppressWarnings(source(config, local = TRUE, echo = FALSE))
  if (is.null(prj)) prj <- "+proj=utm +zone=20 +units=km +datum=WGS84 +no_defs +type=crs"

  out <- list(
    bathy  = suppressWarnings(rast(bathy)),
    land   = suppressWarnings(rast(land)),
    d2land = suppressWarnings(rast(d2land)),
    grad   = suppressWarnings(rast(grad))
  )

  ## Collect u and v file paths across all requested months in calendar order.
  ## terra::rast(character_vector) spanning multiple months creates a single
  ## file-backed SpatRaster with no data loaded into RAM. Layer names encode
  ## the timestamp of each time step, so cross-month continuity is maintained
  ## automatically as long as the day-files are correctly ordered.
  u_files <- character(0)
  v_files <- character(0)

  for (m in month) {
    u_dir <- file.path(fvcom, "u", m)
    v_dir <- file.path(fvcom, "v", m)

    for (d in c(u_dir, v_dir))
      if (!dir.exists(d))
        stop("Directory not found: ", d,
             "\n  Run process_month('", m, "', year = '", year,
             "') in create_envt.R first.")

    mf_u <- sort(list.files(u_dir, pattern = "\\.tif$", full.names = TRUE))
    mf_v <- sort(list.files(v_dir, pattern = "\\.tif$", full.names = TRUE))

    if (length(mf_u) == 0) stop("No u raster files found in: ", u_dir)
    if (length(mf_v) == 0) stop("No v raster files found in: ", v_dir)
    if (length(mf_u) != length(mf_v))
      stop("Unequal number of u (", length(mf_u), ") and v (",
           length(mf_v), ") day-files in ", m, " directories.")

    u_files <- c(u_files, mf_u)
    v_files <- c(v_files, mf_v)
  }

  out[["u"]] <- suppressWarnings(rast(u_files))
  out[["v"]] <- suppressWarnings(rast(v_files))

  ## Derive fvcom.origin and step interval from the first two u layer names
  ## (encoded as "ua_YYYYMMDDTHHMMz"). Stored in data so sim_drifter() and
  ## validate_mpar() never need to read layer names themselves.
  stopifnot("Need at least 2 u layers to derive fvcom_step_secs" =
              terra::nlyr(out[["u"]]) >= 2L)

  parse_lyr_time <- function(nm)
    as.POSIXct(sub("^ua_", "", nm), format = "%Y%m%dT%H%Mz", tz = "UTC")

  u_names <- terra::names(out[["u"]])
  out[["fvcom.origin"]]    <- parse_lyr_time(u_names[1])
  out[["fvcom_step_secs"]] <- as.numeric(
    difftime(parse_lyr_time(u_names[2]),
             parse_lyr_time(u_names[1]), units = "secs")
  )

  out[["month"]] <- month   ## character vector, calendar-ordered
  out[["year"]]  <- year
  out[["prj"]]   <- prj

  return(out)
}
