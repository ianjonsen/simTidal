#' @title Pre-simulation setup
#'
#' @description load required rasters, receiver locations
#'
#' @author Ian Jonsen \email{jonsen@stochastic-qc.org}
#'
#' @param config - path to config.R script containing file.paths for required
#'   & optional environmental layers (see Details)
#' @param month - month name matching the u/v subdirectories produced by
#'   \code{process_month()} in \code{create_envt.R}, e.g. \code{"July"}.
#'   Must be one of \code{"May"}, \code{"June"}, \code{"July"}, \code{"August"}.
#' @param year  - 2-digit year suffix, e.g. \code{"22"}
#'
#' @details Loads u and v as lazy file-backed SpatRasters built directly from
#'   the per-day files written by \code{process_month()}. No large assembled
#'   file is needed or written: \code{terra::rast(file_vector)} provides an
#'   identical interface to a single-file stack but reads only the specific
#'   tiles needed at each \code{extract()} call.
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

  month <- match.arg(month, choices = c("May", "June", "July", "August"))

  suppressWarnings(source(config, local = TRUE, echo = FALSE))
  if (is.null(prj)) prj <- "+proj=utm +zone=20 +units=km +datum=WGS84 +no_defs +type=crs"

  out <- list(
    bathy  = suppressWarnings(rast(bathy)),
    land   = suppressWarnings(rast(land)),
    d2land = suppressWarnings(rast(d2land)),
    grad   = suppressWarnings(rast(grad))
  )

  ## Load u and v as lazy stacks from per-day files produced by process_month().
  ## terra::rast(character_vector) creates a file-backed SpatRaster that spans
  ## all files without loading any data into RAM — equivalent to a single large
  ## file but without the cost of writing one (no 14+ GB assembly step).
  ## Files are sorted numerically by day index (zero-padded: u_July01.tif, ...).
  u_dir <- file.path(fvcom, "u", month)
  v_dir <- file.path(fvcom, "v", month)

  for (d in c(u_dir, v_dir))
    if (!dir.exists(d))
      stop("Directory not found: ", d,
           "\n  Run process_month('", month, "', year = '", year,
           "') in create_envt.R first.")

  u_files <- sort(list.files(u_dir, pattern = "\\.tif$", full.names = TRUE))
  v_files <- sort(list.files(v_dir, pattern = "\\.tif$", full.names = TRUE))

  if (length(u_files) == 0) stop("No u raster files found in: ", u_dir)
  if (length(v_files) == 0) stop("No v raster files found in: ", v_dir)
  if (length(u_files) != length(v_files))
    stop("Unequal number of u (", length(u_files), ") and v (",
         length(v_files), ") day-files in ", month, " directories.")

  out[["u"]] <- suppressWarnings(rast(u_files))
  out[["v"]] <- suppressWarnings(rast(v_files))

  ## Derive fvcom.origin and step interval from the first two u layer names
  ## (encoded as "ua_YYYYMMDDTHHMMz"). Stored in data so sim_drifter() and
  ## validate_mpar() never need to read layer names themselves.
  parse_lyr_time <- function(nm)
    as.POSIXct(sub("^ua_", "", nm), format = "%Y%m%dT%H%Mz", tz = "UTC")
  u_names <- terra::names(out[["u"]])
  out[["fvcom.origin"]]    <- parse_lyr_time(u_names[1])
  out[["fvcom_step_secs"]] <- as.numeric(
    difftime(parse_lyr_time(u_names[2]),
             parse_lyr_time(u_names[1]), units = "secs")
  )

  out[["month"]] <- month
  out[["year"]]  <- year
  out[["prj"]]   <- prj

  return(out)
}
