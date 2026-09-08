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
#'   Must be a subset of \code{c("Jan", "Feb", "Mar", "April", "May",
#'   "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec")}. Multiple months
#'   are automatically sorted into calendar order. E.g. \code{"July"},
#'   \code{c("July", "Aug")}, or \code{c("Nov", "Dec")}.
#'
#'   Note these are the exact directory names \code{process_month()} writes,
#'   not \code{month.abb} or \code{month.name}: "April" not "Apr", "Sept"
#'   not "Sep", "Aug" not "August".
#' @param year  - 4-digit year as a character string, e.g. \code{"2025"}.
#'   Used to locate the year-level subdirectory in the FVCOM current data tree
#'   (\code{{fvcom}/u/{year}/{Month}/}).
#'
#' @details Loads u and v as lazy file-backed SpatRasters built directly from
#'   the per-day files written by \code{process_month()}. When multiple months
#'   are supplied, the day-files from each month are concatenated in calendar
#'   order into a single file-backed stack — no large assembled file is needed
#'   or written. \code{terra::rast(file_vector)} reads only the specific tiles
#'   needed at each \code{extract()} call.
#'
#'   Current data are expected at
#'   \code{{fvcom}/u/{year}/{Month}/u_{Month}DD.tif} (and the same under
#'   \code{v/}), where \code{{fvcom}} is the root directory defined in
#'   \code{config.R} and \code{year} is the 4-digit year passed to this
#'   function. This layout supports multi-year data archives.
#'
#'   \code{fvcom.origin} and \code{fvcom_step_secs} are derived from the first
#'   two u layer names so neither needs to be set manually in \code{mpar}.
#'
#'   \strong{Important:} the returned list contains \pkg{terra} SpatRasters
#'   backed by C++ objects that do \emph{not} survive \code{saveRDS()} /
#'   \code{readRDS()}. Always call \code{sim_setup()} fresh in the same R
#'   session as \code{sim_fish()} or \code{sim_drifter()}. If you need to
#'   cache the object use \code{terra::wrap()} before saving and
#'   \code{terra::unwrap()} after loading.
#'
#' @importFrom terra rast nlyr
#' @export
#'
sim_setup <- function(config = config,
                      month  = "July",
                      year) {

  if (missing(year) || !is.character(year) || length(year) != 1 ||
      !grepl("^[0-9]{4}$", year))
    stop("year must be a 4-digit character string, e.g. \"2025\"")

  valid_months <- c("Jan", "Feb", "Mar", "April", "May", "June",
                    "July", "Aug", "Sept", "Oct", "Nov", "Dec")
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
    u_dir <- file.path(fvcom, "u", year, m)
    v_dir <- file.path(fvcom, "v", year, m)

    for (d in c(u_dir, v_dir))
      if (!dir.exists(d))
        stop("Directory not found: ", d,
             "\n  Expected layout: {fvcom}/u/", year, "/", m, "/",
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
  u_times <- parse_lyr_time(u_names)

  if (anyNA(u_times))
    stop("Could not parse timestamps from ", sum(is.na(u_times)),
         " u layer name(s), e.g. '", u_names[which(is.na(u_times))[1]], "'.\n",
         "  Expected the form 'ua_YYYYMMDDTHHMMz' written by process_month().")

  out[["fvcom.origin"]]    <- u_times[1]
  out[["fvcom_step_secs"]] <- as.numeric(difftime(u_times[2], u_times[1],
                                                  units = "secs"))

  ## ---- Assert a regular time axis -------------------------------------------
  ##
  ## sim_fish() and sim_drifter() locate the FVCOM layer for a given time
  ## arithmetically, as round((t - fvcom.origin) / fvcom_step_secs). That is
  ## only correct if layer position tracks wall-clock time exactly. A missing
  ## day-file, or a day-file with the wrong number of layers, breaks the
  ## mapping for every step after it and advects by the wrong tidal phase with
  ## no error raised. Two such gaps have already occurred (May 2024, Aug 2019),
  ## so this is checked rather than assumed.

  d_secs <- as.numeric(diff(u_times), units = "secs")
  bad    <- which(d_secs != out[["fvcom_step_secs"]])

  if (length(bad)) {
    gap1 <- bad[1]
    stop("FVCOM time axis is not evenly spaced — layer indexing would be wrong.\n",
         "  Expected a constant ", out[["fvcom_step_secs"]] / 60, "-min step.\n",
         "  First break: layer ", gap1, " (", format(u_times[gap1]), ") -> layer ",
         gap1 + 1L, " (", format(u_times[gap1 + 1L]), "), gap = ",
         round(d_secs[gap1] / 3600, 2), " h.\n",
         "  ", length(bad), " break(s) in total across ", length(u_times), " layers.\n",
         "  Re-run process_month() for the affected month(s) in create_envt.R.")
  }

  ## v must carry the same time axis as u
  v_times <- parse_lyr_time(terra::names(out[["v"]]))
  if (length(v_times) != length(u_times) || !all(v_times == u_times))
    stop("u and v time axes differ (", length(u_times), " vs ", length(v_times),
         " layers). Re-run process_month() for the affected month(s).")

  out[["month"]] <- month   ## character vector, calendar-ordered
  out[["year"]]  <- year
  out[["prj"]]   <- prj

  return(out)
}
