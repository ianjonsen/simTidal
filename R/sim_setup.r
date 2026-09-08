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
#' @param axis.tol - tolerance for the FVCOM time-axis check, as a fraction of
#'   one step. Layer \code{k} must sit within \code{axis.tol * step} of
#'   \code{origin + k * step}, otherwise \code{sim_fish()}'s arithmetic layer
#'   lookup would select the wrong layer. Default 0.5 (half a step), which is
#'   exactly the point at which rounding would pick a neighbour. Label jitter
#'   from float32 time storage is typically ~0.3 of a step and passes; a missing
#'   day-file is off by 144 steps and does not.
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
                      year,
                      axis.tol = 0.5) {

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

  ## Derive fvcom.origin and fvcom_step_secs from the u layer names (encoded as
  ## "ua_YYYYMMDDTHHMMz"). Stored in data so sim_fish(), sim_drifter() and
  ## validate_mpar() never need to read layer names themselves.
  stopifnot("Need at least 3 u layers to fit a time axis" =
              terra::nlyr(out[["u"]]) >= 3L)

  parse_lyr_time <- function(nm)
    as.POSIXct(sub("^[uv]a_", "", nm), format = "%Y%m%dT%H%Mz", tz = "UTC")

  u_names <- terra::names(out[["u"]])
  u_times <- parse_lyr_time(u_names)

  if (anyNA(u_times))
    stop("Could not parse timestamps from ", sum(is.na(u_times)),
         " u layer name(s), e.g. '", u_names[which(is.na(u_times))[1]], "'.\n",
         "  Expected the form 'ua_YYYYMMDDTHHMMz' written by process_month().")

  n_lyr <- length(u_times)
  k     <- seq_len(n_lyr) - 1L

  ## ---- Fit the time axis, rather than reading it off two adjacent layers ----
  ##
  ## FVCOM records time as days since 1858-11-17 (MJD). Where that variable is
  ## stored as NC_FLOAT, the 24-bit mantissa quantises modern dates — MJD is
  ## around 60,000, so the representable spacing is 2^-8 d = 5.625 min. A clean
  ## 10-min output axis is therefore *labelled* with timestamps that jitter by
  ## up to +/- 2.8 min, giving minute-steps of 11, 6, 11, 11, 11, 6, 11, 12, ...
  ## in a 9-step cycle. The layers themselves are exactly 10 min apart; only
  ## their labels are not.
  ##
  ## Deriving the step from the first two labels therefore returns 11 min for a
  ## 10-min dataset. Fit the axis across the whole stack instead: the endpoints
  ## carry bounded quantisation error, so the mean step recovers the true
  ## interval to within a fraction of a second over thousands of layers.

  step_raw  <- as.numeric(difftime(u_times[n_lyr], u_times[1L], units = "secs")) /
               (n_lyr - 1L)
  step_secs <- round(step_raw / 60) * 60      ## FVCOM intervals are whole minutes

  if (step_secs <= 0)
    stop("Could not determine the FVCOM step interval: fitted ",
         round(step_raw, 1), " s across ", n_lyr, " layers.\n",
         "  Are the day-files sorted correctly?")

  ## Intercept that centres the residuals, rounded to the nearest whole minute
  ## (the true axis sits on whole minutes; only the labels do not).
  origin_num <- mean(as.numeric(u_times) - k * step_secs)
  origin     <- as.POSIXct(round(origin_num / 60) * 60,
                           origin = "1970-01-01", tz = "UTC")

  out[["fvcom.origin"]]    <- origin
  out[["fvcom_step_secs"]] <- step_secs
  out[["u_times"]]         <- u_times   ## the labels as written, for diagnostics

  ## ---- Verify the assumption sim_fish() actually makes ----------------------
  ##
  ## sim_fish() and sim_drifter() locate a layer as
  ## round((t - fvcom.origin) / fvcom_step_secs). What that needs is not that
  ## consecutive labels differ by a constant, but that layer k sits close enough
  ## to origin + k * step for the rounding to pick the right layer. Label jitter
  ## of a few minutes is harmless; a missing day-file displaces everything after
  ## it by a full day (144 steps at 10 min) and is caught immediately. A day
  ## file with 143 or 145 layers shows up the same way.

  dev   <- as.numeric(u_times) - (as.numeric(origin) + k * step_secs)   ## seconds
  tol   <- axis.tol * step_secs
  worst <- which.max(abs(dev))

  if (abs(dev[worst]) >= tol) {
    first_bad <- which(abs(dev) >= tol)[1]
    stop("FVCOM layer ", first_bad, " is not where its index places it — ",
         "layer indexing would be wrong.\n",
         "  Fitted axis: origin ", format(origin), ", step ", step_secs / 60,
         " min, ", n_lyr, " layers.\n",
         "  Layer ", first_bad, " is labelled ", format(u_times[first_bad]),
         " but index ", first_bad, " places it at ",
         format(origin + (first_bad - 1L) * step_secs), " (off by ",
         round(dev[first_bad] / 60, 1), " min; tolerance +/- ",
         round(tol / 60, 1), " min).\n",
         "  ", sum(abs(dev) >= tol), " of ", n_lyr, " layers are out of place.\n",
         "  This is what a missing or short day-file looks like. Check the day ",
         "files for the month(s) around that date and re-run process_month().")
  }

  if (max(abs(dev)) > 60)
    message(sprintf(
      paste0("  note: layer labels jitter by up to %.1f min around a clean %g-min axis\n",
             "        (float32 MJD precision in the FVCOM time variable; harmless —\n",
             "        the layers are regularly spaced, only their labels are not)."),
      max(abs(dev)) / 60, step_secs / 60))

  ## v must carry the same labels as u
  v_times <- parse_lyr_time(terra::names(out[["v"]]))
  if (anyNA(v_times) || length(v_times) != n_lyr || !all(v_times == u_times))
    stop("u and v time axes differ (", n_lyr, " u layers vs ", length(v_times),
         " v layers). Re-run process_month() for the affected month(s).")

  out[["month"]] <- month   ## character vector, calendar-ordered
  out[["year"]]  <- year
  out[["prj"]]   <- prj

  return(out)
}
