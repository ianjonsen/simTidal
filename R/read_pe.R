## =============================================================================
## read_pe()          — load and merge passing-event detection data
## fish_par_from_pe() — build fish_par from a pair of detection events
## =============================================================================


##' @title Load and merge acoustic passing-event detection data
##'
##' @description Reads a passing-events CSV and a station-position CSV,
##'   builds precise POSIXct timestamps from the component columns
##'   (handling the \code{SS = 60} data-entry artefact via integer arithmetic),
##'   converts station positions from metres to km, and joins the two tables
##'   on station number.
##'
##' @param events   Path to the passing-events CSV, or an already-loaded
##'   \code{data.frame}. Expected columns (order matters; names are flexible):
##'   \code{fish_id}, \code{pe}, \code{stn}, \code{YEAR}, \code{MM},
##'   \code{DD}, \code{HH}, \code{MN}, \code{SS}.
##' @param stations Path to the station-position CSV, or a \code{data.frame}.
##'   Expected columns: \code{stn}, \code{x_m}, \code{y_m}, \code{lat},
##'   \code{lon}, \code{stn_name}. Stations with empty/NA positions are
##'   retained with \code{NA} coordinates and a warning is issued.
##' @param tz       Time zone for the timestamps. Default \code{"UTC"}.
##'   Change to e.g. \code{"America/Halifax"} if detection times are in
##'   local Atlantic time, but note that \code{sim_fish()} compares these
##'   against FVCOM data whose origin is in UTC.
##'
##' @return A tibble with one row per detection event and columns:
##'   \describe{
##'     \item{\code{fish_id}}{Fish identifier (integer).}
##'     \item{\code{pe}}{Passing-event index within fish (integer).}
##'     \item{\code{stn}}{Station number (integer).}
##'     \item{\code{stn_name}}{Station label (character), e.g. "MPS_{07}".}
##'     \item{\code{datetime}}{Detection time (POSIXct, \code{tz}).}
##'     \item{\code{x_km}}{Easting in km, UTM Zone 20N. \code{NA} for
##'       stations with no recorded position.}
##'     \item{\code{y_km}}{Northing in km, UTM Zone 20N.}
##'     \item{\code{lat}}{Latitude (decimal degrees).}
##'     \item{\code{lon}}{Longitude (decimal degrees).}
##'   }
##'
##' @importFrom tibble as_tibble
##' @export
read_pe <- function(events, stations, tz = "UTC") {

  ## ---- Load events ----------------------------------------------------------
  ev <- if (is.character(events)) {
    read.csv(events, header = TRUE, strip.white = TRUE, check.names = FALSE)
  } else {
    as.data.frame(events)
  }

  ## ---- Load stations --------------------------------------------------------
  st <- if (is.character(stations)) {
    read.csv(stations, header = TRUE, strip.white = TRUE, check.names = FALSE,
             na.strings = c("", " ", "NA"))
  } else {
    as.data.frame(stations)
  }

  ## ---- Rename to canonical names --------------------------------------------
  ## Events: 9 columns in order
  if (ncol(ev) != 9)
    stop("events must have exactly 9 columns: ",
         "fish_id, pe, stn, YEAR, MM, DD, HH, MN, SS.\n",
         "  Got ", ncol(ev), " columns.")
  names(ev) <- c("fish_id", "pe", "stn", "YEAR", "MM", "DD", "HH", "MN", "SS")

  ## Stations: 6 columns in order
  if (ncol(st) != 6)
    stop("stations must have exactly 6 columns: ",
         "stn, x_m, y_m, lat, lon, stn_name.\n",
         "  Got ", ncol(st), " columns.")
  names(st) <- c("stn", "x_m", "y_m", "lat", "lon", "stn_name")

  ## ---- Build POSIXct timestamps --------------------------------------------
  ## Add date components as total elapsed seconds since midnight to handle
  ## SS = 60 (a data-entry artefact meaning "the next whole minute"):
  ##   e.g. 17:29:60  ->  midnight + 63000 s  ->  17:30:00
  midnight <- as.POSIXct(
    paste(ev$YEAR,
          sprintf("%02d", as.integer(ev$MM)),
          sprintf("%02d", as.integer(ev$DD)),
          sep = "-"),
    tz = tz
  )
  total_secs  <- as.integer(ev$HH) * 3600L +
                 as.integer(ev$MN) * 60L   +
                 as.integer(ev$SS)
  ev$datetime <- midnight + total_secs

  ## ---- Clean station positions ----------------------------------------------
  st$x_m    <- suppressWarnings(as.numeric(st$x_m))
  st$y_m    <- suppressWarnings(as.numeric(st$y_m))
  st$lat    <- suppressWarnings(as.numeric(st$lat))
  st$lon    <- suppressWarnings(as.numeric(st$lon))
  st$x_km   <- st$x_m / 1000
  st$y_km   <- st$y_m / 1000
  st$stn_name <- trimws(as.character(st$stn_name))

  ## Warn about stations that have no position data
  no_pos <- st$stn[is.na(st$x_km) | is.na(st$y_km)]
  if (length(no_pos) > 0)
    message("Stations with no position data (will have NA coordinates): ",
            paste(sort(no_pos), collapse = ", "))

  ## ---- Merge events and positions ------------------------------------------
  merged <- merge(
    ev[, c("fish_id", "pe", "stn", "datetime")],
    st[, c("stn", "stn_name", "x_km", "y_km", "lat", "lon")],
    by     = "stn",
    all.x  = TRUE,
    sort   = FALSE
  )

  ## Sort by fish, then by event index
  merged <- merged[order(merged$fish_id, merged$pe), ]

  tibble::as_tibble(
    merged[, c("fish_id", "pe", "stn", "stn_name",
               "datetime", "x_km", "y_km", "lat", "lon")]
  )
}


##' @title Build fish_par from a pair of acoustic detection events
##'
##' @description Convenience wrapper around \code{\link{fish_par}} that
##'   extracts \code{start.dt}, \code{end.dt}, \code{start}, and \code{end}
##'   from a passing-event tibble (produced by \code{\link{read_pe}}) for a
##'   specified fish and pair of consecutive (or any) detection events.
##'   All remaining \code{fish_par} parameters (e.g. \code{n_sim},
##'   \code{bearing}, \code{rho}) are passed through \code{...}.
##'
##' @param pe_data  Tibble returned by \code{\link{read_pe}}.
##' @param fish_id  Fish identifier (matches the \code{fish_id} column).
##' @param pe_start Index of the start detection event (\code{pe} column).
##' @param pe_end   Index of the end detection event (\code{pe} column).
##'   Must correspond to a later timestamp than \code{pe_start}.
##' @param ...      Additional arguments forwarded to \code{\link{fish_par}}
##'   (e.g. \code{n_sim}, \code{move}, \code{bearing}, \code{rho},
##'   \code{bl}, \code{fl}, \code{det.range}, \code{advect}).
##'
##' @return Object of class \code{c("fish_par", "mpar")} from
##'   \code{\link{fish_par}}.
##'
##' @examples
##' \dontrun{
##' pe <- read_pe("passing_events2024.csv", "stn_position2024.csv")
##'
##' ## Parameters for dogfish 6, events 1 -> 2
##' mpar <- fish_par_from_pe(pe, fish_id = 6, pe_start = 1, pe_end = 2,
##'                          n_sim   = 300,
##'                          bearing = pi,   ## southward bias
##'                          rho     = 0.6,
##'                          bl      = 2.0,
##'                          fl      = 0.80)
##' out <- sim_fish(data, mpar)
##' }
##'
##' @export
fish_par_from_pe <- function(pe_data, fish_id, pe_start, pe_end, ...) {

  ## Subset to the requested fish
  fd <- pe_data[pe_data$fish_id == fish_id, ]
  if (nrow(fd) == 0)
    stop("No detection events found for fish_id = ", fish_id,
         ".\n  Available fish_ids: ",
         paste(sort(unique(pe_data$fish_id)), collapse = ", "))

  ## Extract start and end rows
  row_s <- fd[fd$pe == pe_start, ]
  row_e <- fd[fd$pe == pe_end,   ]

  if (nrow(row_s) == 0)
    stop("pe_start = ", pe_start, " not found for fish_id = ", fish_id,
         ".\n  Available pe values: ",
         paste(sort(fd$pe), collapse = ", "))
  if (nrow(row_e) == 0)
    stop("pe_end = ", pe_end, " not found for fish_id = ", fish_id,
         ".\n  Available pe values: ",
         paste(sort(fd$pe), collapse = ", "))

  ## Guard against multiple matched rows (shouldn't happen with clean data)
  if (nrow(row_s) > 1)
    stop("Multiple rows matched pe_start = ", pe_start,
         " for fish_id = ", fish_id, ". Check pe_data for duplicates.")
  if (nrow(row_e) > 1)
    stop("Multiple rows matched pe_end = ", pe_end,
         " for fish_id = ", fish_id, ". Check pe_data for duplicates.")

  ## Require position data for both endpoints
  if (is.na(row_s$x_km) || is.na(row_s$y_km))
    stop("Start station (stn = ", row_s$stn, ") has no position data.\n",
         "  Check the stations file.")
  if (is.na(row_e$x_km) || is.na(row_e$y_km))
    stop("End station (stn = ", row_e$stn, ") has no position data.\n",
         "  Check the stations file.")

  fish_par(
    start.dt = row_s$datetime,
    end.dt   = row_e$datetime,
    start    = c(row_s$x_km, row_s$y_km),
    end      = c(row_e$x_km, row_e$y_km),
    ...
  )
}
