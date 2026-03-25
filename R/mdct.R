#' @title Mean Distance at Common Times (MDCT)
#'
#' @description Compares a simulated drifter track from \code{sim_drifter()}
#'   against a real observed track by linearly interpolating the simulated
#'   positions to every observed timestamp, then computing the Euclidean
#'   distance at each moment. Returns a tidy data frame of per-observation
#'   distances plus a scalar summary.
#'
#' @param obs   Data frame of observed drifter positions with columns:
#'   \code{date} (POSIXct), \code{x} and \code{y} (UTM km). Typically
#'   sampled at 1–2 second intervals.
#' @param sim   A \code{simTidal} object from \code{sim_drifter()}, or the
#'   \code{$sim} data frame extracted from one.
#' @param trim  Logical (default \code{TRUE}). If \code{TRUE}, observed
#'   timestamps are clipped to the temporal extent of the simulation so
#'   no extrapolation occurs outside the simulated period.
#' @param smooth_obs  Integer; if > 1, the observed track is smoothed with a
#'   rolling mean of this window width (seconds) before comparison, reducing
#'   GPS noise at 1–2 sec sampling. Default \code{NULL} = no smoothing.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{\code{distances}}{Tibble with columns \code{date}, \code{x_obs},
#'       \code{y_obs}, \code{x_sim}, \code{y_sim}, \code{dist_km} — one row
#'       per observed position within the simulation window.}
#'     \item{\code{summary}}{Named numeric vector: \code{mdct} (mean distance,
#'       km), \code{rmsd} (root mean square distance, km), \code{max_dist}
#'       (maximum distance, km), \code{n} (number of matched observations).}
#'   }
#'
#' @examples
#' \dontrun{
#' ## obs is your real drifter data frame with date (POSIXct), x, y in UTM km
#' ## out is the simTidal object from sim_drifter()
#' result <- mdct(obs, out)
#' print(result$summary)
#' plot(result)
#' }
#'
#' @importFrom dplyr tibble mutate
#' @importFrom stats approx
#' @export

mdct <- function(obs, sim, trim = TRUE, smooth_obs = NULL) {

  ## Accept either a simTidal object or its $sim data frame directly
  if (inherits(sim, "simTidal")) sim <- sim$sim

  ## ---- Input checks --------------------------------------------------------
  required_cols <- c("date", "x", "y")
  if (!all(required_cols %in% names(obs)))
    stop("obs must have columns: date, x, y")
  if (!all(required_cols %in% names(sim)))
    stop("sim must have columns: date, x, y")
  if (!inherits(obs$date, "POSIXct"))
    stop("obs$date must be POSIXct")
  if (!inherits(sim$date, "POSIXct"))
    stop("sim$date must be POSIXct")
browser()
  ## ---- Optional: trim obs to simulation window -----------------------------
  if (trim) {
    obs <- obs[obs$date >= min(sim$date) & obs$date <= max(sim$date), ]
    if (nrow(obs) == 0)
      stop("No observed positions fall within the simulation time window.\n",
           "  Simulation: ", format(min(sim$date)), " to ", format(max(sim$date)), "\n",
           "  Observed:   ", format(min(obs$date)),  " to ", format(max(obs$date)))
  }

  ## ---- Optional: smooth observed track to reduce GPS noise -----------------
  ## Rolling mean over a window of smooth_obs seconds, using the actual
  ## irregular timestamps (not assuming constant rate).
  if (!is.null(smooth_obs) && smooth_obs > 1) {
    obs_t <- as.numeric(obs$date)
    obs$x <- .rolling_mean_time(obs$x, obs_t, smooth_obs)
    obs$y <- .rolling_mean_time(obs$y, obs_t, smooth_obs)
  }

  ## ---- Interpolate simulated track to observation times --------------------
  ## approx() performs linear interpolation. rule = 1 returns NA outside the
  ## sim range — trim = TRUE above ensures this won't happen.
  sim_t <- as.numeric(sim$date)
  obs_t <- as.numeric(obs$date)

  x_sim_interp <- approx(sim_t, sim$x, xout = obs_t, method = "linear", rule = 1)$y
  y_sim_interp <- approx(sim_t, sim$y, xout = obs_t, method = "linear", rule = 1)$y

  ## ---- Compute Euclidean distance at each observation ----------------------
  dist_km <- sqrt((obs$x - x_sim_interp)^2 + (obs$y - y_sim_interp)^2)

  ## Drop any NAs (shouldn't occur after trimming but guard anyway)
  valid   <- !is.na(dist_km)
  dist_km <- dist_km[valid]

  distances <- dplyr::tibble(
    date    = obs$date[valid],
    x_obs   = obs$x[valid],
    y_obs   = obs$y[valid],
    x_sim   = x_sim_interp[valid],
    y_sim   = y_sim_interp[valid],
    dist_km = dist_km
  )

  ## ---- Summary scalars -----------------------------------------------------
  summary_stats <- c(
    mdct     = mean(dist_km),
    rmsd     = sqrt(mean(dist_km^2)),
    max_dist = max(dist_km),
    n        = length(dist_km)
  )

  result <- list(distances = distances, summary = summary_stats)
  class(result) <- "mdct"
  result
}


## -----------------------------------------------------------------------------
## Internal helper: time-aware rolling mean
## Each output value is the mean of all inputs within +/- (window/2) seconds.
## -----------------------------------------------------------------------------
.rolling_mean_time <- function(x, t, window_secs) {
  half <- window_secs / 2
  vapply(seq_along(x), function(i) {
    in_window <- abs(t - t[i]) <= half
    mean(x[in_window], na.rm = TRUE)
  }, numeric(1))
}


## -----------------------------------------------------------------------------
## print method
## -----------------------------------------------------------------------------

#' @exportS3Method print mdct
print.mdct <- function(x, ...) {
  s <- x$summary
  cat("-- MDCT summary --\n")
  cat(sprintf("  Observations matched : %d\n",       as.integer(s["n"])))
  cat(sprintf("  Mean distance (MDCT) : %.4f km  (%.1f m)\n",
              s["mdct"], s["mdct"] * 1000))
  cat(sprintf("  RMSD                 : %.4f km  (%.1f m)\n",
              s["rmsd"], s["rmsd"] * 1000))
  cat(sprintf("  Max distance         : %.4f km  (%.1f m)\n",
              s["max_dist"], s["max_dist"] * 1000))
  invisible(x)
}


## -----------------------------------------------------------------------------
## plot method — three-panel diagnostic
## -----------------------------------------------------------------------------

#' @exportS3Method plot mdct
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_histogram
#'   facet_wrap labs theme_minimal theme element_text scale_colour_manual
#'   scale_x_datetime unit margin
plot.mdct <- function(x, bins = 60, ...) {

  d <- x$distances
  s <- x$summary

  subtitle_str <- sprintf(
    "MDCT = %.1f m  |  RMSD = %.1f m  |  max = %.1f m  |  n = %d",
    s["mdct"] * 1000, s["rmsd"] * 1000, s["max_dist"] * 1000, as.integer(s["n"])
  )

  ## Panel 1: distance over time
  p1 <- ggplot2::ggplot(d, ggplot2::aes(x = date, y = dist_km * 1000)) +
    ggplot2::geom_line(colour = "#2a6496", linewidth = 0.4, alpha = 0.7) +
    ggplot2::geom_hline(yintercept = s["mdct"] * 1000,
                        linetype = "dashed", colour = "#E63946", linewidth = 0.5) +
    ggplot2::labs(x = NULL, y = "Distance (m)",
                  title = "Distance between observed and simulated track over time",
                  subtitle = subtitle_str) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 7, colour = "grey40"))

  ## Panel 2: distance distribution
  p2 <- ggplot2::ggplot(d, ggplot2::aes(x = dist_km * 1000)) +
    ggplot2::geom_histogram(bins = bins, fill = "#2a6496", colour = "white",
                             linewidth = 0.2, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = s["mdct"] * 1000,
                        linetype = "dashed", colour = "#E63946", linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = s["rmsd"] * 1000,
                        linetype = "dotted", colour = "#f4a261", linewidth = 0.6) +
    ggplot2::labs(x = "Distance (m)", y = "Count",
                  title = "Distribution of distances",
                  caption = "Dashed = MDCT (mean)  |  Dotted = RMSD") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 7, colour = "grey40"))

  ## Panel 3: spatial divergence — obs vs sim coloured tracks
  track_df <- rbind(
    data.frame(x = d$x_obs, y = d$y_obs, track = "Observed"),
    data.frame(x = d$x_sim, y = d$y_sim, track = "Simulated")
  )

  p3 <- ggplot2::ggplot(track_df, ggplot2::aes(x = x, y = y, colour = track)) +
    ggplot2::geom_path(linewidth = 0.5, alpha = 0.8) +
    ggplot2::scale_colour_manual(values = c("Observed"  = "#2c7bb6",
                                            "Simulated" = "#E63946")) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Easting (km)", y = "Northing (km)",
                  title = "Track comparison", colour = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")

  ## Arrange with patchwork if available, otherwise print sequentially
  if (requireNamespace("patchwork", quietly = TRUE)) {
    print((p1 / p2) | p3)
  } else {
    print(p1)
    print(p2)
    print(p3)
  }
  invisible(x)
}
