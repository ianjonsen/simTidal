#' @title Compare speed along observed and simulated drifter tracks
#'
#' @description Computes and compares speed along a high-resolution observed
#'   track (1–2 sec) and a 10-minute simulated track. Two comparisons are
#'   provided:
#'   \enumerate{
#'     \item Observed speed aggregated to 10-minute windows, paired with
#'       simulated speed at those same intervals — for direct step-by-step
#'       comparison.
#'     \item Full speed distributions compared via empirical CDF and a
#'       two-sample Kolmogorov-Smirnov test — for asking whether the
#'       simulation reproduces the overall speed regime.
#'   }
#'
#' @param obs   Data frame of observed drifter positions with columns
#'   \code{date} (POSIXct UTC), \code{x} and \code{y} (UTM km).
#' @param sim   A \code{simTidal} object or its \code{$sim} data frame.
#' @param trim  Logical (default \code{TRUE}). Clip observed track to the
#'   simulation time window before computing speeds.
#'
#' @return A list of class \code{speed_comparison} with elements:
#'   \describe{
#'     \item{\code{paired}}{Tibble of 10-min window speeds: \code{date},
#'       \code{speed_obs_mean} (mean observed speed in window, m/s),
#'       \code{speed_obs_sd} (sd of observed speeds in window, m/s),
#'       \code{speed_sim} (simulated step speed, m/s).}
#'     \item{\code{obs_speeds}}{All instantaneous observed speeds (m/s),
#'       one per 1-sec step.}
#'     \item{\code{sim_speeds}}{All simulated step speeds (m/s).}
#'     \item{\code{ks_test}}{Result of \code{ks.test()} comparing the two
#'       speed distributions.}
#'   }
#'
#' @importFrom dplyr tibble
#' @export

compare_speed <- function(obs, sim, trim = TRUE) {

  if (inherits(sim, "simTidal")) sim <- sim$sim

  if (!all(c("date", "x", "y") %in% names(obs)))
    stop("obs must have columns: date, x, y")
  if (!all(c("date", "x", "y") %in% names(sim)))
    stop("sim must have columns: date, x, y")

  if (trim) {
    obs <- obs[obs$date >= min(sim$date) & obs$date <= max(sim$date), ]
    if (nrow(obs) == 0)
      stop("No observed positions fall within the simulation time window.")
  }

  ## ---- Simulated step speeds -----------------------------------------------
  ## Displacement between consecutive 10-min positions, converted to m/s.
  ## step_secs derived from the actual time differences to be robust.
  sim_dt   <- as.numeric(diff(sim$date), units = "secs")   # seconds per step
  sim_dx   <- diff(sim$x) * 1000   # km -> m
  sim_dy   <- diff(sim$y) * 1000
  sim_dist <- sqrt(sim_dx^2 + sim_dy^2)                    # metres per step

  sim_speeds_df <- dplyr::tibble(
    date      = sim$date[-1],          # timestamp at end of each step
    speed_sim = sim_dist / sim_dt      # m/s
  )

  ## ---- Observed instantaneous speeds ---------------------------------------
  ## Step distances between consecutive 1-sec positions, converted to m/s.
  obs   <- obs[order(obs$date), ]
  obs_dt   <- as.numeric(diff(obs$date), units = "secs")
  obs_dx   <- diff(obs$x) * 1000
  obs_dy   <- diff(obs$y) * 1000
  obs_dist <- sqrt(obs_dx^2 + obs_dy^2)
  obs_speed <- obs_dist / obs_dt   # m/s, one value per 1-sec step

  obs_speeds_df <- dplyr::tibble(
    date       = obs$date[-1],
    speed_obs  = obs_speed
  )

  ## ---- Aggregate observed speed to 10-minute windows ----------------------
  ## For each simulation step interval [t_k, t_{k+1}), compute the mean and
  ## sd of all observed instantaneous speeds falling within that window.
  ## This is the only valid basis for a paired comparison given the resolution
  ## mismatch — comparing a distribution summary to a single simulation value.

  n_steps   <- nrow(sim_speeds_df)
  step_mean <- numeric(n_steps)
  step_sd   <- numeric(n_steps)
  step_n    <- integer(n_steps)

  sim_boundaries <- c(sim$date[1], sim_speeds_df$date)  # n_steps + 1 boundaries

  for (k in seq_len(n_steps)) {
    in_window   <- obs_speeds_df$date >= sim_boundaries[k] &
                   obs_speeds_df$date <  sim_boundaries[k + 1]
    window_vals <- obs_speeds_df$speed_obs[in_window]
    step_mean[k] <- if (length(window_vals) > 0) mean(window_vals,  na.rm = TRUE) else NA_real_
    step_sd[k]   <- if (length(window_vals) > 1) sd(window_vals,    na.rm = TRUE) else NA_real_
    step_n[k]    <- length(window_vals)
  }

  paired <- dplyr::tibble(
    date           = sim_speeds_df$date,
    speed_obs_mean = step_mean,
    speed_obs_sd   = step_sd,
    n_obs          = step_n,
    speed_sim      = sim_speeds_df$speed_sim
  )

  ## ---- KS test on full speed distributions ---------------------------------
  ## Compares whether sim and obs speeds are drawn from the same distribution.
  ## Uses all 1-sec observed speeds vs all 10-min sim speeds — the KS test
  ## is valid for samples of unequal size.
  ks <- ks.test(obs_speeds_df$speed_obs, sim_speeds_df$speed_sim)

  result <- list(
    paired     = paired,
    obs_speeds = obs_speeds_df$speed_obs,
    sim_speeds = sim_speeds_df$speed_sim,
    ks_test    = ks
  )
  class(result) <- "speed_comparison"
  result
}


## -----------------------------------------------------------------------------
## print method
## -----------------------------------------------------------------------------

#' @exportS3Method print speed_comparison
print.speed_comparison <- function(x, ...) {
  p <- x$paired[!is.na(x$paired$speed_obs_mean), ]
  cat("-- Speed comparison summary --\n")
  cat(sprintf("  Simulation steps with obs coverage : %d / %d\n",
              nrow(p), nrow(x$paired)))
  cat(sprintf("  Mean observed speed (windowed)     : %.3f m/s\n",
              mean(p$speed_obs_mean)))
  cat(sprintf("  Mean simulated speed               : %.3f m/s\n",
              mean(p$speed_sim)))
  cat(sprintf("  Mean |obs - sim| per step          : %.3f m/s\n",
              mean(abs(p$speed_obs_mean - p$speed_sim), na.rm = TRUE)))
  cat("\n  KS test (observed vs simulated speed distribution):\n")
  cat(sprintf("    D = %.4f,  p = %.4f  %s\n",
              x$ks_test$statistic,
              x$ks_test$p.value,
              if (x$ks_test$p.value < 0.05) "[distributions differ]"
              else "[distributions not distinguishable]"))
  invisible(x)
}


## -----------------------------------------------------------------------------
## plot method — two panels
## -----------------------------------------------------------------------------

#' @exportS3Method plot speed_comparison
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_abline
#'   geom_step scale_colour_manual labs theme_minimal theme element_text
#'   element_blank coord_equal unit margin
plot.speed_comparison <- function(x, ...) {

  p <- x$paired[!is.na(x$paired$speed_obs_mean), ]

  ## Panel 1: paired scatter — simulated vs mean observed per window
  ## Points on the 1:1 line indicate perfect agreement.
  ## Vertical error bars show within-window variability of observed speeds.
  speed_max <- max(c(p$speed_obs_mean + p$speed_obs_sd,
                     p$speed_sim), na.rm = TRUE) * 1.05

  p1 <- ggplot2::ggplot(p, ggplot2::aes(x = speed_sim, y = speed_obs_mean)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = speed_obs_mean - speed_obs_sd,
                   ymax = speed_obs_mean + speed_obs_sd),
      width = 0, colour = "#2a6496", alpha = 0.3, linewidth = 0.3
    ) +
    ggplot2::geom_point(colour = "#2a6496", alpha = 0.6, size = 1.2) +
    ggplot2::coord_equal(xlim = c(0, speed_max), ylim = c(0, speed_max)) +
    ggplot2::labs(
      x        = "Simulated speed (m/s)",
      y        = "Mean observed speed per 10-min window (m/s)",
      title    = "Paired speed comparison",
      subtitle = "Each point = one 10-min step  |  bars = SD of observed speeds in window\nDashed line = 1:1"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 7, colour = "grey40"))

  ## Panel 2: empirical CDF comparison of full speed distributions
  ## Shows where and how much the distributions differ.
  n_obs  <- length(x$obs_speeds)
  n_sim  <- length(x$sim_speeds)

  cdf_df <- rbind(
    data.frame(
      speed = sort(x$obs_speeds),
      cdf   = seq_len(n_obs) / n_obs,
      track = "Observed (1-sec)"
    ),
    data.frame(
      speed = sort(x$sim_speeds),
      cdf   = seq_len(n_sim) / n_sim,
      track = "Simulated (10-min)"
    )
  )

  ks_label <- sprintf("KS: D = %.3f,  p = %.3f",
                      x$ks_test$statistic, x$ks_test$p.value)

  p2 <- ggplot2::ggplot(cdf_df, ggplot2::aes(x = speed, y = cdf,
                                               colour = track)) +
    ggplot2::geom_step(linewidth = 0.6) +
    ggplot2::scale_colour_manual(
      values = c("Observed (1-sec)"   = "#2c7bb6",
                 "Simulated (10-min)" = "#E63946")
    ) +
    ggplot2::labs(
      x        = "Speed (m/s)",
      y        = "Cumulative proportion",
      title    = "Speed distribution comparison",
      subtitle = ks_label,
      colour   = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.subtitle    = ggplot2::element_text(size = 8, colour = "grey30"),
      legend.position  = "bottom"
    )

  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(p1 | p2)
  } else {
    print(p1)
    print(p2)
  }
  invisible(x)
}
