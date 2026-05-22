#' @title Simulate tagged fish movement between acoustic detection events
#'
#' @description Runs an ensemble of \code{n_sim} fish movement simulations
#'   conditioned on two acoustic detection events (known times, uncertain
#'   locations). Each simulation:
#'   \enumerate{
#'     \item Samples a start position uniformly within \code{det.range} of
#'       the start receiver. For \code{move = "crw.bridge"}, an end
#'       position is also sampled within \code{det.range} of the end
#'       receiver.
#'     \item Propagates a movement model — with or without FVCOM tidal
#'       advection — for \code{N} steps, where \code{N} is derived
#'       automatically from the elapsed time between detections.
#'     \item Is flagged as \strong{accepted} if the final position falls
#'       within \code{det.range} of the end receiver at \code{end.dt}.
#'   }
#'
#' @param data  Output list from \code{sim_setup()}, or \code{NULL}.
#'   When \code{mpar$advect = TRUE} (the default), \code{data} must be a
#'   full \code{sim_setup()} object containing u/v current rasters,
#'   \code{fvcom.origin}, and \code{fvcom_step_secs}. When
#'   \code{mpar$advect = FALSE}, \code{data} is optional: pass it to
#'   retain land-avoidance behaviour, or omit it (\code{data = NULL}) to
#'   run unconstrained swimming-only simulations (useful when contemporaneous
#'   current data are unavailable).
#' @param mpar  Parameter object from \code{fish_par()}.
#' @param pb    Show a progress bar (logical, default \code{TRUE}).
#'
#' @return An object of class \code{sim_fish} with elements:
#'   \describe{
#'     \item{\code{sims}}{List of \code{n_sim} tibbles, each with columns
#'       \code{id}, \code{date}, \code{x}, \code{y}, \code{u}, \code{v}.
#'       All simulations are returned, including rejected ones.}
#'     \item{\code{accepted}}{Logical vector of length \code{n_sim};
#'       \code{TRUE} when the simulation ended within \code{det.range}
#'       of the end receiver.}
#'     \item{\code{end_dist}}{Numeric vector (km) — distance from the
#'       final simulated position to the end receiver, for each simulation.
#'       \code{NA} for simulations that stopped early (land or boundary).}
#'     \item{\code{n_accepted}}{Integer; number of accepted simulations.}
#'     \item{\code{acceptance_rate}}{Fraction of simulations accepted.}
#'     \item{\code{params}}{The \code{fish_par} object used.}
#'   }
#'
#' @seealso \code{\link{fish_par}}, \code{\link{fish_par_from_pe}},
#'   \code{\link{read_pe}}, \code{\link{sim_drifter}}, \code{\link{sim_setup}}
#'
#' @examples
#' \dontrun{
#' ## ------------------------------------------------------------------
#' ## Workflow A: coordinates and times supplied directly
#' ## ------------------------------------------------------------------
#' mpar <- fish_par(
#'   n_sim    = 200,
#'   start.dt = as.POSIXct("2024-06-01 12:46:31", tz = "UTC"),
#'   end.dt   = as.POSIXct("2024-06-02 23:11:27", tz = "UTC"),
#'   start    = c(387.487, 5023.234),   ## FORCE_15 receiver (km, UTM Zone 20N)
#'   end      = c(386.807, 5022.452),   ## MPS_09 receiver
#'   bearing  = pi,                     ## southward swimming bias
#'   rho      = 0.6,
#'   bl       = 2.0,
#'   fl       = 0.80
#' )
#'
#' ## With tidal advection (requires sim_setup() output)
#' data <- sim_setup(config, month = c("June", "July"))
#' out  <- sim_fish(data, mpar)
#' print(out)
#'
#' ## Accepted tracks only
#' accepted_tracks <- out$sims[out$accepted]
#'
#' ## ------------------------------------------------------------------
#' ## Workflow B: build mpar directly from passing-event detection files
#' ## ------------------------------------------------------------------
#' pe <- read_pe("passing_events2024.csv", "stn_position2024.csv")
#'
#' mpar <- fish_par_from_pe(
#'   pe, fish_id = 6, pe_start = 1, pe_end = 2,
#'   n_sim   = 300,
#'   bearing = pi,
#'   rho     = 0.6,
#'   bl      = 2.0,
#'   fl      = 0.80
#' )
#'
#' out <- sim_fish(data, mpar)
#'
#' ## ------------------------------------------------------------------
#' ## Workflow C: swimming only — no current data required
#' ## ------------------------------------------------------------------
#' mpar_noflow <- fish_par_from_pe(
#'   pe, fish_id = 6, pe_start = 1, pe_end = 2,
#'   n_sim   = 300,
#'   advect  = FALSE,
#'   bearing = pi,
#'   rho     = 0.6,
#'   bl      = 2.0,
#'   fl      = 0.80
#' )
#'
#' ## data can be omitted entirely when advect = FALSE
#' out_noflow <- sim_fish(NULL, mpar_noflow)
#' }
#'
#' @importFrom terra extract xyFromCell nlyr
#' @importFrom CircStats rwrpcauchy
#' @importFrom dplyr mutate
#' @importFrom tibble tibble
#' @importFrom stats runif rnorm
#' @export

sim_fish <- function(
    data = NULL,
    mpar = NULL,
    pb   = TRUE
) {

  ## ---- Validate mpar --------------------------------------------------------
  if (is.null(mpar))
    stop("mpar is NULL. Create with fish_par().")
  if (!inherits(mpar, "fish_par"))
    stop("mpar must be created with fish_par(), got class: ",
         paste(class(mpar), collapse = ", "))

  ## step_secs is always needed (date-sequence, swimming step size)
  step_secs <- mpar$time.step * 60L
  N    <- mpar$N
  nsim <- mpar$n_sim

  ## ---- FVCOM setup (advect = TRUE only) -------------------------------------
  ##
  ## When advect = FALSE the user may omit data entirely; land-avoidance
  ## still works if data$land is present (move_kernel already guards it).
  ## All FVCOM-specific objects are initialised to safe sentinels here so
  ## the loop body (which is unconditional) can reference them safely.

  start.dt     <- mpar$start.dt   ## exact detection time; may be snapped below
  advect_scale <- step_secs / 1000 ## m/s -> km/step (harmless when advect=FALSE)
  layer_idx    <- NULL             ## populated only when advect=TRUE
  w            <- 0L               ## temporal interpolation weight

  if (mpar$advect) {

    if (is.null(data))
      stop("data is NULL. Pass the output of sim_setup() when advect = TRUE.\n",
           "  To run without tidal advection set advect = FALSE in fish_par().")

    if (is.null(data$fvcom_step_secs))
      stop("data$fvcom_step_secs not found. Re-run sim_setup() to regenerate data.")

    if (step_secs != data$fvcom_step_secs)
      stop(sprintf(
        "mpar$time.step (%g min) does not match FVCOM layer interval (%g min).\n",
        mpar$time.step, data$fvcom_step_secs / 60),
        "  Set time.step = ", data$fvcom_step_secs / 60, " in fish_par()."
      )

    fvcom.origin <- if (!is.null(data$fvcom.origin)) {
      data$fvcom.origin
    } else {
      stop("fvcom.origin not found in data. Re-run sim_setup().")
    }

    n_u_layers  <- nlyr(data$u)
    data_end_dt <- fvcom.origin + n_u_layers * step_secs

    ## Sub-interval handling: align start.dt to the nearest FVCOM layer
    offset_secs <- as.numeric(difftime(mpar$start.dt, fvcom.origin, units = "secs"))
    remainder   <- offset_secs %% step_secs

    if (remainder == 0) {
      w <- 0L
    } else if (mpar$interp) {
      w <- remainder / step_secs
    } else {
      original_dt <- start.dt
      start.dt <- if (remainder < step_secs / 2) {
        start.dt - remainder
      } else {
        start.dt + (step_secs - remainder)
      }
      warning(sprintf(
        "start.dt (%s) is not on a %g-min FVCOM boundary.\n  Snapped to nearest layer: %s\n  Use interp = TRUE to interpolate instead.",
        format(original_dt), mpar$time.step, format(start.dt)
      ))
      w <- 0L
    }

    ## Bounds checks
    if (start.dt < fvcom.origin)
      stop("start.dt (", format(start.dt), ") is earlier than available FVCOM data.\n",
           "  data covers: ", format(fvcom.origin), " to ", format(data_end_dt))
    if (start.dt > data_end_dt)
      stop("start.dt (", format(start.dt), ") is later than available FVCOM data.\n",
           "  data covers: ", format(fvcom.origin), " to ", format(data_end_dt))

    sim_end_dt <- start.dt + N * step_secs
    if (sim_end_dt > data_end_dt)
      stop("Simulation end time (", format(sim_end_dt),
           ") extends beyond available FVCOM data (", format(data_end_dt), ").\n",
           "  Re-run sim_setup() with additional months.")

    fvcom.idx <- as.integer(round(
      as.numeric(difftime(start.dt, fvcom.origin, units = "secs")) / step_secs
    ))
    layer_idx <- seq_len(N - 1L) + fvcom.idx

    if (mpar$interp && w > 0 && max(layer_idx) >= n_u_layers)
      stop("Simulation timeframe requires layer ", max(layer_idx) + 1L,
           " but data$u only has ", n_u_layers, " layers.\n",
           "  Reduce N or re-run sim_setup() with additional months.")

  } ## end if (mpar$advect)

  ## Swimming step size: body-lengths/s -> km/step
  s_swim <- mpar$fl / 1000 * mpar$bl * 60 * mpar$time.step

  ## For bcrw.coa, pre-compute the heading from start receiver toward CoA
  ## (constant across simulations; small jitter in start_xy is negligible).
  ## For crw, a fresh random heading is drawn inside the loop each simulation.
  init_heading_coa <- if (mpar$move == "bcrw.coa") {
    delta <- mpar$coa - mpar$start
    atan2(delta[1], delta[2])
  } else {
    NULL
  }

  ## ---- Ensemble loop --------------------------------------------------------
  tracks   <- vector("list", nsim)
  accepted <- logical(nsim)
  end_dist <- rep(NA_real_, nsim)

  if (pb) {
    cat(sprintf(
      "Running %d fish simulations (N = %d steps, %.1f h)  [advect = %s]...\n",
      nsim, N, N * mpar$time.step / 60, mpar$advect
    ))
    tpb <- txtProgressBar(min = 0, max = nsim, style = 3)
  }

  for (sim_i in seq_len(nsim)) {

    ## -- Sample start position uniformly within det.range ---------------------
    ## sqrt(runif) gives uniform area density (not concentrated at centre).
    angle    <- runif(1, 0, 2 * pi)
    radius   <- sqrt(runif(1)) * mpar$det.range
    start_xy <- mpar$start + c(sin(angle), cos(angle)) * radius

    ## Fall back to receiver location if sampled point is on land
    if (!is.null(data) && !is.null(data$land) &&
        !is.na(extract(data$land, rbind(start_xy))[1, 1])) {
      start_xy <- mpar$start
    }

    ## For Brownian bridge: sample a per-simulation end position uniformly
    ## within det.range of the end receiver (same approach as start).
    sim_end_xy <- NULL
    if (mpar$move == "crw.bridge") {
      e_angle    <- runif(1, 0, 2 * pi)
      e_radius   <- sqrt(runif(1)) * mpar$det.range
      sim_end_xy <- mpar$end + c(sin(e_angle), cos(e_angle)) * e_radius
      ## Fall back to receiver location if sampled point is on land
      if (!is.null(data) && !is.null(data$land) &&
          !is.na(extract(data$land, rbind(sim_end_xy))[1, 1])) {
        sim_end_xy <- mpar$end
      }
    }

    ## State matrix: x | y | heading (radians).
    ## Heading is carried forward so move_kernel can implement correlated
    ## turns. Initial heading depends on move type:
    ##   crw        — random uniform draw (fresh each simulation)
    ##   bcrw       — fixed bias bearing from mpar
    ##   bcrw.coa   — direction from start receiver toward CoA (pre-computed)
    ##   crw.bridge — direction from sampled start toward sampled end
    init_heading <- switch(mpar$move,
      crw        = runif(1, 0, 2 * pi),
      bcrw       = mpar$bearing,
      bcrw.coa   = init_heading_coa,
      crw.bridge = atan2(sim_end_xy[1] - start_xy[1],
                         sim_end_xy[2] - start_xy[2]),
      NA_real_
    )

    xy <- matrix(NA_real_, N, 3)
    xy[1, ] <- c(start_xy, init_heading)

    u <- v <- numeric(N)
    early_stop <- FALSE

    ## -- N-step movement loop -------------------------------------------------
    for (i in 2:N) {

      ## Proposed position from swimming kernel (returns c(x, y, heading))
      ds_i <- move_kernel(data, xy = xy[i - 1L, ], mpar = mpar,
                          s = s_swim, i = i, end_xy = sim_end_xy)

      ## FVCOM advection (skipped entirely when advect = FALSE)
      if (mpar$advect) {
        k <- layer_idx[i - 1L]

        if (w == 0) {
          tmp.u <- extract(data$u[[k]], rbind(xy[i - 1L, 1:2]),
                           method = "simple")[1, 1] * advect_scale
          tmp.v <- extract(data$v[[k]], rbind(xy[i - 1L, 1:2]),
                           method = "simple")[1, 1] * advect_scale
        } else {
          tmp.u <- ((1 - w) * extract(data$u[[k]],      rbind(xy[i - 1L, 1:2]), method = "simple")[1, 1] +
                          w  * extract(data$u[[k + 1L]], rbind(xy[i - 1L, 1:2]), method = "simple")[1, 1]) *
                   advect_scale
          tmp.v <- ((1 - w) * extract(data$v[[k]],      rbind(xy[i - 1L, 1:2]), method = "simple")[1, 1] +
                          w  * extract(data$v[[k + 1L]], rbind(xy[i - 1L, 1:2]), method = "simple")[1, 1]) *
                   advect_scale
        }

        ## Separate multipliers for ebb (u < 0) vs flood (u >= 0) phases
        u[i] <- if (!is.na(tmp.u) && tmp.u < 0) tmp.u * mpar$uvm[1] else tmp.u * mpar$uvm[2]
        v[i] <- if (!is.na(tmp.u) && tmp.u < 0) tmp.v * mpar$uvm[3] else tmp.v * mpar$uvm[4]
      }

      ## Combine swim displacement + advection; carry heading from kernel
      xy[i, 1] <- ds_i[1] + u[i]
      xy[i, 2] <- ds_i[2] + v[i]
      xy[i, 3] <- ds_i[3]

      ## Land / boundary checks.
      ## The land check is skipped when data or data$land is NULL
      ## (move_kernel already applies the same guard internally).
      if (!is.null(data) && !is.null(data$land) &&
          !is.na(extract(data$land, rbind(xy[i, 1:2]))[1, 1])) {
        early_stop <- TRUE
        break
      }
      if (any(is.na(xy[i, 1:2]))) {
        early_stop <- TRUE
        break
      }
    }

    ## -- End-constraint check -------------------------------------------------
    if (early_stop || any(is.na(xy[N, 1:2]))) {
      end_dist[sim_i] <- NA_real_
      accepted[sim_i] <- FALSE
    } else {
      end_dist[sim_i] <- sqrt((xy[N, 1] - mpar$end[1])^2 +
                               (xy[N, 2] - mpar$end[2])^2)
      accepted[sim_i] <- end_dist[sim_i] <= mpar$det.range
    }

    ## -- Store track (trim to last valid row if simulation stopped early) ------
    n_valid <- if (early_stop) max(which(!is.na(xy[, 1])), 0L) else N

    tracks[[sim_i]] <- tibble::tibble(
      id   = sim_i,
      date = seq(start.dt, by = step_secs, length.out = n_valid),
      x    = xy[seq_len(n_valid), 1],
      y    = xy[seq_len(n_valid), 2],
      u    = u[seq_len(n_valid)],
      v    = v[seq_len(n_valid)]
    )

    if (pb) setTxtProgressBar(tpb, sim_i)
  }

  if (pb) close(tpb)

  n_accepted <- sum(accepted)
  message(sprintf("%d / %d simulations accepted (%.1f%%)",
                  n_accepted, nsim, 100 * n_accepted / nsim))

  ## ---- Output ----------------------------------------------------------------
  structure(
    list(
      sims            = tracks,
      accepted        = accepted,
      end_dist        = end_dist,
      n_accepted      = n_accepted,
      acceptance_rate = n_accepted / nsim,
      params          = mpar
    ),
    class = "sim_fish"
  )
}


## -----------------------------------------------------------------------------
## Print method
## -----------------------------------------------------------------------------

#' @exportS3Method print sim_fish
print.sim_fish <- function(x, ...) {
  p <- x$params
  cat("-- sim_fish --\n")
  cat(sprintf("  Method:              %s\n",  p$method))
  cat(sprintf("  Advection:           %s\n",  p$advect))
  cat(sprintf("  Simulations run:     %d\n",  p$n_sim))
  cat(sprintf("  Accepted:            %d  (%.1f%%)\n",
              x$n_accepted, 100 * x$acceptance_rate))
  cat(sprintf("  Detection range:     %.3f km\n", p$det.range))
  cat(sprintf("  N steps:             %d  (%.1f hours at %g-min intervals)\n",
              p$N, p$N * p$time.step / 60, p$time.step))
  cat(sprintf("  start.dt:            %s\n", format(p$start.dt)))
  cat(sprintf("  end.dt:              %s\n", format(p$end.dt)))
  ed <- x$end_dist[!is.na(x$end_dist)]
  if (length(ed) > 0) {
    cat(sprintf("  End distance (all):  mean %.3f km  (range %.3f – %.3f km)\n",
                mean(ed), min(ed), max(ed)))
    if (x$n_accepted > 0) {
      ea <- x$end_dist[x$accepted]
      cat(sprintf("  End distance (acc):  mean %.3f km  (range %.3f – %.3f km)\n",
                  mean(ea), min(ea), max(ea)))
    }
  }
  invisible(x)
}
