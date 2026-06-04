#' @title Simulate tagged fish movement between acoustic detection events
#'
#' @description Runs an ensemble of \code{n_sim} fish movement simulations
#'   conditioned on two acoustic detection events (known times, uncertain
#'   locations). Each simulation:
#'   \enumerate{
#'     \item Samples a start position (and, for \code{move = "crw.bridge"},
#'       also an end position) uniformly within \code{det.range} of the
#'       respective receiver.
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
#'       Accepted tracks are truncated at the detection step; all others
#'       run to step \code{N} (or the last valid step if stopped early).}
#'     \item{\code{accepted}}{Logical vector of length \code{n_sim};
#'       \code{TRUE} when the simulation passed within \code{det.range}
#'       of the end receiver at any step.}
#'     \item{\code{end_dist}}{Numeric vector (km) — distance to the end
#'       receiver at the detection step (accepted) or final step (rejected).
#'       \code{NA} for simulations stopped early by a boundary condition.}
#'     \item{\code{det_step}}{Integer vector — step index at which each
#'       accepted simulation first entered \code{det.range}. \code{NA}
#'       for rejected simulations.}
#'     \item{\code{det_locs}}{Tibble with columns \code{id}, \code{x},
#'       \code{y}, \code{date} — one row per accepted simulation giving
#'       the position and time of detection.}
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
#' data <- sim_setup(config, month = c("June", "July"), year = "2024")
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
#' @importFrom terra extract nlyr
#' @importFrom CircStats rwrpcauchy
#' @importFrom tibble tibble
#' @importFrom stats runif
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

  ## ---- Validate terra objects -----------------------------------------------
  ## terra SpatRasters rely on a C++ object that can become invalid in several
  ## ways — most commonly by saving with saveRDS() and reloading without
  ## terra::wrap() / terra::unwrap(), or if terra was updated without
  ## restarting R (causing a DLL / Rcpp method-pointer mismatch).
  ## Detect either failure early and emit a clear message rather than the
  ## cryptic "NULL value passed as symbol address" error from deep inside terra.
  if (!is.null(data)) {
    for (.nm in intersect(names(data), c("u", "v", "land", "bathy", "d2land", "grad"))) {
      .r <- data[[.nm]]
      if (inherits(.r, "SpatRaster")) {
        .ok <- tryCatch({ terra::nlyr(.r); TRUE }, error = function(e) FALSE)
        if (!.ok)
          stop(
            "data$", .nm, " is an invalid SpatRaster.\n",
            "  Common causes and fixes:\n",
            "  1. data was saved with saveRDS() and reloaded without terra::wrap() /\n",
            "     terra::unwrap() — re-run sim_setup() in the current session.\n",
            "  2. terra was updated without restarting R (DLL / Rcpp pointer mismatch)\n",
            "     — restart R, then re-run sim_setup().\n",
            "  3. Corrupted terra installation — reinstall with install.packages('terra').",
            call. = FALSE
          )
      }
    }
    rm(.nm, .r, .ok)
  }

  step_secs <- mpar$time.step * 60L
  N    <- mpar$N
  nsim <- mpar$n_sim

  ## ---- FVCOM setup (advect = TRUE only) -------------------------------------
  start.dt     <- mpar$start.dt
  advect_scale <- step_secs / 1000   ## m/s -> km/step
  layer_idx    <- NULL
  w            <- 0L

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

  ## Hoist invariant flag (avoids repeated NULL checks inside the step loop)
  has_land <- !is.null(data) && !is.null(data$land)

  ## Pre-compute date sequence (reused for every track)
  date_seq <- seq(start.dt, by = step_secs, length.out = N)

  ## For bcrw.coa: initial heading from start receiver toward CoA
  init_heading_coa <- if (mpar$move == "bcrw.coa") {
    delta <- mpar$coa - mpar$start
    atan2(delta[1], delta[2])
  } else {
    NULL
  }

  ## ---- Pre-allocate arrays --------------------------------------------------
  ##
  ## Step-first architecture: all nsim simulations are advanced one step at
  ## a time, allowing terra::extract() and rwrpcauchy() to be called once
  ## per step rather than once per simulation per step.
  ##
  ## xy_all[sim, step, c(x,y,heading)] — full position history.
  ## u_mat / v_mat store advection displacements for the output tibbles.

  xy_all         <- array(NA_real_,    c(nsim, N, 3L))
  u_mat          <- matrix(0,           nsim, N)
  v_mat          <- matrix(0,           nsim, N)
  active         <- rep(TRUE,           nsim)
  detected       <- rep(FALSE,          nsim)
  early_stop_vec <- rep(FALSE,          nsim)
  det_step       <- rep(NA_integer_,    nsim)
  end_dist       <- rep(NA_real_,       nsim)

  ## ---- Vectorised start-position sampling ----------------------------------
  ang_s <- runif(nsim, 0, 2 * pi)
  rad_s <- sqrt(runif(nsim)) * mpar$det.range
  s_x   <- mpar$start[1] + sin(ang_s) * rad_s
  s_y   <- mpar$start[2] + cos(ang_s) * rad_s

  if (has_land) {
    on_land_s <- !is.na(terra::extract(data$land, cbind(s_x, s_y))[, 1])
    s_x[on_land_s] <- mpar$start[1]
    s_y[on_land_s] <- mpar$start[2]
  }

  ## Bridge: vectorised end-position sampling
  sim_end_x <- sim_end_y <- NULL
  if (mpar$move == "crw.bridge") {
    ang_e     <- runif(nsim, 0, 2 * pi)
    rad_e     <- sqrt(runif(nsim)) * mpar$det.range
    sim_end_x <- mpar$end[1] + sin(ang_e) * rad_e
    sim_end_y <- mpar$end[2] + cos(ang_e) * rad_e
    if (has_land) {
      on_land_e <- !is.na(terra::extract(data$land, cbind(sim_end_x, sim_end_y))[, 1])
      sim_end_x[on_land_e] <- mpar$end[1]
      sim_end_y[on_land_e] <- mpar$end[2]
    }
  }

  ## Vectorised initial headings
  init_heads <- switch(mpar$move,
    crw        = runif(nsim, 0, 2 * pi),
    bcrw       = rep(mpar$bearing, nsim),
    bcrw.coa   = rep(init_heading_coa, nsim),
    crw.bridge = atan2(sim_end_x - s_x, sim_end_y - s_y),
    rep(NA_real_, nsim)
  )

  xy_all[, 1L, 1L] <- s_x
  xy_all[, 1L, 2L] <- s_y
  xy_all[, 1L, 3L] <- init_heads

  ## ---- Progress bar ---------------------------------------------------------
  if (pb) {
    cat(sprintf(
      "Running %d fish simulations (N = %d steps, %.1f h)  [advect = %s]...\n",
      nsim, N, N * mpar$time.step / 60, mpar$advect
    ))
    tpb <- txtProgressBar(min = 1L, max = N, style = 3)
  }

  ## ---- Step-first ensemble loop --------------------------------------------
  ##
  ## At each step, all active simulations are advanced simultaneously:
  ##   - terra::extract() is called once per step (not once per simulation),
  ##     reducing calls from nsim*(N-1) to (N-1) for u, v, and land.
  ##   - rwrpcauchy() is called once per step with a vector of headings.
  ##
  ## Land avoidance: if a proposed position is on land, the fish reverts to
  ## its previous position and continues (reflecting boundary). This replaces
  ## the gradient-based repulsion in move_kernel, which cannot be batched.

  for (i in seq.int(2L, N)) {

    act   <- which(active)
    n_act <- length(act)
    if (n_act == 0L) break

    px <- xy_all[act, i - 1L, 1L]
    py <- xy_all[act, i - 1L, 2L]
    ph <- xy_all[act, i - 1L, 3L]

    ## Step size (handles time-varying bl in bcrw)
    s_i <- if (mpar$move == "bcrw" && length(mpar$bl) > 1L)
      mpar$fl / 1000 * mpar$bl[i] * 60 * mpar$time.step
    else
      s_swim

    ## Angular concentration (handles time-varying rho in bcrw)
    rho_i <- if (length(mpar$rho) > 1L) mpar$rho[i] else mpar$rho

    ## Mean heading for each active simulation
    phi <- switch(mpar$move,
      bcrw = rep(
        if (length(mpar$bearing) > 1L) mpar$bearing[i] else mpar$bearing,
        n_act),
      crw = ph,
      bcrw.coa = {
        d_x <- mpar$coa[1] - px
        d_y <- mpar$coa[2] - py
        psi <- atan2(d_x, d_y)
        atan2(sin(ph) + mpar$nu * sin(psi),
              cos(ph) + mpar$nu * cos(psi))
      },
      crw.bridge = {
        nu_i <- mpar$nu * i / max(mpar$N - i, 1L)
        d_x  <- sim_end_x[act] - px
        d_y  <- sim_end_y[act] - py
        psi  <- atan2(d_x, d_y)
        atan2(sin(ph) + nu_i * sin(psi),
              cos(ph) + nu_i * cos(psi))
      }
    )

    ## Draw headings and compute proposed positions (all active sims at once)
    mu    <- rwrpcauchy(n_act, phi, rho_i)
    new_x <- px + sin(mu) * s_i
    new_y <- py + cos(mu) * s_i

    ## Batched FVCOM advection — 2 extract calls regardless of nsim
    if (mpar$advect) {
      k     <- layer_idx[i - 1L]
      pos_m <- cbind(px, py)

      if (w == 0) {
        u_raw <- terra::extract(data$u[[k]],      pos_m, method = "simple")[, 1]
        v_raw <- terra::extract(data$v[[k]],      pos_m, method = "simple")[, 1]
      } else {
        u_raw <- ((1 - w) * terra::extract(data$u[[k]],        pos_m, method = "simple")[, 1] +
                        w  * terra::extract(data$u[[k + 1L]],  pos_m, method = "simple")[, 1])
        v_raw <- ((1 - w) * terra::extract(data$v[[k]],        pos_m, method = "simple")[, 1] +
                        w  * terra::extract(data$v[[k + 1L]],  pos_m, method = "simple")[, 1])
      }

      u_adj <- ifelse(!is.na(u_raw) & u_raw < 0,
                      u_raw * mpar$uvm[1], u_raw * mpar$uvm[2]) * advect_scale
      v_adj <- ifelse(!is.na(u_raw) & u_raw < 0,
                      v_raw * mpar$uvm[3], v_raw * mpar$uvm[4]) * advect_scale
      u_adj[is.na(u_adj)] <- 0
      v_adj[is.na(v_adj)] <- 0

      new_x         <- new_x + u_adj
      new_y         <- new_y + v_adj
      u_mat[act, i] <- u_adj
      v_mat[act, i] <- v_adj
    }

    ## Land check — 1 batched extract; revert fish that landed on land
    if (has_land) {
      on_land <- !is.na(terra::extract(data$land, cbind(new_x, new_y))[, 1])
      if (any(on_land)) {
        new_x[on_land]          <- px[on_land]
        new_y[on_land]          <- py[on_land]
        mu[on_land]             <- ph[on_land]
        u_mat[act[on_land], i]  <- 0
        v_mat[act[on_land], i]  <- 0
      }
    }

    ## NA check (raster boundary): mark simulation as stopped
    is_na <- is.na(new_x) | is.na(new_y)
    if (any(is_na)) {
      early_stop_vec[act[is_na]] <- TRUE
      active[act[is_na]]         <- FALSE
    }

    ## Store updated positions for all active (non-NA) sims
    still <- !is_na
    xy_all[act[still], i, 1L] <- new_x[still]
    xy_all[act[still], i, 2L] <- new_y[still]
    xy_all[act[still], i, 3L] <- mu[still]

    ## Detection check — record first entry within det.range of end receiver.
    ## Simulations continue running to N steps after detection.
    step_dists <- sqrt((new_x[still] - mpar$end[1])^2 +
                       (new_y[still] - mpar$end[2])^2)
    just_det <- !is.na(step_dists) & step_dists <= mpar$det.range &
                !detected[act[still]]   ## only record first detection
    if (any(just_det)) {
      det_sims           <- act[still][just_det]
      detected[det_sims] <- TRUE
      det_step[det_sims] <- i
      end_dist[det_sims] <- step_dists[just_det]
      ## active[det_sims] is NOT set FALSE — simulations run to N
    }

    if (pb) setTxtProgressBar(tpb, i)
  }

  if (pb) close(tpb)

  ## ---- Final distance for simulations that ran to completion ---------------
  completed <- !detected & !early_stop_vec
  if (any(completed)) {
    end_dist[completed] <- sqrt(
      (xy_all[completed, N, 1L] - mpar$end[1])^2 +
      (xy_all[completed, N, 2L] - mpar$end[2])^2
    )
  }
  accepted <- detected

  n_accepted <- sum(accepted)
  message(sprintf("%d / %d simulations accepted (%.1f%%)",
                  n_accepted, nsim, 100 * n_accepted / nsim))

  ## ---- Extract per-simulation tracks from xy_all ---------------------------
  tracks <- vector("list", nsim)
  for (sim_i in seq_len(nsim)) {
    n_valid <- if (early_stop_vec[sim_i]) {
      max(which(!is.na(xy_all[sim_i, , 1L])), 1L)
    } else {
      N
    }

    tracks[[sim_i]] <- tibble::tibble(
      id   = sim_i,
      date = date_seq[seq_len(n_valid)],
      x    = xy_all[sim_i, seq_len(n_valid), 1L],
      y    = xy_all[sim_i, seq_len(n_valid), 2L],
      u    = u_mat[sim_i,  seq_len(n_valid)],
      v    = v_mat[sim_i,  seq_len(n_valid)]
    )
  }

  ## ---- Detection locations --------------------------------------------------
  acc_idx  <- which(accepted)
  det_locs <- if (length(acc_idx) > 0L) {
    tibble::tibble(
      id   = acc_idx,
      x    = vapply(acc_idx, function(i) tracks[[i]]$x[   det_step[i]], numeric(1)),
      y    = vapply(acc_idx, function(i) tracks[[i]]$y[   det_step[i]], numeric(1)),
      date = do.call(c, lapply(acc_idx, function(i) tracks[[i]]$date[det_step[i]]))
    )
  } else {
    tibble::tibble(id = integer(0), x = numeric(0), y = numeric(0),
                  date = as.POSIXct(character(0), tz = "UTC"))
  }

  ## ---- Output ----------------------------------------------------------------
  structure(
    list(
      sims            = tracks,
      accepted        = accepted,
      end_dist        = end_dist,
      det_step        = det_step,
      det_locs        = det_locs,
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
    cat(sprintf("  End distance (all):  mean %.3f km  (range %.3f - %.3f km)\n",
                mean(ed), min(ed), max(ed)))
    if (x$n_accepted > 0) {
      ea <- x$end_dist[x$accepted]
      cat(sprintf("  End distance (acc):  mean %.3f km  (range %.3f - %.3f km)\n",
                  mean(ea), min(ea), max(ea)))
    }
  }
  invisible(x)
}
