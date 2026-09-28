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
#'       EVERY simulation runs to step \code{N} (or to the last valid step if
#'       it stopped early); detection does not truncate a track. An earlier
#'       version of this note said accepted tracks were cut at the detection
#'       step, which the code has not done for some time.}
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

  ## ---- RNG seed -------------------------------------------------------------
  ##
  ## When mpar$seed is set, the ensemble is exactly reproducible. The caller's
  ## RNG stream is restored on exit so that running a batch of passages gives
  ## the same result for each passage regardless of the order they are run in,
  ## and so that seeding one simulation does not silently reseed whatever the
  ## caller does next.
  if (!is.null(mpar$seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      .old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", .old_seed, envir = globalenv()), add = TRUE)
    } else {
      on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
    }
    set.seed(mpar$seed)
  }

  step_secs <- mpar$time.step * 60L
  ## ---- RUN PAST THE OBSERVED ARRIVAL, TO THE NEXT TURN OF THE TIDE --------
  ##
  ## `N` is the number of steps between the two real detections. Stopping there
  ## makes a LATE arrival unobservable: a fish delayed past N is recorded as
  ## never having crossed rather than as having crossed late, so the delay that
  ## swimming against the flood produces is converted into an absence.
  ##
  ## Measured on the 2023 fish 24 passage under tidal rheotaxis: truncating at
  ## N, the crossing rate falls from 1.000 at 0.25 body lengths per second to
  ## 0.520 at 4. Allowed to run on, it falls only to 0.780, and 26 % of the
  ## fast simulations cross AFTER N. Half of the fast fish were being thrown
  ## away, which is why the bl profile looked flat.
  ##
  ## The run therefore continues past N until each simulation next sees the
  ## tide turn -- the same phase test used for arming -- capped at
  ## `extend.max` steps. A fish that misses its tide gets the chance to come
  ## back on the following one, and no further.
  N_obs <- mpar$N
  ext   <- if (is.null(mpar$extend.max)) 0L else as.integer(mpar$extend.max)
  N     <- N_obs + ext
  turned_after  <- rep(FALSE, mpar$n_sim)   # tide has turned again since N_obs
  post_sign     <- rep(NA_real_, mpar$n_sim)
  post_run      <- rep(0L, mpar$n_sim)
  nsim <- mpar$n_sim

  ## ---- FVCOM setup (advect = TRUE only) -------------------------------------
  start.dt     <- mpar$start.dt
  advect_scale <- step_secs / 1000   ## m/s -> km/step
  layer_idx    <- NULL
  is_flood     <- FALSE              ## tidal phase: set from u at start location & time
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

    ## uses the EXTENDED N, so a run that would overrun the rasters fails here
    ## rather than part way through the extension
    sim_end_dt <- start.dt + N * step_secs
    if (sim_end_dt > data_end_dt)
      stop("Simulation end time (", format(sim_end_dt),
           ") extends beyond available FVCOM data (", format(data_end_dt), ").\n",
           "  Re-run sim_setup() with additional months.")

    fvcom.idx <- as.integer(round(
      as.numeric(difftime(start.dt, fvcom.origin, units = "secs")) / step_secs
    ))
    layer_idx <- seq_len(N - 1L) + fvcom.idx

    ## Determine tidal phase (flood vs ebb) once at the simulation start time.
    ## u is extracted at the start receiver location using the first step's FVCOM
    ## layer; this single scalar governs uvm selection for the entire simulation.
    ##
    ## FLOOD = u > 0: water flowing EASTWARD through Minas Passage, filling
    ## Minas Basin. EBB = u <= 0: water flowing westward, out toward the Bay of
    ## Fundy.
    ##
    ## This test used to read u < 0 and was described as "flood = flow into the
    ## Bay of Fundy", which is self-contradictory: flow into the Bay of Fundy is
    ## the ebb. The consequence was that uvm[1:2], the pair named for the flood,
    ## was applied to westward flow and uvm[3:4] to eastward -- the multipliers
    ## were swapped on every run. The raster data agree with the corrected
    ## reading: at FORCE007 on 11 May 2022 eastward flow peaks at 2.75 m/s
    ## against 2.09 m/s westward, and a shorter, stronger flood is what a
    ## convergent macrotidal basin should show.
    ref_pos  <- matrix(mpar$start, nrow = 1L)
    u_start  <- terra::extract(data$u[[layer_idx[1L]]], ref_pos,
                               method = "simple")[1L, 1L]
    is_flood <- !is.na(u_start) && u_start > 0

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
  ## ---- second detection record: ANY receiver, not just the end one --------
  ##
  ## `det.xy` is every receiver in this passage's own deployment. Station
  ## numbers are NOT unique across deployments -- station 1 sits at three
  ## different positions in the Minas files -- so the caller must pass the set
  ## belonging to this passage's file_tag, not a global list.
  ##
  ## ARMING IS NOT OPTIONAL. A simulation starts uniformly within det.range of
  ## the START receiver, and in this array receivers are often closer together
  ## than det.range -- several pairs are 20 to 200 m apart against a 300 m
  ## range. So at step one almost every simulation is already inside some
  ## receiver's radius, and an unarmed "any receiver" test would accept
  ## everything at step one regardless of how the fish swims. A simulation
  ## therefore becomes eligible only once it has been farther than det.range
  ## from EVERY receiver at least once: the event being recorded is leaving the
  ## array and coming back, which is what the observation is.
  ##
  ## `det_step_end_armed` applies the same arming to the ORIGINAL end-receiver
  ## test, and exists to measure a contaminant rather than to replace anything:
  ## 8 % of kept passages have their two stations closer together than
  ## det.range, so those simulations begin inside the target and are accepted
  ## at step one whatever bl is. `det_step` keeps the unarmed definition so
  ## every result computed to date stays comparable.
  ## `any_x`/`any_y` are WHERE the simulation crossed the array, and
  ## `any_end_dist` is how far that is from the receiver the real fish was
  ## detected at. That distance is the continuous version of the binary
  ## detection record, and it is what movement scenarios can be compared on.
  ##
  ## Two things it is not. It has a floor: the real fish's own crossing point
  ## is known only to within det.range, so a perfect simulation still scores a
  ## few hundred metres. And it is the FIRST armed crossing, not the nearest --
  ## a fish that crosses, leaves, and crosses again nearer the receiver is
  ## scored on the first one, to match how det_step is defined.
  ## ---- ARMING IS A TIDAL-PHASE TEST, NOT A DISTANCE ONE -------------------
  ##
  ## A crossing counts only once the tide has TURNED since the simulation
  ## started. That is the same rule the passages themselves were selected on:
  ## a fish crosses, is carried away, the tide turns, and it is carried back.
  ## If the phase has not changed, the fish has not been anywhere and a
  ## recorded crossing is the fish milling beside the array on its way out.
  ##
  ## The previous version armed on leaving a 300 m radius, which is far too
  ## weak: receivers on this line sit 84 to 391 m apart, so a fish drifting a
  ## few hundred metres off its release point re-enters a neighbour's circle
  ## within two or three steps. 13 % of cells in the 2026-09-15 sweep had their
  ## median crossing inside the first tenth of the passage, which is that
  ## artefact and nothing else.
  ##
  ## Phase is the sign of the along-channel component u at the fish's OWN
  ## position: u < 0 is westward, the ebb. A turn must persist for
  ## `phase.min.steps` consecutive steps, so a single eddy or a numerical
  ## wobble at slack water cannot arm a simulation on its own.
  phase0             <- rep(NA_real_,    nsim)   # sign of u at the first step
  phase_run          <- rep(0L,          nsim)   # consecutive steps at the new sign
  armed              <- rep(FALSE,       nsim)   # TRUE once the tide has turned
  armed_step         <- rep(NA_integer_, nsim)
  pmin <- if (is.null(mpar$phase.min.steps)) 2L else as.integer(mpar$phase.min.steps)

  detected_any       <- rep(FALSE,       nsim)
  det_step_any       <- rep(NA_integer_, nsim)
  any_dist           <- rep(NA_real_,    nsim)
  any_end_dist       <- rep(NA_real_,    nsim)
  any_x              <- rep(NA_real_,    nsim)
  any_y              <- rep(NA_real_,    nsim)
  det_which_any      <- rep(NA_integer_, nsim)

  ## ---- LINE CROSSING, the primary record ----------------------------------
  ##
  ## Proximity to a receiver is not what the model can resolve. The residual
  ## velocity error puts a simulated particle of order a kilometre from truth
  ## over a passage of this length, so asking it to thread within 300 m of one
  ## of 24 points is a lottery: on the 2023 fish 24 passage, 99 % of the
  ## simulations that FAILED the 300 m test came within 2 km of a receiver on
  ## the way back, 66 % within 1 km, median closest approach 851 m, and their
  ## journeys were otherwise identical to the ones that passed.
  ##
  ## So the primary test is whether the step SEGMENT intersects the receiver
  ## polyline -- did the fish cross the line, and when. That is resolvable, and
  ## it is what the observation actually tells us.
  crossed_line       <- rep(FALSE,       nsim)
  line_step          <- rep(NA_integer_, nsim)
  line_x             <- rep(NA_real_,    nsim)
  line_y             <- rep(NA_real_,    nsim)
  line_end_dist      <- rep(NA_real_,    nsim)
  ## closest approach to the line, over the whole armed part of the track, so a
  ## near miss is measurable rather than just absent
  near_line          <- rep(Inf,         nsim)
  near_line_step     <- rep(NA_integer_, nsim)
  detected_end_armed <- rep(FALSE,       nsim)
  det_step_end_armed <- rep(NA_integer_, nsim)
  det_xy <- mpar$det.xy
  if (!is.null(det_xy)) {
    det_xy <- as.matrix(det_xy)
    if (ncol(det_xy) != 2L || !nrow(det_xy) || !all(is.finite(det_xy)))
      stop("det.xy must be a finite two-column matrix of receiver positions")
  }

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
  ##
  ## Every movement model needs one, because the step-1 heading is the value
  ## the wrapped Cauchy is centred on (crw) or the documented fallback when the
  ## model's own rule cannot produce a bearing (the rheotaxis models, at slack
  ## water). Leaving a model out of this switch used to drop it through to
  ## NA_real_, and an NA heading propagates silently: rwrpcauchy() returns NA,
  ## the proposed position is NA, and the simulation is recorded as an early
  ## stop rather than an error. That would have depressed the acceptance rate
  ## of the rheotaxis models only, which is exactly the comparison they exist
  ## to support, so the fall-through now stops rather than guessing.
  init_heads <- switch(mpar$move,
    crw        = runif(nsim, 0, 2 * pi),
    bcrw       = rep(mpar$bearing, nsim),
    bcrw.coa   = rep(init_heading_coa, nsim),
    crw.bridge = atan2(sim_end_x - s_x, sim_end_y - s_y),
    ## Fixed bearing, so the previous heading is never consulted; the value is
    ## set here only so that the stored heading column is never NA.
    west       = rep(atan2(-1, 0), nsim),
    ## The rheotaxis models take their bearing from the current at each step
    ## and fall back on the previous heading only where the flow is too weak to
    ## give one. At step 1 there is no previous heading, so an arbitrary one is
    ## the honest choice, drawn the same way crw draws its first heading.
    rheo.pos   = runif(nsim, 0, 2 * pi),
    rheo.neg   = runif(nsim, 0, 2 * pi),
    rheo.tidal = runif(nsim, 0, 2 * pi),
    stop("no initial heading defined for move = '", mpar$move, "'")
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

    ## ---- Local current, extracted BEFORE the heading is chosen --------------
    ##
    ## The rheotaxis movement models need to know which way the water is going
    ## at each fish's own position before deciding where that fish swims, so
    ## the extraction is hoisted above the heading switch. It used to sit
    ## inside the advection block below. Two extract calls per step either way,
    ## regardless of the number of simulations.
    needs_flow <- mpar$move %in% c("rheo.pos", "rheo.neg", "rheo.tidal")
    u_adj <- v_adj <- NULL

    if (mpar$advect || needs_flow) {
      k       <- layer_idx[i - 1L]
      flood_i <- is_flood            ## scalar: determined once at simulation start time
      pos_m   <- cbind(px, py)

      if (w == 0) {
        u_raw <- terra::extract(data$u[[k]],      pos_m, method = "simple")[, 1]
        v_raw <- terra::extract(data$v[[k]],      pos_m, method = "simple")[, 1]
      } else {
        u_raw <- ((1 - w) * terra::extract(data$u[[k]],        pos_m, method = "simple")[, 1] +
                        w  * terra::extract(data$u[[k + 1L]],  pos_m, method = "simple")[, 1])
        v_raw <- ((1 - w) * terra::extract(data$v[[k]],        pos_m, method = "simple")[, 1] +
                        w  * terra::extract(data$v[[k + 1L]],  pos_m, method = "simple")[, 1])
      }

      ## uvm convention: c(u.flood, v.flood, u.ebb, v.ebb), flood = u > 0.
      ##
      ## phase = "track" (the default and the long-standing behaviour): the
      ## whole simulation uses the phase decided once at start.dt. Simple, but
      ## wrong for any passage that spans a slack -- between 31 and 67 per cent
      ## of them, depending on species -- because those get flood multipliers
      ## applied through an ebb, or the reverse.
      ##
      ## phase = "step": the phase is decided fresh at every step, from the
      ## sign of u at each fish's own position. More nearly correct, but it
      ## changes what uvm MEANS: if uvm was calibrated against drifters under
      ## the whole-track assumption, the fitted values absorbed that assumption
      ## and should be refitted before this is trusted.
      fl_i <- if (identical(mpar$phase, "step")) {
        !is.na(u_raw) & u_raw > 0
      } else {
        rep(isTRUE(flood_i), length(u_raw))
      }

      u_adj <- ifelse(!is.na(u_raw),
                      u_raw * ifelse(fl_i, mpar$uvm[1], mpar$uvm[3]),
                      0) * advect_scale
      v_adj <- ifelse(!is.na(v_raw),
                      v_raw * ifelse(fl_i, mpar$uvm[2], mpar$uvm[4]),
                      0) * advect_scale
    }

    ## ---- Mean heading for each active simulation ----------------------------
    ##
    ## Headings are compass bearings in radians: 0 is north, and a step of
    ## length s is taken as (x + sin(mu) * s, y + cos(mu) * s). So the bearing
    ## a current vector (u, v) is flowing TOWARD is atan2(u, v), and west is
    ## atan2(-1, 0) = -pi/2.
    ##
    ## The four models added below all orient on the UVM-CORRECTED current
    ## (u_adj, v_adj) rather than the raw field, because that is the model's
    ## best estimate of the flow the fish is actually in, and because the two
    ## components carry different multipliers, so the corrected vector points
    ## in a slightly different direction from the raw one.
    phi <- switch(mpar$move,
      bcrw = rep(
        if (length(mpar$bearing) > 1L) mpar$bearing[i] else mpar$bearing,
        n_act),
      crw = ph,

      ## Swim westward, out through Minas Passage toward the Bay of Fundy,
      ## whatever the tide is doing. A fixed bearing, so this is bcrw with the
      ## bearing supplied rather than asked for.
      west = rep(atan2(-1, 0), n_act),

      ## Positive rheotaxis: head into the flow, on the reciprocal of the
      ## bearing the water is travelling toward.
      rheo.pos = {
        if (is.null(u_adj)) stop("rheo.pos requires a current field; ",
                                 "run sim_setup() and leave advect = TRUE")
        h <- atan2(-u_adj, -v_adj)
        ## Slack water has no direction to orient to. Keep the previous
        ## heading rather than drawing an arbitrary one.
        h[!is.finite(h) | (u_adj == 0 & v_adj == 0)] <-
          ph[!is.finite(h) | (u_adj == 0 & v_adj == 0)]
        h
      },

      ## Negative rheotaxis: head downstream, with the flow.
      rheo.neg = {
        if (is.null(u_adj)) stop("rheo.neg requires a current field; ",
                                 "run sim_setup() and leave advect = TRUE")
        h <- atan2(u_adj, v_adj)
        h[!is.finite(h) | (u_adj == 0 & v_adj == 0)] <-
          ph[!is.finite(h) | (u_adj == 0 & v_adj == 0)]
        h
      },

      ## Tidal rheotaxis: go with the ebb, hold against the flood.
      ##
      ## Phase is decided HERE from the sign of u at this fish, at this step,
      ## not from the whole-track is_flood scalar. is_flood is fixed for an
      ## entire simulation, and between 31 and 67 per cent of passages span a
      ## slack, so using it would apply the wrong rule for part of many tracks.
      ##
      ## Stated in terms of DIRECTION so it cannot be knocked over by a
      ## labelling change: u < 0 is westward flow, out of Minas Basin toward
      ## the Bay of Fundy, which is the ebb. Westward flow is ridden, eastward
      ## (flood) flow is opposed.
      rheo.tidal = {
        if (is.null(u_adj)) stop("rheo.tidal requires a current field; ",
                                 "run sim_setup() and leave advect = TRUE")
        with_flow <- u_adj < 0                      # westward: ride it
        h <- ifelse(with_flow, atan2(u_adj, v_adj), atan2(-u_adj, -v_adj))
        h[!is.finite(h) | (u_adj == 0 & v_adj == 0)] <-
          ph[!is.finite(h) | (u_adj == 0 & v_adj == 0)]
        h
      },

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

    ## ---- Apply the advection computed above ---------------------------------
    if (mpar$advect) {
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

    ## ---- the any-receiver record, and the armed end-receiver record --------
    if (!is.null(det_xy)) {
      ix <- act[still]
      sx <- new_x[still]; sy <- new_y[still]
      dmin <- rep(Inf, length(sx)); wmin <- rep(NA_integer_, length(sx))
      for (r in seq_len(nrow(det_xy))) {
        dd <- sqrt((sx - det_xy[r, 1L])^2 + (sy - det_xy[r, 2L])^2)
        up <- !is.na(dd) & dd < dmin
        dmin[up] <- dd[up]; wmin[up] <- r
      }
      ## ---- arm on a sustained change of tidal phase ------------------------
      if (mpar$advect) {
        sg <- sign(u_adj[still])
        first <- is.na(phase0[ix]) & sg != 0
        if (any(first)) phase0[ix[first]] <- sg[first]
        flipped <- !is.na(phase0[ix]) & sg != 0 & sg != phase0[ix]
        phase_run[ix] <- ifelse(flipped, phase_run[ix] + 1L, 0L)
        new_arm <- !armed[ix] & phase_run[ix] >= pmin
        if (any(new_arm)) {
          armed[ix[new_arm]]      <- TRUE
          armed_step[ix[new_arm]] <- i
        }
      } else {
        ## with advection off there is no tide to turn; fall back to the old
        ## distance rule so the machinery still runs for a currents-off test
        new_arm <- !armed[ix] & is.finite(dmin) & dmin > mpar$det.range
        if (any(new_arm)) {
          armed[ix[new_arm]]      <- TRUE
          armed_step[ix[new_arm]] <- i
        }
      }
      hit <- armed[ix] & is.finite(dmin) & dmin <= mpar$det.range & !detected_any[ix]
      if (any(hit)) {
        h <- ix[hit]
        detected_any[h]  <- TRUE
        det_step_any[h]  <- i
        any_dist[h]      <- dmin[hit]
        det_which_any[h] <- wmin[hit]
        any_x[h]         <- sx[hit]
        any_y[h]         <- sy[hit]
        ## step_dists is already the distance to the observed end receiver for
        ## exactly these simulations, so this is free.
        any_end_dist[h]  <- step_dists[hit]
      }
      hit_e <- armed[ix] & !is.na(step_dists) & step_dists <= mpar$det.range &
               !detected_end_armed[ix]
      if (any(hit_e)) {
        detected_end_armed[ix[hit_e]] <- TRUE
        det_step_end_armed[ix[hit_e]] <- i
      }

      ## ---- did the step segment cross the receiver line? -------------------
      if (nrow(det_xy) > 1L && i > 1L) {
        px0 <- xy_all[ix, i - 1L, 1L]; py0 <- xy_all[ix, i - 1L, 2L]
        ok  <- is.finite(px0) & is.finite(py0) & is.finite(sx) & is.finite(sy)
        if (any(ok)) {
          bestt <- rep(NA_real_, length(sx))
          for (r in seq_len(nrow(det_xy) - 1L)) {
            ax <- det_xy[r, 1L];     ay <- det_xy[r, 2L]
            bx <- det_xy[r + 1L, 1L]; by <- det_xy[r + 1L, 2L]
            rx <- sx - px0; ry <- sy - py0
            ssx <- bx - ax;  ssy <- by - ay
            den <- rx * ssy - ry * ssx
            g   <- ok & abs(den) > 1e-12
            tt  <- rep(NA_real_, length(sx)); uu <- tt
            tt[g] <- ((ax - px0[g]) * ssy - (ay - py0[g]) * ssx) / den[g]
            uu[g] <- ((ax - px0[g]) * ry[g] - (ay - py0[g]) * rx[g]) / den[g]
            hitseg <- g & !is.na(tt) & tt >= 0 & tt <= 1 & uu >= 0 & uu <= 1
            upd <- hitseg & (is.na(bestt) | tt < bestt)
            bestt[upd] <- tt[upd]
          }
          xed <- armed[ix] & !crossed_line[ix] & !is.na(bestt)
          if (any(xed)) {
            h  <- ix[xed]; tt <- bestt[xed]
            cx <- px0[xed] + tt * (sx[xed] - px0[xed])
            cy <- py0[xed] + tt * (sy[xed] - py0[xed])
            crossed_line[h]  <- TRUE
            line_step[h]     <- i
            line_x[h]        <- cx
            line_y[h]        <- cy
            line_end_dist[h] <- sqrt((cx - mpar$end[1])^2 + (cy - mpar$end[2])^2)
          }
        }
      }

      ## ---- past N_obs: stop each simulation at the next turn of the tide ---
      ## The extension exists to make a LATE arrival visible, not to run the
      ## fish indefinitely. One further turn of the tide is the fish's next
      ## opportunity to be carried back; after that it is genuinely gone.
      if (ext > 0L && i > N_obs && mpar$advect) {
        sg2  <- sign(u_adj[still])
        newp <- is.na(post_sign[ix]) & sg2 != 0
        if (any(newp)) post_sign[ix[newp]] <- sg2[newp]
        fl2  <- !is.na(post_sign[ix]) & sg2 != 0 & sg2 != post_sign[ix]
        post_run[ix] <- ifelse(fl2, post_run[ix] + 1L, 0L)
        done <- !turned_after[ix] & post_run[ix] >= pmin
        if (any(done)) {
          turned_after[ix[done]] <- TRUE
          active[ix[done]]       <- FALSE
        }
      }

      ## ---- closest approach to the line, armed part of the track only ------
      if (any(armed[ix])) {
        dl <- .dist_to_polyline(sx, sy, det_xy)
        upd <- armed[ix] & is.finite(dl) & dl < near_line[ix]
        if (any(upd)) {
          near_line[ix[upd]]      <- dl[upd]
          near_line_step[ix[upd]] <- i
        }
      }
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
      ## second and third records; all NA when det.xy was not supplied
      detected_any       = detected_any,
      det_step_any       = det_step_any,
      any_dist           = any_dist,
      any_end_dist       = any_end_dist,
      any_x              = any_x,
      any_y              = any_y,
      det_which_any      = det_which_any,
      armed              = armed,
      armed_step         = armed_step,
      phase0             = phase0,
      crossed_line       = crossed_line,
      line_step          = line_step,
      line_x             = line_x,
      line_y             = line_y,
      line_end_dist      = line_end_dist,
      near_line          = ifelse(is.finite(near_line), near_line, NA_real_),
      near_line_step     = near_line_step,
      N_obs              = N_obs,
      ## TRUE when the crossing happened after the fish was actually detected.
      ## Under the old truncated design these simulations were indistinguishable
      ## from ones that never crossed at all.
      crossed_late       = !is.na(line_step) & line_step > N_obs,
      turned_after       = turned_after,
      detected_end_armed = detected_end_armed,
      det_step_end_armed = det_step_end_armed,
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
  cat(sprintf("  Seed:                %s\n",
              if (is.null(p$seed)) "none (not reproducible)" else format(p$seed)))
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
