## =============================================================================
## fish_par() — parameter constructor for sim_fish()
## =============================================================================

##' @title Parameters for tagged fish simulation between detection events
##'
##' @description Constructs and validates parameters for \code{sim_fish()}.
##'   The simulation time span (number of steps, \code{N}) is derived
##'   automatically from the elapsed time between \code{start.dt} and
##'   \code{end.dt}; \code{N} is not a user parameter.
##'
##' @param n_sim     Number of independent simulations (positive integer).
##'   Default 100.
##' @param time.step Step duration in minutes (must match FVCOM output
##'   interval). Default 10.
##' @param start.dt  Precise detection time at the start receiver
##'   (POSIXct, UTC). Required.
##' @param end.dt    Precise detection time at the end receiver
##'   (POSIXct, UTC). Required.
##' @param start     Start receiver coordinates, \code{c(x, y)} in km
##'   (UTM Zone 20N). Required.
##' @param end       End receiver coordinates, \code{c(x, y)} in km.
##'   Required.
##' @param det.range Acoustic detection range in km. Used both to sample
##'   start positions (uniformly within a circle of this radius around
##'   \code{start}) and as the acceptance radius around \code{end} at
##'   \code{end.dt}. Default 0.3 km (300 m).
##' @param move      Movement model controlling active swimming behaviour.
##'   \code{"crw"} — correlated random walk with no directional bias; each
##'   heading is drawn from a Wrapped Cauchy distribution centred on the
##'   \emph{previous} heading, so persistence is governed solely by
##'   \code{rho}. \code{"bcrw"} (default) — biased correlated RW toward a
##'   fixed \code{bearing}. \code{"bcrw.coa"} — biased correlated RW
##'   toward an intermediate Centre of Attraction (\code{coa}).
##'   \code{"crw.bridge"} — Brownian bridge CRW; each simulation samples
##'   both a start position (within \code{det.range} of \code{start}) and
##'   an end position (within \code{det.range} of \code{end}), and
##'   progressively biases heading toward the sampled end as
##'   \code{nu * i / max(N - i, 1)}. Acceptance is still evaluated against
##'   the end \emph{receiver} location, not the sampled end point.
##' @param bearing   For \code{move = "bcrw"}: swimming bias bearing in
##'   radians. Scalar only. Default 0 (north). Ignored for \code{"crw"}
##'   and \code{"bcrw.coa"}.
##' @param coa       For \code{move = "bcrw.coa"}: intermediate Centre of
##'   Attraction, \code{c(x, y)} in km. This is a swimming-bias waypoint,
##'   not the end receiver — the end constraint is handled separately by
##'   rejection sampling. Required when \code{move = "bcrw.coa"};
##'   ignored otherwise.
##' @param nu        Attraction strength (>= 0). For
##'   \code{move = "bcrw.coa"}: constant pull toward \code{coa}; larger
##'   values = stronger attraction. For \code{move = "crw.bridge"}: base
##'   attraction toward the per-simulation sampled end position, scaled
##'   as \code{nu * i / max(N - i, 1)} so the pull grows from near-zero
##'   at the first step to very strong at the final step. Default 1.
##'   Ignored for \code{"crw"} and \code{"bcrw"}.
##' @param rho       Angular concentration (Wrapped Cauchy): 0 = isotropic,
##'   approaching 1 = highly directed. Range [0, 1). Default 0.5.
##' @param bl        Swimming speed in body-lengths per second. Default 1.5.
##' @param fl        Fork length in metres, used to convert \code{bl} to
##'   km/step. Default 0.15.
##' @param noise     Optional positional noise: sd in km added to each
##'   position. \code{NULL} (default) = off. Currently applied for
##'   \code{move = "drifter"} only; included here for future use.
##' @param beta      Land-avoidance potential coefficients, both <= 0.
##'   More negative = stronger repulsion. Default \code{c(-0.15, -0.15)}.
##' @param buffer    Water-search radius (km) when fish grounds on land.
##'   Default 1.
##' @param uvm       Current multiplier. Scalar, 2-element \code{c(u, v)},
##'   or 4-element \code{c(u_neg, u_pos, v_neg, v_pos)} applying separate
##'   multipliers for negative and positive \code{u} current. Default
##'   \code{c(1, 1)}.
##' @param advect    Logical; advect by FVCOM currents. Default \code{TRUE}.
##' @param interp    Logical; if \code{TRUE}, linearly interpolate u and v
##'   between the two FVCOM layers bracketing \code{start.dt} (2 extracts
##'   per step). Default \code{FALSE} (snap to nearest layer).
##' @param method    Simulation method. \code{"rejection"} (default): accept
##'   simulations whose final position falls within \code{det.range} of
##'   \code{end}. Future option \code{"bridge"} will add a temporal
##'   attraction toward \code{end} before filtering.
##'
##' @return Object of class \code{c("fish_par", "mpar")}
##'
##' @examples
##' \dontrun{
##' mpar <- fish_par(
##'   n_sim    = 200,
##'   start.dt = as.POSIXct("2022-07-14 06:00:00", tz = "UTC"),
##'   end.dt   = as.POSIXct("2022-07-14 08:30:00", tz = "UTC"),
##'   start    = c(512.3, 4948.7),
##'   end      = c(487.1, 4942.2),
##'   bearing  = pi    ## southward bias
##' )
##' }
##'
##' @export
fish_par <- function(
    n_sim     = 100,
    time.step = 10,
    start.dt,
    end.dt,
    start     = c(0, 0),
    end       = c(0, 0),
    det.range = 0.3,
    move      = "bcrw",
    bearing   = 0,
    coa       = NULL,
    nu        = 1,
    rho       = 0.5,
    bl        = 1.5,
    fl        = 0.15,
    noise     = NULL,
    beta      = c(-0.15, -0.15),
    buffer    = 1,
    uvm       = c(1, 1),
    advect    = TRUE,
    interp    = FALSE,
    method    = "rejection"
) {

  ## ---- Validation ------------------------------------------------------------
  .check_positive_int(n_sim, "n_sim")

  if (!is.numeric(time.step) || length(time.step) != 1 || time.step <= 0)
    stop("time.step must be a positive number (minutes), e.g. 10")

  .check_posixct(start.dt, "start.dt")
  .check_posixct(end.dt,   "end.dt")
  if (end.dt <= start.dt)
    stop("end.dt must be later than start.dt")

  .check_xy(start, "start")
  .check_xy(end,   "end")

  if (!is.numeric(det.range) || length(det.range) != 1 || det.range <= 0)
    stop("det.range must be a positive scalar (km)")

  move   <- match.arg(move,   c("crw", "bcrw", "bcrw.coa", "crw.bridge"))
  method <- match.arg(method, "rejection")   ## extend to c("rejection", "bridge") for option 3

  if (move == "crw") {
    ## No directional bias: bearing, coa, and nu are not used
    bearing <- NULL
    coa     <- NULL
    nu      <- NULL
  } else if (move == "bcrw") {
    if (!is.numeric(bearing) || length(bearing) != 1)
      stop("bearing must be a numeric scalar (radians)")
    coa <- NULL
    nu  <- NULL
  } else if (move == "bcrw.coa") {
    if (is.null(coa) || !is.numeric(coa) || length(coa) != 2)
      stop("coa must be c(x, y) in km when move = 'bcrw.coa'")
    if (!is.numeric(nu) || length(nu) != 1 || nu < 0)
      stop("nu must be a non-negative scalar")
    bearing <- NULL
  } else {
    ## crw.bridge: time-varying attraction toward per-simulation end position
    if (!is.numeric(nu) || length(nu) != 1 || nu < 0)
      stop("nu must be a non-negative scalar for move = 'crw.bridge'")
    bearing <- NULL
    coa     <- NULL
  }

  .check_scalar_range(rho, "rho", 0, 0.9999)

  if (!is.numeric(bl) || length(bl) != 1 || bl <= 0)
    stop("bl must be a positive scalar (body-lengths/sec)")
  if (!is.numeric(fl) || length(fl) != 1 || fl <= 0)
    stop("fl must be a positive scalar (fork length in metres)")

  .check_beta(beta)

  if (!is.numeric(buffer) || length(buffer) != 1 || buffer <= 0)
    stop("buffer must be a positive scalar (km)")

  .check_uvm(uvm)

  if (!is.logical(advect) || length(advect) != 1)
    stop("advect must be TRUE or FALSE")

  if (!is.null(noise) && (!is.numeric(noise) || length(noise) != 1 || noise < 0))
    stop("noise must be NULL or a non-negative scalar (sd in km)")

  if (!is.logical(interp) || length(interp) != 1)
    stop("interp must be TRUE or FALSE")

  ## ---- Derive N from detection times -----------------------------------------
  step_secs <- time.step * 60
  N <- as.integer(round(
    as.numeric(difftime(end.dt, start.dt, units = "secs")) / step_secs
  ))
  if (N < 2)
    stop("Time span between start.dt and end.dt covers fewer than 2 steps.\n",
         "  span = ", format(round(as.numeric(difftime(end.dt, start.dt, units = "mins")), 1)),
         " min,  time.step = ", time.step, " min")

  ## Expand uvm to 4-element form (u_neg, u_pos, v_neg, v_pos)
  uvm <- switch(as.character(length(uvm)),
    "1" = rep(uvm, 4),
    "2" = c(uvm[1], uvm[1], uvm[2], uvm[2]),
    "4" = uvm,
    stop("uvm must have 1, 2, or 4 elements")
  )

  ## ---- Assemble parameter list -----------------------------------------------
  structure(
    list(
      n_sim     = n_sim,
      N         = N,
      time.step = time.step,
      start.dt  = start.dt,
      end.dt    = end.dt,
      start     = start,
      end       = end,
      det.range = det.range,
      move      = move,
      method    = method,
      bearing   = bearing,
      coa       = coa,
      nu        = nu,
      rho       = rho,
      bl        = bl,
      fl        = fl,
      noise     = noise,
      beta      = beta,
      buffer    = buffer,
      uvm       = uvm,
      advect    = advect,
      interp    = interp,
      land      = FALSE,
      boundary  = FALSE
    ),
    class = c("fish_par", "mpar")
  )
}


## -----------------------------------------------------------------------------
## Print method
## -----------------------------------------------------------------------------

#' @exportS3Method print fish_par
print.fish_par <- function(x, ...) {
  cat("-- fish_par (tagged fish simulation) --\n")
  cat(sprintf("  n_sim:     %d simulations   method: %s\n",
              x$n_sim, x$method))
  cat(sprintf("  N:         %d steps  (%.1f hours at %g-min intervals)\n",
              x$N, x$N * x$time.step / 60, x$time.step))
  cat(sprintf("  start.dt:  %s\n", format(x$start.dt)))
  cat(sprintf("  end.dt:    %s\n", format(x$end.dt)))
  cat(sprintf("  start:     x = %.3f km,  y = %.3f km\n", x$start[1], x$start[2]))
  cat(sprintf("  end:       x = %.3f km,  y = %.3f km\n", x$end[1],   x$end[2]))
  cat(sprintf("  det.range: %.3f km\n", x$det.range))
  cat(sprintf("  move:      %s\n", x$move))
  if (x$move == "crw") {
    cat("  (no directional bias — heading self-correlated via rho)\n")
  } else if (x$move == "bcrw") {
    cat(sprintf("  bearing:   %.4f rad\n", x$bearing))
  } else if (x$move == "bcrw.coa") {
    cat(sprintf("  coa:       x = %.3f km,  y = %.3f km   nu: %.2f\n",
                x$coa[1], x$coa[2], x$nu))
  } else {
    ## crw.bridge
    cat(sprintf("  nu (base): %.2f  (attraction scales as nu * i / max(N - i, 1))\n",
                x$nu))
  }
  cat(sprintf("  rho:       %.3f   bl: %.2f BL/s   fl: %.3f m\n",
              x$rho, x$bl, x$fl))
  cat(sprintf("  advect:    %s    uvm: c(%.2f, %.2f, %.2f, %.2f)    interp: %s\n",
              x$advect, x$uvm[1], x$uvm[2], x$uvm[3], x$uvm[4], x$interp))
  cat(sprintf("  beta:      c(%.3f, %.3f)   buffer: %.2f km\n",
              x$beta[1], x$beta[2], x$buffer))
  invisible(x)
}
