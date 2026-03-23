## =============================================================================
## Parameter constructors for simTidal
##
## Replaces sim_par() with three move-type-specific validated constructors:
##
##   drifter_par()   -- passive advection, move = "drifter"
##   bcrw_par()      -- biased correlated RW toward a bearing, move = "bcrw"
##   bcrw_coa_par()  -- biased correlated RW toward a CoA, move = "bcrw.coa"
##
## Each constructor:
##   - accepts only the parameters relevant to that move type
##   - validates types and ranges at construction time
##   - returns a classed list (S3) with a print method
##   - passes validate_mpar() check at the top of sim_drifter()
##
## Backward compatibility: sim_par() is retained as a thin wrapper that routes
## to the appropriate constructor based on the move argument.
## =============================================================================


## -----------------------------------------------------------------------------
## Internal validators — called by all constructors
## -----------------------------------------------------------------------------

.check_positive_int <- function(x, name) {
  if (!is.null(x) && (!is.numeric(x) || length(x) != 1 || x != round(x) || x <= 0))
    stop(name, " must be a positive integer, got: ", deparse(x))
}

.check_posixct <- function(x, name) {
  if (!inherits(x, "POSIXct"))
    stop(name, " must be a POSIXct datetime, got ", class(x)[1],
         "\n  Use: as.POSIXct('YYYY-MM-DD HH:MM:SS', tz = 'UTC')")
}

.check_xy <- function(x, name) {
  if (!is.numeric(x) || length(x) != 2)
    stop(name, " must be a numeric vector of length 2 (x, y in km)")
}

.check_scalar_range <- function(x, name, lo, hi) {
  if (!is.numeric(x) || length(x) != 1 || x < lo || x > hi)
    stop(name, " must be a single number in [", lo, ", ", hi, "], got: ", deparse(x))
}

.check_beta <- function(x) {
  if (!is.numeric(x) || length(x) != 2 || any(x > 0))
    stop("beta must be a numeric vector of length 2 with both values <= 0,\n",
         "  e.g. c(-0.15, -0.15). More negative = stronger land avoidance.")
}

.check_uvm <- function(x) {
  if (!is.numeric(x) || !length(x) %in% c(1L, 2L))
    stop("uvm must be a numeric scalar or 2-element vector")
  if (any(x < 0))
    stop("uvm values must be >= 0")
}

## Shared land/advection fields validated and returned by all constructors
.common_fields <- function(time.step, start.dt, start,
                            beta, buffer, uvm, advect, noise, interp) {
  if (!is.numeric(time.step) || length(time.step) != 1 || time.step <= 0)
    stop("time.step must be a positive number (minutes), e.g. 10")
  if (!is.null(start.dt)) .check_posixct(start.dt, "start.dt")
  .check_xy(start, "start")
  .check_beta(beta)
  if (!is.numeric(buffer) || length(buffer) != 1 || buffer <= 0)
    stop("buffer must be a positive number (km)")
  .check_uvm(uvm)
  if (!is.logical(advect) || length(advect) != 1)
    stop("advect must be TRUE or FALSE")
  if (!is.null(noise) && (!is.numeric(noise) || length(noise) != 1 || noise < 0))
    stop("noise must be NULL or a non-negative number (sd in km)")
  if (!is.logical(interp) || length(interp) != 1)
    stop("interp must be TRUE or FALSE")

  list(
    time.step = time.step,
    start.dt  = start.dt,
    start     = start,
    beta      = beta,
    buffer    = buffer,
    uvm       = if (length(uvm) == 1) c(uvm, uvm) else uvm,
    advect    = advect,
    noise     = noise,
    interp    = interp,
    land      = FALSE,
    boundary  = FALSE
  )
}


## -----------------------------------------------------------------------------
##' @title Parameters for passive drifter simulation
##'
##' @description Construct and validate parameters for \code{sim_drifter()}
##'   with \code{move = "drifter"}. The drifter is advected entirely by FVCOM
##'   currents; optional Gaussian noise can be added to each position.
##'
##' @param N          Number of simulation time steps (required, positive integer)
##' @param time.step  Duration of each step in minutes (default 10, must match
##'                   FVCOM output interval)
##' @param start.dt   Simulation start datetime (POSIXct, UTC). Defaults to
##'                   \code{NULL}, resolved to the first FVCOM timestamp
##'                   (\code{data$fvcom.origin}) inside \code{sim_drifter()}
##' @param start      Start coordinates, c(x, y) in km (UTM Zone 20N)
##' @param noise      Optional positional noise: sd in km added to each step.
##'                   \code{NULL} (default) = purely advection-driven
##' @param beta       Land-avoidance potential function coefficients, both <= 0.
##'                   More negative = stronger repulsion from land.
##'                   Default \code{c(-0.15, -0.15)}
##' @param buffer     Search radius (km) for water when drifter grounds on land
##' @param uvm        Current multiplier: scalar or c(u_mult, v_mult).
##'                   Default \code{c(1, 1)} = unmodified FVCOM currents
##' @param advect     Logical; advect by FVCOM currents (default \code{TRUE})
##' @param interp     Logical; if \code{TRUE}, linearly interpolate u and v
##'                   between the two FVCOM layers bracketing \code{start.dt}
##'                   (2 extracts per step). If \code{FALSE} (default), snap
##'                   \code{start.dt} to the nearest layer boundary and use
##'                   a single extract per step.
##'
##' @return Object of class \code{c("drifter_par", "mpar")}
##'
##' @examples
##' ## Snap to nearest boundary (default — fastest)
##' mpar <- drifter_par(N = 1440, start = c(512.3, 214.8))
##'
##' ## Linear interpolation between layers
##' mpar <- drifter_par(N = 1440, start = c(512.3, 214.8), interp = TRUE)
##'
##' @export
## -----------------------------------------------------------------------------
drifter_par <- function(
    N,
    time.step = 10,
    start.dt  = NULL,
    start     = c(0, 0),
    noise     = NULL,
    beta      = c(-0.15, -0.15),
    buffer    = 1,
    uvm       = c(1, 1),
    advect    = TRUE,
    interp    = FALSE
) {
  .check_positive_int(N, "N")

  common <- .common_fields(time.step, start.dt, start,
                            beta, buffer, uvm, advect, noise, interp)

  structure(
    c(list(N = N, move = "drifter"), common),
    class = c("drifter_par", "mpar")
  )
}


## -----------------------------------------------------------------------------
##' @title Parameters for biased correlated random walk toward a bearing
##'
##' @description Construct and validate parameters for \code{sim_drifter()}
##'   with \code{move = "bcrw"}. Swimming is biased toward a fixed or
##'   time-varying bearing, with FVCOM advection overlaid.
##'
##' @param N          Number of simulation time steps (required)
##' @param time.step  Step duration in minutes (default 10)
##' @param start.dt   Simulation start datetime (POSIXct, UTC). Defaults to
##'                   \code{NULL}, resolved to the first FVCOM timestamp
##'                   (\code{data$fvcom.origin}) inside \code{sim_drifter()}
##' @param start      Start coordinates, c(x, y) in km
##' @param bearing    Bearing(s) in radians. Scalar for constant direction;
##'                   numeric vector of length N for time-varying
##' @param rho        Angular concentration (Wrapped Cauchy): 0 = isotropic,
##'                   approaching 1 = highly directed. Range [0, 1)
##' @param bl         Swimming speed in body-lengths per second. Scalar or
##'                   vector of length N
##' @param fl         Fork length in metres (used to convert bl to km/step)
##' @param noise      Optional positional noise sd in km. \code{NULL} = off
##' @param beta       Land-avoidance coefficients, both <= 0
##' @param buffer     Water-search radius (km) when grounded
##' @param uvm        Current multiplier: scalar or c(u_mult, v_mult)
##' @param advect     Logical; advect by FVCOM currents (default \code{TRUE})
##' @param interp     Logical; if \code{TRUE}, linearly interpolate u and v
##'                   between the two FVCOM layers bracketing \code{start.dt}.
##'                   Default \code{FALSE}.
##'
##' @return Object of class \code{c("bcrw_par", "mpar")}
##' @export
## -----------------------------------------------------------------------------
bcrw_par <- function(
    N,
    time.step = 10,
    start.dt  = NULL,
    start     = c(0, 0),
    bearing   = 0,
    rho       = 0.5,
    bl        = 1.5,
    fl        = 0.15,
    noise     = NULL,
    beta      = c(-0.15, -0.15),
    buffer    = 1,
    uvm       = c(1, 1),
    advect    = TRUE,
    interp    = FALSE
) {
  .check_positive_int(N, "N")

  if (!is.numeric(bearing))
    stop("bearing must be numeric (radians)")
  if (length(bearing) > 1 && length(bearing) != N)
    stop("bearing must be a scalar or a vector of length N (", N, "), got length ",
         length(bearing))

  .check_scalar_range(rho, "rho", 0, 0.9999)

  if (!is.numeric(bl) || any(bl <= 0))
    stop("bl must be a positive numeric scalar or vector")
  if (length(bl) > 1 && length(bl) != N)
    stop("bl must be a scalar or a vector of length N (", N, "), got length ",
         length(bl))

  if (!is.numeric(fl) || length(fl) != 1 || fl <= 0)
    stop("fl must be a positive scalar (fork length in metres)")

  common <- .common_fields(time.step, start.dt, start,
                            beta, buffer, uvm, advect, noise, interp)

  structure(
    c(list(N = N, move = "bcrw",
           bearing = bearing, rho = rho, bl = bl, fl = fl),
      common),
    class = c("bcrw_par", "mpar")
  )
}


## -----------------------------------------------------------------------------
##' @title Parameters for biased correlated random walk toward a Centre of Attraction
##'
##' @description Construct and validate parameters for \code{sim_drifter()}
##'   with \code{move = "bcrw.coa"}. Swimming is biased toward one or more
##'   sequential Centres of Attraction, with FVCOM advection overlaid.
##'
##' @param N          Number of simulation time steps (required)
##' @param time.step  Step duration in minutes (default 10)
##' @param start.dt   Simulation start datetime (POSIXct, UTC). Defaults to
##'                   \code{NULL}, resolved to the first FVCOM timestamp
##'                   (\code{data$fvcom.origin}) inside \code{sim_drifter()}
##' @param start      Start coordinates, c(x, y) in km
##' @param coa        CoA coordinates. Either c(x, y) for a single CoA, or a
##'                   2-column matrix with one row per sequential CoA
##' @param coa.tol    Tolerance radius (km) around each CoA
##' @param nu         Strength of attraction toward CoA (>= 0)
##' @param rho        Angular concentration (Wrapped Cauchy), range [0, 1)
##' @param bl         Swimming speed in body-lengths per second
##' @param fl         Fork length in metres
##' @param noise      Optional positional noise sd in km. \code{NULL} = off
##' @param beta       Land-avoidance coefficients, both <= 0
##' @param buffer     Water-search radius (km) when grounded
##' @param uvm        Current multiplier: scalar or c(u_mult, v_mult)
##' @param advect     Logical; advect by FVCOM currents (default \code{TRUE})
##' @param interp     Logical; if \code{TRUE}, linearly interpolate u and v
##'                   between the two FVCOM layers bracketing \code{start.dt}.
##'                   Default \code{FALSE}.
##'
##' @return Object of class \code{c("bcrw_coa_par", "mpar")}
##' @export
## -----------------------------------------------------------------------------
bcrw_coa_par <- function(
    N,
    time.step = 10,
    start.dt  = NULL,
    start     = c(0, 0),
    coa       = cbind(c(0, 0)),
    coa.tol   = 1,
    nu        = 1,
    rho       = 0.5,
    bl        = 1.5,
    fl        = 0.15,
    noise     = NULL,
    beta      = c(-0.15, -0.15),
    buffer    = 1,
    uvm       = c(1, 1),
    advect    = TRUE,
    interp    = FALSE
) {
  .check_positive_int(N, "N")

  coa <- if (is.vector(coa) && length(coa) == 2) {
    matrix(coa, nrow = 1, ncol = 2)
  } else if (is.matrix(coa) && ncol(coa) == 2) {
    coa
  } else {
    stop("coa must be c(x, y) or a 2-column matrix of CoA coordinates")
  }

  if (!is.numeric(coa.tol) || length(coa.tol) != 1 || coa.tol <= 0)
    stop("coa.tol must be a positive scalar (km)")
  if (!is.numeric(nu) || length(nu) != 1 || nu < 0)
    stop("nu must be a non-negative scalar")
  .check_scalar_range(rho, "rho", 0, 0.9999)
  if (!is.numeric(bl) || any(bl <= 0))
    stop("bl must be a positive numeric scalar")
  if (!is.numeric(fl) || length(fl) != 1 || fl <= 0)
    stop("fl must be a positive scalar (fork length in metres)")

  common <- .common_fields(time.step, start.dt, start,
                            beta, buffer, uvm, advect, noise, interp)

  structure(
    c(list(N = N, move = "bcrw.coa",
           coa = coa, coa.tol = coa.tol,
           nu = nu, rho = rho, bl = bl, fl = fl),
      common),
    class = c("bcrw_coa_par", "mpar")
  )
}


## -----------------------------------------------------------------------------
## Print methods
## -----------------------------------------------------------------------------

#' @exportS3Method print drifter_par
print.drifter_par <- function(x, ...) {
  cat("-- drifter_par (passive advection) --\n")
  cat(sprintf("  N:         %d steps  (%.1f hours at %g-min intervals)\n",
              x$N, x$N * x$time.step / 60, x$time.step))
  cat(sprintf("  start.dt:  %s\n", if (is.null(x$start.dt)) "<= data$fvcom.origin>" else format(x$start.dt)))
  cat(sprintf("  start:     x = %.2f km,  y = %.2f km\n", x$start[1], x$start[2]))
  cat(sprintf("  advect:    %s    interp: %s\n", x$advect, x$interp))
  if (!is.null(x$noise))
    cat(sprintf("  noise sd:  %.4f km\n", x$noise))
  cat(sprintf("  uvm:       c(%.2f, %.2f)\n", x$uvm[1], x$uvm[2]))
  cat(sprintf("  beta:      c(%.3f, %.3f)  buffer: %.2f km\n",
              x$beta[1], x$beta[2], x$buffer))
  invisible(x)
}

#' @exportS3Method print bcrw_par
print.bcrw_par <- function(x, ...) {
  cat("-- bcrw_par (biased correlated RW toward bearing) --\n")
  cat(sprintf("  N:         %d steps  (%.1f hours at %g-min intervals)\n",
              x$N, x$N * x$time.step / 60, x$time.step))
  cat(sprintf("  start.dt:  %s\n", if (is.null(x$start.dt)) "<= data$fvcom.origin>" else format(x$start.dt)))
  cat(sprintf("  start:     x = %.2f km,  y = %.2f km\n", x$start[1], x$start[2]))
  cat(sprintf("  bearing:   %s\n",
              if (length(x$bearing) == 1) sprintf("%.4f rad", x$bearing)
              else sprintf("vector [length %d]", length(x$bearing))))
  cat(sprintf("  rho:       %.3f   bl: %s BL/s   fl: %.3f m\n",
              x$rho,
              if (length(x$bl) == 1) sprintf("%.2f", x$bl)
              else sprintf("vector [length %d]", length(x$bl)),
              x$fl))
  cat(sprintf("  advect:    %s    uvm: c(%.2f, %.2f)    interp: %s\n",
              x$advect, x$uvm[1], x$uvm[2], x$interp))
  cat(sprintf("  beta:      c(%.3f, %.3f)  buffer: %.2f km\n",
              x$beta[1], x$beta[2], x$buffer))
  invisible(x)
}

#' @exportS3Method print bcrw_coa_par
print.bcrw_coa_par <- function(x, ...) {
  cat("-- bcrw_coa_par (biased correlated RW toward CoA) --\n")
  cat(sprintf("  N:         %d steps  (%.1f hours at %g-min intervals)\n",
              x$N, x$N * x$time.step / 60, x$time.step))
  cat(sprintf("  start.dt:  %s\n", if (is.null(x$start.dt)) "<= data$fvcom.origin>" else format(x$start.dt)))
  cat(sprintf("  start:     x = %.2f km,  y = %.2f km\n", x$start[1], x$start[2]))
  cat(sprintf("  CoA(s):    %d point(s)   tolerance: %.2f km   nu: %.2f\n",
              nrow(x$coa), x$coa.tol, x$nu))
  cat(sprintf("  rho:       %.3f   bl: %.2f BL/s   fl: %.3f m\n",
              x$rho, x$bl, x$fl))
  cat(sprintf("  advect:    %s    uvm: c(%.2f, %.2f)    interp: %s\n",
              x$advect, x$uvm[1], x$uvm[2], x$interp))
  cat(sprintf("  beta:      c(%.3f, %.3f)  buffer: %.2f km\n",
              x$beta[1], x$beta[2], x$buffer))
  invisible(x)
}


## -----------------------------------------------------------------------------
## validate_mpar() — called at the top of sim_drifter() instead of ad-hoc checks
## -----------------------------------------------------------------------------

##' @title Validate an mpar object against a loaded data list
##'
##' @description Checks that the mpar object is the correct S3 class, that
##'   \code{move} matches the constructor used, and that \code{time.step}
##'   matches the FVCOM layer interval cached in \code{data$fvcom_step_secs}
##'   (computed by \code{sim_setup()} from u layer names).
##'   Called internally by \code{sim_drifter()}.
##'
##' @param mpar  An object produced by \code{drifter_par()}, \code{bcrw_par()},
##'              or \code{bcrw_coa_par()}
##' @param data  Output list from \code{sim_setup()}
##' @keywords internal
validate_mpar <- function(mpar, data) {

  if (!inherits(mpar, "mpar"))
    stop("mpar must be produced by drifter_par(), bcrw_par(), or bcrw_coa_par().\n",
         "  Got class: ", paste(class(mpar), collapse = ", "))

  expected_move <- switch(class(mpar)[1],
    drifter_par  = "drifter",
    bcrw_par     = "bcrw",
    bcrw_coa_par = "bcrw.coa",
    stop("Unrecognised mpar class: ", class(mpar)[1])
  )
  if (mpar$move != expected_move)
    stop("mpar$move ('", mpar$move, "') does not match constructor class (",
         class(mpar)[1], ")")

  ## Confirm time.step matches FVCOM layer interval.
  ## fvcom_step_secs is pre-computed in sim_setup() from the first two ua layer
  ## names — avoids reading all 8640 layer names on every sim_drifter() call.
  if (is.null(data$fvcom_step_secs))
    stop("data$fvcom_step_secs not found. Re-run sim_setup() to regenerate data.")
  step_secs <- mpar$time.step * 60L
  if (step_secs != data$fvcom_step_secs)
    stop(sprintf(
      "mpar$time.step (%g min) does not match FVCOM layer interval (%g min).\n",
      mpar$time.step, data$fvcom_step_secs / 60),
      "  Set time.step = ", data$fvcom_step_secs / 60, " in your *_par() call."
    )

  invisible(TRUE)
}


## -----------------------------------------------------------------------------
## sim_par() — backward-compatible wrapper routing to the typed constructors
## -----------------------------------------------------------------------------

##' @title Backward-compatible parameter constructor (deprecated)
##'
##' @description Routes to \code{drifter_par()}, \code{bcrw_par()}, or
##'   \code{bcrw_coa_par()} based on the \code{move} argument. New code should
##'   call those constructors directly for explicit parameter sets and
##'   tab-completion support.
##'
##' @param ... Named simulation parameters. See \code{drifter_par()},
##'   \code{bcrw_par()}, \code{bcrw_coa_par()} for parameter details.
##' @export
sim_par <- function(...) {
  dots <- list(...)
  move <- dots[["move"]] %||% "drifter"
  dots[["move"]] <- NULL

  .Deprecated(
    switch(move,
      "drifter"  = "drifter_par()",
      "bcrw"     = "bcrw_par()",
      "bcrw.coa" = "bcrw_coa_par()"
    ),
    msg = paste0(
      "sim_par() is deprecated. Use the move-type constructor:\n",
      "  drifter_par() | bcrw_par() | bcrw_coa_par()"
    )
  )

  do.call(
    switch(move,
      "drifter"  = drifter_par,
      "bcrw"     = bcrw_par,
      "bcrw.coa" = bcrw_coa_par
    ),
    dots
  )
}

## NULL-coalescing operator (avoids rlang dependency)
`%||%` <- function(x, y) if (is.null(x)) y else x
