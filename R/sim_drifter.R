#' @title simulate drifter advection by surface currents from user-specified start location(s)
#'
#' @description simulates drifter advection
#'
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#'
#' @param id   identifier for simulation run (individual drifter)
#' @param data a list of required data from \code{sim_setup()}
#' @param mpar simulation control parameters from \code{drifter_par()},
#'   \code{bcrw_par()}, or \code{bcrw_coa_par()}
#' @param pb   use progress bar (logical)
#'
#' @details \code{fvcom.origin} and \code{fvcom_step_secs} are read from
#'   \code{data} (set automatically by \code{sim_setup()}), with fallback to
#'   \code{mpar$fvcom.origin} for backward compatibility.
#'
#' @importFrom terra extract xyFromCell nlyr
#' @importFrom CircStats rwrpcauchy
#' @importFrom dplyr %>% mutate lag filter
#' @importFrom tibble as_tibble
#' @importFrom stats runif rbinom rnorm
#' @importFrom lubridate week yday
#' @importFrom stringr str_split
#' @export

sim_drifter <- function(
    id   = 1,
    data = NULL,
    mpar = NULL,
    pb   = TRUE
) {

  ## ---- Validation ----------------------------------------------------------

  if (is.null(data))
    stop("data is NULL. Pass the output of sim_setup().")

  if (is.null(mpar))
    stop("mpar is NULL. Create with e.g.:\n",
         "  mpar <- drifter_par(N = 1440, start = c(x, y))")

  validate_mpar(mpar, data)

  N         <- mpar$N
  step_secs <- mpar$time.step * 60L

  ## advection scaling: m/s -> km/step  (= m/s * step_secs / 1000)
  advect_scale <- step_secs / 1000

  ## fvcom.origin from sim_setup() (preferred) or mpar fallback
  fvcom.origin <- if (!is.null(data$fvcom.origin)) {
    data$fvcom.origin
  } else if (!is.null(mpar$fvcom.origin)) {
    mpar$fvcom.origin
  } else {
    stop("fvcom.origin not found in data or mpar.\n",
         "  Re-run sim_setup() or set mpar$fvcom.origin manually.")
  }

  ## Resolve NULL start.dt to first available FVCOM timestamp
  if (is.null(mpar$start.dt)) {
    mpar$start.dt <- fvcom.origin
    message("start.dt not set — using data$fvcom.origin: ", format(fvcom.origin))
  }

  ## ---- Initialise state ----------------------------------------------------

  xy <- matrix(NA, N, 2)
  xy[1, ] <- mpar$start
  ds <- matrix(NA, N, 2)

  u <- v <- vector("numeric", N)

  n_u_layers <- nlyr(data$u)

  ## ---- Handle sub-interval start.dt ----------------------------------------
  ##
  ## w is the fractional offset of start.dt within its FVCOM 10-min interval.
  ## Because each simulation step advances by exactly step_secs, w is constant
  ## for the entire run — computed once here, used in every loop iteration.
  ##
  ## interp = FALSE (default): snap start.dt to nearest layer boundary (w -> 0),
  ##   warn the user, use 1 extract per component per step.
  ## interp = TRUE: keep start.dt as-is, interpolate linearly between the two
  ##   bracketing layers using w, at the cost of 2 extracts per component per step.

  offset_secs <- as.numeric(difftime(mpar$start.dt, fvcom.origin, units = "secs"))
  remainder   <- offset_secs %% step_secs

  if (remainder == 0) {
    w <- 0
  } else if (mpar$interp) {
    w <- remainder / step_secs    # fractional weight toward next layer; in (0, 1)
  } else {
    ## Save original start.dt for the warning message before modifying
    original_dt <- mpar$start.dt
    if (remainder < step_secs / 2) {
      mpar$start.dt <- mpar$start.dt - remainder
    } else {
      mpar$start.dt <- mpar$start.dt + (step_secs - remainder)
    }
    warning(sprintf(
      "start.dt (%s) is not on a %g-min FVCOM boundary.\n  Snapped to nearest layer: %s\n  Use interp = TRUE to interpolate instead.",
      format(original_dt), mpar$time.step, format(mpar$start.dt)
    ))
    w <- 0
  }

  ## Bounds checks are placed after snapping so that a start.dt snapped
  ## forward cannot silently exceed the available data range.
  if (mpar$start.dt < fvcom.origin)
    stop("Simulation start time is earlier than available FVCOM data\n",
         "  fvcom.origin: ", format(fvcom.origin), "\n",
         "  start.dt:     ", format(mpar$start.dt))

  if (mpar$start.dt > (fvcom.origin + n_u_layers * step_secs))
    stop("Simulation start time is later than available FVCOM data")

  if ((mpar$start.dt + N * step_secs) > (fvcom.origin + n_u_layers * step_secs))
    stop("Simulation timeframe extends beyond the available FVCOM data")

  ## round() before as.integer() guards against floating-point precision issues
  ## (e.g. 599.9999.../600 truncating to 0 instead of 1 after snapping).
  fvcom.idx <- as.integer(round(
    as.numeric(difftime(mpar$start.dt, fvcom.origin, units = "secs")) / step_secs
  ))

  ## Pre-compute layer indices for every simulation step.
  ## For loop step i (2..N), the lower bracketing FVCOM layer is (i-1) + fvcom.idx.
  ## When interp = TRUE, the upper layer is layer_idx + 1.
  ## Bounds check ensures layer_idx + 1 never exceeds nlyr(data$u).
  layer_idx <- seq_len(N - 1L) + fvcom.idx

  if (mpar$interp && w > 0 && max(layer_idx) >= n_u_layers)
    stop("Simulation timeframe requires layer ", max(layer_idx) + 1L,
         " but data$u only has ", n_u_layers, " layers.\n",
         "  Reduce N or choose an earlier start.dt.")

  ## ---- Main loop -----------------------------------------------------------

  for (i in 2:N) {
    if (i == 2 && pb) tpb <- txtProgressBar(min = 2, max = N, style = 3)

    ds[i, ] <- move_kernel(data, xy = xy[i - 1, ], mpar = mpar, s = NULL, i = i)

    if (mpar$advect) {
      k <- layer_idx[i - 1L]
      if (w == 0) {
        ## Snap: single extract per component
        u[i] <- extract(data$u[[k]], rbind(xy[i - 1, ]),
                        method = "simple")[1, 1] * advect_scale * mpar$uvm[1]
        v[i] <- extract(data$v[[k]], rbind(xy[i - 1, ]),
                        method = "simple")[1, 1] * advect_scale * mpar$uvm[2]
      } else {
        ## Interpolate: weighted average of layers k and k+1
        u[i] <- ((1 - w) * extract(data$u[[k]],      rbind(xy[i - 1, ]), method = "simple")[1, 1] +
                       w  * extract(data$u[[k + 1L]], rbind(xy[i - 1, ]), method = "simple")[1, 1]) *
                 advect_scale * mpar$uvm[1]
        v[i] <- ((1 - w) * extract(data$v[[k]],      rbind(xy[i - 1, ]), method = "simple")[1, 1] +
                       w  * extract(data$v[[k + 1L]], rbind(xy[i - 1, ]), method = "simple")[1, 1]) *
                 advect_scale * mpar$uvm[2]
      }
    }

    xy[i, 1:2] <- cbind(ds[i, 1] + u[i], ds[i, 2] + v[i])

    if (!is.na(extract(data$land, rbind(xy[i, ]))[1, 1])) {
      mpar$land <- TRUE
      cat("\n stopping simulation: drifter stuck on land")
      break
    }

    if (any(is.na(xy[i, ]))) {
      mpar$boundary <- TRUE
      cat("\n stopping simulation: drifter hit a boundary")
      break
    }

    if (pb) {
      setTxtProgressBar(tpb, i)
      if (i == N) close(tpb)
    }
  }

  ## ---- Assemble output -----------------------------------------------------

  N <- ifelse(!is.na(which(is.na(xy[, 1]))[1] - 1),
              which(is.na(xy[, 1]))[1] - 1, N)

  sim <- data.frame(
    x  = xy[, 1],
    y  = xy[, 2],
    dx = ds[, 1] - lag(xy[, 1]),
    dy = ds[, 2] - lag(xy[, 2]),
    u  = u,
    v  = v
  )[1:N, ] %>%
    as_tibble()

  if (mpar$land | mpar$boundary)
    sim <- sim %>% filter(!is.na(x) & !is.na(y))

  nsim <- nrow(sim)

  sim <- sim %>%
    mutate(id   = id,
           date = seq(mpar$start.dt, by = step_secs, length.out = nsim)) %>%
    dplyr::select(id, date, everything())

  structure(list(sim = sim, params = mpar), class = "simTidal")
}
