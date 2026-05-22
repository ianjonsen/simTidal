#' @title Movement kernel for fish and drifter simulations
#'
#' @description Internal step function called at each time step by
#'   \code{\link{sim_fish}} and \code{\link{sim_drifter}}. Computes a
#'   proposed new position from the movement process specified in
#'   \code{mpar$move}, applies the land-avoidance potential when
#'   \code{data$land} is present, and returns the accepted position.
#'
#'   Four active swimming processes are supported:
#'   \describe{
#'     \item{\code{"crw"}}{Correlated random walk. The next heading is
#'       drawn from a Wrapped Cauchy distribution centred on the
#'       \emph{previous} heading (\code{xy[3]}). No external directional
#'       bias; step-to-step persistence is controlled entirely by
#'       \code{mpar$rho}.}
#'     \item{\code{"bcrw"}}{Biased correlated random walk toward a fixed
#'       bearing. The Wrapped Cauchy distribution is centred on
#'       \code{mpar$bearing[i]} rather than the previous heading, so the
#'       animal is persistently pulled toward that compass direction.}
#'     \item{\code{"bcrw.coa"}}{Biased correlated random walk toward a
#'       Centre of Attraction. The bias direction is a weighted blend of
#'       the previous heading (\code{xy[3]}) and the instantaneous bearing
#'       toward \code{mpar$coa}, with blend strength \code{mpar$nu}.}
#'     \item{\code{"crw.bridge"}}{Brownian bridge correlated random walk.
#'       The bias direction is a weighted blend of the previous heading and
#'       the bearing toward \code{end_xy} (the per-simulation sampled end
#'       position). Blend weight \code{nu_i = mpar$nu * i / max(N - i, 1)}
#'       grows from near-zero at the first step to very large at the final
#'       step, ensuring the path is progressively drawn toward
#'       \code{end_xy}.}
#'   }
#'   A fifth passive process, \code{"drifter"}, produces no active
#'   displacement; optional Gaussian noise can be added via
#'   \code{mpar$noise}.
#'
#' @param data Environmental data list from \code{\link{sim_setup}}, or
#'   \code{NULL}. Elements \code{data$land} and \code{data$grad} are used
#'   for land-avoidance. The entire land-avoidance block is skipped when
#'   \code{data$land} is \code{NULL}.
#' @param xy   Current state vector \code{c(x, y, heading)} in km and
#'   radians. The heading (element 3) is the animal's current travel
#'   direction and is read by \code{"crw"} (bias = previous heading),
#'   \code{"bcrw.coa"} (blended with CoA bearing), and
#'   \code{"crw.bridge"} (blended with bearing toward \code{end_xy});
#'   it is not read by \code{"bcrw"} or \code{"drifter"}.
#' @param mpar Parameter list produced by \code{\link{fish_par}} or one
#'   of the \code{*_par} constructors. Fields accessed depend on
#'   \code{mpar$move}; see \code{\link{fish_par}} for details.
#' @param s    Swimming step size in km per step, pre-computed by the
#'   calling function as
#'   \code{fl / 1000 * bl * 60 * time.step}.
#'   Ignored when \code{move = "drifter"}.
#' @param i    Current simulation step index (integer). Used by
#'   \code{"bcrw"} to index time-varying \code{mpar$bearing} or
#'   \code{mpar$bl} vectors. Also used by \code{"crw.bridge"} to compute
#'   the time-varying attraction weight \code{nu_i}. Not used by
#'   \code{"crw"}, \code{"bcrw.coa"}, or \code{"drifter"}.
#' @param end_xy Per-simulation sampled end position, \code{c(x, y)} in
#'   km. Required (non-\code{NULL}) when \code{mpar$move = "crw.bridge"};
#'   ignored by all other movement models. Default \code{NULL}.
#'
#' @return For \code{move \%in\% c("crw", "crw.bridge", "bcrw", "bcrw.coa")}: a
#'   length-3 numeric vector \code{c(x_new, y_new, mu)}, where \code{mu}
#'   is the step heading in radians (to be stored and passed back as
#'   \code{xy[3]} on the next call). For \code{move = "drifter"}: a
#'   length-2 vector \code{c(x_new, y_new)}.
#'
#' @seealso \code{\link{sim_fish}}, \code{\link{sim_drifter}},
#'   \code{\link{fish_par}}
#'
#' @importFrom CircStats rwrpcauchy
#' @importFrom stats rnorm
#' @importFrom terra extract xyFromCell
#' @keywords internal

move_kernel <- function(data, xy = NULL, mpar, s, i, end_xy = NULL) {

  switch(mpar$move,

    crw = {
      ## Correlated random walk: next heading drawn around the previous
      ## heading (xy[3]) with concentration rho. No external bias.
      phi    <- xy[3]
      mu     <- rwrpcauchy(1, phi, mpar$rho)
      new.xy <- cbind(xy[1] + sin(mu) * s, xy[2] + cos(mu) * s)
    },

    crw.bridge = {
      ## Brownian bridge CRW: time-varying attraction toward end_xy.
      ## nu_i grows from near-zero at the first step to mpar$nu * N at
      ## the last step, ensuring the path converges toward end_xy.
      if (is.null(end_xy))
        stop("end_xy must be provided for move = 'crw.bridge'")
      delta <- c(end_xy[1] - xy[1], end_xy[2] - xy[2])
      psi   <- atan2(delta[1], delta[2])
      nu_i  <- mpar$nu * i / max(mpar$N - i, 1L)
      phi   <- atan2(sin(xy[3]) + nu_i * sin(psi),
                     cos(xy[3]) + nu_i * cos(psi))
      mu     <- rwrpcauchy(1, phi, mpar$rho)
      new.xy <- cbind(xy[1] + sin(mu) * s, xy[2] + cos(mu) * s)
    },

    bcrw.coa = {
      ## Biased correlated RW toward a Centre of Attraction.
      ## Blend previous heading (xy[3]) with bearing to CoA (psi),
      ## weighted by nu; larger nu = stronger attraction.
      if (all(!is.na(mpar$coa))) {
        delta <- c(mpar$coa[1] - xy[1], mpar$coa[2] - xy[2])
        psi   <- atan2(delta[1], delta[2])
        phi   <- atan2(sin(xy[3]) + mpar$nu * sin(psi),
                       cos(xy[3]) + mpar$nu * cos(psi))
      } else {
        phi <- atan2(sin(xy[3]), cos(xy[3]))
      }
      mu     <- rwrpcauchy(1, phi, mpar$rho)
      new.xy <- cbind(xy[1] + sin(mu) * s, xy[2] + cos(mu) * s)
    },

    bcrw = {
      ## Biased correlated RW toward a fixed bearing.
      phi <- mpar$bearing[i]

      ## rho may be a scalar or a step-indexed vector
      mu <- if (length(mpar$rho) == 1) {
        rwrpcauchy(1, phi, mpar$rho)
      } else {
        rwrpcauchy(1, phi, mpar$rho[i])
      }

      ## bl may be a scalar or a step-indexed vector
      if (length(mpar$bl) > 1)
        s <- mpar$fl / 1000 * mpar$bl[i] * 60 * mpar$time.step

      new.xy <- cbind(xy[1] + sin(mu) * s, xy[2] + cos(mu) * s)
    },

    drifter = {
      ## Passive advection — no active swimming displacement.
      new.xy <- cbind(xy[1], xy[2])
      if (!is.null(mpar$noise))
        new.xy <- new.xy + rnorm(2, 0, mpar$noise)
    }

  ) ## end switch

  ## ---- Land-avoidance potential ----------------------------------------------
  if (!is.null(data$land)) {

    pv    <- c(0, 0)
    pv[1] <- as.numeric(extract(data$grad[[1]], new.xy))
    pv[2] <- as.numeric(extract(data$grad[[2]], new.xy))
    if (any(is.na(pv))) pv <- c(0, 0)

    new2.xy <- new.xy + pv * mpar$beta

    if (!is.na(extract(data$land, new2.xy))) {
      ## First repulsion attempt landed on land — try doubling beta
      pv    <- c(0, 0)
      pv[1] <- as.numeric(extract(data$grad[[1]], new2.xy))
      pv[2] <- as.numeric(extract(data$grad[[2]], new2.xy))
      if (any(is.na(pv))) pv <- c(0, 0)

      new3.xy <- new.xy + pv * (mpar$beta * 2)

      if (!is.na(extract(data$land, new3.xy))) {
        ## Still on land — find nearest water cell within buffer
        message("finding water...")
        cells <- try(
          extract(data$land, rbind(new.xy),
                  buffer = mpar$buffer, cellnumbers = TRUE, df = TRUE),
          silent = TRUE
        )
        if (inherits(cells, "try-error")) {
          new.xy <- xy[1:2]   ## fall back to previous position
        } else {
          idx        <- which(is.na(cells[, 3]))[1]
          cell.water <- cells[idx, 2]
          new.xy     <- xyFromCell(data$land, cell.water) %>% rbind()
        }
      } else {
        new.xy <- new3.xy
      }
    } else {
      new.xy <- new2.xy
    }
  }

  ## ---- Return ---------------------------------------------------------------
  if (mpar$move != "drifter") {
    return(cbind(new.xy[1], new.xy[2], mu))
  } else {
    return(new.xy)
  }
}
