## Internal helpers shared across plotting functions

## Compute the most probable path (MPP) from a list of sim_fish track tibbles.
## Returns a data frame with columns step, x, y using the requested summary
## statistic across all tracks in sims_list at each step index.
## Tracks shorter than max(N) contribute NA beyond their last position, which
## is dropped by the na.rm argument of the summary function.

.compute_mpp <- function(sims_list, method = "median") {

  if (length(sims_list) == 0L)
    return(data.frame(step = integer(0), x = numeric(0), y = numeric(0)))

  N    <- max(vapply(sims_list, nrow, integer(1L)))
  nsim <- length(sims_list)

  x_mat <- matrix(NA_real_, nsim, N)
  y_mat <- matrix(NA_real_, nsim, N)

  for (i in seq_len(nsim)) {
    n <- nrow(sims_list[[i]])
    x_mat[i, seq_len(n)] <- sims_list[[i]]$x
    y_mat[i, seq_len(n)] <- sims_list[[i]]$y
  }

  fn <- if (method == "median")
    function(v) stats::median(v, na.rm = TRUE)
  else
    function(v) mean(v, na.rm = TRUE)

  data.frame(
    step = seq_len(N),
    x    = apply(x_mat, 2L, fn),
    y    = apply(y_mat, 2L, fn)
  )
}


##' Perpendicular distance from points to a polyline, vectorised over points.
##'
##' Used by sim_fish() to record how close a simulation came to the receiver
##' line even when it did not cross it. A near miss is information; absence is
##' not. `xy` is the polyline as a two-column matrix in order along the line.
.dist_to_polyline <- function(px, py, xy) {
  best <- rep(Inf, length(px))
  for (r in seq_len(nrow(xy) - 1L)) {
    ax <- xy[r, 1L]; ay <- xy[r, 2L]
    bx <- xy[r + 1L, 1L]; by <- xy[r + 1L, 2L]
    vx <- bx - ax; vy <- by - ay
    L2 <- vx * vx + vy * vy
    t  <- if (L2 <= 0) rep(0, length(px)) else
      pmin(1, pmax(0, ((px - ax) * vx + (py - ay) * vy) / L2))
    d <- sqrt((px - (ax + t * vx))^2 + (py - (ay + t * vy))^2)
    best <- pmin(best, d, na.rm = TRUE)
  }
  best
}
