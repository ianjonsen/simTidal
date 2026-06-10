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
