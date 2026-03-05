#' @title surface current transport kernel for drifters
#'
#' @description utility function not to be called by user
#'
#' @export
#'
moveDrifter <- function(data, xy = NULL, mpar) {

  ## calculate distance to land
  d2l <- extract(data$land, rbind(xy))

  if (d2l > mpar$pars$buffer) {
    newxy <- cbind(xy[1], xy[2])

  } else {
    ## keep drifter 2 km > buffer off land
    ## find all nearby cells within 5 km & select the one farthest from land
    cells <- extract(data$land, rbind(xy), buffer = mpar$pars$buffer - d2l + 2,
                     cellnumbers = TRUE, df = TRUE)
    cell.max <- cells[cells[, 3] == max(cells[, 3], na.rm = TRUE), 2][1]
    new.xy <- xyFromCell(data$land, cell.max) %>% as.vector()
    newxy <- cbind(new.xy[1], new.xy[2])
  }
  if(!is.null(mpar$noise)) {
    return(newxy + rnorm(2, 0, mpar$noise))
  } else {
    return(newxy)
  }
}
