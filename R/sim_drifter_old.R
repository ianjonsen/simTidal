#' @title simulate drifter advection by surface currents from user-specified start location(s)
#'
#' @description simulates drifter advection
#'
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#'
#' @param id - identifier for simulation run (individual drifter)
#' @param data - a list of required data from \code{sim_setup}
#' @param mpar - simulation control parameters supplied as a list, see details
#' @param pb - use progress bar (logical)
#' @importFrom terra extract xyFromCell nlyr
#' @importFrom CircStats rwrpcauchy
#' @importFrom dplyr %>% mutate lag
#' @importFrom tibble as_tibble
#' @importFrom stats runif rbinom
#' @importFrom lubridate week yday
#' @importFrom stringr str_split
#' @export

sim_drifter_old <- function(
          id=1,
          data = NULL,
          mpar = sim_par(),
          pb = TRUE
  ) {

    N <- mpar$N

    if (is.null(data))
      stop("Can't find output from sim_setup()\n")

    if(mpar$start.dt < mpar$fvcom.origin)
      stop("Simulation start time is earlier than available FVCOM data")

    if(mpar$start.dt > (mpar$fvcom.origin + nlyr(data$u) * 600))
      stop("Simulation start time is later than available FVCOM data")

    if((mpar$start.dt + N * 600) > (mpar$fvcom.origin + nlyr(data$u) * 600))
      stop("Simulation timeframe extends beyond the available FVCOM data")

    if (class(data$land)[1] != "SpatRaster") stop("land must be a SpatRaster")


    if(length(mpar$uvm) == 1) mpar$uvm <- c(1,1)

    ## define location matrix & initialise start position
    ## ds - active swimming displacements
    ## dl - displacements to deflect away from land
    xy <- matrix(NA, N, 2)
    xy[1,] <- cbind(mpar$start)       #cbind(sloc[1], sloc[2])
    ds <- matrix(NA, N, 2)

    ## define other vectors
    reten <- dir <- u <- v <- vector("numeric", N)

    ## determine simulation start date/time offset to start of FVCOM current data
    ##  in 10 min intervals
    fvcom.idx <- as.numeric(difftime(mpar$start.dt, mpar$fvcom.origin, units = "secs")/600)

    ## iterate movement
    for (i in 2:N) {
      if(i==2 && pb)  tpb <- txtProgressBar(min = 2, max = N, style = 3)

      ## Movement kernel
      ds[i, ] <- move_kernel(data, xy = xy[i-1,], mpar = mpar, s = NULL, i = i)

      ### Current Advection
      if (mpar$advect) {
        ## determine envt'l forcing
        ## determine advection due to current, convert from m/s to km/10 min
        u[i] <- as.numeric(extract(data$u[[i-1+fvcom.idx]],
                        rbind(xy[i - 1, ]), method = "simple") * 0.6 * mpar$uvm[1])
        v[i] <- as.numeric(extract(data$v[[i-1+fvcom.idx]],
                        rbind(xy[i - 1, ]), method = "simple") * 0.6 * mpar$uvm[2])
      }

      xy[i, 1:2] <- cbind(ds[i, 1] + u[i],
                          ds[i, 2] + v[i])

      if(!is.na(extract(data$land, rbind(xy[i, ])))) {
        mpar$land <- TRUE
        cat("\n stopping simulation: drifter stuck on land")
        break
      }

      if(any(is.na(xy[i, ]))) {
        mpar$boundary <- TRUE
        cat("\n stopping simulation: drifter hit a boundary")
        break
      }

      if(pb){
        setTxtProgressBar(tpb, i)
        if(i==N) close(tpb)
      }
    }

    N <- ifelse(!is.na(which(is.na(xy[,1]))[1] - 1), which(is.na(xy[,1]))[1] - 1, N)
    X <- data.frame(
      x = xy[, 1],
      y = xy[, 2],
      dx = ds[, 1] - lag(xy[, 1]),
      dy = ds[, 2] - lag(xy[, 2]),
      u = u,
      v = v)[1:N, ]

    sim <- X %>% as_tibble()

    ## remove records after sim is stopped for being stuck on land, etc...
    if(mpar$land | mpar$boundary) {
      sim <- sim %>%
        filter((!is.na(x) & !is.na(y)))
    }

    nsim <- nrow(sim)

    sim <- sim %>%
      mutate(id = id) %>%
      mutate(date = seq(mpar$start.dt, by = 600, length.out = nsim)) %>%
      dplyr::select(id, date, everything())

    param <- mpar
    out <- list(sim = sim, params = param)
    class(out) <- "simTidal"

    return(out)

  }
