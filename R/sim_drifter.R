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
#' @importFrom raster extract xyFromCell
#' @importFrom CircStats rwrpcauchy
#' @importFrom dplyr %>% mutate lag
#' @importFrom tibble as_tibble
#' @importFrom stats runif rbinom
#' @importFrom lubridate week yday
#' @importFrom stringr str_split
#' @export

sim_drifter <- function(id=1,
           data = NULL,
           mpar = sim_par(),
           pb = TRUE
  ) {
    
    if (is.null(data))
      stop("Can't find output from sim_setup()\n")
    if (class(data$land)[1] != "RasterLayer") stop("distance2land must be a RasterLayer")
    
    N <- mpar$pars$N
    ## define location matrix & initialise start position
    ## ds - active swimming displacements
    ## dl - displacements to deflect away from land
    xy <- matrix(NA, N, 2)
    xy[1,] <- cbind(mpar$pars$start)       #cbind(sloc[1], sloc[2])
    ds <- matrix(NA, N, 2)
    
    ## define other vectors
    reten <- dir <- u <- v <- vector("numeric", N)
    
    ## iterate movement
    for (i in 2:N) {
      if(i==2 && pb)  tpb <- txtProgressBar(min = 2, max = N, style = 3)
      
      ## Movement kernel
      ds[i,] <- moveDrifter(data, xy = xy[i-1,], mpar = mpar)
      
      ### Current Advection
      if (mpar$advect) {
        ## determine envt'l forcing
        ## determine advection due to current, convert from m/s to km/h
        u[i] <- extract(data$u[[yday(mpar$pars$start.dt + i * 3600)]],
                        rbind(xy[i - 1, ]), method = "simple") * 3.6 * mpar$par$uvm
        v[i] <- extract(data$v[[yday(mpar$pars$start.dt + i * 3600)]],
                        rbind(xy[i - 1, ]), method = "simple") * 3.6 * mpar$par$uvm
        
        ## turn off advection in sobi.box b/c too challenging to get through w currents...
      } else if(!mpar$advect | all(xy[1] >= data$sobi.box[1],
                                   xy[1] <= data$sobi.box[2],
                                   xy[2] >= data$sobi.box[3],
                                   xy[2] <= data$sobi.box[4])) {
        u[i] <- v[i] <- 0
      }
      
      xy[i, 1:2] <- cbind(ds[i, 1] + u[i],
                          ds[i, 2] + v[i])
      
      if((extract(data$land, rbind(xy[i, ])) == 0 |
          is.na(extract(data$land, rbind(xy[i, ]))))  & any(!is.na(xy[i,]))) {
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
      mutate(date = seq(mpar$pars$start.dt, by = 3600, length.out = nsim)) %>%
      select(id, date, everything())
    
    param <- mpar
    out <- list(sim = sim, params = param)
    class(out) <- "simsmolt"
    
    return(out)
    
  }
