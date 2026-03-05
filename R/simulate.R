#' @title simulate smolt migration 
#' 
#' @description simulates smolt migration with or without: 1) temperature-dependent growth; 2) temperature-dependent migration; 3) current advection, using inputs from sim_setup. Movements are simulated as either a biased random walk (brw) or a random walk (rw)
#' 
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#' 
#' @param id - identifier for simulation run (individual animal)
#' @param N - number of 1-h time steps to simulate
#' @param data - a list of required data from \code{presim}
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
#' 
simulate <-
  function(id=1, 
           N = 2760,
           data = NULL,
           mpar = sim_par(),
           pb = TRUE
  ) {
    
    if (is.null(data))
      stop("Can't find output from sim_setup()\n")
    if (class(data$land)[1] != "RasterLayer") stop("distance2land must be a RasterLayer")
      
    ## define location matrix & initialise start position
    ## ds - active swimming displacements
    ## dl - displacements to deflect away from land
    xy <- matrix(NA, N, 2)
    xy[1,] <- cbind(mpar$pars$start)       #cbind(sloc[1], sloc[2])
    ds <- matrix(NA, N, 2)
    ds[1,] <- c(NA, NA)
    
    ## define other vectors
    reten <- dir <- surv <- m <- u <- v <- vector("numeric", N)
    reten[1] <- 1
    surv[1] <- 1
    ## movement type
    move <- mpar$move
    ## present migration direction
    md <- 1
    dir[1] <- mpar$pars$mdir[md]
    r12 <- 1 ## start out using mpar$pars$rho[r12]
    
    if(mpar$growth) {
      s <- ts <- w <- fl <- vector("numeric", N)
     # move_dir <- rep(NA, N)
      w[1] <- mpar$pars$w0
      fl[1] <- (w[1] / 8987.9) ^ (1 / 2.9639)
      s[1] <- fl[1] * mpar$pars$b * 3.6 ## initial swim speed (fl * b body-lengths / s) - in km/h
    } else {
      w <- rep(w[1], N)
      fl <- rep((w[1] / 8987.9) ^ (1 / 2.9639), N)
      s <- fl * mpar$pars$b * 3.6
    }
    
    ## what is start week in env data
    d1 <- as.numeric(str_split(names(data$ts)[1], "d", simplify = TRUE)[,2]) - 1

    ## iterate movement
    for (i in 2:N) {
      if(i==2 && pb)  tpb <- txtProgressBar(min = 2, max = N, style = 3)
      
      ### Apply Energetics
      if (mpar$growth) {
        ## extract Temperature
        ts[i - 1] <- extract(data$ts[[(yday(mpar$pars$start.dt + i * 3600) - d1)]], 
                             rbind(xy[i - 1,])) - 273
        if (is.na(ts[i - 1])) {
          ## calc mean Temp within 2 km buffer of location @ time i-1
          ts[i - 1] <-
            extract(data$ts[[(yday(mpar$pars$start.dt + i * 3600) - d1)]],
                    rbind(xy[i - 1,]),
                    buffer = 2,
                    df = TRUE)[, 2] %>%
            mean(., na.rm = TRUE) - 273
        }
        
        ## calculate growth in current time step based on water temp at previous location, etc...
        if (ts[i - 1] <= 0 & !is.na(ts[i - 1])) {
          cat("\n stopping simulation: smolt has entered water <= 0 deg C")
          break
        } else if (ts[i - 1] > 0 & !is.na(ts[i - 1])) {
          w[i] <- growth(w[i - 1], ts[i - 1], s[i - 1])
          
        } else if (is.na(ts[i - 1])) {
          cat("\n stopping simulation: NA value for temperature")
          break
        }
          
        ## convert W to fL - based on Byron et al 2014 (fig A1 in S1)
        fl[i] <- (w[i] / 8987.9) ^ (1 / 2.9639)
        ## smolts can't shrink their length... (this would unrealistically affect swim speed & energetics), so
        ## if change in weight implies reduction in forklength, then stick with last forklength  - mimics loss of condition
        if (fl[i] < fl[i - 1]) {
          fl[i] <- fl[i - 1]
        }
        
        ## determine size-based average move step for current time step
        ## assume avg swim speed of b bL/s
        s[i] <- fl[i] * mpar$pars$b * 3.6 ## forklength * b m/s converted to km/h
      } 
      
      ## Migration direction scenarios
      if (mpar$scenario == 2 & xy[i-1, 2] <= 1000) {
        ## stop migration if arrived on Grand Bank  
        move <- "rw"
        mpar$pars$b <- 1
        mpar$pars$uvm <- 0.25
        
      } else if(mpar$scenario == 3) { 
        ## 1) travel E from NS/S NL and turn N at random pt & at random rate after passing 
        ##            Avalon Penninsula;
        ## 2) reverse BRW migration direction if SST <= min growth Temp for 3 h
        
        if(xy[i-1,1] < mpar$pars$NFline & xy[i-1,2] < 850) {
          dir[i] <- dir[i-1]
        } else if(xy[i-1,1] >=  mpar$pars$NFline & xy[i-1,2] < 850) {
          ## change direction bias gradually once around SE NF, 
          ##   first to 0 N and then to mdir once N of 850
          ##   this should stop smolts from banging into St John's
          dir[i] <- dir[i - 1] - mpar$pars$turn / 180 * pi
          dir[i] <- ifelse(dir[i] < -0.1745, -0.1745, dir[i]) 
        } else if (xy[i - 1, 2] >= 850) {
          md <- 2
          dir[i] <- mpar$pars$mdir[md]
        }
        
      } else if(mpar$scenario == 4) {
        ## 1) smolts travel S from NB around NS
        ## 2) smolts travel E from NS/S NL and turn N at random pt & 
        ##           at random rate after passing Avalon Penninsula;
        ## 3) smolts reverse BRW migration direction if SST <= min growth Temp for 3 h
        
        # head S from StJohn river to S NS
        if(xy[i-1,2] > 290) {
          ## reset land buffer to 1 km to ensure migration out of BoF
          buff <- mpar$pars$buffer
          mpar$pars$buffer <- 1
          dir[i] <- dir[i-1]
        } else if (xy[i-1,2] < 290 & dir[i] > mpar$pars$mdir[2]) {
          ## start turning toward new migration heading (to Grand Banks)
          md <- 2
          dir[i] <- dir[i-1] - mpar$pars$turn/180*pi
          ## re-establish original land buffer once out of BoF
          mpar$pars$buffer <- buff
          
        } else if (xy[i-1,2] < 290 & dir[i] < mpar$pars$mdir[2]) {
          ## stop the turn once new migration heading is reached
            md <- 2
            dir[i] <- mpar$pars$mdir[md]
        }
        if (xy[i-1,1] >= mpar$pars$NFline) {
          ## change direction bias gradually once around SE NF, 
          ##   first to 0 N and then to mdir once N of 950
          ##   this should stop smolts from banging into St John's
          if(xy[i-1,2] < 850 & dir[i] > -10/180*pi) {
            dir[i] <- dir[i-1] - mpar$pars$turn/180*pi
          } else if(xy[i-1,2] >= 850 & dir[i] > mpar$pars$mdir[3]) {
            md <- 3
            dir[i] <- dir[i-1] - mpar$pars$turn/180*pi
          } else if(xy[i-1,2] >= 850 & dir[i] < mpar$pars$mdir[3]) {
            md <- 3
            dir[i] <- mpar$pars$mdir[md]
          }
        }
        
      } else if(mpar$scenario == 5) {
        ## 1) smolts depart Campbellton River and head N with migration influenced by SST
        ## Campbellton River, NL
        if(xy[i-1,2] < 1200) {
          dir[i] <- dir[i-1]
          if(dir[i-1] != mpar$pars$mdir[2]) {
            r12 <- 1
          }
        } else if(xy[i-1,2] >= 1200) {
          md <- 2
          dir[i] <- mpar$pars$mdir[md]
          r12 <- 2
        } else {
          dir[i] <- dir[i-1]
        }
      }
      
      ## Temperature-dependent alteration of migration direction - 
      ##      reverse course if water too cold.This over-rides migration 
      ##      scenario in the short-term
      if (mpar$temp) {
        ## reverse migration direction from mpar$pars$mdir to 
        ##    mpar$pars$mdir - pi, if smolt in SST <= min growth C for 3 h
        if (i > 3) { 
          ## movement direction influenced by Temp experienced over past 3 h
          ##  implies 0 or -ve growth, for current mass(w[i]) and speed (s[i])
          g.rng <- growth(w[i], seq(6, 20, l = 100), s[i])
          ts.mig <- seq(6, 20, l=100)[which(g.rng >= w[i])] %>% min()
          dir[i] <- ifelse(all(ts[i - 1:3] <= ts.mig), 
                           (mpar$pars$mdir[md] - pi) %% pi, 
                           dir[i])
          dir[i] <- ifelse(abs(dir[i] - mpar$pars$mdir[md]) > 45 & 
                             all(ts[i - 1:3] > ts.mig), 
                           mpar$pars$mdir[md], 
                           dir[i])
          
          ## can use this for longer runs (ie. Kelts w 440 d tags)
          # ## if smolt in optimal T range for growth then switch from biased RW to simple RW
          # ts.rng <- seq(6, 20, l=100)[which(g.rng >= quantile(g.rng, mpar$pars$ts.q))]  %>% range()
          # b.tmp <- mpar$pars$b
          # if(ts[i-1] >= ts.rng[1] & ts[i-1] <= ts.rng[2]) {
          #   move <- "rw"
          #   mpar$pars$b <- 1
          # } else {
          #   move <- mpar$move
          #   mpar$pars$b <- b.tmp
          # }
            
        } 
      
      } 
        
      ## Movement kernel
      ## First check if smolt is in SoBI & within 25 km of land, if so then move toward Lab Sea after which reg movement rules can be applied
      d2l <- extract(data$land, rbind(xy[i-1,]))
        if(all(xy[i-1,1] >= data$sobi.box[1], 
               xy[i-1,1] <= data$sobi.box[2], 
               xy[i-1,2] >= data$sobi.box[3], 
               xy[i-1,2] <= data$sobi.box[4]) & d2l < 15) {

          phi <- rwrpcauchy(1, ifelse(mpar$pars$psi< 0.5, 0.28, 0.35) * pi, 0.8)
          ds[i, 1] <- xy[i-1, 1] + s[i] * sin(phi)
          ds[i, 2] <- xy[i-1, 2] + s[i] * cos(phi)
          
        } else {
          
        ds[i, ] <- switch(move,
                          brw = {
                            
                            brw(n=1, 
                                      data, 
                                      xy = xy[i-1,], 
                                      coa = mpar$pars$coa, 
                                      dir = dir[i],
                                      buffer = mpar$pars$buffer, 
                                      rho = mpar$pars$rho[r12],
                                      a = mpar$pars$a,
                                      b = s[i], 
                                      shelf = mpar$shelf,
                                      beta = mpar$pars$beta)
                          },
                          rw = {

                            rw(n=1,
                                data,
                                xy = xy[i-1,],
                                buffer = mpar$pars$buffer,
                                rho = mpar$pars$rho[2],
                                a = mpar$pars$a,
                                b = s[i],
                                shelf = mpar$shelf,
                                beta = mpar$pars$beta)
                          },
                          drift = {
                            ds[i, ] <- rbind(xy[i-1, 1:2])
                          })
        }
      
      ### Current Advection
      if (mpar$advect) {
        ## determine envt'l forcing
        ## determine advection due to current, convert from m/s to km/h
        u[i] <- extract(data$u[[(yday(mpar$pars$start.dt + i * 3600) - d1)]], 
                        rbind(xy[i - 1, ])) * 3.6 * mpar$par$uvm
        v[i] <- extract(data$v[[(yday(mpar$pars$start.dt + i * 3600) - d1)]], 
                        rbind(xy[i - 1, ])) * 3.6 * mpar$par$uvm
 
        
      } else if(!mpar$advect) {
        u[i] <- v[i] <- 0
      }
      
      xy[i, 1:2] <- cbind(ds[i, 1] + u[i], 
                          ds[i, 2] + v[i])
      
      if((extract(data$land, rbind(xy[i, ])) == 0 | is.na(extract(data$land, rbind(xy[i, ]))))  & any(!is.na(xy[i,]))) {
        mpar$land <- TRUE
        cat("\n stopping simulation: stuck on land")
        break
      }
      
      if(any(is.na(xy[i, ]))) {
        mpar$boundary <- TRUE
        cat("\n stopping simulation: hit a boundary")
        break
      } 
      ## determine survival
      if(!is.na(mpar$pars$surv)) {
        surv[i] <-
          rbinom(1, 1, mpar$pars$surv ^ (1 / 24)) # rescales daily survival to hourly
        if (surv[i] == 0) {
          cat("\n smolt is dead")
          break
        }
      }
      
      ## determine tag retention
      if(!is.na(mpar$pars$reten) & i <= mpar$pars$Dreten*24) {
        reten[i] <- rbinom(1, 1, mpar$pars$reten ^ (1 / 24))
        if(reten[i] == 0) {
          cat("\n tag expulsion")
          break
        }
      } else if(!is.na(mpar$pars$reten) & i > mpar$pars$Dreten*24) {
        reten[i] <- 1
      }
      
      if(pb){
        setTxtProgressBar(tpb, i)
        if(i==N) close(tpb)
      }
    }
  
    N <- ifelse(!is.na(which(is.na(xy[,1]))[1] - 1), which(is.na(xy[,1]))[1] - 1, N)
    if(mpar$growth) {
    X <-
      data.frame(
        x = xy[, 1],
        y = xy[, 2],
        dx = ds[, 1] - lag(xy[, 1]),
        dy = ds[, 2] - lag(xy[, 2]),
        u = u,
        v = v,
        ts = ts,
        w = w,
        fl = fl,
        s = s, 
        surv = surv,
        reten = reten,
        dir = dir
      )[1:N, ]
    } else if(!mpar$growth) {
      X <-
        data.frame(
          x = xy[, 1],
          y = xy[, 2],
          dx = ds[, 1] - lag(xy[, 1]),
          dy = ds[, 2] - lag(xy[, 2]),
          u = u,
          v = v, 
          ts = ifelse(mpar$temp, ts, NA),
          surv = surv,
          reten = reten
        )[1:N, ] 
    } 
    if(sum(is.na(X$ts)) == nrow(X)) {
      X <- X %>% select(-ts)
    } else {
      X$ts[nrow(X)] <- ifelse(X$ts[nrow(X)] == 0, NA, X$ts[nrow(X)])
    }
    
    sim <- X %>% as_tibble() 
    
    ## remove records after sim is stopped for being stuck on land, etc...
    if(mpar$land | mpar$boundary) {
      if(mpar$growth) {
        sim <- sim %>%
          filter((!is.na(x) & !is.na(y) & w != 0 & fl != 0 & s != 0) | surv != 1 | reten != 1)
      } else {
        sim <- sim %>%
          filter((!is.na(x) & !is.na(y)) | surv != 1 | reten != 1)
      }
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
