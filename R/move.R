#' @title random walk movement kernel for Campbellton River kelts
#' 
#' @description utility function not to be called by user
#' 
#' @importFrom CircStats rwrpcauchy
#' @importFrom stats rweibull
#' @importFrom raster extract xyFromCell
#' @export
#' 
moveKcam <- function(n = 1, data, xy = NULL, coa = NULL, dir = NULL, buffer = NULL, rho, a, b, shelf, beta) {
  
  ## Migration direction scenarios
  ## 1) kelts depart Campbellton River and head N with migration influenced by SST
  ## Campbellton River, NL
  if (xy[i - 1, 2] < 1200) {
    dir[i] <- dir[i - 1]
    if (dir[i - 1] != mpar$pars$mdir[2]) {
      r12 <- 1
    }
  } else if (xy[i - 1, 2] >= 1200) {
    md <- 2
    dir[i] <- mpar$pars$mdir[md]
    r12 <- 2
  } else {
    dir[i] <- dir[i - 1]
  }
  
  
  ## Temperature-dependent alteration of migration direction - 
  ##      reverse course if water too cold. This over-rides migration 
  ##      scenario in the short-term

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