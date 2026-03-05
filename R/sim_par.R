##' \code{sim_par} defines the simulation parameters & control scenarios used by \code{simulate}.
##'
##' The movement process used predominantly in the simulation is
##' selected by the \code{move} argument.  Additional
##' parameters include: \code{temp} is movement temperature-dependent; \code{advect} do ocean
##' currents influence movement; \code{growth} do smolts grow in temp-dependent fashion;
##' \code{start.date};
##' \code{start} (location); \code{coa} centre of attraction (can be NULL);
##' \code{mdir} directional bias (a 2-element vector);
##' \code{rho} concentration of step directions (for wrapped-Cauchy, a 2-element vector); ...
##'
##' @title Control Values for \code{simulate}.
##' @param move "brw" a biased RW; "rw" simple RW; "drift" no active swimming
##' @param temp logical
##' @param advect logical
##' @param noise numeric; add small random component to u,v vectors - only applies when mode = "drift"
##' @param growth logical
##' @param shelf logical; should smolts be constrained to stay in water > - 1000 m depth (continental shelf)
##' @param migs migration scenario 1 or 2 (stop migration upon arrival to Grand Banks)
##' @param land keep track of sim rep hitting land (TRUE)
##' @param boundary keep track of sim rep hitting sim boundary (TRUE)
##' @param ... additional simulation control parameters
##' @return Returns a list with components
##'   \item{\code{move}}{the main move process}
##'   \item{\code{temp}}{temperature-dependent movements}
##'   \item{\code{advect}}{ocean-current-dependent movements}
##'   \item{\code{noise}}{add random noise to current u,v values}
##'   \item{\code{growth}}{temperature-dependent growth}
##'   \item{\code{shelf}}{movements constrained to stay on shelf}
##'   \item{\code{migs}}{migration strategy}
##'   \item{\code{land}}{record if smolt gets "stuck" on land, thereby ending simulation}
##'   \item{\code{boundary}}{record if smolt hits simulation boundary, thereby ending simulation}
##'   \item{\code{pars}}{list of additional, required control parameters}
##' @export

sim_par <-
  function(move = c("brw","rw","drift"),
           temp = TRUE,
           advect = TRUE,
           noise = NULL,
           growth = TRUE,
           shelf = TRUE,
           scenario = "rs",
           land = FALSE,
           boundary = FALSE,
           ...) {

    move <- match.arg(move)

    dots <- list(...)

    pars <- list(
      N = 1440,
      start.dt = ISOdatetime(2021,07,09,00,00,00, tz = "UTC"),
      start = c(7305, 1425),
      coa = NULL,
      mdir = c(75,-45)/180*pi, # bias direction for N migration (S migration is mdir - pi)
      NFline.x = runif(1, 7700, 7900), # location on x-axis at which smolt turns N to Lab Shelf
      NFline.y = runif(1, 1950, 2125),
      pN = 0.75,
      rho = c(0.6, 0.4), # directional persistence for brw [1] and rw [2]
      turn = 2.5, # rate at which smolts turn N after rounding NF
      ntries = 1,
      ts.q = 0.75,
      psi = 0.9,
      uvm = 1, # magnitude of current vectors: if uvm < 1 current strength is down-scaled
      buffer = 5,
      b = 2,
      a = 0,
      w0 = 185,
      surv = 0.9936, ## daily survival rate
      reten = 0.845^(1/60), ## daily V7/8 tag retention rate (in first ~ 60 d - Brundson et al 2019 ICES JMarSci 76:7)
      Dreten = 60, ## number of days within which tags can be expulsed
      pdrf = c(5, -0.02), # = p(0.5) @ 250 m  + < 0.01 @ 500 m   [c(4.865, -0.0139)  (~ consistent w HFX line V9 @ high power)]
      beta = c(-2, -2) # potential fn params to keep smolts on shelf
    )

    ## overide default control pars
    pars[names(dots)] <- dots

    list(move = move,
         temp = temp,
         advect = advect,
         noise = noise,
         growth = growth,
         shelf = shelf,
         scenario = scenario,
         land = land,
         boundary = boundary,
         pars = pars)
  }
