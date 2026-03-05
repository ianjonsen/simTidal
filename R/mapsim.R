##' Produce a simple map of simulated track(s) with or without receivers
##'
##' @title mapsim
##' @param x a fitted object of class simsmolt
##' @param data data object created by sim_setup
##' @param xlim plot x limits
##' @param ylim plot y limits
##' @param res downsampling factor to product a plot faster (default = 5, full-resolution = 0)
##' @param rec should receiver locations be displayed (logical)
##' @param alpha translucence for smolt track(s)
##' @param lwd width of smolt track(s)
##' @param col colour for receiver locations
##' @param size of smolt track end point(s)

##'
##' @importFrom ggplot2 ggplot coord_fixed aes coord_sf geom_sf guide_colourbar
##' @importFrom ggplot2 scale_fill_gradientn stat_smooth geom_line guides
##' @importFrom ggplot2 geom_path theme_minimal geom_point element_blank
##' @importFrom ggplot2 aes_string theme ylab xlab guides guide_legend xlim ylim
##' @importFrom raster extent crop nlayers crs projectRaster
##' @importFrom dplyr "%>%"
##' @importFrom patchwork wrap_plots area
##' @importFrom sf st_transform st_as_sf st_crop st_bbox st_coordinates st_geometry<-
##' @importFrom stars geom_stars st_as_stars st_contour
##' @export

mapsim <- function(x,
                   data = NULL,
                   xlim = NULL,
                   ylim = NULL,
                   res = 2,
                   rec = TRUE,
                   track = TRUE,
                   last = TRUE,
                   det = TRUE,
                   dt = NULL,
                   ol = TRUE,
                   hl = "orange",
                   tcol = "salmon",
                   alpha = 0.5,
                   lwd = 0.25,
                   reccol = "blue",
                   detcol = "red",
                   pal = "Blues 3",
                   prj = "+proj=laea +lat_0=41 +lon_0=-71 +units=km +datum=WGS84",
                   ...) {

  if(!is.null(dt) & ol == FALSE) {
    track <- FALSE
    last <- FALSE
  } else if(!is.null(dt) & ol == TRUE) {
    track <- TRUE
    last <- last
  }

  if(!is.null(dt) & !inherits(dt, "POSIXt")) stop("dt must be a POSIX class date")

  ## process simulated tracks
  if (!is.na(class(x)[2]) && (class(x)[2] == "tbl_df")) {

    ## handle multiple tracks
    compl <- sapply(x$rep, function(.) !inherits(., "try-error"))
    if(sum(!compl) > 0) cat(sprintf("dropping %i failed runs", sum(!compl)))
    x <- x[compl, ]

    if(any(sapply(x$rep, function(.)
      "detect" %in% names(.)))) {
      detect <-
        lapply(x$rep, function(.)
          .$detect) %>% do.call(rbind, .)
    } else {
      detect <- NULL
    }

    sim <- lapply(x$rep, function(.) .$sim) %>% do.call(rbind, .)
    sim.last <- lapply(x$rep, function(.) .$sim[nrow(.$sim), ]) %>% do.call(rbind, .)
    Nsim <- nrow(x)

  } else if(is.na(class(x)[2])){

    ## handle single track
    sim <- x$sim
    sim.last <- x$sim[nrow(x$sim), ]
    if("detect" %in% names(x)) detect <- x$detect
    else detect <- NULL
    Nsim <- 1
  }

  ## project sim locs using default crs
  sim <- sim %>%
    st_as_sf(., coords = c("x","y"), crs = crs(data$bathy)) %>%
    st_transform(., crs = prj)
  xy <- st_coordinates(sim) %>% as.data.frame() %>% rename(x=X, y=Y)
  sim <- bind_cols(sim, xy)
  st_geometry(sim) <- NULL

  ## project bathy
  bathy <- projectRaster(data$bathy, crs = prj)

  if (is.null(xlim))
    xlim <- c(extent(bathy)[1], extent(bathy)[2])
  if (is.null(ylim))
    ylim <- c(extent(bathy)[3], extent(bathy)[4])

  ## convert raster to stars
    bathy <- stars::st_as_stars(bathy)
#    bathy.c <- stars::st_contour(bathy, contour_lines = TRUE, breaks = c(-1000, -999.99)) # to get continental shelf without many holes

  ## project other spatial data to default prj
  esrf <- st_as_sf(data$esrf, coords = c("x","y"), crs = crs(data$bathy)) %>%
    st_transform(., crs = prj)
  recLocs <- st_as_sf(data$recLocs, coords = c("x","y"), crs = crs(data$bathy)) %>%
    st_transform(., crs = prj)

  ## generate plot
  m <- ggplot() +
    suppressWarnings(geom_stars(data = bathy, downsample = res)) +
    # geom_sf(
    #   data = bathy.c,
    #   col = "white",
    #   alpha = 0.75,
    #   lwd = 0.4
    # ) +
    scale_fill_gradientn(colours = hcl.colors(n=100, pal), 
                         na.value = grey(0.8), guide = "none") +
    theme_minimal()

  ## ESRF Oil & Gas polygon
  m <- m + geom_sf(data = esrf, col = "snow2", fill = NA, lwd = 1, alpha = 0.35)

  if(rec) {
    m <-
      m + geom_sf(
        data = recLocs,
        colour = reccol,
        size = 0.5
      )
  }

  if (track) {
    m <- m + geom_path(
      data = sim,
      aes(x, y, group = id, col = ts),
#      colour = tcol,
      alpha = alpha,
      size = lwd
    ) +
      scale_colour_gradientn(colours = hcl.colors(n=100, "Temp"),
                             na.value = "black")

    if (!is.null(detect) & det) {
      if (nrow(detect) > 0) {
        m <- m + geom_point(
          data = detect,
          aes(recv_x, recv_y),
          colour = detcol,
          size = 0.8
        )
      }
    }
  }

  ## plot overlay after track
  if(!is.null(dt)) {
    if(length(dt) == 1) {
      m <- m + geom_path(
        data = subset(sim, date > dt & date <= dt + 86400),
        aes(x, y, group = id),
        colour = hl,
        size = lwd+0.3
      )
    } else if(length(dt == 2)) {
      m <- m + geom_path(
        data = subset(sim, date > dt[1] & date <= dt[2]),
        aes(x, y, group = id),
        colour = hl,
        size = lwd+0.3
      )
    }

  }

  if(last){
    m <- m + geom_point(data = sim.last %>% filter(surv==1 & reten==1),
                 aes(x, y),
                 colour = "dodgerblue",
                 size = 0.5,
                 alpha = 1) +
      geom_point(data = sim.last %>% filter(surv==0 | reten==0),
                 aes(x, y),
                 colour = "black",
                 size = 0.5,
                 alpha = 1)
  }

  if(!is.null(prj)) {
    m <- m +
      coord_sf(crs = prj) +
      xlim(xlim[1], xlim[2]) +
      ylim(ylim[1], ylim[2])
  }


  m <- m + theme(
    axis.title = element_blank()
  ) +
    guides(colour = guide_colourbar(title = "Temp"))

  suppressWarnings(m)
}

