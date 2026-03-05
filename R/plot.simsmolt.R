##' Map simulated track
##'
##' @title plot
##' @param x a fitted object of class simsmolt
##' @param data data object created by sim_setup
##' @param xlim plot x limits
##' @param ylim plot y limits
##' @param raster select a raster (from data) to use as background (bathy, uv, or NULL)
##' @param res downsampling factor to product a plot faster (default = 5, full-resolution = 0)
##' @param rec should receiver locations be displayed (logical)
##' @param track should smolt track(s) be displayed (logical)
##' @param last display position at end of simulation only (logical; if T then track forced to be F)
##' @param alpha translucence for smolt track(s)
##' @param lwd width of smolt track(s)
##' @param col colour for receiver locations
##' @param size of smolt track end point(s)

##'
##' @importFrom ggplot2 ggplot coord_fixed aes coord_sf geom_sf
##' @importFrom ggplot2 scale_fill_gradientn stat_smooth geom_line
##' @importFrom ggplot2 geom_path theme_minimal geom_point scale_colour_manual
##' @importFrom ggplot2 aes_string theme ylab xlab guides guide_legend
##' @importFrom raster extent crop nlayers
##' @importFrom dplyr "%>%" mutate
##' @importFrom patchwork wrap_plots area
##' @importFrom sf st_transform st_as_sf st_crop st_bbox
##' @importFrom stars geom_stars st_as_stars st_contour
##' @importFrom wesanderson wes_palette
##' @method plot simsmolt
##' @export

plot.simsmolt <- function(x, data = NULL, xlim = NULL, ylim = NULL, 
                          raster = "ts", res = 5, esrf = TRUE, layer = NULL, 
                          rec = FALSE, track = TRUE, last = FALSE,
                          ts.plots = FALSE, 
                          alpha = 0.5, lwd = 0.2, size = 0.2, 
                          col = c("salmon", "red", "dodgerblue"), 
                          pal = "Light Grays", maxp = 5e4,
                          landfill = grey(0.6),
                          crs = "+proj=laea +lat_0=41 +lon_0=-71 +units=km +datum=WGS84",
                          ...) {
  
  ## process simulated tracks
  if (!is.na(class(x)[2]) && (class(x)[2] == "tbl_df")) {
    
    ## handle multiple tracks
    compl <- sapply(x$rep, function(.) !inherits(., "try-error"))
    if(sum(!compl) > 0) cat(sprintf("dropping %i failed runs", sum(!compl)))
    x <- x[compl, ]
    
    detect <-
      lapply(x$rep, function(.)
        .$detect) %>% do.call(rbind, .)
    Nsim <- nrow(x)
    sim <- lapply(x$rep, function(.) .$sim) %>% do.call(rbind, .)
    sim.last <- lapply(x$rep, function(.) .$sim[nrow(.$sim), ]) %>% do.call(rbind, .)
    hg.ids <- which(sim.last$fl > quantile(sim.last$fl, 0.85))
    lg.ids <- which(sim.last$fl < quantile(sim.last$fl, 0.15))
    sim <- sim %>% 
      mutate(growth = ifelse(id %in% hg.ids, "high", ifelse(id %in% lg.ids, "low", "intermediate")))
    
  } else if(is.na(class(x)[2])){
    Nsim <- 1
    ## handle single track
    sim <- x$sim
    sim.last <- x$sim[nrow(x$sim), ] 
    detect <- x$detect
    
  }
  
  wespal <- wes_palette("Darjeeling2", type = "discrete")
  
  if(!is.null(data)) {
    bathy <- stars::st_as_stars(data$bathy)
    bathy.c <- stars::st_contour(bathy, contour_lines = TRUE, breaks = c(-1000, -999.99)) # to get continental shelf without many holes
  }
    
    if (!is.null(raster)) {
      switch(raster,
             uv = {
               if (nlayers(data$u) > 1) {
                 u <- data$u[[ifelse(!is.null(layer), layer, nlayers(data$u))]]
                 v <- data$v[[ifelse(!is.null(layer), layer, nlayers(data$u))]]
                 ras <- sqrt(u ^ 2 + v ^ 2) * 3.6
               } else {
                 ras <- sqrt(data$u ^ 2 + data$v ^ 2) * 3.6
               }
               
             },
             ts = {
               if (nlayers(data$ts) > 1) {
                 ras <- data$ts[[ifelse(!is.null(layer), layer, nlayers(data$u))]] - 273
               } else {
                 ras <- data$ts - 273
               }
             })
      if (inherits(ras, "RasterLayer")) {
        ras <- stars::st_as_stars(ras)
      }
    } else {
      ras <- bathy
    }
  
  if (is.null(xlim))
    xlim <- c(extent(data$bathy)[1], extent(data$bathy)[2])
  if (is.null(ylim))
    ylim <- c(extent(data$bathy)[3], extent(data$bathy)[4])
  
  ## get coastline
#  coast <- sf::st_as_sf(rworldmap::countriesLow) %>%
#    st_transform(crs = crs)
  
  ## generate plot
    m <- ggplot() +
      geom_stars(data = ras, downsample = res) + 
      geom_sf(
        data = bathy.c,
        col = "white",
        alpha = 0.25,
        lwd = 0.4
      ) +
      scale_fill_gradientn(colours = hcl.colors(n=100, pal), na.value = "transparent") +
      theme_minimal() + 
      guides(fill = guide_legend(title = ifelse(raster == "ts", "T ºC", "||uv|| km/h")))
  
      if(esrf) {
        m <- m + geom_sf(data = data$esrfPoly, colour="snow2", fill=NA, lwd=1)
      }
    
  ## add land  
#    m <- m + geom_sf(data = coast, 
#            fill = landfill,
#            lwd = 0
#            )
  
  if(track & !last) {
    if(Nsim > 1) {
      m <- m + geom_path(data = sim,
                         aes(x, y, group=id, colour = growth),
                         alpha = alpha, 
                         size = lwd)
    } else {
      m <- m + geom_path(data = sim,
                         aes(x, y, group=id),
                         colour = "salmon",
                         alpha = alpha, 
                         size = lwd)
    }
   
    m <- m + geom_point(data = sim.last, 
               aes(x, y),
               colour = col[2],
               size = size,
               alpha = 1) +
    scale_colour_manual(values = wespal[c(4,3,2)])
    
  } else if(last) {
     m <- m + geom_point(data = sim.last, 
                         aes(x, y),
                         colour = col[2],
                         size = size + 0.5,
                         alpha = 1)
   }
  
  if(rec) {
      m <-
        m + geom_point(
          data = data$recLocs,
          aes(x, y),
          colour = col[3],
          size = 0.4
        ) 
    }
    
    
  if (rec && !is.null(detect) && nrow(detect) > 0) {
    m <-
      m + geom_point(data = detect, 
                     aes(recv_x, recv_y), 
                     colour = col[3],
                     size = 0.6)
    
  }
  if(!is.null(crs)) {
    m <- m + 
      coord_sf(crs = crs) +
      xlim(xlim[1], xlim[2]) +
      ylim(ylim[1], ylim[2])
  }
    
  
  m <- m + theme(
    axis.title = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(1, 0.8)
  ) +
    guides(colour = "none")
  
  if (ts.plots) {
    ## generate time-series plots
    env24 <- lapply(x$rep, function(.) {
      .$sim[seq(1, nrow(.$sim), by = 24),]
    }) %>%
      do.call(rbind, .)
    
    env24 <- env24 %>%
      mutate(growth = ifelse(
        id %in% hg.ids,
        "high",
        ifelse(id %in% lg.ids, "low", "intermediate")
      ))
    
    ## Mass gain/loss daily time-series
    wp <- ggplot(env24, aes(date, w, group = id, colour = growth)) +
      geom_line(lwd = 0.5,
                alpha = 0.6) +
      theme_minimal() +
      labs(title = "Mass (g)") +
      theme(axis.title = element_blank()) +
      scale_colour_manual(values = wespal[c(4, 3, 2)])
    
    ## Fork-length daily time-series
    fp <- ggplot(env24, aes(date, fl, group = id, colour = growth)) +
      geom_line(lwd = 0.5,
                alpha = 0.6) +
      theme_minimal() +
      labs(title = "Fork-length (m)") +
      theme(axis.title = element_blank()) +
      scale_colour_manual(values = wespal[c(4, 3, 2)])
    
    ## Temperature daily time-series
    tp <- ggplot(env24, aes(date, ts, group = id, colour = growth)) +
      geom_line(lwd = 0.5,
                alpha = 0.6) +
      theme_minimal() +
      labs(title = "Temperature (ºC)") +
      theme(axis.title = element_blank()) +
      scale_colour_manual(values = wespal[c(4, 3, 2)])
    
    
    ## daily displacement time-series
    sim24 <- lapply(x$rep, function(.) {
      .$sim %>%
        mutate(yday = yday(date)) %>%
        group_by(id, yday) %>%
        summarise(
          "swimming" = sqrt(sum(dx) ^ 2 + sum(dy) ^ 2),
          "current" = sqrt(sum(u) ^ 2 + sum(v) ^ 2),
          "total" = sqrt(sum(dx + u) ^ 2 + sum(dy + v) ^ 2),
          .groups = "drop"
        )
    }) %>%
      do.call(rbind, .) %>%
      mutate(date = as.POSIXct(as.Date("2018-01-01") + yday)) %>%
      pivot_longer(., cols = 3:5, names_to = "disp")
    
    ## combined displacement plots
    dp <- ggplot(sim24, aes(date, value, group = id, col = disp)) +
      geom_point(size = 0.2, alpha = 0.3) +
      geom_smooth(
        aes(group = disp),
        method = "gam",
        formula = y ~ s(x),
        alpha = 0.3
      ) +
      theme_minimal() +
      scale_colour_manual(values = wespal[c(3, 4, 2)],
                          name = element_blank()) +
      labs(title = "Daily displacement (km)") +
      theme(axis.title = element_blank())
    
    layout <- c(
      area(
        t = 1,
        l = 1,
        b = 4,
        r = 4
      ),
      area(
        t = 1,
        l = 5,
        b = 2,
        r = 6
      ),
      area(
        t = 3,
        l = 5,
        b = 4,
        r = 6
      ),
      area(
        t = 5,
        l = 5,
        b = 6,
        r = 6
      ),
      area(
        t = 5,
        l = 1,
        b = 6,
        r = 4
      )
    )
    
    wrap_plots(m, wp, fp, tp, dp, design = layout)
  }
  
  if(!ts.plots) return(m)
  
}
