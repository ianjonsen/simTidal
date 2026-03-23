##' Map simulated drifter track
##'
##' @title map_drifter
##' @param x a fitted object of class simTidal
##' @param data data object created by sim_setup
##' @param xlim plot x limits
##' @param ylim plot y limits
##' @param res downsampling factor to product a plot faster (default = 5, full-resolution = 0)

##'
##' @importFrom ggplot2 ggplot aes coord_sf geom_point geom_path
##' @importFrom ggplot2 theme ylab xlab
##' @importFrom stars geom_stars st_as_stars
##' @export

map_drifter <- function(x,
                       data,
                       xlim = c(350,430),
                       ylim = c(4990,5030),
                       raster,
                       res = 1) {

  m <- ggplot() +
    geom_stars(data = st_as_stars(data$land)) +
    coord_sf(xlim=xlim,
             ylim=ylim) +
    theme(legend.position = "none") +
    geom_path(data = x$sim,
              aes(x,y),
              col = "orange",
              linewidth = 0.5) +
    xlab("Easting (km)") +
    ylab("Northing (km)") +
    geom_point(data = x$sim[1,],
               aes(x,y),
               shape = 17,
               colour = "blue") +
    geom_point(data = x$sim[nrow(x$sim),],
               aes(x,y),
               shape = 17,
               colour = "firebrick")

  m
}
