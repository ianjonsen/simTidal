#' @title Plot a simulated drifter or fish track on a bathymetric background
#'
#' @description Renders a \code{simTidal} object as a ggplot2 map. Bathymetry
#'   is shown as a muted fill behind the land mask, with the simulated track
#'   overlaid. The map extent is automatically fitted to the track with a
#'   small configurable margin.
#'
#' @param x          A \code{simTidal} object returned by \code{sim_drifter()}
#' @param data       The data list from \code{sim_setup()}, used to access the
#'   \code{bathy} and \code{land} SpatRasters
#' @param margin     Margin added around the track extent in km (default 5)
#' @param pt.size    Size of the start/end point markers (default 2.5)
#' @param ln.col     Colour of the track line (default \code{"#E63946"})
#' @param ln.size    Line width of the track (default 0.6)
#' @param bathy.col  Two-element character vector giving the low and high colours
#'   for the bathymetry fill scale (default a muted blue range)
#' @param bathy.alpha Alpha transparency for the bathymetry layer (default 0.6)
#'
#' @return A \code{ggplot} object that can be further modified or saved with
#'   \code{ggsave()}
#'
#' @importFrom ggplot2 ggplot aes geom_raster geom_tile geom_path geom_point
#'   scale_fill_gradientn scale_fill_manual coord_equal theme_minimal theme
#'   labs element_text element_blank element_rect margin unit guide_colourbar
#'   element_line
#' @importFrom ggnewscale new_scale_fill
#' @importFrom terra crop ext as.data.frame
#' @export

plot_sim <- function(x,
                     data,
                     margin      = 5,
                     pt.size     = 2.5,
                     ln.col      = "#E63946",
                     ln.size     = 0.6,
                     bathy.col   = c("#2a6496","#d0e8f2"),
                     bathy.alpha = 0.6) {

  if (!inherits(x, "simTidal"))
    stop("x must be a simTidal object returned by sim_drifter()")
  if (is.null(data$bathy) || is.null(data$land))
    stop("data must contain $bathy and $land SpatRasters from sim_setup()")

  sim <- x$sim

  ## ---- Compute track extent + margin ----------------------------------------
  xlim     <- range(sim$x, na.rm = TRUE) + c(-margin, margin)
  ylim     <- range(sim$y, na.rm = TRUE) + c(-margin, margin)
  crop_ext <- ext(xlim[1], xlim[2], ylim[1], ylim[2])

  ## ---- Crop rasters to plot window ------------------------------------------
  bathy_crop <- crop(data$bathy, crop_ext)
  land_crop  <- crop(data$land,  crop_ext)

  ## Bathymetry: ocean cells only (h < 0); negate so deeper = larger value
  bathy_df <- as.data.frame(bathy_crop, xy = TRUE, na.rm = TRUE)
  names(bathy_df)[3] <- "depth"
  bathy_df$depth <- bathy_df$depth * -1
  bathy_df <- bathy_df[bathy_df$depth <= 0, ]
  #bathy_df$depth <- -bathy_df$depth

  ## Land: any non-NA cell
  land_df <- as.data.frame(land_crop, xy = TRUE, na.rm = FALSE)
  names(land_df)[3] <- "val"
  land_df <- land_df[!is.na(land_df$val), ]

  ## ---- Start / end markers --------------------------------------------------
  start_pt <- sim[1, ]
  end_pt   <- sim[nrow(sim), ]

  ## ---- Assemble plot --------------------------------------------------------
  ggplot() +

    ## Muted bathymetry (ocean only) — low-saturation blues
    geom_raster(
      data = bathy_df,
      aes(x = x, y = y, fill = depth),
      alpha = bathy.alpha
    ) +
    scale_fill_gradientn(
      colours = bathy.col,
      name    = "Depth (m)",
      guide   = guide_colourbar(
        barwidth    = unit(0.4, "cm"),
        barheight   = unit(3.0, "cm"),
        title.theme = element_text(size = 7),
        label.theme = element_text(size = 6)
      )
    ) +

    ## Land mask — neutral warm grey sits quietly behind the track
    ggnewscale::new_scale_fill() +
    geom_tile(
      data = land_df,
      aes(x = x, y = y, fill = "Land"),
      colour = NA
    ) +
    scale_fill_manual(
      values = c("Land" = "#c8c0a8"),
      name   = NULL,
      guide  = "none"
    ) +

    ## Simulated track
    geom_path(
      data      = sim,
      aes(x = x, y = y),
      colour    = ln.col,
      linewidth = ln.size,
      lineend   = "round",
      linejoin  = "round"
    ) +

    ## Start point — open circle, white fill
    geom_point(
      data   = start_pt,
      aes(x = x, y = y),
      shape  = 21,
      fill   = "white",
      colour = ln.col,
      size   = pt.size,
      stroke = 0.8
    ) +

    ## End point — filled circle
    geom_point(
      data   = end_pt,
      aes(x = x, y = y),
      shape  = 21,
      fill   = ln.col,
      colour = "white",
      size   = pt.size,
      stroke = 0.8
    ) +

    ## Equal-aspect crop to track extent
    coord_equal(
      xlim   = xlim,
      ylim   = ylim,
      expand = FALSE
    ) +

    labs(
      x        = "Easting (km, UTM 20N)",
      y        = "Northing (km, UTM 20N)",
      title    = paste0("Simulated track  \u2014  id: ", sim$id[1]),
      subtitle = paste0(
        format(sim$date[1],        "%Y-%m-%d %H:%M"), "  \u2192  ",
        format(sim$date[nrow(sim)], "%Y-%m-%d %H:%M"), " UTC",
        "  |  ", nrow(sim), " steps"
      )
    ) +

    theme_minimal(base_size = 10) +
    theme(
      plot.title       = element_text(size = 10, face = "bold",
                           margin = margin(b = 2)),
      plot.subtitle    = element_text(size = 7, colour = "grey40",
                           margin = margin(b = 6)),
      axis.text        = element_text(size = 7, colour = "grey40"),
      axis.title       = element_text(size = 8, colour = "grey30"),
      panel.grid.major = element_line(colour = "white", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#eef3f7", colour = NA),
      legend.position  = "right",
      legend.key.size  = unit(0.5, "cm"),
      plot.margin      = margin(6, 6, 6, 6)
    )
}
