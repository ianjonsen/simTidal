#' @title Plot simulated drifter or fish tracks on a bathymetric background
#'
#' @description Renders a \code{simTidal} or \code{sim_fish} object as a
#'   ggplot2 map. Bathymetry is shown as a muted fill behind the land mask,
#'   with the simulated track(s) overlaid. The map extent is automatically
#'   fitted to the track(s) with a small configurable margin.
#'
#'   For \code{sim_fish} objects, all tracks are drawn in \code{ln.col} as a
#'   background layer, with detected (accepted) tracks overlaid in
#'   \code{ln.col.det}. Detection locations are annotated as small markers
#'   when \code{show.det = TRUE}. An inset ECDF showing the distribution of
#'   detection times (as a fraction of the total transit \code{N}) is added
#'   to the lower-right corner when \code{show.inset = TRUE}.
#'
#' @param x           A \code{simTidal} object returned by
#'   \code{sim_drifter()}, or a \code{sim_fish} object returned by
#'   \code{sim_fish()}.
#' @param data        The data list from \code{sim_setup()}, used to access
#'   the \code{bathy} and \code{land} SpatRasters.
#' @param accepted_only Logical; for \code{sim_fish} objects only. If
#'   \code{FALSE} (default), all tracks are drawn in \code{ln.col} with
#'   detected tracks overlaid in \code{ln.col.det}. If \code{TRUE}, only
#'   detected tracks are shown. Ignored for drifter objects.
#' @param show.det    Logical; for \code{sim_fish} objects only. If
#'   \code{TRUE} (default), a small marker is placed at the detection
#'   location of each accepted track. Ignored for drifter objects.
#' @param show.inset  Logical; for \code{sim_fish} objects with at least one
#'   accepted simulation. If \code{TRUE} (default), a compact ECDF inset is
#'   placed in the lower-right corner showing detection times as a fraction
#'   of total transit steps (\code{det_step / N}). Ignored for drifter objects.
#' @param margin      Margin added around the track extent in km (default 5).
#' @param pt.size     Size of the start/end receiver point markers (default
#'   2.5). Detection markers are scaled to \code{pt.size * 0.65}.
#' @param ln.col      Colour of all tracks (background layer) and receiver
#'   markers (default \code{"#E63946"}, red).
#' @param ln.col.det  Colour of detected (accepted) tracks and detection
#'   markers (default \code{"#1D6FA4"}, blue). Ignored for drifter objects.
#' @param ln.size     Line width (default 0.6).
#' @param ln.alpha    Alpha transparency for the background (all-tracks) layer.
#'   Default \code{NULL} sets it automatically: 0.5 for fish ensembles, 1 for
#'   a single drifter track.
#' @param ln.alpha.det Alpha transparency for the detected (blue) track layer
#'   (default 1). Ignored for drifter objects.
#' @param bathy.col   Two-element character vector giving the low and high
#'   colours for the bathymetry fill scale (default a muted blue range).
#' @param bathy.alpha Alpha transparency for the bathymetry layer (default
#'   0.6).
#'
#' @return A \code{ggplot} object that can be further modified or saved with
#'   \code{ggsave()}.
#'
#' @importFrom ggplot2 ggplot aes geom_raster geom_tile geom_path geom_point
#'   geom_step scale_fill_gradientn scale_fill_manual scale_x_continuous
#'   scale_y_continuous coord_equal coord_cartesian theme_minimal theme
#'   annotation_custom ggplotGrob labs element_text element_blank element_rect
#'   element_line margin unit guide_colourbar
#' @importFrom ggnewscale new_scale_fill
#' @importFrom dplyr bind_rows
#' @importFrom terra crop ext as.data.frame
#' @export

plot_sim <- function(x,
                     data,
                     accepted_only = FALSE,
                     show.det      = TRUE,
                     show.inset    = TRUE,
                     margin        = 5,
                     pt.size       = 2.5,
                     ln.col        = "#E63946",
                     ln.col.det    = "#1D6FA4",
                     ln.size       = 0.4,
                     ln.alpha      = NULL,
                     ln.alpha.det  = 1,
                     bathy.col     = c("#2a6496","#d0e8f2"),
                     bathy.alpha   = 0.6) {

  if (!inherits(x, c("simTidal", "sim_fish")))
    stop("x must be a simTidal object (sim_drifter()) or sim_fish object (sim_fish())")
  if (is.null(data$bathy) || is.null(data$land))
    stop("data must contain $bathy and $land SpatRasters from sim_setup()")

  is_fish <- inherits(x, "sim_fish")

  ## ---- Extract track data and start/end markers -----------------------------
  if (is_fish) {

    acc_idx <- which(x$accepted)
    rej_idx <- which(!x$accepted)

    if (accepted_only) {
      if (length(acc_idx) == 0L)
        stop("No accepted simulations to plot. Use accepted_only = FALSE to plot all.")
      sim_all <- NULL
      sim_acc <- dplyr::bind_rows(x$sims[acc_idx])
      sim_ext <- sim_acc
    } else {
      sim_all <- dplyr::bind_rows(x$sims)
      sim_acc <- if (length(acc_idx) > 0L) dplyr::bind_rows(x$sims[acc_idx]) else NULL
      sim_ext <- sim_all
    }

    start_pt <- data.frame(x = x$params$start[1], y = x$params$start[2])
    end_pt   <- data.frame(x = x$params$end[1],   y = x$params$end[2])
    if (is.null(ln.alpha)) ln.alpha <- 0.5

  } else {

    sim_all  <- NULL
    sim_acc  <- x$sim
    sim_ext  <- x$sim
    start_pt <- x$sim[1, ]
    end_pt   <- x$sim[nrow(x$sim), ]
    if (is.null(ln.alpha)) ln.alpha <- 1

  }

  ## ---- Compute track extent + margin ----------------------------------------
  xlim     <- range(sim_ext$x, na.rm = TRUE) + c(-margin, margin)
  ylim     <- range(sim_ext$y, na.rm = TRUE) + c(-margin, margin)
  crop_ext <- ext(xlim[1], xlim[2], ylim[1], ylim[2])

  ## ---- Crop rasters to plot window ------------------------------------------
  bathy_crop <- crop(data$bathy, crop_ext)
  land_crop  <- crop(data$land,  crop_ext)

  bathy_df <- as.data.frame(bathy_crop, xy = TRUE, na.rm = TRUE)
  names(bathy_df)[3] <- "depth"
  bathy_df$depth <- bathy_df$depth * -1
  bathy_df <- bathy_df[bathy_df$depth <= 0, ]

  land_df <- as.data.frame(land_crop, xy = TRUE, na.rm = FALSE)
  names(land_df)[3] <- "val"
  land_df <- land_df[!is.na(land_df$val), ]

  ## ---- Build title and subtitle ---------------------------------------------
  if (is_fish) {
    n_acc    <- length(acc_idx)
    n_rej    <- length(rej_idx)
    title    <- paste0("sim_fish  —  move: ", x$params$move,
                       "  |  n = ", x$params$n_sim)
    subtitle <- paste0(
      format(x$params$start.dt, "%Y-%m-%d %H:%M"), "  ->  ",
      format(x$params$end.dt,   "%Y-%m-%d %H:%M"), " UTC",
      "  |  ", x$params$N, " steps",
      if (accepted_only) {
        paste0("  |  ", n_acc, " track", if (n_acc != 1L) "s",
               " accepted (", round(100 * x$acceptance_rate, 1), "%)")
      } else {
        paste0("  |  ", n_acc, " accepted  /  ", n_rej, " rejected",
               "  (", round(100 * x$acceptance_rate, 1), "%)")
      }
    )
  } else {
    title    <- paste0("Simulated track  —  id: ", sim_acc$id[1])
    subtitle <- paste0(
      format(sim_acc$date[1],             "%Y-%m-%d %H:%M"), "  ->  ",
      format(sim_acc$date[nrow(sim_acc)], "%Y-%m-%d %H:%M"), " UTC",
      "  |  ", nrow(sim_acc), " steps"
    )
  }

  ## ---- Assemble main plot ---------------------------------------------------
  p <- ggplot() +

    geom_raster(
      data  = bathy_df,
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

    ggnewscale::new_scale_fill() +
    geom_tile(
      data   = land_df,
      aes(x = x, y = y, fill = "Land"),
      colour = NA
    ) +
    scale_fill_manual(
      values = c("Land" = "#c8c0a8"),
      name   = NULL,
      guide  = "none"
    ) +

    ## All tracks — background layer in ln.col (skipped when accepted_only)
    { if (!is.null(sim_all))
        geom_path(
          data      = sim_all,
          aes(x = x, y = y, group = factor(id)),
          colour    = ln.col,
          linewidth = ln.size,
          alpha     = ln.alpha,
          lineend   = "round",
          linejoin  = "round"
        )
      else NULL } +

    ## Detected / accepted tracks — foreground in ln.col.det (or drifter track)
    { if (!is.null(sim_acc))
        geom_path(
          data      = sim_acc,
          aes(x = x, y = y,
              group = if (is_fish) factor(id) else NULL),
          colour    = if (is_fish) ln.col.det else ln.col,
          linewidth = ln.size,
          alpha     = if (is_fish) ln.alpha.det else ln.alpha,
          lineend   = "round",
          linejoin  = "round"
        )
      else NULL } +

    ## Detection location markers
    { if (is_fish && show.det && nrow(x$det_locs) > 0L)
        geom_point(
          data   = x$det_locs,
          aes(x = x, y = y),
          shape  = 21,
          fill   = ln.col.det,
          colour = "white",
          size   = pt.size * 0.65,
          stroke = 0.7
        )
      else NULL } +

    ## Start receiver
    geom_point(
      data   = start_pt,
      aes(x = x, y = y),
      shape  = 21,
      fill   = "white",
      colour = ln.col,
      size   = pt.size,
      stroke = 0.8
    ) +

    ## End receiver
    geom_point(
      data   = end_pt,
      aes(x = x, y = y),
      shape  = 21,
      fill   = ln.col,
      colour = "white",
      size   = pt.size,
      stroke = 0.8
    ) +

    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +

    labs(
      x        = "Easting (km, UTM 20N)",
      y        = "Northing (km, UTM 20N)",
      title    = title,
      subtitle = subtitle
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

  ## ---- Detection-time inset (sim_fish only) ---------------------------------
  if (is_fish && show.inset && length(acc_idx) > 0L) {

    det_frac   <- x$det_step[acc_idx] / x$params$N
    det_sorted <- sort(det_frac)
    n_det      <- length(det_sorted)
    ecdf_df    <- data.frame(
      frac = c(0, det_sorted),
      prop = c(0, seq_len(n_det) / n_det)
    )

    inset_p <- ggplot(ecdf_df, aes(x = frac, y = prop)) +
      geom_step(colour = ln.col.det, linewidth = 0.6) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
      scale_y_continuous(breaks = c(0, 1),       labels = c("0", "1")) +
      labs(x = "step / N", y = "prop.") +
      theme_minimal(base_size = 6) +
      theme(
        plot.background  = element_rect(fill     = "white",
                                        colour   = "grey75",
                                        linewidth = 0.4),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        axis.text        = element_text(size = 5, colour = "grey40"),
        axis.title.x     = element_text(size = 5.5, colour = "grey30",
                             margin = margin(t = 2)),
        axis.title.y     = element_text(size = 5.5, colour = "grey30",
                             margin = margin(r = 2)),
        plot.margin      = margin(3, 4, 3, 3)
      )

    iw <- (xlim[2] - xlim[1]) * 0.36
    ih <- (ylim[2] - ylim[1]) * 0.28

    p <- p + annotation_custom(
      grob = ggplotGrob(inset_p),
      xmin = xlim[2] - iw,
      xmax = xlim[2],
      ymin = ylim[1],
      ymax = ylim[1] + ih
    )
  }

  p
}
