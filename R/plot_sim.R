#' @title Plot simulated drifter or fish tracks on a bathymetric background
#'
#' @description Renders a \code{simTidal} or \code{sim_fish} object as a
#'   ggplot2 map. Bathymetry is shown as a muted fill behind the land mask,
#'   with the simulated track(s) overlaid. The map extent is automatically
#'   fitted to the track(s) with a small configurable margin.
#'
#'   For \code{sim_fish} objects, all tracks are drawn in \code{ln.col} as a
#'   background layer, with detected (accepted) tracks overlaid in
#'   \code{ln.col.det}. Direction-of-travel arrows are added at 1/3 and 2/3
#'   of each track when \code{show.arrows = TRUE}. Detection locations are
#'   annotated as small markers when \code{show.det = TRUE}. A separate ECDF
#'   panel below the map shows the distribution of detection times; when
#'   \code{tidal.period} is set, a secondary x-axis and cycle-boundary lines
#'   show how many tidal cycles the simulated transit covered.
#'
#' @param x             A \code{simTidal} object returned by
#'   \code{sim_drifter()}, or a \code{sim_fish} object returned by
#'   \code{sim_fish()}.
#' @param data          The data list from \code{sim_setup()}, used to access
#'   the \code{bathy} and \code{land} SpatRasters.
#' @param accepted_only Logical; for \code{sim_fish} objects only. If
#'   \code{FALSE} (default), all tracks are drawn in \code{ln.col} with
#'   detected tracks overlaid in \code{ln.col.det}. If \code{TRUE}, only
#'   detected tracks are shown. Ignored for drifter objects.
#' @param show.det      Logical; if \code{TRUE} (default), a small marker is
#'   placed at the detection location of each accepted track.
#' @param show.inset    Logical; if \code{TRUE} (default), a ECDF panel is
#'   placed below the map showing the distribution of detection times.
#' @param show.arrows   Logical; if \code{TRUE} (default), direction-of-travel
#'   arrows are added at 1/3 and 2/3 along each simulated track, coloured to
#'   match the track layer. Ignored for drifter objects.
#' @param tidal.period  Tidal period in hours used to annotate the ECDF panel
#'   with a secondary x-axis and cycle-boundary lines. Default \code{12.4}
#'   (M2 semi-diurnal). Set to \code{NULL} to disable tidal cycle annotation.
#' @param margin        Margin added around the track extent in km (default 5).
#' @param pt.size       Size of the start/end receiver point markers (default
#'   2.5). Detection markers are scaled to \code{pt.size * 0.65}, track
#'   endpoints to \code{pt.size * 0.35}.
#' @param ln.col        Colour of all tracks (background layer) and receiver
#'   markers (default \code{"#E63946"}, red).
#' @param ln.col.det    Colour of detected (accepted) tracks and detection
#'   markers (default \code{"#1D6FA4"}, blue).
#' @param ln.size       Line width (default 0.6).
#' @param ln.alpha      Alpha transparency for the background track layer.
#'   Default \code{NULL} auto-sets to 0.5 for fish, 1 for drifter.
#' @param ln.alpha.det  Alpha transparency for the detected track layer
#'   (default 1).
#' @param bathy.col     Two-element character vector for the bathymetry fill
#'   scale (default muted blue range).
#' @param bathy.alpha   Alpha transparency for the bathymetry layer (default
#'   0.6).
#'
#' @return A \code{ggplot} or \code{patchwork} object saveable with
#'   \code{ggsave()}.
#'
#' @importFrom ggplot2 ggplot aes geom_raster geom_tile geom_path geom_point
#'   geom_segment geom_step geom_hline geom_vline
#'   scale_fill_gradientn scale_fill_manual scale_x_continuous scale_y_continuous
#'   sec_axis coord_equal coord_cartesian theme_minimal theme
#'   labs element_text element_blank element_rect element_line
#'   margin unit arrow
#' @importFrom ggnewscale new_scale_fill
#' @importFrom dplyr bind_rows
#' @importFrom patchwork wrap_plots
#' @importFrom terra crop ext as.data.frame
#' @export

plot_sim <- function(x,
                     data,
                     accepted_only  = FALSE,
                     show.det       = TRUE,
                     show.inset     = TRUE,
                     show.real.det  = TRUE,
                     show.arrows    = TRUE,
                     mpp            = "median",
                     mpp.sims       = "accepted",
                     mpp.col        = "white",
                     mpp.lwd        = 1,
                     tidal.period   = 12.4,
                     margin         = 5,
                     pt.size       = 2.5,
                     ln.col        = "firebrick",
                     ln.col.det    = "dodgerblue",
                     ln.size       = 0.6,
                     ln.alpha      = NULL,
                     ln.alpha.det  = 1,
                     bathy.col     = c("#2a6496","#d0e8f2"),
                     bathy.alpha   = 0.6) {

  if (!inherits(x, c("simTidal", "sim_fish")))
    stop("x must be a simTidal object (sim_drifter()) or sim_fish object (sim_fish())")
  if (is.null(data$bathy) || is.null(data$land))
    stop("data must contain $bathy and $land SpatRasters from sim_setup()")

  is_fish <- inherits(x, "sim_fish")

  ## ---- Extract track data, endpoints, and direction arrows ------------------
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

    ## Endpoint (last position) of each track
    ends_all <- if (!accepted_only)
      dplyr::bind_rows(lapply(x$sims, tail, 1L))
    else NULL
    ends_acc <- if (length(acc_idx) > 0L)
      dplyr::bind_rows(lapply(x$sims[acc_idx], tail, 1L))
    else NULL

    ## Direction arrows at 1/3 and 2/3 along each track
    if (show.arrows) {
      .mk_arrows <- function(sims_list) {
        dplyr::bind_rows(lapply(sims_list, function(tr) {
          n <- nrow(tr)
          if (n < 3L) return(NULL)
          pts <- unique(pmax(1L, pmin(n - 1L, round(c(1/3, 2/3) * n))))
          data.frame(x    = tr$x[pts],
                     y    = tr$y[pts],
                     xend = tr$x[pts + 1L],
                     yend = tr$y[pts + 1L])
        }))
      }
      arrow_df_all <- if (!accepted_only) .mk_arrows(x$sims) else NULL
      arrow_df_acc <- if (length(acc_idx) > 0L) .mk_arrows(x$sims[acc_idx]) else NULL
    } else {
      arrow_df_all <- arrow_df_acc <- NULL
    }

    ## Most probable path
    mpp_df <- if (!is.null(mpp)) {
      mpp <- match.arg(mpp, c("median", "mean"))
      mpp.sims <- match.arg(mpp.sims, c("accepted", "all"))
      src <- if (mpp.sims == "accepted" && length(acc_idx) > 0L) {
        x$sims[acc_idx]
      } else {
        if (mpp.sims == "accepted" && length(acc_idx) == 0L)
          message("mpp.sims = 'accepted' but no accepted sims; using all.")
        x$sims
      }
      .compute_mpp(src, method = mpp)
    } else NULL

    start_pt <- data.frame(x = x$params$start[1], y = x$params$start[2])
    end_pt   <- data.frame(x = x$params$end[1],   y = x$params$end[2])
    if (is.null(ln.alpha)) ln.alpha <- 0.5

  } else {

    sim_all      <- NULL
    sim_acc      <- x$sim
    sim_ext      <- x$sim
    ends_all     <- NULL
    ends_acc     <- NULL
    arrow_df_all <- NULL
    arrow_df_acc <- NULL
    start_pt     <- x$sim[1, ]
    end_pt       <- x$sim[nrow(x$sim), ]
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
    n_acc <- length(acc_idx)
    n_rej <- length(rej_idx)
    title    <- paste0("sim_fish  —  move: ", x$params$move,
                       "  |  n = ", x$params$n_sim)
    subtitle <- paste0(
      format(x$params$start.dt, "%Y-%m-%d %H:%M"), "  ->  ",
      format(x$params$end.dt,   "%Y-%m-%d %H:%M"), " UTC",
      "  |  ", x$params$N, " steps",
      if (!is.null(tidal.period) && tidal.period > 0) {
        n_cyc <- x$params$N * x$params$time.step / 60 / tidal.period
        sprintf("  |  %.2g tidal cycles", n_cyc)
      },
      if (accepted_only)
        paste0("  |  ", n_acc, " track", if (n_acc != 1L) "s",
               " accepted (", round(100 * x$acceptance_rate, 1), "%)")
      else
        paste0("  |  ", n_acc, " accepted  /  ", n_rej, " rejected",
               "  (", round(100 * x$acceptance_rate, 1), "%)")
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
      guide   = "none"
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

    ## All tracks — background layer in ln.col
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

    ## Detected / accepted tracks — foreground in ln.col.det (or drifter)
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

    ## Direction arrows — all tracks (red)
    { if (!is.null(arrow_df_all) && nrow(arrow_df_all) > 0L)
        geom_segment(
          data      = arrow_df_all,
          aes(x = x, y = y, xend = xend, yend = yend),
          colour    = ln.col,
          linewidth = ln.size * 0.7,
          alpha     = ln.alpha,
          arrow     = arrow(length = unit(0.1, "cm"), type = "open")
        )
      else NULL } +

    ## Direction arrows — accepted tracks (blue)
    { if (!is.null(arrow_df_acc) && nrow(arrow_df_acc) > 0L)
        geom_segment(
          data      = arrow_df_acc,
          aes(x = x, y = y, xend = xend, yend = yend),
          colour    = ln.col.det,
          linewidth = ln.size * 0.7,
          alpha     = ln.alpha.det,
          arrow     = arrow(length = unit(0.1, "cm"), type = "open")
        )
      else NULL } +

    ## Track endpoints — all tracks (red), accepted on top (blue)
    { if (!is.null(ends_all))
        geom_point(
          data   = ends_all,
          aes(x = x, y = y),
          shape  = 16,
          colour = ln.col,
          size   = pt.size * 0.35,
          alpha  = ln.alpha
        )
      else NULL } +
    { if (!is.null(ends_acc))
        geom_point(
          data   = ends_acc,
          aes(x = x, y = y),
          shape  = 16,
          colour = ln.col.det,
          size   = pt.size * 0.35,
          alpha  = ln.alpha.det
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

    ## Most probable path — topmost map layer (above all symbols)
    { if (is_fish && !is.null(mpp_df))
        geom_path(
          data      = mpp_df,
          aes(x = x, y = y),
          colour    = mpp.col,
          linewidth = mpp.lwd,
          lineend   = "round",
          linejoin  = "round"
        )
      else NULL } +
    { if (is_fish && !is.null(mpp_df)) {
        mpp_end <- mpp_df[nrow(mpp_df), , drop = FALSE]
        geom_point(
          data   = mpp_end,
          aes(x = x, y = y),
          shape  = 23,
          fill   = mpp.col,
          colour = "black",
          size   = pt.size,
          stroke = 0.8
        )
      } else NULL } +

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

  ## ---- Detection-time ECDF panel (sim_fish only) ---------------------------
  if (is_fish && show.inset && length(acc_idx) > 0L) {

    det_frac   <- x$det_step[acc_idx] / x$params$N
    det_sorted <- sort(det_frac)
    n_det      <- length(det_sorted)
    ecdf_df    <- data.frame(
      frac = c(0, det_sorted),
      prop = c(0, seq_len(n_det) / n_det)
    )

    ## Tidal cycle info for secondary axis and cycle-boundary lines
    total_hrs  <- x$params$N * x$params$time.step / 60
    use_tidal  <- !is.null(tidal.period) && is.numeric(tidal.period) &&
                  tidal.period > 0 && total_hrs > 0

    if (use_tidal) {
      n_cycles      <- total_hrs / tidal.period
      cycle_bounds  <- seq_len(floor(n_cycles)) / n_cycles  ## fractions of N
      x_scale <- scale_x_continuous(
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        labels = c("0", "0.25", "0.5", "0.75", "1"),
        sec.axis = sec_axis(
          ~ . * n_cycles,
          name   = sprintf("Tidal cycles (%.4g h period)", tidal.period),
          breaks = seq(0, ceiling(n_cycles), by = 1)
        )
      )
    } else {
      cycle_bounds <- numeric(0)
      x_scale <- scale_x_continuous(
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        labels = c("0", "0.25", "0.5", "0.75", "1")
      )
    }

    ## Real-world detection fraction: end.dt expressed as a fraction of N steps.
    ## N = round(elapsed / step_secs) so end.dt may not sit exactly at 1.0.
    real_det_frac <- if (show.real.det && is_fish) {
      elapsed_secs <- as.numeric(
        difftime(x$params$end.dt, x$params$start.dt, units = "secs"))
      min(1, elapsed_secs / (x$params$N * x$params$time.step * 60))
    } else NULL

    ecdf_p <- ggplot(ecdf_df, aes(x = frac, y = prop)) +
      { if (length(cycle_bounds) > 0L)
          geom_vline(xintercept = cycle_bounds, linetype = "dotted",
                     colour = "grey55", linewidth = 0.4)
        else NULL } +
      geom_hline(yintercept = 0.5, linetype = "dashed",
                 colour = "grey70", linewidth = 0.3) +
      { if (!is.null(real_det_frac))
          geom_vline(xintercept = real_det_frac,
                     colour = ln.col.det, linewidth = 0.7, linetype = "solid")
        else NULL } +
      geom_step(colour = ln.col.det, linewidth = 0.7) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      x_scale +
      scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
      labs(x = "Detection time (step / N)", y = "Cumulative\nproportion") +
      theme_minimal(base_size = 9) +
      theme(
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.text        = element_text(size = 7, colour = "grey40"),
        axis.title       = element_text(size = 8, colour = "grey30"),
        plot.margin      = margin(2, 6, 4, 6)
      )

    return(patchwork::wrap_plots(p, ecdf_p, ncol = 1, heights = c(4, 1)))
  }

  p
}
