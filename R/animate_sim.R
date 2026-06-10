#' @title Animate simulated fish track ensemble with fading wake
#'
#' @description Renders a \code{sim_fish} object as an animated GIF (or other
#'   format) using \pkg{gganimate}. At each frame only the most recent
#'   \code{wake_steps} path segments per track are drawn, with older segments
#'   fading linearly to transparent (comet-tail effect). Rejected tracks appear
#'   in \code{ln.col} and accepted tracks in \code{ln.col.det}. Detection
#'   markers pulse large at the moment of detection and decay to normal size
#'   over \code{pulse_frames} steps, then persist for the rest of the
#'   animation. When \code{show.ecdf = TRUE} an ECDF inset is overlaid on the
#'   map, growing step-by-step as detections accumulate.
#'
#' @param x       A \code{sim_fish} object returned by \code{sim_fish()}.
#' @param data    The data list from \code{sim_setup()}, used to access the
#'   \code{bathy} and \code{land} SpatRasters.
#' @param accepted_only Logical; animate only accepted (detected) tracks.
#'   Default \code{FALSE}.
#' @param show.det    Logical; show detection markers with pulse effect.
#'   Default \code{TRUE}.
#' @param show.ecdf   Logical; overlay a cumulative detection ECDF inset on
#'   the map. Default \code{TRUE}. Silently ignored when there are no
#'   accepted simulations.
#' @param ecdf.pos    Corner of the map for the ECDF inset. One of
#'   \code{"bottomright"} (default), \code{"bottomleft"}, \code{"topright"},
#'   or \code{"topleft"}.
#' @param wake_steps  Number of past steps kept in the wake (default 20).
#'   Current-head segment is drawn at full opacity; each earlier step fades
#'   linearly, reaching zero at \code{wake_steps} behind the head. Capped at
#'   \code{N - 1}.
#' @param pulse_frames Number of animation frames over which the detection
#'   symbol decays from peak size to normal size (default 6). At 10 fps the
#'   default gives a ~0.6 s pulse.
#' @param fps     Frames per second (default 10).
#' @param nframes Number of frames to render. Default \code{NULL} uses
#'   \code{mpar$N} (one frame per simulation step). Reduce to speed up
#'   rendering at the cost of temporal resolution.
#' @param width   Output width in pixels (default 800).
#' @param height  Output height in pixels (default 600).
#' @param renderer A \pkg{gganimate} renderer. Default
#'   \code{gganimate::gifski_renderer()} for GIF (requires \pkg{gifski}).
#'   Use \code{gganimate::av_renderer()} for MP4 (requires \pkg{av}).
#' @param margin  Margin around the track extent in km (default 5).
#' @param pt.size Size of start/end receiver markers (default 2.5). Detection
#'   markers pulse from \code{pt.size * 2.5} at detection down to
#'   \code{pt.size * 0.65} at rest.
#' @param ln.col      Colour for rejected tracks (default \code{"#E63946"}).
#' @param ln.col.det  Colour for accepted tracks and detection markers
#'   (default \code{"#1D6FA4"}).
#' @param ln.size     Line width (default 0.6).
#' @param ln.alpha    Peak alpha for rejected track wake (default 0.5).
#' @param ln.alpha.det Peak alpha for accepted track wake (default 1.0).
#' @param bathy.col   Two-element colour vector for bathymetry fill scale.
#' @param bathy.alpha Alpha for the bathymetry layer (default 0.6).
#'
#' @return A \code{magick-image} animation object. Display with \code{print()},
#'   save with \code{gganimate::anim_save("out.gif", animation = .)}.
#'
#' @details
#'   \strong{Wake:} each segment (step \eqn{t-1 \to t}) is replicated in the
#'   data frame at frames \eqn{t, t+1, \ldots, t+\text{wake\_steps}-1} with
#'   alpha \eqn{\alpha_{\text{base}} \times (1 - a/\text{wake\_steps})} where
#'   \eqn{a} is the age (0 = head). \code{transition_manual} is used so static
#'   layers (bathymetry, land, receivers) are never affected.
#'
#'   \strong{ECDF inset:} drawn entirely in map-km coordinates so it is
#'   compatible with \code{transition_manual} (patchwork cannot be used with
#'   gganimate). At frame \eqn{f} the ECDF shows all detections at steps
#'   \eqn{\le f}. The x-axis spans 0 to \eqn{N} (simulation steps) and the
#'   y-axis spans 0 to 1 (cumulative proportion of total accepted).
#'
#'   \strong{Performance:} the wake data frame has approximately
#'   \eqn{n_\text{sim} \times N \times \text{wake\_steps}} rows. For
#'   large ensembles reduce \code{nframes} or subsample sims beforehand.
#'
#' @seealso \code{\link{plot_sim}}, \code{\link{sim_fish}}
#'
#' @importFrom ggplot2 ggplot aes geom_raster geom_tile geom_segment geom_step
#'   geom_point geom_text annotate
#'   scale_fill_gradientn scale_fill_manual scale_alpha_identity scale_size_identity
#'   coord_equal theme_minimal theme labs
#'   element_text element_blank element_rect element_line margin unit
#' @importFrom ggnewscale new_scale_fill
#' @importFrom dplyr bind_rows
#' @importFrom terra crop ext as.data.frame
#' @importFrom gganimate transition_manual animate gifski_renderer
#' @export

animate_sim <- function(x, data,
                        accepted_only  = FALSE,
                        show.det       = TRUE,
                        show.ecdf      = TRUE,
                        show.real.det  = TRUE,
                        ecdf.pos       = "bottomright",
                        mpp            = "median",
                        mpp.sims       = "accepted",
                        mpp.col        = "white",
                        mpp.lwd        = 1,
                        wake_steps    = 20L,
                        pulse_frames  = 6L,
                        fps           = 10,
                        nframes       = NULL,
                        width         = 800,
                        height        = 600,
                        renderer      = gganimate::gifski_renderer(),
                        margin        = 5,
                        pt.size       = 2.5,
                        ln.col        = "firebrick",
                        ln.col.det    = "dodgerblue",
                        ln.size       = 0.6,
                        ln.alpha      = 0.5,
                        ln.alpha.det  = 1.0,
                        bathy.col     = c("#2a6496", "#d0e8f2"),
                        bathy.alpha   = 0.6) {

  if (!inherits(x, "sim_fish"))
    stop("x must be a sim_fish object from sim_fish()", call. = FALSE)
  if (is.null(data$bathy) || is.null(data$land))
    stop("data must contain $bathy and $land SpatRasters from sim_setup()", call. = FALSE)

  N            <- x$params$N
  nsim         <- x$params$n_sim
  acc_idx      <- which(x$accepted)
  n_acc_total  <- length(acc_idx)
  wake_steps   <- min(as.integer(wake_steps),  N - 1L)
  pulse_frames <- max(1L, as.integer(pulse_frames))
  ecdf.pos     <- match.arg(ecdf.pos, c("bottomright","bottomleft","topright","topleft"))

  ## ---- Select simulations to animate -----------------------------------------
  sim_ids <- if (accepted_only) {
    if (n_acc_total == 0L)
      stop("No accepted simulations. Use accepted_only = FALSE.", call. = FALSE)
    acc_idx
  } else {
    seq_len(nsim)
  }

  ## ---- Wake segment data frame -----------------------------------------------
  ## Each segment (step t-1 -> t) is replicated for frames t .. t+wake_steps-1.
  ## seg_alpha = base_alpha * (1 - age / wake_steps), age = 0 at head.
  est_rows <- length(sim_ids) * (N - 1L) * wake_steps
  if (est_rows > 500000L)
    message(sprintf(
      "Building wake data: ~%.1fM rows (nsim = %d, N = %d, wake = %d). This may take a moment.",
      est_rows / 1e6, length(sim_ids), N, wake_steps))

  wake_df <- dplyr::bind_rows(lapply(sim_ids, function(i) {
    tr <- x$sims[[i]]
    n  <- nrow(tr)
    if (n < 2L) return(NULL)

    steps    <- seq(2L, n)
    ns       <- length(steps)
    step_rep <- rep(steps, each = wake_steps)
    age_rep  <- rep(seq(0L, wake_steps - 1L), times = ns)
    frame_v  <- step_rep + age_rep
    valid    <- frame_v <= N

    base_a <- if (x$accepted[i]) ln.alpha.det else ln.alpha

    data.frame(
      sim       = i,
      frame     = frame_v[valid],
      x         = tr$x[(step_rep - 1L)[valid]],
      y         = tr$y[(step_rep - 1L)[valid]],
      xend      = tr$x[step_rep[valid]],
      yend      = tr$y[step_rep[valid]],
      accepted  = x$accepted[i],
      seg_alpha = base_a * pmax(0, 1 - age_rep[valid] / wake_steps),
      stringsAsFactors = FALSE
    )
  }))

  wake_rej <- wake_df[!wake_df$accepted, , drop = FALSE]
  wake_acc <- wake_df[ wake_df$accepted, , drop = FALSE]

  ## ---- Detection pulse data frame --------------------------------------------
  ## Detection markers appear at det_step (large), decay to normal size over
  ## pulse_frames steps, then persist at normal size for the rest of the run.
  det_pulse_df <- NULL
  if (show.det && n_acc_total > 0L) {
    acc_in_set <- intersect(acc_idx, sim_ids)
    if (length(acc_in_set) > 0L) {
      pulse_list <- lapply(acc_in_set, function(sim_i) {
        ds      <- x$det_step[sim_i]
        if (is.na(ds)) return(NULL)
        loc_row <- which(x$det_locs$id == sim_i)[1L]

        frames_v  <- seq(ds, N)
        ages_v    <- frames_v - ds          ## 0 at detection, grows after

        normal_sz <- pt.size * 0.65
        pulse_sz  <- pt.size * 2.5

        ## Linear decay from pulse_sz to normal_sz, then constant
        size_v <- ifelse(
          ages_v < pulse_frames,
          pulse_sz + (normal_sz - pulse_sz) * ages_v / pulse_frames,
          normal_sz
        )

        data.frame(
          sim      = sim_i,
          frame    = frames_v,
          x        = x$det_locs$x[loc_row],
          y        = x$det_locs$y[loc_row],
          det_size = size_v,
          stringsAsFactors = FALSE
        )
      })
      det_pulse_df <- do.call(rbind, Filter(Negate(is.null), pulse_list))
    }
  }

  ## ---- Most probable path ----------------------------------------------------
  ## Compute at each step; then expand so frame f contains steps 1:f, giving
  ## transition_manual(cumulative=FALSE) a growing path at each frame.
  ## mpp_end_df holds the current head position per frame for the endpoint symbol.
  mpp_path    <- NULL
  mpp_anim_df <- NULL
  mpp_end_df  <- NULL
  if (!is.null(mpp)) {
    mpp      <- match.arg(mpp, c("median", "mean"))
    mpp.sims <- match.arg(mpp.sims, c("accepted", "all"))
    src <- if (mpp.sims == "accepted" && length(acc_idx) > 0L) {
      x$sims[intersect(acc_idx, sim_ids)]
    } else {
      if (mpp.sims == "accepted" && length(acc_idx) == 0L)
        message("mpp.sims = 'accepted' but no accepted sims; using all.")
      x$sims[sim_ids]
    }
    mpp_path    <- .compute_mpp(src, method = mpp)
    mpp_anim_df <- dplyr::bind_rows(lapply(seq_len(N), function(f) {
      n_pts <- min(f, nrow(mpp_path))
      data.frame(frame = f,
                 x     = mpp_path$x[seq_len(n_pts)],
                 y     = mpp_path$y[seq_len(n_pts)])
    }))
    mpp_end_df <- data.frame(
      frame = seq_len(N),
      x     = mpp_path$x[seq_len(N)],
      y     = mpp_path$y[seq_len(N)]
    )
  }

  ## ---- Raster crop -----------------------------------------------------------
  xlim     <- range(c(wake_df$x, wake_df$xend), na.rm = TRUE) + c(-margin, margin)
  ylim     <- range(c(wake_df$y, wake_df$yend), na.rm = TRUE) + c(-margin, margin)
  crop_ext <- terra::ext(xlim[1], xlim[2], ylim[1], ylim[2])

  bathy_crop <- terra::crop(data$bathy, crop_ext)
  land_crop  <- terra::crop(data$land,  crop_ext)

  bathy_df <- terra::as.data.frame(bathy_crop, xy = TRUE, na.rm = TRUE)
  names(bathy_df)[3] <- "depth"
  bathy_df$depth <- bathy_df$depth * -1
  bathy_df <- bathy_df[bathy_df$depth <= 0, ]

  land_df <- terra::as.data.frame(land_crop, xy = TRUE, na.rm = FALSE)
  names(land_df)[3] <- "val"
  land_df <- land_df[!is.na(land_df$val), ]

  ## ---- Receiver positions (static — no 'frame' column) -----------------------
  start_pt <- data.frame(x = x$params$start[1], y = x$params$start[2])
  end_pt   <- data.frame(x = x$params$end[1],   y = x$params$end[2])

  ## ---- ECDF inset layout and data --------------------------------------------
  ## The ECDF is drawn in map-km coordinates (compatible with transition_manual).
  ## to_ix / to_iy scale the unit square [0,1]x[0,1] into the inset box.
  do_ecdf <- show.ecdf && n_acc_total > 0L

  ecdf_df      <- NULL
  count_df     <- NULL
  ecdf_statics <- list()   ## annotate calls assembled below

  if (do_ecdf) {
    inset_w <- diff(xlim) * 0.30
    inset_h <- diff(ylim) * 0.20
    pad     <- margin * 0.7

    ix0 <- if (ecdf.pos %in% c("bottomright","topright")) xlim[2]-pad-inset_w else xlim[1]+pad
    ix1 <- ix0 + inset_w
    iy0 <- if (ecdf.pos %in% c("bottomleft","bottomright")) ylim[1]+pad else ylim[2]-pad-inset_h
    iy1 <- iy0 + inset_h

    to_ix <- function(v) ix0 + v * (ix1 - ix0)
    to_iy <- function(v) iy0 + v * (iy1 - iy0)

    ## Sorted detection steps for all accepted sims in sim_ids
    all_det_steps <- sort(x$det_step[intersect(acc_idx, sim_ids)])
    all_det_steps <- all_det_steps[!is.na(all_det_steps)]

    ## Cumulative detection count at each frame (vectorised)
    det_counts <- findInterval(seq_len(N), all_det_steps)

    ## Per-frame ECDF step data (x=fraction of N, y=cumulative proportion)
    ecdf_df <- dplyr::bind_rows(lapply(seq_len(N), function(f) {
      k <- det_counts[f]
      if (k == 0L)
        return(data.frame(frame = f, ex = to_ix(c(0, 1)), ey = to_iy(c(0, 0))))
      fracs <- all_det_steps[seq_len(k)] / N
      props <- seq_len(k) / n_acc_total
      data.frame(frame = f,
                 ex    = to_ix(c(0, fracs, 1)),
                 ey    = to_iy(c(0, props, props[k])))
    }))

    ## Per-frame count label: "k / n_acc" in top-right of inset
    count_df <- data.frame(
      frame = seq_len(N),
      ex    = to_ix(0.96),
      ey    = to_iy(0.88),
      label = sprintf("%d / %d", det_counts, n_acc_total),
      stringsAsFactors = FALSE
    )

    ## Static annotation elements (no 'frame' col → visible every frame) ------
    ## Text sizes are in ggplot2 units (~mm); keep small relative to the inset.
    tx <- (ix0 + ix1) / 2
    ty <- iy1 - inset_h * 0.06

    ecdf_statics <- list(
      ## White semi-transparent background
      ggplot2::annotate("rect",
        xmin = ix0, xmax = ix1, ymin = iy0, ymax = iy1,
        fill = "white", colour = "grey45", alpha = 0.88, linewidth = 0.35),

      ## 50% reference line
      ggplot2::annotate("segment",
        x = to_ix(0), xend = to_ix(1),
        y = to_iy(0.5), yend = to_iy(0.5),
        linetype = "dashed", colour = "grey60", linewidth = 0.25),

      ## Title inside box
      ggplot2::annotate("text",
        x = tx, y = ty,
        label = "Cumulative detections",
        hjust = 0.5, size = 2.2, colour = "grey25"),

      ## Corner tick labels (x=0 and x=N, y=0 and y=1)
      ggplot2::annotate("text",
        x = to_ix(0.03), y = to_iy(0.06),
        label = "0", hjust = 0, size = 1.9, colour = "grey40"),
      ggplot2::annotate("text",
        x = to_ix(0.97), y = to_iy(0.06),
        label = "N", hjust = 1, size = 1.9, colour = "grey40"),
      ggplot2::annotate("text",
        x = to_ix(0.03), y = to_iy(0.94),
        label = "1", hjust = 0, size = 1.9, colour = "grey40")
    )

    ## Real-world detection line: end.dt expressed as fraction of N steps.
    ## Appended to ecdf_statics so it is static (visible every frame).
    if (show.real.det) {
      elapsed_secs  <- as.numeric(
        difftime(x$params$end.dt, x$params$start.dt, units = "secs"))
      real_det_frac <- min(1, elapsed_secs / (N * x$params$time.step * 60))
      ecdf_statics  <- c(ecdf_statics, list(
        ggplot2::annotate("segment",
          x     = to_ix(real_det_frac), xend = to_ix(real_det_frac),
          y     = to_iy(0),             yend = to_iy(1),
          colour = ln.col.det, linewidth = 0.7, linetype = "solid")
      ))
    }
  }

  ## ---- Labels ----------------------------------------------------------------
  n_acc      <- n_acc_total
  pct_acc    <- round(100 * x$acceptance_rate, 1)
  step_hrs   <- x$params$time.step / 60

  title_str <- paste0(
    "sim_fish  |  move: ", x$params$move,
    "  |  n = ", nsim,
    "  |  ", n_acc, " accepted (", pct_acc, "%)",
    "  |  wake: ", wake_steps, " steps"
  )
  sub_str <- paste0(
    format(x$params$start.dt, "%Y-%m-%d %H:%M"), " UTC",
    "  |  Step {frame} / ", N,
    "  |  {sprintf('%.1f', (as.integer(frame) - 1L) * ", step_hrs, ")} h elapsed"
  )

  ## ---- Assemble ggplot -------------------------------------------------------
  p <- ggplot2::ggplot() +

    ## Bathymetry (static)
    ggplot2::geom_raster(
      data  = bathy_df,
      ggplot2::aes(x = x, y = y, fill = depth),
      alpha = bathy.alpha
    ) +
    ggplot2::scale_fill_gradientn(
      colours = bathy.col,
      name    = "Depth (m)",
      guide   = "none"
    ) +

    ## Land mask (static)
    ggnewscale::new_scale_fill() +
    ggplot2::geom_tile(
      data   = land_df,
      ggplot2::aes(x = x, y = y, fill = "Land"),
      colour = NA
    ) +
    ggplot2::scale_fill_manual(
      values = c("Land" = "#c8c0a8"),
      name   = NULL, guide = "none"
    ) +

    ## Rejected track wake (red, fading)
    { if (nrow(wake_rej) > 0L)
        ggplot2::geom_segment(
          data = wake_rej,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                       alpha = seg_alpha),
          colour    = ln.col,
          linewidth = ln.size,
          lineend   = "round"
        ) else NULL } +

    ## Accepted track wake (blue, fading)
    { if (nrow(wake_acc) > 0L)
        ggplot2::geom_segment(
          data = wake_acc,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                       alpha = seg_alpha),
          colour    = ln.col.det,
          linewidth = ln.size,
          lineend   = "round"
        ) else NULL } +

    ## Alpha identity: seg_alpha values used directly (no scale mapping)
    ggplot2::scale_alpha_identity(guide = "none") +

    ## Detection pulse markers (animated — size pulses at det_step)
    { if (!is.null(det_pulse_df) && nrow(det_pulse_df) > 0L)
        ggplot2::geom_point(
          data = det_pulse_df,
          ggplot2::aes(x = x, y = y, size = det_size),
          shape  = 21,
          fill   = ln.col.det,
          colour = "white",
          stroke = 0.7
        ) else NULL } +

    ## Size identity: det_size values used directly
    ggplot2::scale_size_identity(guide = "none") +

    ## ECDF inset: static annotations (box, reference line, corner labels)
    ecdf_statics +

    ## ECDF step line (animated — grows as detections accumulate)
    { if (!is.null(ecdf_df))
        ggplot2::geom_step(
          data      = ecdf_df,
          ggplot2::aes(x = ex, y = ey),
          colour    = ln.col.det,
          linewidth = 0.65,
          direction = "hv"
        ) else NULL } +

    ## ECDF count label (animated — "k / n_acc" at top-right of inset)
    { if (!is.null(count_df))
        ggplot2::geom_text(
          data  = count_df,
          ggplot2::aes(x = ex, y = ey, label = label),
          hjust = 1, size = 2.0, colour = ln.col.det
        ) else NULL } +

    ## Start receiver (static)
    ggplot2::geom_point(
      data   = start_pt,
      ggplot2::aes(x = x, y = y),
      shape  = 21, fill = "white", colour = ln.col,
      size   = pt.size, stroke = 0.8
    ) +

    ## End receiver (static)
    ggplot2::geom_point(
      data   = end_pt,
      ggplot2::aes(x = x, y = y),
      shape  = 21, fill = ln.col, colour = "white",
      size   = pt.size, stroke = 0.8
    ) +

    ## Most probable path — topmost map layer (above all symbols)
    { if (!is.null(mpp_anim_df))
        ggplot2::geom_path(
          data      = mpp_anim_df,
          ggplot2::aes(x = x, y = y),
          colour    = mpp.col,
          linewidth = mpp.lwd,
          lineend   = "round",
          linejoin  = "round"
        ) else NULL } +
    { if (!is.null(mpp_end_df))
        ggplot2::geom_point(
          data   = mpp_end_df,
          ggplot2::aes(x = x, y = y),
          shape  = 23,
          fill   = mpp.col,
          colour = "black",
          size   = pt.size,
          stroke = 0.8
        ) else NULL } +

    ggplot2::coord_equal(xlim = xlim, ylim = ylim, expand = FALSE, clip = "off") +

    ggplot2::labs(
      x        = "Easting (km, UTM 20N)",
      y        = "Northing (km, UTM 20N)",
      title    = title_str,
      subtitle = sub_str
    ) +

    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 10, face = "bold",
                           margin = ggplot2::margin(b = 2)),
      plot.subtitle    = ggplot2::element_text(size = 7, colour = "grey40",
                           margin = ggplot2::margin(b = 6)),
      axis.text        = ggplot2::element_text(size = 7, colour = "grey40"),
      axis.title       = ggplot2::element_text(size = 8, colour = "grey30"),
      panel.grid.major = ggplot2::element_line(colour = "white", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "#eef3f7", colour = NA),
      legend.position  = "right",
      legend.key.size  = ggplot2::unit(0.5, "cm"),
      plot.margin      = ggplot2::margin(6, 6, 6, 6)
    ) +

    ## transition_manual(frame): at animation frame f, show rows where frame == f.
    ## Layers without a 'frame' column (bathy, land, receivers, ecdf_statics)
    ## are treated as static and shown in every frame.
    gganimate::transition_manual(frame, cumulative = FALSE)

  ## ---- Render ----------------------------------------------------------------
  n_frames <- if (!is.null(nframes)) nframes else N

  gganimate::animate(
    p,
    nframes  = n_frames,
    fps      = fps,
    width    = width,
    height   = height,
    renderer = renderer
  )
}
