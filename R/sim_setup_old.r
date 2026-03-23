#' @title Pre-simulation setup
#'
#' @description load required rasters, receiver locations
#'
#'
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#'
#' @param config - path to config.R script containing file.paths for required & optional environmental layers (see Details)
#' @importFrom terra rast
#' @export
#'
sim_setup_old <-
  function(config = config) {

    suppressWarnings(source(config, local = TRUE, echo=FALSE))
    if(is.null(prj)) prj <- "+proj=utm +zone=20 +units=km +datum=WGS84 +no_defs +type=crs"

    out <- list(
      bathy = suppressWarnings(rast(bathy)),
      land = suppressWarnings(rast(land)),
      d2land = suppressWarnings(rast(d2land)),
      grad = suppressWarnings(rast(grad))
    )

     out[["u"]] <- suppressWarnings(rast(file.path(fvcom, "u.tif")))
     out[["v"]] <- suppressWarnings(rast(file.path(fvcom, "v.tif")))
    # out[["ts"]] <- suppressWarnings(stack(file.path(riops, "riops_doy_t.grd")))

    # out[["recLocs"]] <- esrf_rec
    # out[["recPoly"]] <- recPoly_sf

    # out[["sobi.box"]] <- c(980,1030,1230,1275)
    # out[["esrfPoly"]] <- readRDS(file.path(polygon.file, "NLpoly.RDS"))


    out[["prj"]] <- prj

    return(out)
  }
