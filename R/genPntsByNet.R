#' Generate points based on distance constraints of a fish net
#'
#' This function generate points by a fish net (center points). The fish net is created
#' based on a input raster providing the bounding box and a given cell size.
#'
#' @param x object of class sf
#' @param r raster layer providing the extent of the fish net
#' @param dist distance between points
#' @param do.PointsOnSurface if TRUE than points on surface are generated for missed objects. Default: TRUE
#' @param quiet no outputs in console. Default: TRUE
#' @param save.path full path to store the result. Default: NULL
#' @return
#' list of 1. points of interest, 2. all points
#'
#'
#' @keywords structured points, distance constraint, fish net
#'
#'
#' @export
#'
genPntsByNet <- function(x, r, dist, do.PointsOnSurface = TRUE, quiet = TRUE, save.path = NULL)
{
  # browser()

  # get start time of process
  process.time.start <- proc.time()

  # https://gis.stackexchange.com/questions/225157/generate-rectangular-fishnet-or-vector-grid-cells-shapefile-in-r
  if(!quiet) cat("generation of bounding box\n")
  x.bbx <- sf::st_bbox(r) %>% sf::st_as_sfc(.) %>% sf::st_sf(Id = 1, geometry = .)

  cat("generation of fish net points\n")
  fnet.pnts <- sf::st_make_grid(x = x.bbx, cellsize = dist, what = "centers") %>%
    sf::st_sf(Id = 1:length(.), geometry = .)

  sf::st_crs(fnet.pnts) <- suppressWarnings(sf::st_crs(x))

  cat("intersection to input \n")
  fnet.pnts.sel <- suppressWarnings(sf::st_intersection(x = fnet.pnts, y = x))
  fnet.pnts.sel$Processed <- "Pixel"
  fnet.pnts.sel <- fnet.pnts.sel[, names(fnet.pnts.sel) %in% c("Id", "Processed", "geometry")]

  # check if some features do not contain points
  checkPts <- sf::st_intersects(x = x, y =  fnet.pnts.sel)
  checkPts.fail <- lapply(X = 1:length(checkPts), FUN = function(i, l){
    if(length(l[[i]]) == 0){return(i)} else {NULL}
  }, l = checkPts) %>% unlist(.)

  if(length(checkPts.fail) > 0){
    warning("Some feature does not contain points\n")

    if(do.PointsOnSurface)
    {
      warning("Points on surface is used!\n")

      ptsOnSurface <- suppressWarnings(sf::st_point_on_surface(x = x[checkPts.fail, ]))
      ptsOnSurface$Processed <- "PtsOnSrfce"

      fnet.pnts.sel <- rbind(fnet.pnts.sel, ptsOnSurface)
    }

  } # end of if check

  fnet.pnts.sel$Id <- c(1:nrow(fnet.pnts.sel))

  if(!is.null(save.path)) sf::st_write(obj = fnet.pnts.sel, dsn = save.path)

  # get time of process
  process.time.run <- proc.time() - process.time.start
  if(!quiet) cat(paste0("------ Run of genRandomPntsByRaster: " , round(x = process.time.run["elapsed"][[1]]/60, digits = 3), " Minutes ------\n"))

  return(list(Pnts_Net = fnet.pnts.sel, Pnts_all = fnet.pnts))
} # end of function genPntsByNet
