#' Generate random points based on distance constraints
#'
#' This function generate random points inside polygons based on distance constraints.
#' The function is based on spatstat::rSSI.
#'
#' @param x object of class sf
#' @param r rsater layer object
#' @param field field of x object for rasterisation
#' @param env.grass setting for GRASS GIS initiation
#' @param do.PointsOnSurface if TRUE than points on surface are generated for missed objects. Default: TRUE
#' @param dist minimum distance between points
#' @param save.path full path to store the result. Default: NULL
#' @param seed number to determine the randomness. Default: 12345
#' @param quiet no outputs in console. Default: TRUE
#' @param ... additional parameters for spatstat::rSSI-function
#' @return
#' Randomly sampled points of class sf
#'
#'
#' @keywords random points, distance constraint
#'
#'
#' @export
#'
genRandomPntsByGrd <- function(x, r, field, env.grass, do.PointsOnSurface = TRUE, dist, save.path = NULL, quiet = TRUE, seed = 12345)
{
  # browser()

  # get start time of process
  process.time.start <- proc.time()

  # get projection
  x.proj <- sf::st_crs(x = x)
  r.res <- raster::res(r) %>% unique(.) %>% .[1]

  # rasterize object
  r.rasterize <- raster::rasterize(x = x, y = r, field = field)

  listRndmPts <- lapply(X = 1:3, FUN = function(i, sf.obj, r.rasterize, field, dist, seed){ # nrow(x)
    # browser()

    cat(i, " from ", nrow(sf.obj), "\n")

    # get sf-object
    sf.i <- sf.obj[i, ]

    # crop raster
    r.i <- raster::crop(x = r.rasterize, y = sf.i)

    # rasterize polygon
    # r.i <- raster::rasterize(x = sf.i, y = r.crop, field = field)
    r.i.name <- paste0("tmp_rCrp_", i)
    r.i.path <- file.path(tempdir(), paste0(r.i.name, ".tif"))
    raster::writeRaster(x = r.i, filename = r.i.path, overwrite = TRUE, NAflag = -99999)

    # init GRASS GIS and write raster into GRASS GIS
    invisible(link2GI::linkGRASS7(x = r.i, default_GRASS7 = env.grass, quiet = TRUE))

    # print(rgrass7::parseGRASS(cmd = "r.in.gdal"))
    # print(rgrass7::parseGRASS(cmd = "r.random.cells"))
    rgrass7::execGRASS(cmd = "r.random.cells", flags = c('overwrite', 'quiet'), Sys_show.output.on.console = FALSE, parameters = list(
      output = paste0("rndmCells_", i), distance = dist, seed = seed))

    # print(rgrass7::parseGRASS(cmd = "r.to.vect"))
    rgrass7::execGRASS(cmd = "r.to.vect", flags = c('overwrite', 'quiet'), Sys_show.output.on.console = FALSE, parameters = list(
      input = paste0("rndmCells_", i), output = paste0("rndmPts_", i), type = 'point'))

    # print(rgrass7::parseGRASS(cmd = "v.out.ogr"))
    out.shp <- file.path(tempdir(), paste0("rndmPts_", i, ".shp"))
    rgrass7::execGRASS(cmd = "v.out.ogr", flags = c('overwrite', 'quiet'), Sys_show.output.on.console = FALSE, parameters = list(
      input = paste0("rndmPts_", i), output = out.shp, format = 'ESRI_Shapefile'))

    sf.tmp <- sf::st_read(out.shp, quiet = TRUE)

  }, sf.obj = x, r.rasterize = r.rasterize, field = field, dist = dist, seed = seed)

  # get points and intersect with input feature
  rndmPts <- do.call(rbind, listRndmPts)
  rndmPts <- sf::st_set_crs(x = rndmPts, value = x.proj)

  # rndmPts <- sf::st_intersection(x = rndmPts, y = x)
  rndmPts$Id <- c(1:nrow(rndmPts))
  rndmPts$Processed <- "Pixel"
  rndmPts <- rndmPts[, grep(pattern = "Id|Processed|geometry", x = names(rndmPts))]


  # check if some features do not contain points
  checkPts <- sf::st_intersects(x = sf::st_buffer(x = x, dist = (r.res/2)), y = rndmPts)
  checkPts.fail <- lapply(X = 1:length(checkPts), FUN = function(i, l){
    if(length(l[[i]]) == 0){return(i)} else {NULL}
  }, l = checkPts) %>% unlist(.)

  if(length(checkPts.fail) > 0){
    warning("Some feature does not contain points: Points on surface is used!\n")

    if(do.PointsOnSurface)
    {
      ptsOnSurface <- sf::st_point_on_surface(x = x[checkPts.fail, ])
      ptsOnSurface$Processed <- "PtsOnSrfce"

      rndmPts <- rbind(rndmPts, ptsOnSurface)
    }

  }

  rndmPts$Id <- c(1:nrow(rndmPts))

  if(!is.null(save.path)) sf::st_write(obj = rndmPt, dsn = save.path)


  # get time of process
  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) cat(paste0("------ Run of genRandomPntsByRaster: " , round(x = process.time.run["elapsed"][[1]]/60, digits = 3), " Minutes ------\n"))

  return(rndmPts)
} # end of function genRandomPntsByRaster
