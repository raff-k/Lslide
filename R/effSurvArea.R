#' Calculation of effective surveyed area
#'
#' This function calculate the effective surveyed area based on viewshed analysis.
#' A valid GRASS GIS session must be initiated before.
#'
#' @param elev \linkS4class{RasterLayer} containing elevation values
#' @param pts \linkS4class{sf} containing the points from which the effective surveyed area is estimated.
#' @param maxdist viewshed distance to analysie. Default: 1000 m
#' @param ID vector containing names for temporary written data. Vector length must be the same as number of input points. Default: NULL.
#' @param do.extend extend output to extent of input elev. Default: FALSE
#' @param path.save path to store files. Default: tempdir()
#' @param path.temp path to store tempory files. Default: tempdir()
#' @param method method for extracting elevation values for points. Default: 'bilinear'. Use 'simple' for neares neighbor.
#' @param memory memory used in viewshed computation. GRASS GIS r.viewshed. Default: 4096.
#' @param NAflag replace value for NA (NULL data). Default: -99999.
#' @param return.geom If TRUE geometry is returned by functtion. Default: TRUE.
#' @param quiet do not show comments and state. Default: TRUE
#' @param show.output.on.console show GRASS GIS console output. Default: FALSE
#'
#' @return
#' list of \linkS4class{RasterLayer} of viewshed parameters
#'
#'
#' @note
#' \itemize{
#'   \item function is taken from Bornaetxea, T., Rossi, M., Marchesini, I., & Alvioli, M. (2018). Effective surveyed area and its role in statistical landslide susceptibility assessments. Natural Hazards and Earth System Sciences, 18(9), 2455-2469.
#' }
#'
#'
#' @keywords effective surveyed area, GRASS GIS, viewshed
#'
#'
#' @export
effSurvArea <- function(elev, pts, ID = NULL, maxdist = 1000, do.extend = FALSE, path.save = tempdir(), path.temp = tempdir(), method = 'bilinear',
                        memory = 4096, NAflag = -99999, return.geom = TRUE, quiet = TRUE, show.output.on.console = FALSE)
{

  # get start time of process
  process.time.start <- proc.time()

  if(!is.null(ID) && length(ID) != nrow(pts))
  {
    stop("Length of ID must be the same as the number of input points!")
  }

  # get dimension and extent of elevation input
  dim.x <- raster::xres(elev)
  dim.y <- raster::yres(elev)
  dim.max <- max(c(dim.x, dim.y), na.rm = TRUE)
  extent.elev <- raster::extent(elev)
  r.extent <- extent.elev %>% matrix(.) %>% c(.) %>% .[c(1, 3, 2, 4)] # order: xmin ymin xmax ymax

  # re-arrange maxdist based on resolution for region, maxdist is calculated from cell center
  maxdist.x <- ceiling(maxdist/dim.x) * dim.x + (dim.x/2)
  maxdist.y <- ceiling(maxdist/dim.y) * dim.y + (dim.y/2)



  ## mask elevation
  if(!quiet) cat("... mask elevation with buffered points (using maxdist) \n")
  pts.buf <- sf::st_buffer(x = pts, dist = (maxdist+2*dim.max)) %>% # adding double cell size to avoid border effects
    sf::st_union(.) %>%
    sf::st_cast(x = ., to = "POLYGON", warn = FALSE) %>%
    sf::st_sf(ID = 1:length(.), geometry = .)

  elev.mask <- raster::mask(x = elev, mask = pts.buf)


  ## write raster to temp path
  path.elev <- file.path(tempdir(), "tmp_elev.tif")
  raster::writeRaster(x = elev.mask, filename = path.elev, NAflag = NAflag, overwrite = TRUE)

  ## load elevation into GRASS GIS
  # print(parseGRASS("r.in.gdal"))
  rgrass7::execGRASS("r.in.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    input = path.elev, output = "dem"))

  ## save region
  # print(parseGRASS("g.region"))
  rgrass7::execGRASS("g.region", flags = c("quiet", "overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    save = "saved_region"))

  ## calculate slope and aspect
  if(!quiet) cat("... calculate slope and aspect \n")
  # print(parseGRASS("r.slope.aspect"))
  rgrass7::execGRASS("r.slope.aspect", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    elevation = "dem",  slope = "slope", aspect = "aspect"))


  # Process points ----------------------------------------
  if(!quiet) cat("... extract point elevation using method ", method, "\n")
  ## ... creating a list of the available points in the input layer
  # If 'bilinear' the returned values are interpolated from the values of the four nearest raster cells
  pts.elev <- raster::extract(x = elev, y = pts, method = method)
  pts.coord <- sf::st_coordinates(x = pts) # real coordinates
  pts.cells <- raster::cellFromXY(object = elev, xy = pts.coord) # get cell in which point is falling
  pts.NA <- which(is.na(pts.cells)) # get position of points outside of extent
  pts.xy <- raster::xyFromCell(object = elev, cell = pts.cells) # get cell coordiantes in which a point is falling
  pts.xy[pts.NA, ] <- NA # set default value to NA



  ## Process elevation ----------------------------------------
  if(!quiet) cat("... process views on elevation \n")
  ## ... evaluation of the azimuth layer
  # print(parseGRASS("r.mapcalc"))
  rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    expression = "azimuth = (450-aspect) - int( (450-aspect) / 360) * 360"))

  ## ... evaluation of the layer of the vertical component of the versor perpendicular to the terrain slope
  rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    expression = "c_dem = cos(slope)"))

  ## ... evaluation of the layer of the north component of the versor perpendicular to the terrain slope
  rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    expression = "b_dem = sin(slope)*cos(azimuth)"))

  ## ... evaluation of the layer of the east component of the versor perpendicular to the terrain slope
  rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    expression = "a_dem = sin(slope)*sin(azimuth)"))


  ## ... creating some empty (0) raster layers
  # if(!quiet) cat("... create empty rasters \n")
  # xxtemp <- raster::raster(ext = raster::extent(elev), crs = raster::crs(elev), res = raster::res(elev), vals = 0)
  # xxtempNA <- raster::raster(ext = raster::extent(elev), crs = raster::crs(elev), res = raster::res(elev), vals = NA)


  if(!quiet) cat("... START CALCULATION \n")
  # STARTING LOOP ----------------------------------------
  results <- lapply(X = 1:nrow(pts), FUN = function(i, pts, ID, pts.coord, pts.xy, pts.elev, r.extent, maxdist, maxdist.x, maxdist.y, dim.max, show.output.on.console, quiet, path.temp, memory)
    # for(i in 1:nrow(pts))
  {
   # browser()

    ## remove files in GRASS GIS session
    # print(parseGRASS("g.remove"))
    rgrass7::execGRASS("g.remove", flags = c("f", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      type = "raster", name = c("angolo_vista", "view", "view90", "view_complete", "a_view", "b_view", "c_view", "angle", "dist_rescaled", "distance")))


    if(!quiet) cat("... ... running point: ", i, " of ", nrow(pts), " points\n")

    coords.i <- pts.coord[i,]
    obselev.i <- pts.elev[i]
    xy.cell.i <- pts.xy[i, ] # is used to build the region


    if(is.na(obselev.i))
    {
      warning("Point ", i, " skipped due to NA in elevation!\n")
      # next
      return(NULL)
    }

    ## using point coordinate result in a shifted raster
    # order: xmin ymin xmax ymax
    bbox.i <- c(xy.cell.i[1] - maxdist.x, xy.cell.i[2] - maxdist.y, xy.cell.i[1] + maxdist.x, xy.cell.i[2] + maxdist.y)
    if(bbox.i[1] < r.extent[1]){bbox.i[1] <-  r.extent[1]}
    if(bbox.i[2] < r.extent[2]){bbox.i[2] <-  r.extent[2]}
    if(bbox.i[3] > r.extent[3]){bbox.i[3] <-  r.extent[3]}
    if(bbox.i[4] > r.extent[4]){bbox.i[4] <-  r.extent[4]}

    bbox.i <- bbox.i %>% as.character(.)


    # running visibility analysis ----------------------
    # print(parseGRASS("g.region"))
    rgrass7::execGRASS("g.region", flags = c("a"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      e = bbox.i[[3]], w = bbox.i[[1]], s = bbox.i[[2]], n = bbox.i[[4]]))

    # print(parseGRASS("r.viewshed"))
    rgrass7::execGRASS("r.viewshed", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "dem", output = "view", coordinates = coords.i, max_distance = maxdist, memory = memory))

    # since r.viewshed set the cell of the output visibility layer to 180 under the point, this cell is set to 0.01
    rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      expression = "view_complete = if(view == 180, 0.01, view)"))

    # estimating the layer of the horizontal angle between point and each visible cell (angle of the horizontal line of sight)
    expr <- paste0("angolo_vista =
                   if( y()>", coords.i[[2]], " && x()>", coords.i[[1]], ", atan((", coords.i[[1]], "-x())/(", coords.i[[2]], "-y())),
                   if( y()<", coords.i[[2]], " && x()>", coords.i[[1]], ", 180+atan((", coords.i[[1]], "-x())/(", coords.i[[2]], "-y())),
                   if( y()<", coords.i[[2]], " && x()<", coords.i[[1]], ", 180+atan((", coords.i[[1]], "-x())/(", coords.i[[2]], "-y())),
                   if( y()>", coords.i[[2]], " && x()<", coords.i[[1]], ", 360+atan((", coords.i[[1]], "-x())/(", coords.i[[2]], "-y())))      )      )    )")

    expr <- gsub(pattern = "\n", replacement = "", x = expr)

    rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      expression = expr))


    # estimating the layer of the vertical angle between point and each visible cell  (angle of the vertical line of sight)
    rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      expression = "view90 = view_complete - 90"))

    # evaluate the vertical component of the versor oriented along the line of sight
    rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      expression = "c_view = sin(view90)"))

    # evaluate the northern component of the versor oriented along the line of sight
    rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      expression = "b_view = cos(view90)*cos(angolo_vista)"))

    # evaluate the eastern component of the versor oriented along the line of sight
    rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      expression = "a_view = cos(view90)*sin(angolo_vista)"))


    # estimate the three-dimensional distance between the point and each visible cell
    expr_distance <- paste0("distance = pow(pow(abs(y()-", coords.i[[2]], "),2)+pow(abs(x()-", coords.i[[1]], "),2)+pow(abs(dem-(", obselev.i, "+1.75)),2),0.5)")

    rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      expression = expr_distance))

    # estimating the layer of the angle between the versor of the terrain and the line of sight
    rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      expression = "angle = acos((a_view*a_dem+b_view*b_dem+c_view*c_dem)/(sqrt(a_view*a_view+b_view*b_view+c_view*c_view)*sqrt(a_dem*a_dem+b_dem*b_dem+c_dem*c_dem)))"))

    # evaluating the layer of the distance scaled by the cosine of the angle
    rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      expression = "dist_rescaled = if(angle>91,(distance/(-cos(angle))),null())"))


    # setting all the null cells to zero
    # print(parseGRASS("r.null"))
    rgrass7::execGRASS("r.null", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      map = "angle", null = 0))

    ## write out data
    # print(parseGRASS("r.out.gdal"))
    if(!is.null(ID))
    {
      path.angle <- file.path(path.temp, paste0("tmp_angle_", ID[i], ".tif"))
    } else {
      path.angle <- file.path(path.temp, paste0("tmp_angle_", i, ".tif"))
    }

    rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "angle", output = path.angle))


    # print(parseGRASS("r.out.gdal"))
    if(!is.null(ID))
    {
      path.dist <- file.path(path.temp, paste0("tmp_dist_", ID[i], ".tif"))
    } else {
      path.dist <- file.path(path.temp, paste0("tmp_dist_", i, ".tif"))
    }
    rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "dist_rescaled", output = path.dist))

    # coming back to the original working region
    rgrass7::execGRASS("g.region", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      region = "saved_region"))

    ## load data into R
    r.angle <- raster::raster(path.angle)
    r.dist <- raster::raster(path.dist)


    ## remove files in GRASS GIS session
    # print(parseGRASS("g.remove"))
    rgrass7::execGRASS("g.remove", flags = c("f", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      type = "raster", name = c("angolo_vista", "view", "view90", "view_complete", "a_view", "b_view", "c_view", "angle", "dist_rescaled", "distance")))


    return(list(angle = r.angle, dist = r.dist))

  }, pts = pts, pts.coord = pts.coord, ID = ID, pts.xy = pts.xy, pts.elev = pts.elev, maxdist = maxdist, maxdist.x = maxdist.x, maxdist.y = maxdist.y, dim.max = dim.max,
  show.output.on.console = show.output.on.console, quiet = quiet, memory = memory, path.temp = path.temp, r.extent = r.extent) # end of loop


  ## .... check data
  results.angle <- lapply(X = results, FUN = function(x) x$angle) # %>%
    # append(xxtemp, .)

  results.angle <- results.angle[!vapply(results.angle, is.null, logical(1))] # ... remove possible NULLS



  results.dist <- lapply(X = results, FUN = function(x) x$dist)  # %>%
   # append(xxtempNA, .)
  results.dist <- results.dist[!vapply(results.dist, is.null, logical(1))] # ... remove possible NULLS


  ## CALCULATE EFFECTIVE SURVEYED AREA -----------------------------
  if(!quiet) cat("... calculate effective surveyed area by mosaicing \n")

  # updating the output layer of the best angle of view among all the points in the path
  if(!quiet) cat("... view angles \n")
  l.xxtemp <- results.angle
  l.xxtemp$fun <- max
  l.xxtemp$na.rm <- TRUE
  l.xxtemp$filename <- file.path(path.temp, "tmp_viewangles.tif")
  xxtemp <- do.call(what = raster::mosaic, args = l.xxtemp)

  # updating the output layer of the number of points from which a cell is visible
  if(!quiet) cat("... number of views \n")
  l.xxtemp2 <- results.angle
  l.xxtemp2$fun <- function(x, na.rm){length(which(x > 0))}
  l.xxtemp2$filename <- file.path(path.temp, "tmp_numberofviews.tif")
  xxtemp2 <- do.call(what = raster::mosaic, args = l.xxtemp2)


  # updating the output layer of the category of the point who has the higher angles with the considered cell
  if(!quiet) cat("... point of view \n")
  l.xxtemp3 <- results.angle
  l.xxtemp3$fun <- function(x, na.rm){ifelse(length(which.max(x)) == 0 || max(x, na.rm = na.rm) == 0, NA, which.max(x))}
  l.xxtemp3$filename <- file.path(path.temp, "tmp_pointofview.tif")
  xxtemp3 <- do.call(what = raster::mosaic, args = l.xxtemp3)


  # updating the output layer of the rescaled distance
  if(!quiet) cat("... distance \n")
  l.xxtemp4 <- results.dist
  l.xxtemp4$fun <- min
  l.xxtemp4$na.rm <- TRUE
  l.xxtemp4$filename <- file.path(path.temp, "tmp_distance.tif")
  xxtemp4 <- do.call(what = raster::mosaic, args = l.xxtemp4)


  # set 0 to NA
  xxtemp[xxtemp == 0] <- NA
  xxtemp2[xxtemp2 == 0] <- NA
  xxtemp3[xxtemp3 == 0] <- NA

  if(do.extend)
  {
    xxtemp <- raster::extend(x = xxtemp, y = extent.elev, value = NA)
    if(!identical(raster::extent(xxtemp), extent.elev)){raster::extent(xxtemp) <- extent.elev}

    xxtemp2 <- raster::extend(x = xxtemp2, y = elev, value = NA)
    if(!identical(raster::extent(xxtemp2), extent.elev)){raster::extent(xxtemp2) <- extent.elev}

    xxtemp3 <- raster::extend(x = xxtemp3, y = elev, value = NA)
    if(!identical(raster::extent(xxtemp3), extent.elev)){raster::extent(xxtemp3) <- extent.elev}

    xxtemp4 <- raster::extend(x = xxtemp4, y = elev, value = NA)
    if(!identical(raster::extent(xxtemp4), extent.elev)){raster::extent(xxtemp4) <- extent.elev}
  }


  # save rasters
  raster::writeRaster(x = xxtemp, filename = file.path(path.save, "viewangles.tif"), NAflag = NAflag, overwrite = TRUE)
  raster::writeRaster(x = xxtemp2, filename = file.path(path.save, "numberofviews.tif"), NAflag = NAflag, overwrite = TRUE)
  raster::writeRaster(x = xxtemp3, filename =  file.path(path.save, "pointofview.tif"), NAflag = NAflag, overwrite = TRUE)
  raster::writeRaster(x = xxtemp4, filename = file.path(path.save, "distance.tif"), NAflag = NAflag, overwrite = TRUE)

  # get time of process
  process.time.run <- proc.time() - process.time.start
  if(!quiet) cat(paste0("------ Run of effSurvArea: " , round(x = process.time.run["elapsed"][[1]]/60, digits = 4), " Minutes \n"))


  if(return.geom)
  {
    # return data
    return(list(viewangles = xxtemp,
                numberofviews = xxtemp2,
                pointofview = xxtemp3,
                distance = xxtemp4))
  }

} # end of function effSurvArea



























# effSurvArea <- function(elev, pts, maxdist = 1000, path.save = tempdir(), path.temp = tempdir(), method = 'bilinear', cores = 1,
#                         memory = 4096, NAflag = -99999, return.geom = TRUE, quiet = TRUE, show.output.on.console = FALSE)
# {
#
#   # get start time of process
#   process.time.start <- proc.time()
#
#   dim.x <- raster::xres(elev)
#   dim.y <- raster::yres(elev)
#   dim.max <- max(c(dim.x, dim.y), na.rm = TRUE)
#
#
#   ## mask elevation
#   if(!quiet) cat("... mask elevation with buffered points (using maxdist) \n")
#   pts.buf <- sf::st_buffer(x = pts, dist = (maxdist+2*dim.max)) %>% # adding double cell size to avoid border effects
#     sf::st_union(.) %>%
#     sf::st_cast(x = ., to = "POLYGON", warn = FALSE) %>%
#     sf::st_sf(ID = 1:length(.), geometry = .)
#
#   elev.mask <- raster::mask(x = elev, mask = pts.buf)
#
#
#   ## write raster to temp path
#   path.elev <- file.path(tempdir(), "tmp_elev.tif")
#   raster::writeRaster(x = elev.mask, filename = path.elev, NAflag = NAflag, overwrite = TRUE)
#
#   ## load elevation into GRASS GIS
#   # print(parseGRASS("r.in.gdal"))
#   rgrass7::execGRASS("r.in.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#     input = path.elev, output = "dem"))
#
#   ## save region
#   # print(parseGRASS("g.region"))
#   rgrass7::execGRASS("g.region", flags = c("quiet", "overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#     save = "saved_region"))
#
#   ## calculate slope and aspect
#   if(!quiet) cat("... calculate slope and aspect \n")
#   # print(parseGRASS("r.slope.aspect"))
#   rgrass7::execGRASS("r.slope.aspect", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#     elevation = "dem",  slope = "slope", aspect = "aspect"))
#
#   # # print(parseGRASS("r.out.gdal"))
#   # rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   #   input = "aspect", output = path.aspect))
#   #
#   # rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   #   input = "slope", output = path.slope))
#
#   # Process points ----------------------------------------
#   if(!quiet) cat("... extract point elevation using method ", method, "\n")
#   ## ... creating a list of the available points in the input layer
#   # If 'bilinear' the returned values are interpolated from the values of the four nearest raster cells
#   pts.elev <- raster::extract(x = elev, y = pts, method = method)
#   pts.coord <- sf::st_coordinates(x = pts)
#
#
#   ## Process elevation ----------------------------------------
#   if(!quiet) cat("... process views on elevation \n")
#   ## ... evaluation of the azimuth layer
#   # print(parseGRASS("r.mapcalc"))
#   rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#     expression = "azimuth = (450-aspect) - int( (450-aspect) / 360) * 360"))
#
#   ## ... evaluation of the layer of the vertical component of the versor perpendicular to the terrain slope
#   rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#     expression = "c_dem = cos(slope)"))
#
#   ## ... evaluation of the layer of the north component of the versor perpendicular to the terrain slope
#   rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#     expression = "b_dem = sin(slope)*cos(azimuth)"))
#
#   ## ... evaluation of the layer of the east component of the versor perpendicular to the terrain slope
#   rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#     expression = "a_dem = sin(slope)*sin(azimuth)"))
#
#
#   ## Create empty output raster ----------------------------------------
#   # creating some empty (0) raster layers
#   # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   #   expression = "xxtemp = 0"))
#   #
#   # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   #   expression = "xxtemp2 = 0"))
#   #
#   # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   #   expression = "xxtemp3 = 0"))
#   #
#   # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   #   expression = "xxtemp4 = 0"))
#
#   ## ... creating some empty (0) raster layers
#   if(!quiet) cat("... create empty rasters \n")
#   xxtemp <- raster::raster(ext = raster::extent(elev), crs = raster::crs(elev), res = raster::res(elev), vals = 0)
#   # xxtemp2 <- raster::raster(ext = extent(elev), crs = crs(elev), res = res(elev), vals = 0)
#   # xxtemp3 <- raster::raster(ext = extent(elev), crs = crs(elev), res = res(elev), vals = 0)
#   # xxtemp4 <- raster::raster(ext = extent(elev), crs = crs(elev), res = res(elev), vals = 0)
#
#
#   if(!quiet) cat("... START CALCULATION \n")
#   # STARTING LOOP ----------------------------------------
#   results <- lapply(X = 1:nrow(pts), FUN = function(i, pts, pts.coord, pts.elev, maxdist, dim.max, show.output.on.console, quiet, path.temp, memory)
#     # for(i in 1:nrow(pts))
#   {
#     # browser()
#
#     if(!quiet) cat("... ... running point: ", i, " of ", nrow(pts), " points\n")
#
#     coords.i <- pts.coord[i,]
#     obselev.i <- pts.elev[i]
#
#     if(is.na(obselev.i))
#     {
#       warning("Point ", i, " skipped due to NA in elevation!\n")
#       # next
#       return(NULL)
#     }
#
#     bbox.i <- (sf::st_bbox(obj = pts[i,]) + c(-(maxdist+dim.max), -(maxdist+dim.max), (maxdist+dim.max), (maxdist+dim.max))) %>%
#       as.numeric(.) %>% as.character(.)
#
#     # running visibility analysis ----------------------
#     # print(parseGRASS("g.region"))
#     rgrass7::execGRASS("g.region", flags = c("a"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       e = bbox.i[[3]], w = bbox.i[[1]], s = bbox.i[[2]], n = bbox.i[[4]]))
#
#     # print(parseGRASS("r.viewshed"))
#     rgrass7::execGRASS("r.viewshed", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       input = "dem", output = "view", coordinates = coords.i, max_distance = maxdist, memory = memory))
#
#     # coming back to the original working region
#     rgrass7::execGRASS("g.region", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       region = "saved_region"))
#
#     # since r.viewshed set the cell of the output visibility layer to 180 under the point, this cell is set to 0.01
#     rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       expression = "view_complete = if(view == 180, 0.01, view)"))
#
#     # estimating the layer of the horizontal angle between point and each visible cell (angle of the horizontal line of sight)
#     expr <- paste0("angolo_vista =
#                    if( y()>", coords.i[[2]], " && x()>", coords.i[[1]], ", atan((", coords.i[[1]], "-x())/(", coords.i[[2]], "-y())),
#                    if( y()<", coords.i[[2]], " && x()>", coords.i[[1]], ", 180+atan((", coords.i[[1]], "-x())/(", coords.i[[2]], "-y())),
#                    if( y()<", coords.i[[2]], " && x()<", coords.i[[1]], ", 180+atan((", coords.i[[1]], "-x())/(", coords.i[[2]], "-y())),
#                    if( y()>", coords.i[[2]], " && x()<", coords.i[[1]], ", 360+atan((", coords.i[[1]], "-x())/(", coords.i[[2]], "-y())))      )      )    )")
#
#     expr <- gsub(pattern = "\n", replacement = "", x = expr)
#
#     rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       expression = expr))
#
#
#     # estimating the layer of the vertical angle between point and each visible cell  (angle of the vertical line of sight)
#     rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       expression = "view90 = view_complete - 90"))
#
#     # evaluate the vertical component of the versor oriented along the line of sight
#     rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       expression = "c_view = sin(view90)"))
#
#     # evaluate the northern component of the versor oriented along the line of sight
#     rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       expression = "b_view = cos(view90)*cos(angolo_vista)"))
#
#     # evaluate the eastern component of the versor oriented along the line of sight
#     rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       expression = "a_view = cos(view90)*sin(angolo_vista)"))
#
#
#     # estimate the three-dimensional distance between the point and each visible cell
#     expr_distance <- paste0("distance = pow(pow(abs(y()-", coords.i[[2]], "),2)+pow(abs(x()-", coords.i[[1]], "),2)+pow(abs(dem-(", obselev.i, "+1.75)),2),0.5)")
#
#     rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       expression = expr_distance))
#
#     # estimating the layer of the angle between the versor of the terrain and the line of sight
#     rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       expression = "angle = acos((a_view*a_dem+b_view*b_dem+c_view*c_dem)/(sqrt(a_view*a_view+b_view*b_view+c_view*c_view)*sqrt(a_dem*a_dem+b_dem*b_dem+c_dem*c_dem)))"))
#
#     # evaluating the layer of the distance scaled by the cosine of the angle
#     rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       expression = "dist_rescaled = if(angle>91,(distance/(-cos(angle))),null())"))
#
#
#     # setting all the null cells to zero
#     # print(parseGRASS("r.null"))
#     rgrass7::execGRASS("r.null", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       map = "angle", null = 0))
#
#     ## write out data
#     # print(parseGRASS("r.out.gdal"))
#     path.angle <- file.path(path.temp, paste0("tmp_angle_", i, ".tif"))
#     rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       input = "angle", output = path.angle))
#
#
#     # print(parseGRASS("r.out.gdal"))
#     path.dist <- file.path(path.temp, paste0("tmp_dist_", i, ".tif"))
#     rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       input = "dist_rescaled", output = path.dist))
#
#     ## THROW ERROR IN WINDOWS AND GRASS GIS 7.4.x
#     # # updating the output layer of the rescaled distance
#     # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#     #   expression = "xxtemp4_out = if(isnull(dist_rescaled),xxtemp4, if(xxtemp4 != 0,min(dist_rescaled,xxtemp4),dist_rescaled))"))
#     #
#     # # updating the output layer of the category of the point who has the higher angles with the considered cell
#     # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#     #   expression = paste0("xxtemp3_out = if(angle==0, xxtemp3, if(angle<xxtemp,xxtemp3,", i, ") ) ")))
#     #
#     # # updating the output layer of the number of points from which a cell is visible
#     # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#     #   expression = paste0("xxtemp2_out = if(angle==0,xxtemp2,xxtemp2+1)")))
#     #
#     # # updating the output layer of the best angle of view among all the points in the path
#     # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#     #   expression = paste0("xxtemp_out = max(xxtemp,angle)")))
#
#
#     ## load data into R
#     r.angle <- raster::raster(path.angle)
#     r.dist <- raster::raster(path.dist)
#
#
#     ## remove files in GRASS GIS session
#     # print(parseGRASS("g.remove"))
#     rgrass7::execGRASS("g.remove", flags = c("f", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#       type = "raster", name = c("angolo_vista", "view", "view90", "view_complete", "a_view", "b_view", "c_view", "angle", "dist_rescaled", "distance")))
#
#     # ## updating temp files
#     # # updating the output layer of the rescaled distance
#     # xxtemp4 <- raster::overlay(x = xxtemp4, y = r.dist, fun = function(x, y){ifelse(is.na(y), x,
#     #                                                                                 ifelse(x != 0, min(y, x, na.rm = TRUE), y))})
#     #
#     # # updating the output layer of the category of the point who has the higher angles with the considered cell
#     # xxtemp3 <- raster::overlay(x = xxtemp3, y = r.angle, fun = function(x, y){ifelse(y == 0, x,
#     #                                                                               ifelse(y < x, x, i))})
#     #
#     # # updating the output layer of the number of points from which a cell is visible
#     # xxtemp2 <- raster::overlay(x = xxtemp2, y = r.angle, fun = function(x, y){ifelse(y == 0, x, x + 1)})
#     #
#     # # updating the output layer of the best angle of view among all the points in the path
#     # xxtemp <- raster::overlay(x = xxtemp, y = r.angle, fun = function(x, y) {ifelse(x >= y, x, y)})
#
#     return(list(angle = r.angle, dist = r.dist))
#
#   }, pts = pts, pts.coord = pts.coord, pts.elev = pts.elev, maxdist = maxdist, dim.max = dim.max,
#   show.output.on.console = show.output.on.console, quiet = quiet, memory = memory, path.temp = path.temp) # end of loop
#
#   ## ... remove possible NULLS
#   results <- results[!vapply(results, is.null, logical(1))]
#
#
#   ## .... stack data
#   if(!quiet) cat("... stack calculation results \n")
#   results.angle <- lapply(X = results, FUN = function(x) x$angle) %>%
#     append(xxtemp, .)
#
#   results.dist <- lapply(X = results, FUN = function(x) x$dist) # %>%
#   # append(xxtemp, .)
#
#   stack.angle <- raster::stack(x = results.angle)
#   stack.dist <- raster::stack(x = results.dist)
#
#
#   ## CALCULATE EFFECTIVE SURVEYED AREA -----------------------------
#   ## updating temp files
#   if(!quiet) cat("... calculate effective surveyed area \n")
#
#   if(cores > 1)
#   {
#     if(!quiet) cat("... ... using parallel mode \n")
#
#     # raster::beginCluster(n = cores)
#     cl <- parallel::makeCluster(cores)
#
#     # updating the output layer of the best angle of view among all the points in the path
#     xxtemp <- raster::clusterR(cl = cl, x = stack.angle, fun = calc, args = list(max, na.rm = TRUE),
#                                filename = file.path(path.temp, "xxtemp.tif"), overwrite = TRUE, NAflag = NAflag)
#
#     # updating the output layer of the number of points from which a cell is visible
#     xxtemp2 <- raster::clusterR(cl = cl, x = stack.angle, fun = calc, args = list(function(x){length(which(x > 0))}),
#                                 filename = file.path(path.temp, "xxtemp.tif"), overwrite = TRUE, NAflag = NAflag)
#
#     # updating the output layer of the category of the point who has the higher angles with the considered cell
#     xxtemp3 <- raster::clusterR(cl = cl, x = stack.angle, fun = calc, args = list(function(x){ifelse(length(unique(x)) == 1 && x == 0, 0, which.max(x)-1)}),
#                                 filename = file.path(path.temp, "xxtemp.tif"), overwrite = TRUE, NAflag = NAflag)
#
#     # updating the output layer of the rescaled distance
#     xxtemp4 <- raster::clusterR(cl = cl, x = stack.dist, fun = calc, args = list(min, na.rm = TRUE),
#                                 filename = file.path(path.temp, "xxtemp.tif"), overwrite = TRUE, NAflag = NAflag)
#
#
#     parallel::stopCluster(cl)
#
#   } else {
#
#     if(!quiet) cat("... ... using sequential mode \n")
#
#     # updating the output layer of the best angle of view among all the points in the path
#     xxtemp <- raster::calc(x = stack.angle, max, na.rm = TRUE)
#
#     # updating the output layer of the number of points from which a cell is visible
#     xxtemp2 <- raster::calc(x = stack.angle, fun = function(x){length(which(x > 0))})
#
#     # updating the output layer of the category of the point who has the higher angles with the considered cell
#     xxtemp3 <- raster::calc(x = stack.angle, fun = function(x){ifelse(length(unique(x)) == 1 && x == 0, 0, which.max(x)-1)})
#
#     # updating the output layer of the rescaled distance
#     xxtemp4 <- suppressWarnings(raster::calc(x = stack.dist, fun = min, na.rm = TRUE))
#   }
#
#   # set 0 to NA
#   xxtemp[xxtemp == 0] <- NA
#   xxtemp2[xxtemp2 == 0] <- NA
#   xxtemp3[xxtemp3 == 0] <- NA
#
#
#   # save rasters
#   raster::writeRaster(x = xxtemp, filename = file.path(path.save, "viewangles.tif"), NAflag = NAflag, overwrite = TRUE)
#   raster::writeRaster(x = xxtemp2, filename = file.path(path.save, "numberofviews.tif"), NAflag = NAflag, overwrite = TRUE)
#   raster::writeRaster(x = xxtemp3, filename =  file.path(path.save, "pointofview.tif"), NAflag = NAflag, overwrite = TRUE)
#   raster::writeRaster(x = xxtemp4, filename = file.path(path.save, "distance.tif"), NAflag = NAflag, overwrite = TRUE)
#
#   # get time of process
#   process.time.run <- proc.time() - process.time.start
#   if(!quiet) cat(paste0("------ Run of effSurvArea: " , round(x = process.time.run["elapsed"][[1]]/60, digits = 4), " Minutes \n"))
#
#
#   if(return.geom)
#   {
#     # return data
#     return(list(viewangles = xxtemp,
#                 numberofviews = xxtemp2,
#                 pointofview = xxtemp3,
#                 distance = xxtemp4))
#   }
#
# } # end of function effSurvArea


