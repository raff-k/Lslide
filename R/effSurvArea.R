pacman::p_load(raster, sf, rgrass7, link2GI, RSAGA)

path.test <- "d:/Users/qo23hel/Desktop/Temp/GRASS TOOL"
path.test <- "D:/Users/rAVer/Desktop/GRASS TOOL"
elev <- raster::raster("L:/EASICLIM/Geom/Parameters/10m/dtm_10.tif")
elev.path <- elev@file@name
elev <- raster::raster(file.path(path.test, "raster_pnt_1.tif"))
elev.path <- elev@file@name
pnt <- sf::st_read("L:/EASICLIM/Geom/Landslides/Lslide_pts_TEST.shp")
maxdist <- 1000
memory <- 2400
method = 'bilinear'
show.output.on.console <- TRUE
path.aspect <- file.path(path.test, "aspect.tif")
path.slope <- file.path(path.test, "slope.tif")
pnt.1 <- pnt[1, ]
pnt.1 <- sf::st_read(file.path(path.test, "pnt_1.shp"))
pnt.1_2 <- pnt[1:2, ]
pnt.1_2 <- sf::st_read(file.path(path.test, "pnt_1_2.shp"))
pnt <- pnt.1_2
NAflag <- -99999

bbox <- (sf::st_bbox(pnt.1) + c(-1500, -1500, 1500, 1500)) %>%
          .[c(1,3,2,4)] %>%
          raster::extent(.)
coords.string <- paste0(sf::st_coordinates(pnt.1), collapse = ",")
coords <- sf::st_coordinates(pnt.1)
obselev <- pnt.1$ALS_EASICLI
obselev <- raster::extract(x = elev, y = pnt)
i <- 1
pref <- "test"
sf::st_write(obj = pnt.1, dsn = file.path(path.test, "pnt_1.shp"))
sf::st_write(obj = pnt.1_2[2,], dsn = file.path(path.test, "pnt_2.shp"))
sf::st_write(obj = pnt.1_2, dsn = file.path(path.test, "pnt_1_2.shp"))
sf::st_coordinates(pnt.1)


env.rsaga <- RSAGA::rsaga.env(path = "C:/OSGeo4W64/apps/saga-6.3.0")

## subset elev to this point
elev.subP1 <- raster::crop(x = elev, y = bbox)
raster::writeRaster(x = elev.subP1, filename =  file.path(path.test, "raster_pnt_1.tif"), NAflag = -99999, overwrite = TRUE)




## init grass gis
link2GI::linkGRASS7(x = elev, default_GRASS7 = c("C:\\OSGeo4W64", "grass-7.4.0", "osgeo4W"),
                    gisdbase = path.test, location = "Tool", quiet = FALSE)



effSurvArea <- function(elev, pts, maxdist = 1000, path.save = tempdir(), method = 'bilinear',
                        memory = 4096, NAflag = -99999, return.geom = TRUE, quiet = TRUE, show.output.on.console = FALSE)
{

  # get start time of process
  process.time.start <- proc.time()

  dim.x <- raster::xres(elev)
  dim.y <- raster::yres(elev)
  dim.max <- max(c(dim.x, dim.y), na.rm = TRUE)


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

  # # print(parseGRASS("r.out.gdal"))
  # rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  #   input = "aspect", output = path.aspect))
  #
  # rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  #   input = "slope", output = path.slope))

  # Process points ----------------------------------------
  if(!quiet) cat("... extract point elevation using method ", method, "\n")
  ## ... creating a list of the available points in the input layer
  # If 'bilinear' the returned values are interpolated from the values of the four nearest raster cells
  pts.elev <- raster::extract(x = elev, y = pts, method = method)
  pts.coord <- sf::st_coordinates(x = pts)


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


  ## Create empty output raster ----------------------------------------
  # creating some empty (0) raster layers
  # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  #   expression = "xxtemp = 0"))
  #
  # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  #   expression = "xxtemp2 = 0"))
  #
  # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  #   expression = "xxtemp3 = 0"))
  #
  # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  #   expression = "xxtemp4 = 0"))

  ## ... creating some empty (0) raster layers
  if(!quiet) cat("... create empty rasters \n")
  xxtemp <- raster::raster(ext = extent(elev), crs = crs(elev), res = res(elev), vals = 0)
  # xxtemp2 <- raster::raster(ext = extent(elev), crs = crs(elev), res = res(elev), vals = 0)
  # xxtemp3 <- raster::raster(ext = extent(elev), crs = crs(elev), res = res(elev), vals = 0)
  # xxtemp4 <- raster::raster(ext = extent(elev), crs = crs(elev), res = res(elev), vals = 0)


  if(!quiet) cat("... START CALCULATION \n")
  # STARTING LOOP ----------------------------------------
  results <- lapply(X = 1:nrow(pts), FUN = function(i, pts, pts.coord, pts.elev, maxdist, dim.max, show.output.on.console, quiet, memory)
  # for(i in 1:nrow(pts))
  {
    if(!quiet) cat("... ... running point: ", i, " of ", nrow(pts), " points\n")

    coords.i <- pts.coord[i,]
    obselev.i <- pts.elev[i]

    if(is.na(obselev.i))
    {
      warning("Point ", i, " skipped due to NA in elevation!\n")
      next
    }

    bbox.i <- (sf::st_bbox(obj = pts[i,]) + c(-(maxdist+dim.max), -(maxdist+dim.max), (maxdist+dim.max), (maxdist+dim.max))) %>%
                as.numeric(.) %>% as.character(.)

    # running visibility analysis ----------------------
    # print(parseGRASS("g.region"))
    rgrass7::execGRASS("g.region", flags = c("a"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      e = bbox.i[[3]], w = bbox.i[[1]], s = bbox.i[[2]], n = bbox.i[[4]]))

    # print(parseGRASS("r.viewshed"))
    rgrass7::execGRASS("r.viewshed", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "dem", output = "view", coordinates = coords.i, max_distance = maxdist, memory = memory))

    # coming back to the original working region
    rgrass7::execGRASS("g.region", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      region = "saved_region"))

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
    path.angle <- file.path(tempdir(), paste0("tmp_angle_", i, ".tif"))
    rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "angle", output = path.angle))


    # print(parseGRASS("r.out.gdal"))
    path.dist <- file.path(tempdir(), paste0("tmp_dist_", i, ".tif"))
    rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "dist_rescaled", output = path.dist))

    ## THROW ERROR IN WINDOWS AND GRASS GIS 7.4.x
    # # updating the output layer of the rescaled distance
    # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    #   expression = "xxtemp4_out = if(isnull(dist_rescaled),xxtemp4, if(xxtemp4 != 0,min(dist_rescaled,xxtemp4),dist_rescaled))"))
    #
    # # updating the output layer of the category of the point who has the higher angles with the considered cell
    # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    #   expression = paste0("xxtemp3_out = if(angle==0, xxtemp3, if(angle<xxtemp,xxtemp3,", i, ") ) ")))
    #
    # # updating the output layer of the number of points from which a cell is visible
    # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    #   expression = paste0("xxtemp2_out = if(angle==0,xxtemp2,xxtemp2+1)")))
    #
    # # updating the output layer of the best angle of view among all the points in the path
    # rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    #   expression = paste0("xxtemp_out = max(xxtemp,angle)")))


    ## load data into R
    r.angle <- raster::raster(path.angle)
    r.dist <- raster::raster(path.dist)


    ## remove files in GRASS GIS session
    # print(parseGRASS("g.remove"))
    rgrass7::execGRASS("g.remove", flags = c("f", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      type = "raster", name = c("angolo_vista", "view", "view90", "view_complete", "a_view", "b_view", "c_view", "angle", "dist_rescaled", "distance")))

    # ## updating temp files
    # # updating the output layer of the rescaled distance
    # xxtemp4 <- raster::overlay(x = xxtemp4, y = r.dist, fun = function(x, y){ifelse(is.na(y), x,
    #                                                                                 ifelse(x != 0, min(y, x, na.rm = TRUE), y))})
    #
    # # updating the output layer of the category of the point who has the higher angles with the considered cell
    # xxtemp3 <- raster::overlay(x = xxtemp3, y = r.angle, fun = function(x, y){ifelse(y == 0, x,
    #                                                                               ifelse(y < x, x, i))})
    #
    # # updating the output layer of the number of points from which a cell is visible
    # xxtemp2 <- raster::overlay(x = xxtemp2, y = r.angle, fun = function(x, y){ifelse(y == 0, x, x + 1)})
    #
    # # updating the output layer of the best angle of view among all the points in the path
    # xxtemp <- raster::overlay(x = xxtemp, y = r.angle, fun = function(x, y) {ifelse(x >= y, x, y)})

    return(list(angle = r.angle, dist = r.dist))

  }, pts = pts, pts.coord = pts.coord, pts.elev = pts.elev, maxdist = maxdist, dim.max = dim.max,
     show.output.on.console = show.output.on.console, quiet = quiet, memory = memory) # end of loop


  ## .... stack data
  if(!quiet) cat("... stack calculation results \n")
  results.angle <- lapply(X = results, FUN = function(x) x$angle) %>%
                   append(xxtemp, .)

  results.dist <- lapply(X = results, FUN = function(x) x$dist) # %>%
                   # append(xxtemp, .)

  stack.angle <- raster::stack(x = results.angle)
  stack.dist <- raster::stack(x = results.dist)



  ## CALCULATE EFFECTIVE SURVEYED AREA -----------------------------
  ## updating temp files
  if(!quiet) cat("... calculate effective surveyed area \n")
  # updating the output layer of the best angle of view among all the points in the path
  xxtemp <- raster::calc(x = stack.angle, max, na.rm = TRUE)

  # updating the output layer of the number of points from which a cell is visible
  xxtemp2 <- raster::calc(x = stack.angle, fun = function(x){length(which(x > 0))})

  # updating the output layer of the category of the point who has the higher angles with the considered cell
  xxtemp3 <- raster::calc(x = stack.angle, fun = function(x){ifelse(length(unique(x)) == 1 && x == 0, 0, which.max(x)-1)})

  # updating the output layer of the rescaled distance
  xxtemp4 <- suppressWarnings(raster::calc(x = stack.dist, fun = min, na.rm = TRUE))


  # set 0 to NA
  xxtemp[xxtemp == 0] <- NA
  xxtemp2[xxtemp2 == 0] <- NA
  xxtemp3[xxtemp3 == 0] <- NA


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




# TESTING FUNCTION -------------------------------
debug(effSurvArea)
undebug(effSurvArea)

# ONE POINT
path.test1 <- file.path(path.test, "Test1")

result.pnt.500 <- effSurvArea(elev = elev, pts = pnt.1_2, maxdist = 500, path.save = path.test1, quiet = FALSE, show.output.on.console = TRUE)
plot(result.pnt.500$distance)

result.pnt.1000 <- effSurvArea(elev = elev, pts = pnt.1_2, maxdist = 1200, path.save = path.test1, quiet = FALSE, show.output.on.console = TRUE)
plot(result.pnt.1000$distance)
plot(result.pnt.1000$viewangles)
plot(result.pnt.1000$pointofview)
plot(result.pnt.1000$numberofviews)


# creating the output layer
# rgrass7::execGRASS("r.null", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   map = "xxtemp_out", setnull = "0"))
#
# rgrass7::execGRASS("r.null", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   map = "xxtemp2_out", setnull = "0"))
#
# rgrass7::execGRASS("r.null", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   map = "xxtemp3_out", setnull = "0"))
#
# rgrass7::execGRASS("r.null", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   map = "xxtemp4_out", setnull = "0"))


# print(parseGRASS("g.copy"))
# rgrass7::execGRASS("g.copy", flags = c("quiet", "overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   raster = paste0("xxtemp_out,", pref,"_viewangles")))
#
# rgrass7::execGRASS("g.copy", flags = c("quiet", "overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   raster = paste0("xxtemp2_out,", pref,"_numberofviews")))
#
# rgrass7::execGRASS("g.copy", flags = c("quiet", "overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   raster = paste0("xxtemp3_out,", pref,"_pointofviews")))
#
# rgrass7::execGRASS("g.copy", flags = c("quiet", "overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   raster = paste0("xxtemp4_out,", pref,"_distance_rescaled")))



# path.viewangles <- "d:/Users/qo23hel/Desktop/Temp/GRASS TOOL/ESA_viewangles.tif"
# path.numberofviews <- "d:/Users/qo23hel/Desktop/Temp/GRASS TOOL/ESA_numberofviews.tif"
# path.pointofview <- "d:/Users/qo23hel/Desktop/Temp/GRASS TOOL/ESA_pointofview.tif"
# path.distance <-  "d:/Users/qo23hel/Desktop/Temp/GRASS TOOL/ESA_distance.tif"
#
# path.viewangles <- "D:/Users/rAVer/Desktop/GRASS TOOL/ESA_viewangles.tif"
# path.numberofviews <- "D:/Users/rAVer/Desktop/GRASS TOOL/ESA_numberofviews.tif"
# path.pointofview <- "D:/Users/rAVer/Desktop/GRASS TOOL/ESA_pointofview.tif"
# path.distance <-  "D:/Users/rAVer/Desktop/GRASS TOOL/ESA_distance.tif"
#
# # print(parseGRASS("r.out.gdal"))
# rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   input = paste0(pref, "_viewangles"), output = path.viewangles))
#
# rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   input = paste0(pref, "_numberofviews"), output = path.numberofviews))
#
# rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   input = paste0(pref, "_pointofviews"), output = path.pointofview))
#
# rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   input = paste0(pref, "_distance_rescaled"), output = path.distance))


# write data
#evaluation of the azimuth layer
# grass.mapcalc("azimuth = (450-aspect) - int( (450-aspect) / 360) * 360", overwrite=True, quiet=True)
# #evaluation of the layer of the vertical component of the versor perpendicular to the terrain slope
# grass.mapcalc("c_dem = cos(slope)", overwrite=True, quiet=True)
# #evaluation of the layer of the north component of the versor perpendicular to the terrain slope
# grass.mapcalc("b_dem = sin(slope)*cos(azimuth)", overwrite=True, quiet=True)
# #evaluation of the layer of the east component of the versor perpendicular to the terrain slope
# grass.mapcalc("a_dem = sin(slope)*sin(azimuth)", overwrite=True, quiet=True)

# r.aspect <- raster::raster(path.aspect)
# r.slope <- raster::raster(path.slope) %>%
#            raster::calc(x =., fun = function(x){x*pi/180})
#
# r.azimuth <- raster::calc(x = r.aspect, fun = function(x){(450-x) - as.integer((450-x)/360)* 360}) %>% # as.integer
#              raster::calc(x =., fun = function(x){x*pi/180})
#
# r.c_dem <-  raster::calc(x = r.slope, fun = function(x){cos(x)})
# r.b_dem <-  raster::overlay(x = r.slope, y = r.azimuth, fun = function(x, y){sin(x) * cos(y)})
# r.a_dem <-  raster::overlay(x = r.slope, y = r.azimuth, fun = function(x, y){sin(x) * sin(y)})
#
# raster::writeRaster(x = r.azimuth, filename = file.path(path.test, "azimuth.tif"), overwrite = TRUE, NAflag = -99999)
# raster::writeRaster(x = r.c_dem, filename = file.path(path.test, "c_dem.tif"), overwrite = TRUE, NAflag = -99999)
# raster::writeRaster(x = r.b_dem, filename = file.path(path.test, "b_dem.tif"), overwrite = TRUE, NAflag = -99999)
# raster::writeRaster(x = r.a_dem, filename = file.path(path.test, "a_dem.tif"), overwrite = TRUE, NAflag = -99999)


# https://sourceforge.net/p/saga-gis/discussion/790705/thread/338ba707/


# if( y()>520210 && x()>569369.1, atan((569369.1-x())/(5202104-y())),
#     if( y()<5202104 && x()>569369.1, 180+atan((569369.1-x())/(5202104-y())),
#         if( y()<5202104 && x()<569369.1, 180+atan((569369.1-x())/(5202104-y())),
#             if( y()>5202104 && x()<569369.1, 360+atan((569369.1-x())/(5202104-y()))))))


# pow(pow(abs(y()-5202104),2)+pow(abs(x()-569369.1),2)+pow(abs( raster_pnt_1@PERMANENT   -(341.39+1.75)),2),0.5)




