pacman::p_load(raster, sf, rgrass7, link2GI)

path.test <- "d:/Users/qo23hel/Desktop/Temp/GRASS TOOL"
elev <- raster::raster("L:/EASICLIM/Geom/Parameters/10m/dtm_10.tif")
elev.path <- elev@file@name
pnt <- sf::st_read("L:/EASICLIM/Geom/Landslides/Lslide_pts_TEST.shp")
maxdist <- 1000
show.output.on.console <- TRUE
path.aspect <- file.path(path.test, "aspect.tif")
path.slope <- file.path(path.test, "slope.tif")
pnt.1 <- pnt[1, ]
pnt.1_2 <- pnt[1:2, ]
bbox <- (sf::st_bbox(pnt.1) + c(-1500, -1500, 1500, 1500)) %>%
          .[c(1,3,2,4)] %>%
          raster::extent(.)
coords.string <- paste0(sf::st_coordinates(pnt.1), collapse = ",")
coords <- sf::st_coordinates(pnt.1)
obselev <- pnt.1$ALS_EASICLI
i <- 1
pref <- "test"
sf::st_write(obj = pnt.1, dsn = file.path(path.test, "pnt_1.shp"))
sf::st_write(obj = pnt.1_2, dsn = file.path(path.test, "pnt_1_2.shp"))
sf::st_coordinates(pnt.1)


## subset elev to this point
elev.subP1 <- raster::crop(x = elev, y = bbox)
raster::writeRaster(x = elev.subP1, filename =  file.path(path.test, "raster_pnt_1.tif"), NAflag = -99999, overwrite = TRUE)




## init grass gis
link2GI::linkGRASS7(x = elev, default_GRASS7 = c("C:\\OSGeo4W64", "grass-7.4.1", "osgeo4W"), quiet = FALSE)

# print(parseGRASS("r.in.gdal"))
rgrass7::execGRASS("r.in.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  input = elev.path, output = "dem"))

# # print(parseGRASS("g.region"))
# rgrass7::execGRASS("g.region", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
#   save = "saved_region"))


# print(parseGRASS("r.slope.aspect"))
rgrass7::execGRASS("r.slope.aspect", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  elevation = "dem",  slope = "slope", aspect = "aspect"))

# print(parseGRASS("r.out.gdal"))
rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  input = "aspect", output = path.aspect))

rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  input = "slope", output = path.slope))


# print(parseGRASS("r.mapcalc"))
rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = "azimuth = (450-aspect) - int( (450-aspect) / 360) * 360"))

rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = "c_dem = cos(slope)"))

rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = "b_dem = sin(slope)*cos(azimuth)"))

rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = "a_dem = sin(slope)*sin(azimuth)"))

# creating some empty (0) raster layers
rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = "xxtemp = 0"))

rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = "xxtemp2 = 0"))

rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = "xxtemp3 = 0"))

rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = "xxtemp4 = 0"))


# running visibility analysis ----------------------
# print(parseGRASS("r.viewshed"))
rgrass7::execGRASS("r.viewshed", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  input = "dem", output = "view", coordinates = coords, max_distance = maxdist, memory = 5000))

rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = "view_complete = if(view == 180, 0.01, view)"))


rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = paste0("angolo_vista = if( y() > ", coords[[2]], " && x() > ", coords[[1]], ", atan((", coords[[1]], "-x())/(", coords[[2]], "-y())), \
    if(y() < ", coords[[2]], " && x() > ", coords[[1]], ", 180 + atan((", coords[[1]], " -x())/(", coords[[2]], "-y())), \
        if(y()< ", coords[[2]], " && x()<", coords[[1]], ", 180+atan((", coords[[1]], "-x())/(", coords[[2]], "-y())), \
            if(y() > ", coords[[2]], " && x()<", coords[[1]], ", 360+atan((", coords[[1]], "-x())/(", coords[[2]], "-y()))) ) ) )")))

expr <- paste0("angolo_vista =
            if( y()>", coords[[2]], " && x()>", coords[[1]], ", atan((", coords[[1]], "-x())/(", coords[[2]], "-y())),
  if( y()<", coords[[2]], " && x()>", coords[[1]], ", 180+atan((", coords[[1]], "-x())/(", coords[[2]], "-y())),
  if( y()<", coords[[2]], " && x()<", coords[[1]], ", 180+atan((", coords[[1]], "-x())/(", coords[[2]], "-y())),
  if( y()>", coords[[2]], " && x()<", coords[[1]], ", 360+atan((", coords[[1]], "-x())/(", coords[[2]], "-y())))      )      )    )")

expr <- gsub(pattern = "\n", replacement = "", x = expr)

rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = expr))

rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = "view90 = view_complete - 90"))

rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = "c_view = sin(view90)"))

rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = "b_view = cos(view90)*cos(angolo_vista)"))

rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = "a_view = cos(view90)*sin(angolo_vista)"))


expr_distance <- paste0("distance = pow(pow(abs(y()-", coords[[2]], "),2)+pow(abs(x()-", coords[[1]], "),2)+pow(abs(dem-(", obselev, "+1.75)),2),0.5)")

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

# updating the output layer of the rescaled distance
rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = "xxtemp4_out = if(isnull(dist_rescaled),xxtemp4, if(xxtemp4 != 0,min(dist_rescaled,xxtemp4),dist_rescaled))"))

# updating the output layer of the category of the point who has the higher angles with the considered cell
rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = paste0("xxtemp3_out = if(angle==0, xxtemp3, if(angle<xxtemp,xxtemp3,", i, ") ) ")))

# updating the output layer of the number of points from which a cell is visible
rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = paste0("xxtemp2_out = if(angle==0,xxtemp2,xxtemp2+1)")))

# updating the output layer of the best angle of view among all the points in the path
rgrass7::execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  expression = paste0("xxtemp_out = max(xxtemp,angle)")))






# creating the output layer
rgrass7::execGRASS("r.null", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  map = "xxtemp_out", setnull = "0"))

rgrass7::execGRASS("r.null", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  map = "xxtemp2_out", setnull = "0"))

rgrass7::execGRASS("r.null", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  map = "xxtemp3_out", setnull = "0"))

rgrass7::execGRASS("r.null", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  map = "xxtemp4_out", setnull = "0"))


# print(parseGRASS("g.copy"))
rgrass7::execGRASS("g.copy", flags = c("quiet", "overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  raster = paste0("xxtemp_out,", pref,"_viewangles")))

rgrass7::execGRASS("g.copy", flags = c("quiet", "overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  raster = paste0("xxtemp2_out,", pref,"_numberofviews")))

rgrass7::execGRASS("g.copy", flags = c("quiet", "overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  raster = paste0("xxtemp3_out,", pref,"_pointofviews")))

rgrass7::execGRASS("g.copy", flags = c("quiet", "overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  raster = paste0("xxtemp4_out,", pref,"_distance_rescaled")))



path.viewangles <- "d:/Users/qo23hel/Desktop/Temp/GRASS TOOL/ESA_viewangles.tif"
path.numberofviews <- "d:/Users/qo23hel/Desktop/Temp/GRASS TOOL/ESA_numberofviews.tif"
path.pointofview <- "d:/Users/qo23hel/Desktop/Temp/GRASS TOOL/ESA_pointofview.tif"
path.distance <-  "d:/Users/qo23hel/Desktop/Temp/GRASS TOOL/ESA_distance.tif"

# print(parseGRASS("r.out.gdal"))
rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  input = paste0(pref, "_viewangles"), output = path.viewangles))

rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  input = paste0(pref, "_numberofviews"), output = path.numberofviews))

rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  input = paste0(pref, "_pointofviews"), output = path.pointofview))

rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
  input = paste0(pref, "_distance_rescaled"), output = path.distance))


# write data
#evaluation of the azimuth layer
# grass.mapcalc("azimuth = (450-aspect) - int( (450-aspect) / 360) * 360", overwrite=True, quiet=True)
# #evaluation of the layer of the vertical component of the versor perpendicular to the terrain slope
# grass.mapcalc("c_dem = cos(slope)", overwrite=True, quiet=True)
# #evaluation of the layer of the north component of the versor perpendicular to the terrain slope
# grass.mapcalc("b_dem = sin(slope)*cos(azimuth)", overwrite=True, quiet=True)
# #evaluation of the layer of the east component of the versor perpendicular to the terrain slope
# grass.mapcalc("a_dem = sin(slope)*sin(azimuth)", overwrite=True, quiet=True)

r.aspect <- raster::raster(path.aspect)
r.slope <- raster::raster(path.slope) %>%
           raster::calc(x =., fun = function(x){x*pi/180})

r.azimuth <- raster::calc(x = r.aspect, fun = function(x){(450-x) - as.integer((450-x)/360)* 360}) %>% # as.integer
             raster::calc(x =., fun = function(x){x*pi/180})

r.c_dem <-  raster::calc(x = r.slope, fun = function(x){cos(x)})
r.b_dem <-  raster::overlay(x = r.slope, y = r.azimuth, fun = function(x, y){sin(x) * cos(y)})
r.a_dem <-  raster::overlay(x = r.slope, y = r.azimuth, fun = function(x, y){sin(x) * sin(y)})

raster::writeRaster(x = r.azimuth, filename = file.path(path.test, "azimuth.tif"), overwrite = TRUE, NAflag = -99999)
raster::writeRaster(x = r.c_dem, filename = file.path(path.test, "c_dem.tif"), overwrite = TRUE, NAflag = -99999)
raster::writeRaster(x = r.b_dem, filename = file.path(path.test, "b_dem.tif"), overwrite = TRUE, NAflag = -99999)
raster::writeRaster(x = r.a_dem, filename = file.path(path.test, "a_dem.tif"), overwrite = TRUE, NAflag = -99999)


# https://sourceforge.net/p/saga-gis/discussion/790705/thread/338ba707/


# if( y()>520210 && x()>569369.1, atan((569369.1-x())/(5202104-y())),
#     if( y()<5202104 && x()>569369.1, 180+atan((569369.1-x())/(5202104-y())),
#         if( y()<5202104 && x()<569369.1, 180+atan((569369.1-x())/(5202104-y())),
#             if( y()>5202104 && x()<569369.1, 360+atan((569369.1-x())/(5202104-y()))))))


# pow(pow(abs(y()-5202104),2)+pow(abs(x()-569369.1),2)+pow(abs( raster_pnt_1@PERMANENT   -(341.39+1.75)),2),0.5)
