#' Initialisation of rgrass7 in Windows
#'
#' This funcion helps to initialise rgrass7 using windows.
#'
#' @param x \linkS4class{RasterLayer}  for region and projection. Default: NULL
#' @param osgeo4w.root root of OSGeo installation. Default: "C:\\OSGEO4W64"
#' @param grass.version GRASS GIS version. Default: "grass-7.2.2"
#' @param initGRASS.path initiation path of GRASS GIS through rgrass7. Default: "C:/OSGeo4W64/apps/grass/grass-7.2.2"
#' @param set.SysEnv set windows environmental path in R session. Default: TRUE. If FALSE no path is set.
#' @param use.link2GI use of link2GI package of initialisation. Default: FALSE
#' @param link2GI.defaultGrass initialisation using link2GI package. Default: c("C:/OSGeo4W64", "grass-7.2.2", "OSGeo4W64"))
#' @param quiet show output on console. Default: TRUE
#'
#' @keywords rgrass7, windows
#'
#'
#' @export
#'
initGRASSWin <- function(x = NULL, osgeo4w.root = "C:\\OSGEO4W64", grass.version = "grass-7.2.2",
                         initGRASS.path = "C:/OSGeo4W64/apps/grass/grass-7.2.2", set.SysEnv = TRUE,
                         use.link2GI = FALSE, link2GI.defaultGrass = c("C:/OSGeo4W64", "grass-7.2.2", "OSGeo4W64"), quiet = TRUE)
{
  if(is.null(x))
  {
    stop("RasterLayer as input is missing.")
  }


  if(set.SysEnv)
  {
    Sys.setenv(OSGEO4W_ROOT=osgeo4w.root)

    # define GISBASE
    grass.gis.base<-paste0(osgeo4w.root,"\\apps\\grass\\", grass.version)
    Sys.setenv(GISBASE=grass.gis.base)

    Sys.setenv(GRASS_PYTHON=paste0(Sys.getenv("OSGEO4W_ROOT"),"\\bin\\python.exe"))
    Sys.setenv(PYTHONHOME=paste0(Sys.getenv("OSGEO4W_ROOT"),"\\apps\\Python27"))
    Sys.setenv(PYTHONPATH=paste0(Sys.getenv("OSGEO4W_ROOT"),"\\apps\\grass\\", grass.version, "\\etc\\python"))
    Sys.setenv(GRASS_PROJSHARE=paste0(Sys.getenv("OSGEO4W_ROOT"),"\\share\\proj"))
    Sys.setenv(PROJ_LIB=paste0(Sys.getenv("OSGEO4W_ROOT"),"\\share\\proj"))
    Sys.setenv(GDAL_DATA=paste0(Sys.getenv("OSGEO4W_ROOT"),"\\share\\gdal"))
    Sys.setenv(GEOTIFF_CSV=paste0(Sys.getenv("OSGEO4W_ROOT"),"\\share\\epsg_csv"))
    Sys.setenv(FONTCONFIG_FILE=paste0(Sys.getenv("OSGEO4W_ROOT"),"\\etc\\fonts.conf"))


    Sys.setenv(PATH=paste0(grass.gis.base,";",
                           "C:\\OSGEO4~1\\apps\\Python27\\lib\\site-packages\\numpy\\core",";",
                           paste0("C:\\", basename(osgeo4w.root), "\\apps\\grass\\", grass.version, "\\bin",";"),
                           paste0("C:\\", basename(osgeo4w.root), "\\apps\\grass\\", grass.version, "\\lib",";"),
                           paste0("C:\\", basename(osgeo4w.root), "\\apps\\grass\\", grass.version, "\\etc",";"),
                           paste0("C:\\", basename(osgeo4w.root), "\\apps\\grass\\", grass.version, "\\etc\\python",";"),
                           paste0("C:\\", basename(osgeo4w.root), "\\apps\\Python27\\Scripts",";"),
                           paste0("C:\\", basename(osgeo4w.root), "\\bin",";"),
                           paste0("C:\\", basename(osgeo4w.root), "\\apps",";"),
                           paste0(Sys.getenv("WINDIR"),"/WBem"),";",
                           Sys.getenv("PATH")))
  }



  if(use.link2GI)
    {
      tryCatch({invisible(link2GI::linkGRASS7(x = elevation, defaultGrass))},
               error = function(e) {
                 cat("link2GI::linkGRASS7 throwed an Error... Trying grass7::initGrass() now...\n")

                 if(quiet)
                 {
                   tryCatch({invisible(rgrass7::initGRASS(initGRASS.path, home=tempdir(), override = TRUE,
                                                          SG = as(x, 'SpatialGridDataFrame')))},
                            error = function(err){
                              stop("Something wrong with the initialisation of GRASS GIS: ", err)
                            })
                 } else {
                   tryCatch({rgrass7::initGRASS(initGRASS.path, home=tempdir(), override = TRUE, mapset = "PERMANENT",
                                                          SG = as(x, 'SpatialGridDataFrame'))},
                            error = function(err){
                              stop("Something wrong with the initialisation of GRASS GIS: ", err)
                            })
                 }

               })
  } else {
    if(quiet)
    {
      tryCatch({invisible(rgrass7::initGRASS(initGRASS.path, home=tempdir(), override = TRUE, mapset = "PERMANENT"
                                             ))}, # SG = as(x, 'SpatialGridDataFrame')
               error = function(err){
                 stop("Something wrong with the initialisation of GRASS GIS: ", err)
               })
    } else {
      tryCatch({rgrass7::initGRASS(initGRASS.path, home=tempdir(), override = TRUE, mapset = "PERMANENT"
                                   )}, # SG = as(x, 'SpatialGridDataFrame')
               error = function(err){
                 stop("Something wrong with the initialisation of GRASS GIS: ", err)
               })
    }


    # get set projection and region
    resolution <- raster::res(x)[1]
    proj4 <- as.character(x@crs)
    ymax <- x@extent@ymax
    ymin <- x@extent@ymin
    xmax <- x@extent@xmax
    xmin <- x@extent@xmin

    # assign GRASS projection according to data set
    rgrass7::execGRASS('g.proj', flags = c('c','quiet'),  proj4 = proj4)

    # assign GRASS extent
    rgrass7::execGRASS('g.region', flags = c('quiet','d'),
                         n = as.character(ymax),
                         s = as.character(ymin),
                         e = as.character(xmax),
                         w = as.character(xmin),
                         res = as.character(resolution))

  if(quiet==FALSE){rgrass7::gmeta()}

 }
} # end function initGRASSWin


# EXAMPLE --------------------
# library(RQGIS, rgrass7)
# initGRASSWin(x = RQGIS::dem, quiet = FALSE)

