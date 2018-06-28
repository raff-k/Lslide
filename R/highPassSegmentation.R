#' Contrast-Filter-Segmentation
#'
#' This function performs a contrast filter segmentation on an input \linkS4class{RasterLayer}.
#' In the firest stept the raster is converted to a gray-scale matrix, followed by filter functions of the \code{\link[EBImage] package.
#' Then, the filtering is applied as mask on the segmentation-input. Afterwards, a segmentation is performed
#' using GRASS GIS region growing. Finally, the segments are transformed to \linkS4class{SpatialPolygonsDataFrame}
#' and returned as function output.
#' In the future this funcion will be included in the segmentation-function.
#'
#' @param input.filter \linkS4class{RasterLayer} on which the filtering is applied
#' @param input.segmentation \linkS4class{RasterLayer} for segmentation. Default: input.filter
#' @param offset Offset of values in filter matrix. Default: 0.06
#' @param makeBrush.size matrix size for smoothing input.filter. Default: 25
#' @param makeBrush.shape form of smoothin matrix. Default: "disc"
#' @param NA.val.in replace value for no data (NA) in \emph{input.filter}. Default: 0
#' @param fill.holes polygon segmentation output contains holes. Default: FALSE
#' @param clump.thresh area thresholds for removing clumbs - measured in cell sizes. Default: NULL
#' @param clump.directions directions for "clump-searching". Default: 8
#' @param env.rsaga environment of SAGA GIS. Default: RSAGA::rsaga.env()
#' @param CFR.buf buffer applied on contrast-filter-results, measured in cell-sizes. Default: "1"
#' @param output.path path for output files. Default: tempdir()
#' @param writeCFRaster write contrast-filter results. Default: FALSE
#' @param writeRaster.NAflag NoData values for NA values in raster. Default: -99999
#' @param freeMemory free memory after contrast filtering. Default: TRUE
#' @param show.output.on.console show output of tools on console. Default: FALSE
#' @param Fast.Representativeness.LevelOfGeneralisation level of generalisation for computing seed points. Default: "2.0"
#' @param quiet do not show comments and process time in console. Default: TRUE
#' @param Grass.Segmentation.Threshold threshold used of merging segments. Default: 0.24
#' @param Grass.Segmentation.Minsize minsize of a segment. Default: 0
#' @param Grass.Segmentation.Memory memory used in computation. Default: 1024
#' @param defaultGrass GRASS GIS environmet for temporal setup, see package link2GI. Default: c("C:/OSGeo4W64", "grass-7.2.2", "OSGeo4W64")
#' @param Segments.Poly output location for segmented objects (shapefile). Default: paste0(tempdir(), "/", "outSegPoly.shp")
#' @param Segments.Grid output location for segmented grid. Default: paste0(tempdir(), "/", "outSegGrid.sgrd")
#'
#' @return
#' \linkS4class{SpatialPolygonsDataFrame} containing the segments/objects
#'
#'
#' @note
#' \itemize{
#'   \item \emph{r} is stretched to gray-scale (0-255) during process by linear normalization
#' }
#'
#'
#' @keywords EBImage, raster, gray-scale-image, matrix, conversion, seeded-region growing, SAGA
#'
#'
#' @export
# input filter and input segmentation should have same extent and projection

highPassSegmentation <- function(input.filter, input.segmentation = input.filter, offset = 0.06, makeBrush.size = 25, makeBrush.shape = "disc", NA.val.in = 0, clump.thresh = NULL, clump.directions = 4, env.rsaga = RSAGA::rsaga.env(),
                                       CFR.buf = 1, output.path = tempdir(), writeCFRaster = FALSE, writeRaster.NAflag = -99999, freeMemory = TRUE, show.output.on.console = FALSE, quiet = TRUE, morph.Closing = TRUE, closing.size = 3, closing.shape = "box",
                                       Grass.Segmentation.Threshold = 0.24, Grass.Segmentation.Minsize = 1, Grass.Segmentation.Memory = 1024, Segments.Poly =  file.path(tempdir(), "outSegPoly.shp"),   Segments.Grid = file.path(tempdir(), "outSegGrid.sgrd"),
                                       defaultGrass = c("C:/OSGeo4W64", "grass-7.2.2", "OSGeo4W64"), load.output = FALSE, fill.holes = FALSE, ...)
{

  # browser()

  if(quiet == FALSE) process.time.start <- proc.time()
  if(quiet == FALSE) cat("Start High-Pass-Segmentation ...\n")

  # Start high-pass filter -----------------
  hipass <- Lslide::highPassThresholding()


  # mask segmentation input to filter
  if(quiet == FALSE) cat("... clip input for segmentation based on filter results\n")
  input.segmentation <- raster::overlay(input.segmentation, hipass, fun = function(x, y){return(ifelse(y > 0, x, NA))})



  # Start segmentation -----------------
  if(quiet == FALSE) cat("... GRASS: Start segmentation\n")


  ## initialize variables for segmentation in GRASS GIS
  input.segmentation.path <- file.path(tempdir(), "tmp_segInp.tif")
  raster::writeRaster(x = input.segmentation, filename = input.segmentation.path, overwrite = TRUE, NAflag = writeRaster.NAflag)

  # print(rgrass7::parseGRASS("r.in.gdal"))
  rgrass7::execGRASS("r.in.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    input = input.segmentation.path, output = "inputGrass"))



  # create an imagery group
  rgrass7::execGRASS("i.group", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    group = "GRASS.Segmentation.Group", input = "inputGrass"))

  # browser()
  # start segmentation
  # print(parseGRASS("i.segment"))
  rgrass7::execGRASS("i.segment", group = "GRASS.Segmentation.Group", flags = c("overwrite"), output = "output.segmentation.GRASS", threshold = Grass.Segmentation.Threshold,
                     memory = Grass.Segmentation.Memory, minsize = Grass.Segmentation.Minsize, Sys_show.output.on.console = show.output.on.console)


  # get data from GRASS
  # print(parseGRASS("r.out.gdal"))
  if(tools::file_ext(Segments.Grid) == "sgrd")
  {
   #  browser()

   Segments.Grid.tmp <- file.path(tempdir(), paste0(tools::file_path_sans_ext(basename(Segments.Grid)), ".tif"))


    rgrass7::execGRASS("r.out.gdal", flags = c("overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "output.segmentation.GRASS", type = "Int32", output =  Segments.Grid.tmp, format = "GTiff",
      nodata =  writeRaster.NAflag))


    RSAGAUsage <- RSAGA::rsaga.get.usage(lib="io_gdal", module = 1, env = env.rsaga,  show = FALSE)
    formatSAGA <- gsub("\\D", "", grep('SAGA GIS Binary', RSAGAUsage, value = TRUE))

    RSAGA::rsaga.geoprocessor(lib = "io_gdal", module = 1, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
      GRIDS = Segments.Grid.tmp, FILE =  paste0(tools::file_path_sans_ext(Segments.Grid), ".sdat"), FORMAT =  formatSAGA))


  } else {

    rgrass7::execGRASS("r.out.gdal", flags = c("overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "output.segmentation.GRASS", type = "Int32", output =  paste0(tools::file_path_sans_ext(Segments.Grid), ".tif"), format = "GTiff",
      nodata =  writeRaster.NAflag))
    }

  # vectorising grid classes
  if(quiet == FALSE) cat("... GRASS: Vectorising Grid Classes\n")
  # shapes_grid
  # rsaga.get.modules("shapes_grid", env = env.rsaga)
  # rsaga.get.usage("shapes_grid", 6, env = env.rsaga)
  # RSAGA::rsaga.geoprocessor(lib = "shapes_grid", module = 6, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
  #   GRID = paste0(tools::file_path_sans_ext(Segments.Grid), ".tif"), POLYGONS =  Segments.Poly))
  # print(parseGRASS("r.to.vect"))
  rgrass7::execGRASS("r.to.vect", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    input = "output.segmentation.GRASS", output = "Segments_Poly", type = "area"))


  # print(parseGRASS("v.out.ogr"))
  # what to do with holes?
  if(fill.holes)
  {
    rgrass7::execGRASS("v.out.ogr", flags = c("overwrite", "quiet", "c"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "Segments_Poly", output = Segments.Poly, type = "area", format = "ESRI_Shapefile"))
  } else {
    rgrass7::execGRASS("v.out.ogr", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "Segments_Poly", output = Segments.Poly, type = "area", format = "ESRI_Shapefile"))
  }



  if(load.output == TRUE)
  {
    if(quiet == FALSE) cat("... Loading Segments: Using sf::() and conversion to SpatialPolygonsDataFrame\n")
    outpoly <- sf::st_read(Segments.Poly, quiet = !show.output.on.console)
    outpoly <- as(outpoly, "Spatial")

    if(is.na(sp::proj4string(outpoly)))
    {
      sp::proj4string(outpoly) <- raster::projection(input.filter)
    }

    if(quiet == FALSE) cat(paste0("------ Run of contrastFilterSegmentation: " , (proc.time() - process.time.start)["elapsed"][[1]]/60, " Minutes ------\n"))

    return(outpoly)
  }

  if(quiet == FALSE) cat(paste0("------ Run of contrastFilterSegmentation: " , (proc.time() - process.time.start)["elapsed"][[1]]/60, " Minutes ------\n"))

} # end function highPassSegmentation





# EXAMPLE -----------------------------------
# library("RQGIS", "reticulate")
# open_app()
# pacman::p_load("EBImage", "raster", "data.table", "igraph", "RQGIS", "reticulate", "RSAGA", "rgrass7", "link2GI")
#
# env.rsaga <- RSAGA::rsaga.env(path = "C:/Program Files (x86)/SAGA-GIS")
#
# setwd("E:/Masterarbeit/Data/Output/Segmentation/EBImage")
# # save.image("cFS_run.RDat")
# load("cFS_run.RDat")
# slope <- raster("Input/OP14_slp_01m.tif")
#
# dem <- RQGIS::dem
# slope <- RQGIS::run_qgis(alg  = "grass7:r.slope.aspect", elevation = dem, slope = "slope.tif", load_output = TRUE)
#
#
#
# # ------ Run of contrastFilterSegmentation: 9.63516666666667 Minutes ------
# segments <- contrastFilterSegmentation(input.filter = slope, makeBrush.size = 25, offset = 0.07, env.rsaga = env.rsaga, show.output.on.console = FALSE, quiet = FALSE,
#                                        Grass.Segmentation.Threshold =  0.02, Grass.Segmentation.Minsize = 50, clump.directions = 4, clump.thresh = 50, freeMemory = TRUE,
#                                        Segments.Poly = paste0(getwd(), "/Output/cFS_Mb25Of007_ClD4T50_SgT002Sz50.shp"), Segments.Grid = paste0(getwd(), "/Output/cFS_Mb25Of007_ClD4T50_SgT002Sz50.sgrd"))
#
#


# # testing:
# save.image("cFS_Debugging_b4Segmentaion.RData")
# load("cFS_Debugging_b4Segmentaion.RData")
# input.filter = slope
# makeBrush.size = 25
# env.rsaga = env.rsaga
# show.output.on.console = FALSE
# quiet = FALSE
# input.segmentation = input.filter
# offset = 0.07
# makeBrush.shape = "disc"
# NA.val.in = 0
# clump.thresh = 50
# clump.directions = 4
# CFR.buf = "1"
# output.path = tempdir()
# writeCFRaster = FALSE
# writeRaster.NAflag = -99999
# freeMemory = TRUE
# # Saga.Segmentation.Sig.2 = "125"
# # Saga.Segmentation.Leafsize = 1024
# Grass.Segmentation.Threshold = 0.02
# Grass.Segmentation.Minsize = 50
# Grass.Segmentation.Memory = 1024
# Segments.Poly =  paste0(tempdir(),  "outSegPoly.shp")
# Segments.Grid = paste0(tempdir(), "outSegGrid.sgrd")
# defaultGrass = c("C:/OSGeo4W64", "grass-7.2.2", "OSGeo4W64")

# debug
# debug.InSeedsMasked <- raster::raster(paste0(tools::file_path_sans_ext(Output.Seeds.Masked), ".sdat"))
# raster::writeRaster(x = debug.InSeedsMasked, filename = paste0(getwd(), "/Debug/seedsMaksed.tif"))
# raster::writeRaster(x = input.segmentation, filename = paste0(getwd(), "/Debug/inputSegmentation.tif"))
#
# debug.SegmentsGrid <- raster::raster(paste0(tools::file_path_sans_ext(Segments.Grid), ".sdat"))
# raster::writeRaster(x = debug.SegmentsGrid, filename = paste0(getwd(), "/Debug/outputSegmentation_SAGA.tif"))
#
#
#
# input.segmentation.NA <- raster::writeRaster(input.segmentation, paste0(getwd(),"/Debug/inputSegmentationTest.tif"), NAflag = -9999)
# input.segmentation.NA <- raster::raster(input.segmentation.NA)
# segmentation.GRASS <-  RQGIS::run_qgis(alg = "grass7:i.segment", input = input.segmentation.NA, threshold = "0.05", minsize = "50", iterations = "20",
#                                        output = "tmpGrassSegmt.tif", load_output = TRUE, show_output_paths = TRUE)

