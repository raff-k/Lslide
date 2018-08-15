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
#' @param env.rsaga environment of SAGA GIS. Default: RSAGA::rsaga.env()
#' @param writeRaster.NAflag NoData values for NA values in raster. Default: -99999
#' @param show.output.on.console show output of tools on console. Default: FALSE
#' @param Fast.Representativeness.LevelOfGeneralisation level of generalisation for computing seed points. Default: "2.0"
#' @param quiet do not show comments and process time in console. Default: TRUE
#' @param Grass.Segmentation.Threshold threshold used of merging segments. Default: 0.24
#' @param Grass.Segmentation.Minsize minsize of a segment. Default: 0
#' @param Grass.Segmentation.Memory memory used in computation. Default: 1024
#' @param Segments.Poly output location for segmented objects (shapefile). Default: paste0(tempdir(), "/", "outSegPoly.shp")
#' @param Segments.Grid output location for segmented grid. Default: paste0(tempdir(), "/", "outSegGrid.sgrd")
#' @param estimateScaleParameter must be only be TRUE when scale parameter function is used. Default: FALSE
#' @param ... parameters for function hiPassThresh
#'
#' @return
#' \linkS4class{SpatialPolygonsDataFrame} containing the segments/objects
#'
#'
#' @keywords SAGA GIS, raster, high-pass filtering, matrix, conversion, seeded-region growing, GRASS GIS
#'
#'
#' @export
# input filter and input segmentation should have same extent and projection

hiPassSegmentation <- function(input.filter, input.segmentation = input.filter, result.filter = NULL, scale.factor, threshold, env.rsaga = RSAGA::rsaga.env(),
                                       writeRaster.NAflag = -99999, show.output.on.console = FALSE, quiet = TRUE,
                                       Grass.Segmentation.Threshold = 0.24, Grass.Segmentation.Minsize = 1, Grass.Segmentation.Memory = 1024, Segments.Poly =  file.path(tempdir(), "outSegPoly.shp"),   Segments.Grid = file.path(tempdir(), "outSegGrid.sgrd"),
                                       load.output = FALSE, fill.holes = FALSE, estimateScaleParameter = FALSE, ...)
{

  # browser()

  if(quiet == FALSE) process.time.start <- proc.time()
  if(quiet == FALSE) cat("Start High-Pass-Segmentation ...\n")

  # Start high-pass filter -----------------
  if(!is.null(result.filter) && class(result.filter) == "RasterLayer")
  {
    hipass <- result.filter
    if(quiet == FALSE) cat('... input of "result.filter" is used as high-pass filter ...\n')
  } else {
    hipass <- Lslide::hiPassThresh(x = input.filter, scale.factor, threshold, NoData = writeRaster.NAflag ,
                                   env.rsaga = env.rsaga, show.output.on.console = show.output.on.console, quiet = quiet,...)
  }


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
  rgrass7::execGRASS("i.segment", group = "GRASS.Segmentation.Group", flags = c("overwrite"), output = "output_seg", threshold = Grass.Segmentation.Threshold,
                     memory = Grass.Segmentation.Memory, minsize = Grass.Segmentation.Minsize, Sys_show.output.on.console = show.output.on.console)


  # get data from GRASS
  # print(parseGRASS("r.out.gdal"))
  if(tools::file_ext(Segments.Grid) == "sgrd")
  {
   #  browser()

   Segments.Grid.tmp <- file.path(tempdir(), paste0(tools::file_path_sans_ext(basename(Segments.Grid)), ".tif"))


    rgrass7::execGRASS("r.out.gdal", flags = c("overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "output_seg", type = "Int32", output =  Segments.Grid.tmp, format = "GTiff",
      nodata =  writeRaster.NAflag))


    RSAGAUsage <- RSAGA::rsaga.get.usage(lib="io_gdal", module = 1, env = env.rsaga,  show = FALSE)
    formatSAGA <- gsub("\\D", "", grep('SAGA GIS Binary', RSAGAUsage, value = TRUE))

    RSAGA::rsaga.geoprocessor(lib = "io_gdal", module = 1, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
      GRIDS = Segments.Grid.tmp, FILE =  paste0(tools::file_path_sans_ext(Segments.Grid), ".sdat"), FORMAT =  formatSAGA))


  } else {

    rgrass7::execGRASS("r.out.gdal", flags = c("overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "output_seg", type = "Int32", output =  paste0(tools::file_path_sans_ext(Segments.Grid), ".tif"), format = "GTiff",
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
    input = "output_seg", output = "Segments_Poly", type = "area"))


  # print(parseGRASS("v.out.ogr"))
  # what to do with holes?
  if(fill.holes)
  {
    rgrass7::execGRASS("v.out.ogr", flags = c("overwrite", "quiet", "c"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "Segments_Poly", output = Segments.Poly, type = "area", format = "ESRI_Shapefile"))

    if(estimateScaleParameter)
    {
      rgrass7::execGRASS("v.in.ogr", flags = c("overwrite", "quiet", "c"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        input = Segments.Poly, output = "Segments_Poly"))
    }

  } else {

    if(!estimateScaleParameter)
    {
      rgrass7::execGRASS("v.out.ogr", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        input = "Segments_Poly", output = Segments.Poly, type = "area", format = "ESRI_Shapefile"))
    }
  }



  if(load.output == TRUE && !estimateScaleParameter)
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


