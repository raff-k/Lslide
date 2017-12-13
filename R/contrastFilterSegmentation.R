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

contrastFilterSegmentation <- function(input.filter, input.segmentation = input.filter, offset = 0.06, makeBrush.size = 25, makeBrush.shape = "disc", NA.val.in = 0, clump.thresh = NULL, clump.directions = 4, env.rsaga = RSAGA::rsaga.env(),
                                       CFR.buf = 1, output.path = tempdir(), writeCFRaster = FALSE, writeRaster.NAflag = -99999, freeMemory = TRUE, show.output.on.console = FALSE, quiet = TRUE, morph.Closing = TRUE, closing.size = 3, closing.shape = "box",
                                       Grass.Segmentation.Threshold = 0.24, Grass.Segmentation.Minsize = 0, Grass.Segmentation.Memory = 1024, Segments.Poly =  paste0(tempdir(), "/", "outSegPoly.shp"),   Segments.Grid = paste0(tempdir(), "/", "outSegGrid.sgrd"),
                                       defaultGrass = c("C:/OSGeo4W64", "grass-7.2.2", "OSGeo4W64"), load.output = FALSE, fill.holes = FALSE, ...)
{

  # browser()

  if(quiet == FALSE) process.time.start <- proc.time()
  if(quiet == FALSE) cat("Start Contrast-Filter-Segmentation ...\n")


  # creation of gray scale image
  if(quiet == FALSE) cat("... convert input raster to grey-scale-matrix\n")
  gsm <- Lslide::convR2GSM(r = input.filter, NA.val.in = NA.val.in)


  # Start contrast filter -----------------
  # create brush and blurr image for contrast filter - using the fast 2D FFT convolution product.
  if(quiet == FALSE) cat("... apply 2d FFT convolution filter on matrix\n")
  brush <- EBImage::makeBrush(size = makeBrush.size, shape = makeBrush.shape, ...)
  brush  <- brush/sum(brush)
  gsm.output <- EBImage::filter2(x = gsm, filter = brush, boundary = "replicate")


  # conversion from matrix to raster
  if(quiet == FALSE) cat("... conversion of filter matrix and grey-scale-matrix to raster\n")
  gsm <- raster::raster(gsm)
  extent(gsm) <- extent(input.filter)
  # projection(gsm) <- projection(input.filter)

  gsm.output <- raster::raster(gsm.output)
  extent(gsm.output) <- extent(input.filter)
  # projection(gsm.output) <- projection(input.filter)


  # selection based on offset
  if(quiet == FALSE) cat("... selection based on offset\n")
  gsm.output <- raster::overlay(gsm, gsm.output , fun = function(x, y){return(ifelse(x > (y + offset), 1, 0))})


  # removal of clumps
  if(!is.null(clump.thresh) && clump.thresh > 0)
  {
    # require(igraph)

    if(quiet == FALSE) cat("... removal of clumbs based on threshold: ", clump.thresh, "\n")
    if(quiet == FALSE) cat("... ... calculate clumps\n")

    clumps <- raster::clump(gsm.output, directions = clump.directions, ...) # need some time

    # calculate pixel frequency for each clumpID
    clumpFreq <- data.table::as.data.table(freq(clumps))

    # clumpID to be excluded from output raster
    excludeVal <- clumpFreq[count <= clump.thresh, value]

    # removal
    if(quiet == FALSE) cat("... ... remove clumps under the threshold area\n")
    gsm.output <- raster::overlay(gsm.output, clumps, fun = function(x, y)
                                        {
                                            x[y %in% excludeVal] <- 0
                                            return(x)
                                            })

    rm(excludeVal, clumpFreq, clumps)

  } # end if doing clump removal

  # raster::writeRaster(x = gsm.output, filename = file.path(getwd(), "Segmentation", "Closing", "gsmOut.tif"), overwrite=TRUE)


  if(morph.Closing)
  {
    gsm.output <- Lslide::convR2GSM(r = gsm.output, NA.val.in = NA)
    if(quiet == FALSE) cat("... morphological filter: closing\n")
    gsm.output <- EBImage::closing(x = gsm.output, kern = EBImage::makeBrush(closing.size, shape = closing.shape))
    gsm.output <- raster::raster(gsm.output)
    extent(gsm.output) <- extent(input.filter)
  }


  # project raster
  projection(gsm.output) <- projection(input.filter)

  if(writeCFRaster == TRUE)
  {
    if(quiet == FALSE) cat("... write contrast filter raster\n")

    if(!is.null(clump.thresh) && clump.thresh > 0)
    {
      raster::writeRaster(x = gsm.output, filename = paste0(output.path, '/', "CFRaster_", sub("[.]", "" ,as.character(offset)), "Ofs_", makeBrush.size, "Bs_", clump.thresh, "Thrs",  ".tif"), overwrite = TRUE)

    } else{
      raster::writeRaster(x = gsm.output, filename = paste0(output.path, '/', "CFRaster_", sub("[.]", "" ,as.character(offset)), "Ofs_", makeBrush.size, "Bs", ".tif"), overwrite = TRUE)

    }

  } # end if writeCFRaster


  if(freeMemory == TRUE)
  {
    if(quiet == FALSE) cat("... free memory\n")
    rm(gsm)
    invisible(gc())
  }


  gsm.output <- raster::calc(gsm.output, fun = function(x){return(ifelse(x == 0, NA, x))})

  if(!is.null(CFR.buf) && CFR.buf > 0)
  {
    # Start buffering contrast filter results -----------------
    if(quiet == FALSE) cat("... buffer contrast filter results\n")

    # RQGIS::find_algorithms(search_term = "buf")
    # RQGIS::get_args_man("saga:gridbuffer")
    # RQGIS::get_usage("saga:gridbuffer")
    # BUFFERTYPE: 1 - [1] Cell value
    # gsm.output <- RQGIS::run_qgis(alg = "saga:gridbuffer", FEATURES = gsm.output, DIST = CFR.buf, BUFFERTYPE = "1", BUFFER = "tmpBuf.tif",
    #                                  load_output = TRUE, show_output_paths = FALSE)

    # rsaga.get.modules("grid_tools", env = env.rsaga)
    # rsaga.get.usage("grid_tools", 8, env = env.rsaga)

    # RSAGA::rsaga.geoprocessor(lib = "io_gdal", module = 1, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
    #   GRIDS = Segments.Grid.tmp, FILE =  paste0(tools::file_path_sans_ext(Segments.Grid), ".sdat"), FORMAT =  "43"))


    # raster::writeRaster(x = raster::raster(Segments.Grid.tmp), filename = paste0(tools::file_path_sans_ext(Segments.Grid), ".sdat"),
    #                     overwrite = TRUE)

   #  browser()
#
#     gsm.output.path <- paste0(tempdir(), "/", "tmpGsmOut.sgrd")
#     gsm.buf.path <- paste0(tempdir(),"/", "tmpBuf.sgrd")
#     # gsm.buf.path <- paste0(tempdir(),"/", "tmpBuf.tif")
    # gsm.output.path <- paste0(tempdir(), "/", "tmpGsmOut.tif")
    # browser()
    # raster::writeRaster(x = gsm.output, filename = gsm.output.path, overwrite = TRUE, NAflag = writeRaster.NAflag)
    # rgdal::writeGDAL(dataset = as(gsm.output, "SpatialGridDataFrame"), fname = paste0(tools::file_path_sans_ext(gsm.output.path), ".sdat"),
    #                   drivername = "SAGA")


    # RSAGA::rsaga.geoprocessor(lib="grid_tools", module = 8, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
    #   FEATURES = gsm.output.path, BUFFER = gsm.buf.path, DIST = CFR.buf, BUFFERTYPE = "1"))

    # RSAGA::rsaga.geoprocessor(lib="grid_tools", module = 8, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
    #   FEATURES = gsm.output.path, BUFFER = gsm.buf.path, DISTANCE = CFR.buf, TYPE = "1"))

    gsm.buf.path <- paste0(tempdir(),"/", "tmpBuf.tif")

    rgrass7::writeRAST(x = as(gsm.output, "SpatialGridDataFrame"), vname = "tmpGsm", zcol = names(gsm.output),
                       overwrite = TRUE, flags = 'quiet')

    # print(rgrass7::parseGRASS("r.buffer"))
    rgrass7::execGRASS("r.buffer", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "tmpGsm", output = "tmpGsmBuf", distances = CFR.buf))


    rgrass7::execGRASS("r.out.gdal", flags = c("overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "tmpGsmBuf", type = "Int32", output = gsm.buf.path, format = "GTiff",
      nodata =  writeRaster.NAflag))


    gsm.output <- raster::raster(gsm.buf.path)
  }



  # mask segmentation input to filter
  if(quiet == FALSE) cat("... clip input for segmentation based on filter results\n")
  input.segmentation <- raster::overlay(input.segmentation, gsm.output, fun = function(x, y){return(ifelse(y > 0, x, NA))})



  # Start segmentation -----------------
  if(quiet == FALSE) cat("... GRASS: Start segmentation\n")


  ## initialize variables for segmentation in GRASS GIS



  # tryCatch({invisible(link2GI::linkGRASS7(x = input.segmentation, defaultGrass = defaultGrass))},
  #          error = function(e) {
  #   stop("Something wrong with the initialisation of GRASS GIS: ", e)
  # })


  # load data into GRASS
  rgrass7::writeRAST(as(input.segmentation, 'SpatialGridDataFrame'), "inputGrass",
                     zcol = "layer", useGDAL = TRUE, flags = c("overwrite"))


  # create an imagery group
  rgrass7::execGRASS("i.group", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    group = "GRASS.Segmentation.Group", input = "inputGrass"))

  # start segmentation
  # print(parseGRASS("i.segment"))
  rgrass7::execGRASS("i.segment", group = "GRASS.Segmentation.Group", flags = c("overwrite"), output = "output.segmentation.GRASS", threshold = Grass.Segmentation.Threshold,
                     memory = Grass.Segmentation.Memory, minsize = Grass.Segmentation.Minsize, Sys_show.output.on.console = show.output.on.console)


  # get data from GRASS
  # print(parseGRASS("r.out.gdal"))
  if(tools::file_ext(Segments.Grid) == "sgrd")
  {
   Segments.Grid.tmp <- paste0(tempdir(), "/", tools::file_path_sans_ext(basename(Segments.Grid)), ".tif")

    rgrass7::execGRASS("r.out.gdal", flags = c("overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "output.segmentation.GRASS", type = "Int32", output =  Segments.Grid.tmp, format = "GTiff",
      nodata =  writeRaster.NAflag))

    # rgrass7::execGRASS("r.out.gdal", flags = c("overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    #   input = "output.segmentation.GRASS", type = "Int32", output =  paste0(tools::file_path_sans_ext(Segments.Grid), ".sdat"), format = "SAGA",
    #   nodata =  writeRaster.NAflag))



    RSAGA::rsaga.geoprocessor(lib = "io_gdal", module = 1, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
      GRIDS = Segments.Grid.tmp, FILE =  paste0(tools::file_path_sans_ext(Segments.Grid), ".sdat"), FORMAT =  "43"))

#
#     raster::writeRaster(x = raster::raster(Segments.Grid.tmp), filename = paste0(tools::file_path_sans_ext(Segments.Grid), ".sdat"),
#                         overwrite = TRUE)

  } else {

    rgrass7::execGRASS("r.out.gdal", flags = c("overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "output.segmentation.GRASS", type = "Int32", output =  paste0(tools::file_path_sans_ext(Segments.Grid), ".tif"), format = "GTiff",
      nodata =  writeRaster.NAflag))
    }


  # Output.Seeds <- paste0(tempdir(),"/", "tmpSeeds.sgrd")
  # Output.Seeds.Masked <- paste0(tempdir(),"/", "tmpSeedsMasked.sgrd")
  # write data to temp
  # Input.Grid <- paste0(tempdir(),"/tmpInputSeg.sgrd")
  # raster::writeRaster(x = input.segmentation, filename = Input.Grid, overwrite = TRUE, NAflag = writeRaster.NAflag) # filtered and clipped input!

  # GRASS
  # RQGIS::find_algorithms("segment")
  # RQGIS::get_args_man("grass7:i.group")
  # RQGIS::get_usage("grass7:i.group")
  # input.segmentation.group <- RQGIS::run_qgis(alg = "grass7:i.group", input = base::list(input.segmentation),
  #                                             load_output = TRUE, show_output_paths = TRUE)
  # segmentation.GRASS <-  RQGIS::run_qgis(alg = "grass7:i.segment", input = input.segmentation, threshold = "0.05", minsize = "50", memory = "5000", iterations = "20",
  #                                       output = "tmpGrassSegmt.tif", load_output = TRUE, show_output_paths = TRUE)



  # TO DO!!!!
  # Lslide::segmentation(Tool = "SAGA", Segments.Grid = segmentation.first.grid, Segments.Poly = segmentation.first, Input.Grid = c(slope, slope.edge.vigra),
  #              Seed.Method = "Fast Representativeness", Fast.Representativeness.LevelOfGeneralisation = 1.25,  Saga.Segmentation.Method = "0", Saga.Segmentation.Sig.2 = "125",
  #              Generalisation.Flac = FALSE, Saga.Segmentation.Leafsize = 1024, NoData = TRUE, Mask = slope)

  # rsaga.get.modules("statistics_grid", env = env.rsaga)
  # rsaga.get.usage("statistics_grid", 0, env = env.rsaga)
  # SAGA SEEDS NOT WORKING!
  # rsaga.geoprocessor(lib="statistics_grid", module = 0, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
  #   INPUT = Input.Grid, RESULT = paste0(tempdir(),"/", "tmpResultFR"), RESULT_LOD = paste0(tempdir(),"/", "tmpLod"), SEEDS = Output.Seeds,
  #   LOD = Fast.Representativeness.LevelOfGeneralisation))


  # to do mask seeds
  # if(quiet == FALSE) cat("... SAGA: Masking of Seed Points in NoData area\n")
  # grid tools - grid masking
  # rsaga.get.modules("grid_tools", env = env.rsaga)
  # rsaga.get.usage("grid_tools", 24, env = env.rsaga)
  # rsaga.geoprocessor(lib = "grid_tools", module = 24, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
  #   GRID = Output.Seeds, MASK = Input.Grid, MASKED = Output.Seeds.Masked))


  # seeded region growing
  # if(quiet == FALSE) cat("... SAGA: Seeded Region Growing\n")
  # rsaga.get.modules("imagery_segmentation", env = env.rsaga)
  # rsaga.get.usage("imagery_segmentation", 3, env = env.rsaga)
   # SRG <- rsaga.geoprocessor(lib = "imagery_segmentation", module = 3, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
   #    SEEDS = Output.Seeds.Masked, FEATURES = Input.Grid, SIG_2 = Saga.Segmentation.Sig.2, SEGMENTS = Segments.Grid,
   #       LEAFSIZE = Saga.Segmentation.Leafsize))

   # if(any(grepl('Error', SRG)))
   # {
   #   stop("Error in Seeded Region Growing (SAGA). Check input and control number of seed points.")
   # }

  # to do mask region growing
  # if(quiet == FALSE) cat("... SAGA: Masking of Segments in NoData area\n")
  # grid tools - grid masking
  # rsaga.get.modules("grid_tools", env = env.rsaga)
  # rsaga.get.usage("grid_tools", 24, env = env.rsaga)
  # rsaga.geoprocessor(lib = "grid_tools", module = 24, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
  #   GRID = Segments.Grid, MASK = Input.Grid, MASKED = Segments.Grid))


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
      proj4string(outpoly) <- projection(input.filter)
    }

    if(quiet == FALSE) cat(paste0("------ Run of contrastFilterSegmentation: " , (proc.time() - process.time.start)["elapsed"][[1]]/60, " Minutes ------\n"))

    return(outpoly)
  }

  if(quiet == FALSE) cat(paste0("------ Run of contrastFilterSegmentation: " , (proc.time() - process.time.start)["elapsed"][[1]]/60, " Minutes ------\n"))

} # end function contrastFilterSegmentation





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

