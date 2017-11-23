#' Contrast-Filter-Segmentation
#'
#' This function performs a contrast filter segmentation on an input \linkS4class{RasterLayer}.
#' In the firest stept the raster is converted to a gray-scale matrix, followed by filter functions of the \code{\link[EBImage] package.
#' Then, the filtering is applied as mask on the segmentation-input. Afterwards, a segmentation is performeds
#' as combination from fast-representativeness and seeded region growing. Finally, the segments are transformed
#' to \linkS4class{SpatialPolygonsDataFrame} and returned as function output.
#' In the future this funcion will be included in the segmentation-function.
#'
#' @param input.filter \linkS4class{RasterLayer} on which the filtering is applied
#' @param input.segmentation \linkS4class{RasterLayer} for segmentation. Default: input.filter
#' @param offset Offset of values in filter matrix. Default: 0.06
#' @param makeBrush.size matrix size for smoothing input.filter. Default: 25
#' @param makeBrush.shape form of smoothin matrix. Default: "disc"
#' @param NA.val.in replace value for no data (NA) in \emph{input.filter}. Default: 0
#' @param clump.thresh area thresholds for removing clumbs - measured in cell sizes. Default: NULL
#' @param clump.directions directions for "clump-searching". Default: 8
#' @param env.rsaga environment of SAGA GIS. Default: RSAGA::rsaga.env()
#' @param CFR.buf buffer applied on contrast-filter-results, measured in cell-sizes. Default: "1"
#' @param output.path path for output files. Default: tempdir()
#' @param writeCFRaster write contrast-filter results. Default: FALSE
#' @param writeRaster.NAflag NoData values for NA values in raster. Default: -99999
#' @param freeMemory free memory after contrast filtering. Default: TRUE
#' @param SOOC show output of tools on console. Default: FALSE
#' @param Fast.Representativeness.LevelOfGeneralisation level of generalisation for computing seed points. Default: "2.0"
#' @param quiet do not show comments and process time in console. Default: TRUE
#' @param Saga.Segmentation.Sig.2 variance of feature space in seeded region growing. Default: "125"
#' @param Saga.Segmentation.Leafsize memory for computing SAGA GIS algorithm. Default: 1024
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

segments <- contrastFilterSegmentation(input.filter = slope, makeBrush.size = 11, env.rsaga = env.rsaga, SOOC = FALSE, quiet = FALSE,
                                       Fast.Representativeness.LevelOfGeneralisation = "3", clump.thresh = 10)

# input filter and input segmentation should have same extent and projection

contrastFilterSegmentation <- function(input.filter, input.segmentation = input.filter, offset = 0.06, makeBrush.size = 25, makeBrush.shape = "disc", NA.val.in = 0, clump.thresh = NULL, clump.directions = 8, env.rsaga = RSAGA::rsaga.env(),
                                       CFR.buf = "1", output.path = tempdir(), writeCFRaster = FALSE, writeRaster.NAflag = -99999, freeMemory = TRUE, SOOC = FALSE, Fast.Representativeness.LevelOfGeneralisation = "2.0", quiet = TRUE,
                                       Saga.Segmentation.Sig.2 = "125", Saga.Segmentation.Leafsize = 1024, Segments.Poly =  paste0(tempdir(), "/", "outSegPoly.shp"),   Segments.Grid = paste0(tempdir(), "/", "outSegGrid.sgrd"), ...)
{

  if(quiet == FALSE) process.time.start <- proc.time()
  if(quiet == FALSE) cat("Start Contrast-Filter-Segmentation ...\n")


  # creation of gray scale image
  if(quiet == FALSE) cat("... convert input raster to grey-scale-matrix\n")
  gsm <- convR2GSM(r = input.filter, NA.val.in = NA.val.in)


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
  # projection(gsm) <- projection(iinput.filter)

  gsm.output <- raster::raster(gsm.output)
  extent(gsm.output) <- extent(input.filter)
  # projection(gsm.output) <- projection(input.filter)


  # selection based on offset
  if(quiet == FALSE) cat("... selection based on offset\n")
  gsm.output <- raster::overlay(gsm, gsm.output , fun = function(x, y){return(ifelse(x > (y + offset), 1, 0))})


  # removal of clumps
  if(!is.null(clump.thresh) && clump.thresh > 0)
  {
    if(quiet == FALSE) cat("... removal of clumbs based on threshold: ", clump.thresh, "\n")
    if(quiet == FALSE) cat("... ... calculate clumps\n")
    clumps <- raster::clump(gsm.output, directions = clump.directions, ...) # need some time

    # calculate pixel frequency for each clumpID
    clumpFreq <- as.data.table(freq(clumps))

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
    gsm.output.path <- paste0(tempdir(), "/", "tmpGsmOut.sgrd")
    gsm.buf.path <- paste0(tempdir(),"/", "tmpBuf.sgrd")
    writeRaster(x = gsm.output, filename = gsm.output.path, overwrite = TRUE, NAflag = writeRaster.NAflag)

    rsaga.geoprocessor(lib="grid_tools", module = 8, env = env.rsaga, show.output.on.console = SOOC, param = list(
      FEATURES = gsm.output.path, BUFFER = gsm.buf.path, DIST = CFR.buf, BUFFERTYPE = "1"))

    gsm.output <- raster::raster(gsm.output.path)
  }

  # mask segmentation input to filter
  if(quiet == FALSE) cat("... clip input for segmentation based on filter results\n")
  input.segmentation <- raster::overlay(input.segmentation, gsm.output, fun = function(x, y){return(ifelse(y > 0, x, NA))})

  ## initialize variables for segmentation in SAGA GIS
  Output.Seeds <- paste0(tempdir(),"/", "tmpSeeds.sgrd")
  Output.Seeds.Masked <- paste0(tempdir(),"/", "tmpSeedsMasked.sgrd")
  # write data to temp
  Input.Grid <- paste0(tempdir(),"/tmpInputSeg.sgrd")
  raster::writeRaster(x = input.segmentation, filename = Input.Grid, overwrite = TRUE, NAflag = writeRaster.NAflag) # filtered and clipped input!


  # Start segmentation -----------------
  if(quiet == FALSE) cat("... SAGA: Fast Representativeness - for optimal Seed Points\n")
  # TO DO!!!!
  # Lslide::segmentation(Tool = "SAGA", Segments.Grid = segmentation.first.grid, Segments.Poly = segmentation.first, Input.Grid = c(slope, slope.edge.vigra),
  #              Seed.Method = "Fast Representativeness", Fast.Representativeness.LevelOfGeneralisation = 1.25,  Saga.Segmentation.Method = "0", Saga.Segmentation.Sig.2 = "125",
  #              Generalisation.Flac = FALSE, Saga.Segmentation.Leafsize = 1024, NoData = TRUE, Mask = slope)

  # rsaga.get.modules("statistics_grid", env = env.rsaga)
  # rsaga.get.usage("statistics_grid", 0, env = env.rsaga)
  rsaga.geoprocessor(lib="statistics_grid", module = 0, env = env.rsaga, show.output.on.console = SOOC, param = list(
    INPUT = Input.Grid, RESULT = paste0(tempdir(),"/", "tmpResultFR"), RESULT_LOD = paste0(tempdir(),"/", "tmpLod"), SEEDS = Output.Seeds,
    LOD = Fast.Representativeness.LevelOfGeneralisation))


  # to do mask seeds
  if(quiet == FALSE) cat("... SAGA: Masking of Seed Points in NoData area\n")
  # grid tools - grid masking
  # rsaga.get.modules("grid_tools", env = env.rsaga)
  # rsaga.get.usage("grid_tools", 24, env = env.rsaga)
  rsaga.geoprocessor(lib = "grid_tools", module = 24, env = env.rsaga, show.output.on.console = SOOC, param = list(
    GRID = Output.Seeds, MASK = Input.Grid, MASKED = Output.Seeds.Masked))



  # seeded region growing
  if(quiet == FALSE) cat("... SAGA: Seeded Region Growing\n")
  # rsaga.get.modules("imagery_segmentation", env = env.rsaga)
  # rsaga.get.usage("imagery_segmentation", 3, env = env.rsaga)
 SRG <- rsaga.geoprocessor(lib = "imagery_segmentation", module = 3, env = env.rsaga, show.output.on.console = SOOC, param = list(
    SEEDS = Output.Seeds.Masked, FEATURES = Input.Grid, SIG_2 = Saga.Segmentation.Sig.2, SEGMENTS = Segments.Grid,
       LEAFSIZE = Saga.Segmentation.Leafsize))

 if(any(grepl('Error', SRG)))
 {
   stop("Error in Seeded Region Growing (SAGA). Check input and control number of seed points.")
 }

  # to do mask region growing
  if(quiet == FALSE) cat("... SAGA: Masking of Segments in NoData area\n")
  # grid tools - grid masking
  # rsaga.get.modules("grid_tools", env = env.rsaga)
  # rsaga.get.usage("grid_tools", 24, env = env.rsaga)
  rsaga.geoprocessor(lib = "grid_tools", module = 24, env = env.rsaga, show.output.on.console = SOOC, param = list(
    GRID = Segments.Grid, MASK = Input.Grid, MASKED = Segments.Grid))


  # vectorising grid classes
  if(quiet == FALSE) cat("... SAGA: Vectorising Grid Classes\n")
  # shapes_grid
  # rsaga.get.modules("shapes_grid", env = env.rsaga)
  # rsaga.get.usage("shapes_grid", 6, env = env.rsaga)
  rsaga.geoprocessor(lib = "shapes_grid", module = 6, env = env.rsaga, show.output.on.console = SOOC, param = list(
    GRID = Segments.Grid, POLYGONS =  Segments.Poly))



  if(quiet == FALSE) cat("... Loading Segments: Using sf::() and conversion to SpatialPolygonsDataFrame\n")
  outpoly <- sf::st_read(Segments.Poly, quiet = !SOOC)
  outpoly <- as(outpoly, "Spatial")

  if(is.na(sp::proj4string(outpoly)))
  {
    proj4string(outpoly) <- projection(input.filter)
  }

  if(quiet == FALSE) cat(paste0("------ Run of contrastFilterSegmentation: " , (proc.time() - process.time.start)["elapsed"][[1]]/60, " Minutes ------"))

  return(outpoly)
} # end function contrastFilterSegmentation





# EXAMPLE -----------------------------------
# library("RQGIS", "reticulate")
# open_app()
# pacman::p_load("EBImage", "raster", "data.table", "igraph", "RQGIS", "reticulate", "RSAGA")
#
#
#
# library("raster", "data.table")
# library("RQGIS")
# open_app()
#
# pacman::p_load("raster", "data.table")
# library("RQGIS")
# open_app()
#
# dem <- RQGIS::dem
#
# slope <- RQGIS::run_qgis(alg  = "grass7:r.slope.aspect", elevation = dem, slope = "slope.tif", load_output = TRUE)
#
#
#
# setwd("E:/Masterarbeit/Data/Output/Segmentation/EBImage")
# # save.image("EBImageSegmentation.RData")
# load("EBImageSegmentation.RData")
#
# # input <- raster("slopeT.tif")
# input <- raster("Input/OP14_slp_01m.tif")
# env.rsaga <- RSAGA::rsaga.env(path = "C:/Program Files (x86)/SAGA-GIS")
#
#
# # testing:
# input.filter = slope
# makeBrush.size = 10
# env.rsaga = env.rsaga
# SOOC = TRUE
# quiet = FALSE
# Fast.Representativeness.LevelOfGeneralisation = "10"
# input.segmentation = input.filter
# offset = 0.06
# makeBrush.shape = "disc"
# NA.val.in = 0
# clump.thresh = NULL
# clump.directions = 8
# CFR.buf = "1"
# output.path = tempdir()
# writeCFRaster = FALSE
# writeRaster.NAflag = -99999
# freeMemory = TRUE
# Saga.Segmentation.Sig.2 = "125"
# Saga.Segmentation.Leafsize = 1024
# Segments.Poly =  paste0(tempdir(), "/", "outSegPoly.shp")
# Segments.Grid = paste0(tempdir(), "/", "outSegGrid.sgrd")


