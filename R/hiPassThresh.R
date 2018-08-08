#' Using high-pass filter and thresholds to filter image
#'
#' This function performs a high-pass filter followed by a thresholding procedure in order to
#' select specific cluster of cells.
#'
#'
#' @param x raster object or input path of file
#' @param scale.factor
#' @param threshold
#' @param path.output = NULL
#' @param do.sieve = TRUE
#' @param sieve.mode mode of cell-neighborhood,  Neumann [0]: the four horizontally and vertically neighboured cells; Moore [1]: all eight adjacent cells. Default: "0"
#' @param sieve.thresh minimum number of cells in a group of adjacent cells (minimum 2). Default: 4
#' @param do.shrink.expand regions with valid data in the input grid can be shrinked or expanded by a certain amount (radius). Default: TRUE
#' @param expand maximum expansion distance (radius). Default: 4 [cells]
#' @param do.buf = FALSE
#' @param buf.size buffer size in cell's value. Default: 1
#' @param do.use.temp.HPF already existing high-pass filter files are used. Default: FALSE
#' @param do.use.temp.Thresh already existing thresholding files are used. Default: FALSE
#' @param path.save output path for saving objects. Default: tempdir()
#' @param NoData no data value. Default: -99999
#' @param show.output.on.console show output on console. Default: FALSE
#' @param env.rsaga environment of SAGA GIS. As default, the environment will be automatically estimated. Default: NULL
#' @param quiet no outputs in console. Default: TRUE
#' @return
#' raster object
#'
#'
#' @note
#' \itemize{
#'   \item SAGA GIS must be installed.
#'
#' }
#'
#'
#' @keywords high-pass filter, SAGA GIS
#'
#'
#' @export
#'
hiPassThresh <- function(x, scale.factor, threshold, path.output = NULL, do.sieve = TRUE, sieve.mode = "0", sieve.thresh = 4, do.shrink.expand = TRUE, expand = 4,
                                 do.buf = FALSE, buf.size = 1, do.use.temp.HPF = FALSE, do.use.temp.Thresh = FALSE, path.save = tempdir(), NoData = -99999, env.rsaga = NULL, show.output.on.console = FALSE, quiet = TRUE, ...)
{

  # get start time of process
  process.time.start <- proc.time()

  if(!is.null(path.output) && (tools::file_ext(path.output) != ".tif" || tools::file_ext(path.output) != ".tiff"))
  {
    stop('File path of "path.output" is not correct. Must be of format .tif or .tiff!')
  }

  if(is.null(env.rsaga))
  {
    tryCatch({
    env.rsaga <- RSAGA::rsaga.env()
    }, error = function(err) {
      stop("Could not initialize SAGA GIS!")
    })
  }

  if(class(x) == "RasterLayer")
  {
    path.x <- file.path(tempdir(), "input.sgrd")
    raster::writeRaster(x = x, filename = paste0(tools::file_path_sans_ext(path.x), ".sdat"), overwrite = TRUE, NAflag = NoData)
  } else {
    path.x <- x
  }

  if(sieve.thresh < 2)
  {
    stop('Parameter "sieve.thresh" must be greater or equal to 2!')
  }


  ## ... performing high-pass filter -------------------
  if(quiet == FALSE) cat("... high-pass filter with scale: ", scale.factor, "\n")

  scale.txt <- gsub(pattern = "\\.", replacement = "", x = as.character(scale.factor))
  path.hipass <- file.path(path.save, paste0("hipass_", scale.txt, ".sgrd"))

  if(do.use.temp.HPF && file.exists(path.hipass))
  {
    cat("Tempfile is used! \n")

  } else {
    # RSAGA::rsaga.get.usage(lib = "grid_filter", module = 11, env = env.rsaga)
    RSAGA::rsaga.geoprocessor(lib = "grid_filter", module = 11, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
      GRID = path.x, HIPASS = path.hipass, SCALE = scale.factor))
  }


  ## ... check raster for continuation ----------------------
  # thresh.txt <- gsub(pattern = "\\.", replacement = "", x = as.character(threshold))
  thresh.txt <- gsub(pattern = "\\.", replacement = "_", x = as.character(threshold))
  path.hipassThresh <- file.path(path.save, paste0("hipass_", scale.txt, "_", thresh.txt.save, ".sgrd"))

  r.check <- raster::raster(paste0(tools::file_path_sans_ext(path.hipass), ".sdat"))
  r.check.max <- suppressWarnings(max(raster::values(r.check), na.rm = TRUE))

  if(is.na(r.check.max) || r.check.max == -Inf || r.check.max < threshold)
  {
    raster::values(r.check) <- NA
    hipass <- r.check
    names(hipass) <- basename(tools::file_path_sans_ext(path.hipassThresh))
    raster::writeRaster(x = hipass, filename = paste0(tools::file_path_sans_ext(path.hipassThresh), ".sdat"), overwrite = TRUE, NAflag = NoData)

    # re-save raster
    RSAGAUsage <- RSAGA::rsaga.get.usage(lib="io_gdal", module = 1, env = env.rsaga,  show = FALSE)
    formatSAGA <- gsub("\\D", "", grep('SAGA GIS Binary', RSAGAUsage, value = TRUE))
    # RSAGA::rsaga.get.usage(lib = "io_gdal", module = 1, env = env.rsaga)
    RSAGA::rsaga.geoprocessor(lib = "io_gdal", module = 1, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
      GRIDS = path.hipassThresh, FILE = paste0(tools::file_path_sans_ext(path.hipassThresh), ".sdat"), FORMAT = formatSAGA, SET_NODATA = "1", NODATA = NoData))

  } else {
    ## ... thresholding ----------------------
    # RSAGA::rsaga.get.usage(lib = "grid_calculus", module = 1, env = env.rsaga)
    if(quiet == FALSE) cat("... thresholding high-pass filter with threshold: ", threshold, "\n")
    formula.grdFil <- paste0("gt(a,", threshold, ")")

    if(!file.exists(path.hipassThresh) || !do.use.temp.Thresh)
    {
      RSAGA::rsaga.geoprocessor(lib = "grid_calculus", module = 1, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
        GRIDS = path.hipass, FORMULA = formula.grdFil, FNAME = "1", RESULT = path.hipassThresh))
    }



    ## ... sieving and clump ----------------------
    if(do.sieve && sieve.thresh > 2) # 2 is minumum
    {
      if(quiet == FALSE) cat("... removal of clumbs based on threshold: ", sieve.thresh, "\n")
      if(quiet == FALSE) cat("... ... sieving\n")

      # RSAGA::rsaga.get.usage(lib = "grid_filter", module = 15, env = env.rsaga)
      # MODE: [0] Neumann, [1] Moore
      if(!file.exists(path.hipassThresh) || !do.use.temp.Thresh)
      {
        RSAGA::rsaga.geoprocessor(lib = "grid_filter", module = 15, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
          INPUT = path.hipassThresh, OUTPUT = path.hipassThresh, MODE = sieve.mode, THRESHOLD = sieve.thresh, ALL = "1"))
      }
    }


    ## ... shrink and expand ----------------------
    if(do.shrink.expand)
    {
      if(quiet == FALSE) cat("... shrink and expand\n")
      # RSAGA::rsaga.get.usage(lib = "grid_tools", module = 28, env = env.rsaga)
      # OPERATION: [3] expand and shrink
      # CIRCLE: [1] Circle
      # EXPAND: [3] majority
      if(!file.exists(path.hipassThresh) || !do.use.temp.Thresh)
      {
        RSAGA::rsaga.geoprocessor(lib = "grid_tools", module = 28, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
          INPUT = path.hipassThresh, RESULT = path.hipassThresh, OPERATION = "3", CIRCLE = "1", RADIUS = expand, EXPAND = "3"))
      }
    }


    ## ... set 0 to NO DATA ----------------------
    # RSAGA::rsaga.get.usage(lib = "grid_tools", module = 15, env = env.rsaga)
    if(!file.exists(path.hipassThresh) || !do.use.temp.Thresh)
    {
      RSAGA::rsaga.geoprocessor(lib = "grid_tools", module = 15, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
        INPUT = path.hipassThresh, RESULT = path.hipassThresh, METHOD = "0", OLD = 0, NEW = NoData, SOPERATOR = "0",
        RESULT_NODATA_CHOICE = "1", RESULT_NODATA_VALUE = NoData))
    }


    if(do.buf && buf.size > 0)
    {
      # Start buffering contrast filter results -----------------
      if(quiet == FALSE) cat("... buffer contrast filter results\n")

      ## ... buffer data
      # RSAGA::rsaga.get.usage(lib = "grid_tools", module = 8, env = env.rsaga)
      # [1] cell's value
      if(!file.exists(path.hipassThresh) || !do.use.temp.Thresh)
      {
        RSAGA::rsaga.geoprocessor(lib="grid_tools", module = 8, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
          FEATURES = path.hipassThresh, BUFFER = path.hipassThresh, DISTANCE = buf.size, TYPE = "1"))
      }
    }


    ## ... load data
    hipass <- raster::raster(paste0(tools::file_path_sans_ext(path.hipassThresh), ".sdat"))

  } # end of if-else check


  if(!is.null(path.output))
  {
    raster::writeRaster(x = hipass, filename = path.output, overwrite = TRUE, NAflag = NoData)
  }

  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) cat("------ Run of hiPassThresh: " , round(x = process.time.run["elapsed"][[1]]/60, digits = 3), " Minutes ------\n")


  return(hipass)

} # end of function highPassThresholding
