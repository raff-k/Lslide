#' Geomorphometry & Gradient Metrics Toolbox
#'
#' This toolbox contains different function for calculating geomorphometric indices or gradient models.
#' The toolbox is based on the Geomorphometry & Gradient Metrics Toolbox from Jeffrey Evans.
#' Actually, it supports the computation of 'Linear Aspect', 'Mean Slope', 'Site Exposure Index',
#' 'Landform' (concavity/convexity landform index (Bolstads variant)) and 'Roughness Index'.
#' 
#' @param elevation \linkS4class{RasterLayer} with elevation data 
#' @param slope \linkS4class{RasterLayer} with slope values. As default the slope will be computet using \emph{RQGIS::run_qgis(alg = "grass7:r.slope.aspect")}. Default: NULL
#' @param aspect \linkS4class{RasterLayer} with aspect values. As default the aspect will be computet using \emph{RQGIS::run_qgis(alg = "grass7:r.slope.aspect")}. Default: NULL
#' @param tool vector containg the name of tools. Actually, the following are supported: "Linear Aspect", "Landfrom", "Mean Slope",  "SEI", "RI". Default: c("Linear Aspect", "Landfrom", "Mean Slope",  "SEI", "RI")
#' @param params list of lists containing function parameter for functions in tool. There are basically three parameters defining the (1) "size" of moving window, (2) the filename, and (3) a logical value if the raster shall be written (writeRaster).  I.e. for "Linear Aspect" - list(size = 3, filename = "LinAspect.tif", writeRaster = TRUE). Default: see notice
#' @param zScale multiplicative factor to convert elevation units to horizontal units. Default: "1.0"
#' @param output.path path to output folder. Default: tempdir()
#' @param quiet no outputs in console. Default: TRUE
#'
#' @return
#' List of \linkS4class{RasterLayer} containing computed geomorphometric indices or gradient models
#'
#' @note
#' \itemize{
#'   \item GRASS GIS algorithm (grass7:r.slope.aspect) is similar to ArcGIS slope/aspect
#'   \item default for params: list(list(size = 3, filename = "LinAspect.tif", writeRaster = TRUE), list(size = 3, filename = "Landform.tif", writeRaster = TRUE), list(size = 3, filename = "MeanSlp.tif", writeRaster = TRUE), list(size = NA, filename = "SEI.tif", writeRaster = TRUE), list(size = 3, filename.RI = "RI.tif", filename.RIw = "RIw.tif", writeRaster = TRUE))
#' }
#'
#'
#' @keywords Geomorphometry & Gradient Metrics Toolbox, Jeffrey Evans, Linear Aspect, Mean Slope, Site Exposure Index, Landform, Roughness Index
#' 
#'
#' @export
#'
GeomorphometryGradientMetricsToolbox <- function(elevation, slope = NULL, aspect = NULL, tool = c("Linear Aspect", "Landfrom", "Mean Slope",  "SEI", "RI"),
                                                 params = params, zScale = "1.0", output.path = tempdir(), quiet = TRUE)
{

  # TO DO: PARALLELIZING, POSSIBILITY TO ADD ANOTHER ALGORYTHM IN RQGIS FOR EXAMPLE ONE OF SAGA GIS

  # get start time of process
  if(quiet == FALSE) process.time.start <- proc.time()


  # Initialize default parameters --------------------
  params <- list(list(size = 3, filename = "LinAspect.tif", writeRaster = TRUE), # Linear Aspect
                         list(size = 3, filename = "Landform.tif", writeRaster = TRUE), # Landfrom
                         list(size = 3, filename = "MeanSlp.tif", writeRaster = TRUE), # Mean Slope
                         list(size = NA, filename = "SEI.tif", writeRaster = TRUE), # SEI
                         list(size = 3, filename.RI = "RI.tif", filename.RIw = "RIw.tif", writeRaster = TRUE)) # RI


  # Check slope or aspect input --------------------
  if(is.null(slope) | is.null(aspect))
  {

    if(quiet == FALSE) cat("Calculate slope and/or aspect parameter ...\n")

    if(zScale == "1.0")
    {
      # RQGIS::get_args_man(alg = "grass7:r.slope.aspect")
      slope.aspect <- RQGIS::run_qgis(alg = "grass7:r.slope.aspect", elevation = elevation, show_output_paths = FALSE,
                                      slope = paste0(output.path , "/","slope.tif"), aspect = paste0(output.path , "/", "aspect.tif"), load_output = TRUE)
    } else {
      slope.aspect <- RQGIS::run_qgis(alg = "grass7:r.slope.aspect", elevation = elevation, zscale = zScale, show_output_paths = FALSE,
                                      slope = paste0(output.path , "/" ,"slope.tif"), aspect = paste0(output.path , "/", "aspect.tif"), load_output = TRUE)
    }

    # get slope from QGIS result
    if(is.null(slope)) slope <- slope.aspect[[1]]

    if(is.null(aspect))
    {
      # get aspect from QGIS result
      aspect <- slope.aspect[[2]]

      # change direction, so that North is 0 Degree
      aspect <- raster::calc(aspect, fun = function(x) { return((450 - x) %% 360)})
    }

  } # end check slope aspect input


  # Check tools to call functions  --------------------
  # supported at the moment: c("Linear Aspect", "Landfrom", "Mean Slope",  "SEI", "RI")
  loc.LinAsp <- match("Linear Aspect", tool)
  names(loc.LinAsp) <- "LinearAspect"

  loc.Landform <- match("Landfrom", tool)
  names(loc.Landform) <- "Landform"

  loc.MeanSlp <- match("Mean Slope", tool)
  names(loc.MeanSlp) <- "MeanSlope"

  loc.SEI <- match("SEI", tool)
  names(loc.SEI) <- "SEI"

  loc.RI <- match("RI", tool)
  names(loc.RI) <- "RI"

  iter <- c(loc.LinAsp, loc.Landform, loc.MeanSlp, loc.SEI, loc.RI)
  iter <- iter[!is.na(iter)] # remove NAs



  # Get params for functions --------------------
  params.LinAsp <- params[[loc.LinAsp]]
  params.Landform <- params[[loc.Landform]]
  params.MeanSlp <- params[[loc.MeanSlp]]
  params.SEI <- params[[loc.SEI]]
  params.RI <- params[[loc.RI]]



  # FUNCTIONS ---------------------------------
  # Linear Aspect ----------------------------
  LinearAspect <- function(aspect.Fun = aspect, size = params.LinAsp$size, filename = params.LinAsp$filename, writeRaster = params.LinAsp$writeRaster)
  {
    if(quiet == FALSE) process.time.start.LinearAspect <- proc.time()

    if(quiet == FALSE) cat("Running Linear Aspect ...\n")

    # remove negative values and convert to radiant
    tmp.aspect <- raster::calc(aspect.Fun, fun = function(x){ifelse(x < 0, NA, (450.0 - x)/57.296)}) # convert to radiant

    if(quiet == FALSE) cat("... calculate cos and sin from aspect\n")
    # calculate sinus and cosinus from aspect
    tmp.sin <- raster::calc(tmp.aspect, sin)
    tmp.cos <- raster::calc(tmp.aspect, cos)

    # calculate sums in focal statistics for sin and cos
    tmp.sin <- raster::focal(tmp.sin, w = matrix(1, size, size), fun = sum, na.rm = TRUE)
    tmp.cos <- raster::focal(tmp.cos, w = matrix(1, size, size), fun = sum, na.rm = TRUE)


    # start final calculations
    tmp.Mod <- raster::overlay(tmp.sin, tmp.cos, fun = function(x, y){return((((450-(atan2(x, y) * 57.296)) * 100) %% 36000)/100)})

    if(quiet == FALSE) cat("... final calculation\n")
    outRaster.linearAspect <- raster::overlay(tmp.sin, tmp.cos, tmp.Mod, fun = function(x, y, z) {ifelse((x == 0) & (y == 0), -1, z)})

    if(writeRaster == TRUE)
    {
      # write data
      if(quiet == FALSE) cat("... write raster\n")
      raster::writeRaster(x = outRaster.linearAspect, filename = paste0(output.path, '/', filename), overwrite = TRUE)
    }

    if(quiet == FALSE) cat(paste0("------ Run of Linear Aspect: " , (proc.time() - process.time.start.LinearAspect)["elapsed"][[1]]/60, " Minutes ------\n"))
    if(quiet == FALSE) cat("-------------------------\n")

    names(outRaster.linearAspect) <- "Linear Aspect"
    return(outRaster.linearAspect)
  } # end function Linear Aspect




  # Landfrom ------------------------------------
  Landform <- function(elevation.Fun = elevation, size = params.Landform$size, filename = params.Landform$filename, writeRaster = params.Landform$writeRaster)
  {
    if(quiet == FALSE) process.time.start.Landform <- proc.time()

    if(quiet == FALSE) cat("Running Landform: concavity/convexity landform index (Bolstads variant) ...\n")


    # calculate focal statistics: mean
    mean.Tmp <- raster::focal(elevation.Fun, w = matrix(1, size, size), fun = mean, na.rm = TRUE)


    # calculate landform
    if(quiet == FALSE) cat("... calculate landform\n")
    outRaster.landfrom <- raster::overlay(elevation.Fun, mean.Tmp, fun = function(x, y) {return(10000 * ((x - y)/1000/36.2))})

    if(writeRaster == TRUE)
    {
      # write data
      if(quiet == FALSE) cat("... write raster\n")
      raster::writeRaster(x = outRaster.landfrom, filename = paste0(output.path, '/', filename), overwrite = TRUE)
    }

    if(quiet == FALSE) cat(paste0("------ Run of Landform: " , (proc.time() - process.time.start.Landform)["elapsed"][[1]]/60, " Minutes ------\n"))
    if(quiet == FALSE) cat("-------------------------\n")

    names(outRaster.landfrom) <- "Landform"
    return(outRaster.landfrom)

  } # end function Landform





  # Mean Slope --------------------------
  MeanSlope <- function(slope.Fun = slope, size = params.MeanSlp$size, filename = params.MeanSlp$filename, writeRaster = params.MeanSlp$writeRaster)
  {

    if(quiet == FALSE) process.time.start.MeanSlope <- proc.time()

    if(quiet == FALSE) cat("Running Mean Slope ...\n")

    if(quiet == FALSE) cat("... calculate mean slope\n")
    outRaster.meanSlope <- raster::focal(slope.Fun, w = matrix(1, size, size), fun = mean, na.rm = TRUE)


    if(writeRaster == TRUE)
    {
      if(quiet == FALSE) cat("... write raster\n")
      raster::writeRaster(x = outRaster.meanSlope, filename = paste0(output.path, '/', filename), overwrite = TRUE)
    }

    if(quiet == FALSE) cat(paste0("------ Run of Mean Slope: " , (proc.time() - process.time.start.MeanSlope)["elapsed"][[1]]/60, " Minutes ------\n"))
    if(quiet == FALSE) cat("-------------------------\n")

    names(outRaster.meanSlope) <- "MeanSlope"
    return(outRaster.meanSlope)

  } # end function mean slope




  # Site Exposure Index ---------------------------
  SEI <- function(slope.Fun = slope, aspect.Fun = aspect, size = params.SEI$size, filename = params.SEI$filename, writeRaster = params.SEI$writeRaster)
  {

    if(quiet == FALSE) process.time.start.SEI <- proc.time()

    if(quiet == FALSE) cat("Running Site Exposure Index ...\n")

    tmp.cosResult <- raster::calc(aspect.Fun, fun = function(x) { return(cos((3.142 * (x - 180))/180))})

    if(quiet == FALSE) cat("... calculate Site Exposure Index\n")
    outraster.SEI <- raster::overlay(slope.Fun, tmp.cosResult, fun = function(x, y) {return(x * y)})

    if(writeRaster == TRUE)
    {
      # write raster
      if(quiet == FALSE) cat("... write raster\n")
      raster::writeRaster(x = outraster.SEI, filename = paste0(output.path, '/', filename), overwrite = TRUE)
    }

    if(quiet == FALSE) cat(paste0("------ Run of Site Exposure Index: " , (proc.time() - process.time.start.SEI)["elapsed"][[1]]/60, " Minutes ------\n"))
    if(quiet == FALSE) cat("-------------------------\n")

    names(outraster.SEI) <- "SEI"
    return(outraster.SEI)

  } # end of function SEI





  # Roughness Index ---------------------------
  RI <- function(elevation.Fun = elevation, size = params.RI$size, filename.RI = params.RI$filename.RI, filename.RIw = params.RI$filename.RIw, writeRaster = params.RI$writeRaster)
  {

    if(quiet == FALSE) process.time.start.RI <- proc.time()

    if(quiet == FALSE) cat("Running Roughness Index as standard deviation of residual topography (Cavalli et al. 2008) ...\n")

    # moving window mean of input
    if(quiet == FALSE) cat("... moving window mean\n")
    tmp.mean <- raster::focal(elevation.Fun, w = matrix(1, size, size), fun = mean, na.rm = TRUE)

    # difference beteen original and smoothed: residual topography
    if(quiet == FALSE) cat("... calculate residual topography\n")
    tmp.dif <- raster::overlay(elevation.Fun, tmp.mean, fun = function(x, y) {return(x - y)}, na.rm = TRUE)

    # roughness index
    # too slow...?
    # outraster.RI <- raster::focal(tmp.dif, w = matrix(1, size, size), fun = function(x) {return(sd(x, na.rm = TRUE))})
    # outraster.RI <- raster::focal(tmp.dif, w = matrix(1, size, size), fun = function(x) {return( sqrt(sum(abs(x-mean(x))^2)/size) )})


    if(writeRaster == TRUE)
    {
      # RQGIS::get_args_man(alg = "grass7:r.neighbors")
      if(quiet == FALSE) cat("... calculate Roughness Index\n")
      outraster.RI  <- RQGIS::run_qgis(alg = "grass7:r.neighbors", size = as.character(size), output = paste0(output.path, '/', filename.RI),
                                       input = tmp.dif , method = "stddev", load_output = TRUE, show_output_paths = FALSE)

      # weighted
      # RI.max <- maxValue(outraster.RI[!is.na(outraster.RI)])
      RI.max <- cellStats(outraster.RI,'max')

      if(quiet == FALSE) cat("... calculate normalized Roughness Index\n")
      outraster.RI.w <- raster::calc(outraster.RI, fun = function(x){return(1 - (x/RI.max))})

      # write weighted raster
      if(quiet == FALSE) cat("... write raster\n")
      raster::writeRaster(x = outraster.RI.w, filename = paste0(output.path, '/', filename.RIw), overwrite = TRUE)

    } else {

      if(quiet == FALSE) cat("... calculate Roughness Index\n")
      outraster.RI  <- RQGIS::run_qgis(alg = "grass7:r.neighbors", size = as.character(size),
                                       input = tmp.dif , method = "stddev", load_output = TRUE)

      RI.max <- cellStats(outraster.RI,'max')
      if(quiet == FALSE) cat("... calculate normalized Roughness Index\n")
      outraster.RI.w <- raster::calc(outraster.RI, fun = function(x){return(1 - (x/RI.max))})

    } # enf if write == TRUE

    if(quiet == FALSE) cat(paste0("------ Run of Roughness Index: " , (proc.time() - process.time.start.RI)["elapsed"][[1]]/60, " Minutes ------\n"))
    if(quiet == FALSE) cat("-------------------------\n")

    names(outraster.RI) <- "RI"
    names(outraster.RI.w) <- "RIw"
    return(list(outraster.RI, outraster.RI.w))

  } # end of function RI



  # CALL FUNCTIONS ------------------------------
  # get iteration for function call
  # TO DO: PARALLELIZING!
  list.i <- list()
  for(i in 1:length(iter))
  {
    out.i <- tryCatch({do.call(names(iter)[i], args = list())
    }, error = function(e)
    {
      warning(paste0("##### Something went wrong in: ", names(iter)[i],  ". Check params and input files #####"))
      if(quiet == FALSE) cat("-------------------------\n")
      return(NA)
    })

    list.i[[i]] <- out.i
  }
  names(list.i) <- names(iter) # name list items


  if(quiet == FALSE) cat(paste0("------ Run of Geomorphometry Gradient Metrics Toolbox: " , (proc.time() - process.time.start)["elapsed"][[1]]/60, " Minutes ------\n"))

  return(list.i)

} # end function GeomorphometryGradientMetricsToolbox






















# # EXAMPLE --------------------------
# pacman::p_load(RQGIS, raster)
#
# data(dem)
# dem <- dem
#
# setwd("D:/Users/rAVer/Desktop/LSlide New")
#
# dtm <- raster::raster("Input/dtm.tif")
# # dtm <- raster::raster("Input/DTM_Po_30m.tif")
#
# elevation <- dem
#
# set_env()
# open_app()
# geomorph <- GeomorphometryGradientMetricsToolbox(elevation = dem, quiet = FALSE)
# # size <- "3"
# size <- 3
# out.path.linearAspect <- "Output/LinearAspect.tif"
# out.path.landform <- "Output/Landform.tif"
# out.path.meanSlope <- "Output/Mean Slope.tif"
# out.path.SEI <- "Output/SEI.tif"
# out.path.RI <- "Output/RI.tif"
# out.path.RIw <- "Output/RIw.tif"


# # BACKUP --------------------------
# # Linear Aspect ----------------------------
# # (Elapsed Time: 3 minutes 30 seconds) ArcGIS
# if(modul == "Linear Aspect")
# {
#   cat("Running Linear Aspect ......")
#
#   # remove negative values and convert to radiant
#   tmp.aspect <- raster::calc(aspect, fun = function(x){ifelse(x < 0, NA, (450.0 - x)/57.296)}) # convert to radiant
#
#
#   # calculate sinus and cosinus from aspect
#   tmp.sin <- raster::calc(tmp.aspect, sin)
#   tmp.cos <- raster::calc(tmp.aspect, cos)
#
#   # RQGIS::get_args_man(alg = "grass7:r.neighbors")
#
#   # tmp5 <- RQGIS::run_qgis(alg = "grass7:r.neighbors", input = tmp.sin, output = "tempSin.tif", size = size,
#   #                         method = "sum", load_output = TRUE)
#
#   # tmp6 <- RQGIS::run_qgis(alg = "grass7:r.neighbors", input = tmp.cos, output = "tempCos.tif", size = size,
#   #                         method = "sum", load_output = TRUE)
#
#
#   # calculate sums in focal statistics
#   tmp.sin <- raster::focal(tmp.sin, w = matrix(1, size, size), fun = sum, na.rm = TRUE)
#   # tmp.sin <- raster::focal(tmp.sin, fun = sum, na.rm = TRUE, w = raster::focalWeight(x = tmp.sin, type = 'rectangle', d = size))
#   # raster::writeRaster(x = tmp.sin , filename = "Output/LinearAspect_FocalSumSin.tif", overwrite = TRUE)
#
#   tmp.cos <- raster::focal(tmp.cos, w = matrix(1, size, size), fun = sum, na.rm = TRUE)
#   # tmp.cos <- raster::focal(tmp.cos, w = raster::focalWeight(x = tmp.cos , type = 'rectangle', d = size))
#   # raster::writeRaster(x = tmp.cos , filename = "Output/LinearAspect_FocalSumCos.tif", overwrite = TRUE)
#
#
#   # start final calculations
#   tmp.Mod <- raster::overlay(tmp.sin, tmp.cos, fun = function(x, y){return((((450-(atan2(x, y) * 57.296)) * 100) %% 36000)/100)})
#   # raster::writeRaster(x = tmpMod , filename = "Output/LinearAspect_Mod.tif", overwrite = TRUE)
#
#   outRaster.linearAspect <- raster::overlay(tmp.sin, tmp.cos, tmp.Mod, fun = function(x, y, z) {ifelse((x == 0) & (y == 0), -1, z)})
#
#   # write data
#   raster::writeRaster(x = outRaster.linearAspect, filename = out.path.linearAspect, overwrite = TRUE)
#
#   # remove temporary data
#   rm(tmp.sin, tmp.cos, tmpMod, tmp.aspect)
#   gc()
#
#   # return outrast
# }
#
#
#
#
#
#
# # Landfrom ------------------------------------
# # ArcGIS (Elapsed Time: 50,89 seconds)
# if(modul == "Landfrom")
# {
#   cat("Landform: concavity/convexity landform index (Bolstads variant)")
#
#
#   # fw <- focalWeight(elevation, 1, type = 'rectangle')
#   # fw[fw>0] <- 1
#   # mean.Tmp <- raster::focal(elevation, w = fw, fun = mean, na.rm = TRUE)
#
#   # calculate focal statistics: mean
#   mean.Tmp <- raster::focal(elevation, w = matrix(1, size, size), fun = mean, na.rm = TRUE)
#
#   # calculate landform
#   outRaster.landfrom <- raster::overlay(elevation, mean.Tmp, fun = function(x, y) {return(10000 * ((x - y)/1000/36.2))})
#
#   # write data
#   raster::writeRaster(x = outRaster.landfrom, filename = out.path.landform, overwrite = TRUE)
#
#   # remove temporary data
#   rm(mean.Tmp)
#   gc()
# }
#
#
#
#
# # Mean Slope --------------------------
# if(modul == "Mean Slope")
# {
#   cat("Running Mean Slope ......")
#
#   outRaster.meanSlope <- raster::focal(slope, w = matrix(1, size, size), fun = mean, na.rm = TRUE)
#   raster::writeRaster(x = outRaster.meanSlope, filename = out.path.meanSlope, overwrite = TRUE)
# }
#
#
#
#
#
#
# # Site Exposure Index ---------------------------
# # spRefType = "Geographic" default: "", inUnits = "Meter" or "Feet" or ""
# # ArcGIS Elapsed Time: 1 minutes 31 seconds
# if(modul == "SEI")
# {
#   # Meter | Foot | Geographic
#
#   cat("Site Exposure Index")
#
#
#   # zFactor <- getZFactor(r = dtm, spRefType = "", inUnits = "Meter")
#   #
#   # # RQGIS::get_args_man(alg = "grass7:r.slope.aspect")
#   # slope.aspectZ <- RQGIS::run_qgis(alg = "grass7:r.slope.aspect", elevation = dtm, zscale = zFactor,
#   #                                 slope = "slope.tif", aspect = "aspect.tif", load_output = TRUE)
#
#   # get slope from result
#   # slopeZ <- slope.aspectZ[[1]]
#   # aspectZ <- slope.aspectZ[[2]]
#
#   # caculate cosinus raster
#   tmp.cosResult <- raster::calc(aspect, fun = function(x) { return(cos((3.142 * (x - 180))/180))})
#
#   outraster.SEI <- raster::overlay(slope, tmp.cosResult, fun = function(x, y) {return(x * y)})
#   raster::writeRaster(x = outraster.SEI, filename = out.path.SEI, overwrite = TRUE)
#
#   rm(tmp.cosResult)
#   gc()
#
# }
#
#
#
#
#
# # Roughness Index ---------------------------
# if(modul == "Roughness Index")
# {
#
#   cat("Roughness Index as standard deviation of residual topography (Cavalli et al. 2008)")
#
#   # moving window mean of input
#   tmp.mean <- raster::focal(elevation, w = matrix(1, size, size), fun = mean, na.rm = TRUE)
#
#   # difference beteen original and smoothed: residual topography
#   tmp.dif <- raster::overlay(elevation, tmp.mean, fun = function(x, y) {return(x - y)}, na.rm = TRUE)
#
#   # roughness index
#   # too slow...?
#   # outraster.RI <- raster::focal(tmp.dif, w = matrix(1, size, size), fun = function(x) {return(sd(x, na.rm = TRUE))})
#   # outraster.RI <- raster::focal(tmp.dif, w = matrix(1, size, size), fun = function(x) {return( sqrt(sum(abs(x-mean(x))^2)/size) )})
#
#   # RQGIS::get_args_man(alg = "grass7:r.neighbors")
#   outraster.RI  <- RQGIS::run_qgis(alg = "grass7:r.neighbors", size = as.character(size), output = out.path.RI,
#                                   input = tmp.dif , method = "stddev", load_output = TRUE)
#
#   # weighted
#   # RI.max <- maxValue(outraster.RI[!is.na(outraster.RI)])
#   RI.max <- cellStats(outraster.RI,'max')
#   outraster.RI.w <- raster::calc(outraster.RI, fun = function(x){return(1 - (x/RI.max))})
#
#   raster::writeRaster(x = outraster.RI.w, filename = out.path.RIw, overwrite = TRUE)
#
# }
#
#
#
#
#
#
#
#
#
#
#
#
#
# # geomorph routines ---------------------------------
#
#
# getZFactor <- function(r, spRefType = "", inUnits = "Meter")
# {
#
#   r.extent <- extent(r)
#
#   if(spRefType == "Geographic")
#   {
#     medianLat <- getMidLat(r, r.extent, spRefType = spRefType)
#     rLat <- medianLat * pi /180 # Convert angle x from degrees to radians
#
#     meters2Degree <- 111412.84 * cos(rLat) - 93.5* cos(3*rLat)
#
#     if(inUnits == "Meter")
#     {
#       # using 1 deg long is 111.412 kms at equator
#       numerator <- meters2Degree
#     } else {
#       #Using 1 deg long is 69.172 miles at equator
#       numerator <- meters2Degree*3.28084
#
#     } # end forif(inZUnits == "Meter")
#
#     zFactor <- 1/numerator
#
#   } else { # spRefType is not geographic
#
#     # inUnits <- inSpRef.linearUnitName
#     #
#     # if(inUnits == "Meter")
#     # {
#     #   #Since they don't equal this means xyUnits == Meter and zUnits == Feet
#     #   zFactor <- 0.3048
#     #
#     # } else if(inUnits == "Feet"){
#     #
#     #   #Means xyUnits == Feet and  zUnits == Meters (not common)
#     #   zFactor <- 3.28084
#     #
#     # } else {
#
#       zFactor <- 1
#
#     # } # end for if(inZUnits == "Meter")
#
#   } # end for if(spRefType == "Geographic")
#
#   return(zFactor)
# }
#
#
#
# getMidLat <- function(r, rExt, spRefType = "")
# {
#   if(spRefType == "Geographic")
#   {
#     rYMax = rExt[4] # YMax
#     rYMin = rExt[3] # YMin
#     medianLat = abs((rYMax - rYMin)/2 + rYMin)
#
#   } else {
#
#     # use as WGS84: EPSG(4326)
#     # proj4string(prjExt) <- sp::proj4string(CRS("+init=epsg:4326"))
#     prjExt <- as(rExt, "SpatialPolygons")
#     sp::proj4string(prjExt) <- raster::projection(dtm)
#     prjExt <- sp::spTransform(prjExt, sp::CRS("+init=epsg:4326"))
#
#     prjExt <- raster::extent(prjExt)
#
#     rYMax = prjExt[4] # YMax
#     rYMin = prjExt[3] # YMin
#     medianLat = abs((rYMax- rYMin)/2 + rYMin)
#
#   }
#
#   return(medianLat)
#
# }

