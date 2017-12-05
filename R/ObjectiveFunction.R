#' Objective Function
#'
#' The Objective Function by ESPINDOLA et al. (2006) finds optimal
#' scale parameters by the comparision of the sum of the normalized
#' local variance and the normalized Moran's I of different
#' segmented images.
#'
#' @param Tool open-source software to compute segmentation analysis. GRASS or SAGA
#' @param Segments.Poly ...
#' @param Seed.Method ""
#' @param env ...
#' @param Grass.Segmentation.Minsize = NA,
#' @param Grass.Segmentation.Threshold = NA
#' @param Scale.Input.Grid input grid for computing segmentation scale parameters
#' @param Scale.Input.Grid.Cell.Size cell size of input grid. Default: "1"
#' @param Scale.Statistic.Min.Size min size of area/polygon which is included in statistics (usefull for SAGA GIS segmentations). Default: "0"
#' @param Objective.Function.save save estimations of function. Default: FALSE
#' @param Objective.Function.save.path path where estimations are stored. Default: input path of \emph{segment.poly}
#' @param Objective.Function.save.name name for file. Default: ""
#' @param Objective.Function.MoransI output text file with Morans'I. Default: temp
#' @param Scales containing scale parameters for loop-segmentation
#' @param Count amount of cells (Grid Statistics). Default:"1"
#' @param Min minimum value (Grid Statistics). Default:"0"
#' @param Max maximum value (Grid Statistics). Default:"0"
#' @param Range range of values (Grid Statistics). Default:"0"
#' @param Sum sum of values (Grid Statistics). Default:"0"
#' @param Mean mean of values (Grid Statistics). Default:"1"
#' @param Var variance of values (Grid Statistics). Default:"1"
#' @param Stddev standard deviation (Grid Statistics). Default:"0"
#' @param Quantile qunatile (Grid Statistics). Default:"0"
#' @param Objective.Function.Mean.Segmentation.Grid - output grid with mean values for segments. Default: temp
#' @param Objective.Function.Count.Grid output grid with amount of cells in segments. Default: temp
#' @param Grass.Objective.Function.Method determining on what parameter the objective function (~loop) is performed. Default: "Threshold"
#' @param other... see \code{\link{segmentation}}
#'
#'
#' @note
#' \itemize{
#'   \item ESPINDOLA, G.M., CAMARA, G., REIS, I.A., BINS, L.S. & A.M. MONTEIRO (2006):
#'         Parameter selection for region-growing image segmentation algorithms using spatial autocorrelation.
#'          - International Journal of Remote Sensing 27, 14, 3035-3040.
#' }
#'
#' @keywords objective function, USPO
#'
#'
#' @export
Objective.Function <- function(Tool, Scale.Input.Grid, Scale.Input.Grid.Cell.Size = "1", Scale.Statistic.Min.Size = "0", Objective.Function.save = FALSE, Objective.Function.save.path = NULL, Count = "1", Min = "0", Max = "0", Range = "0", Sum = "0", Mean = "1", Var = "1", Stddev = "0", Objective.Function.save.name = "",
                               Objective.Function.Mean.Segmentation.Grid = paste0(tmp.path, "Objective.Function.Mean.Segmentation.Grid.sgrd"), Objective.Function.Count.Grid = paste0(tmp.path, "Objective.Function.Count.Grid.sgrd"), Objective.Function.MoransI = paste0(tmp.path, "Objective.Function.MoransI.txt"),
                               Quantile = 0, Scales, Grass.Objective.Function.Method = "Threshold", Segments.Poly, Grass.Segmentation.Minsize = NA, Grass.Segmentation.Threshold = NA, Seed.Method = "", env = RSAGA::rsaga.env(), show.output.on.console = FALSE, ...)
{

  # browser()

  # get start time of process
  process.time.start.scale <- proc.time()

  # Data preparation -----------------------------------------------------------------------
  # create data table for Objective.Function
  if(Tool == "SAGA")
  {
    df.Objective.Function  <- data.frame("Scale Parameter" = double(), "Intrasegment Variance" = double(), "Normalized Intrasegment Variance" = double(),
                                         "Morans I" = double(), "Normalized Morans I" = double(), "Objective Function" = double())
  }

  if(Tool == "GRASS")
  {
    df.Objective.Function  <- data.frame("Threshold" = double(), "Minsize" = double(), "Intrasegment Variance" = double(), "Normalized Intrasegment Variance" = double(),
                                         "Morans I" = double(), "Normalized Morans I" = double(), "Objective Function" = double())
  }


  if(Tool == "GRASS Superpixels SLIC")
  {
    df.Objective.Function  <- data.frame("Superpixels" = double(), "Compactness" = double(), "Intrasegment Variance" = double(), "Normalized Intrasegment Variance" = double(),
                                         "Morans I" = double(), "Normalized Morans I" = double(), "Objective Function" = double())
  }


  Objective.Function.save.path.default <- dirname(Segments.Poly)
  print("Calculate Objective Function")
  # Start loop through scales -----------------------------------------------------------------------
  for(i in Scales)
  {

    multiplier <- 1

    if(Seed.Method == "Fast Representativeness")
    {
      multiplier <- 100
    }

    if(Grass.Objective.Function.Method == "Compactness")
    {
      multiplier <- 100
    }


    print(paste0("Level of Generalisation|Threshold|Minsize|... : ", i))
    segments.poly <- paste0(tools::file_path_sans_ext(Segments.Poly) , (i*multiplier),".", tools::file_ext(Segments.Poly))
    # segments.scale.statistic <- paste0(tools::file_path_sans_ext(segments.poly) , "_scaleStat.", tools::file_ext(segments.poly))
    segments.scale.statistic <- segments.poly # update after changed dbf header

    # Segmentation -----------------------------------------------------------------------
    # perform segmentation
    if(Tool == "SAGA")
    {
      segmentation(Tool = Tool, Segments.Poly = segments.poly, Seed.Method = Seed.Method, env = env,
                   Fast.Representativeness.LevelOfGeneralisation = i, Seed.Generation.Scale = as.character(i),
                   show.output.on.console = show.output.on.console, ...)
    } else if(Tool == "GRASS")
    {
      if(Grass.Objective.Function.Method == "Threshold")
      {
        # browser()
        segmentation(Tool = Tool, Segments.Poly = segments.poly, Seed.Method = Seed.Method, env = env,
                     Grass.Segmentation.Threshold = as.character(i), Grass.Segmentation.Minsize = Grass.Segmentation.Minsize,
                     show.output.on.console = show.output.on.console, ...)
      }

      if(Grass.Objective.Function.Method == "Minsize")
      {
        segmentation(Tool = Tool, Segments.Poly = segments.poly, Seed.Method = Seed.Method, env = env,
                     Grass.Segmentation.Minsize = i, show.output.on.console = show.output.on.console,
                     Grass.Segmentation.Threshold = Grass.Segmentation.Threshold,...)
      }

      if(Grass.Objective.Function.Method == "Seeds")
      {
        segmentation(Tool = Tool, Segments.Poly = segments.poly, Seed.Method = Seed.Method, env = env,
                     Fast.Representativeness.LevelOfGeneralisation = i, Seed.Generation.Scale = as.character(i),
                     show.output.on.console = show.output.on.console,
                     Grass.Segmentation.Minsize = Grass.Segmentation.Minsize, Grass.Segmentation.Threshold = Grass.Segmentation.Threshold, ...)

      }

    } else if(Tool == "GRASS Superpixels SLIC")
    {
      if(Grass.Objective.Function.Method == "Compactness")
      {
        segmentation(Tool = Tool, Segments.Poly = segments.poly, Seed.Method = Seed.Method, env = env,
                     Grass.SLIC.Compactness = i, show.output.on.console = show.output.on.console, ...)
      }

      if(Grass.Objective.Function.Method == "Superpixels")
      {
        segmentation(Tool = Tool, Segments.Poly = segments.poly, Seed.Method = Seed.Method, env = env,
                     Grass.SLIC.Superpixels = i, show.output.on.console = show.output.on.console, ...)
      }

    }

    # browser()

    # Grid Statistics -----------------------------------------------------------------------
    # calculate grid statistics
    # rsaga.get.modules("shapes_grid", env = env)
    # rsaga.get.usage("shapes_grid", 2, env = env)
    # method: [0] standard | [1] shape wise, supports overlapping polygons

    # if(tools::file_ext(Scale.Input.Grid) == "sgrd")
    # {
    #   Scale.Input.Grid.tmp <- Scale.Input.Grid
    # } else {
    #
    #   Scale.Input.Grid.tmp  <- paste0(tempdir(), "/", "tmpScaleGrid.sgrd")
    #   raster::writeRaster(x = raster::raster(Scale.Input.Grid), filename = Scale.Input.Grid.tmp, overwrite = TRUE)
    #
    #   # RSAGA::rsaga.geoprocessor(lib = "io_gdal", module = 1, env = env, show.output.on.console = show.output.on.console, param = list(
    #   #   GRIDS = Segments.Grid.tmp.path, FILE =  paste0(tools::file_path_sans_ext(Segments.Grid.tmp.path), ".sdat"), FORMAT =  "42")
    # }

   # browser()

    # print("Calculate Grid Statistics")
    RSAGA::rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = show.output.on.console, param = list(
      GRIDS = Scale.Input.Grid, POLYGONS = segments.poly, METHOD = "0", NAMING = 1, RESULT = segments.scale.statistic,
      COUNT = Count, MIN = Min, MAX = Max, RANGE = Range, SUM = Sum, MEAN = Mean, VAR = Var, STDDEV = Stddev, QUANTILE = Quantile))

    # read statistic dbf
    print("Extraction for Objective Function Parameter")
    # get start time of process
    process.time.start.extraction <- proc.time()

    stat <- suppressMessages(read.dbf(paste0(tools::file_path_sans_ext(segments.scale.statistic), ".dbf")))    # last: variance | last - 1: mean | last - 2: cell count
    colnames(stat$dbf) <- c("cat", "value", "Cell_Count", "Mean", "Variance")
    write.dbf(stat, paste0(tools::file_path_sans_ext(segments.scale.statistic), ".dbf")) # write dbf with better header


    # Calculation of intrasegment Variance and Morans I -----------------------------------------------------------------------
    ### calculation of intrasegment variance
    print("... Calculation of Intrasegment Variance")
    stat$dbf <- stat$dbf[stat$dbf$Cell_Count > (as.numeric(Scale.Statistic.Min.Size) * as.numeric(Scale.Input.Grid.Cell.Size)),]
    stat.intrasegment.variance <- weighted.mean(stat$dbf[[ncol(stat$dbf)]], stat$dbf[[ncol(stat$dbf)-2]] * as.numeric(Scale.Input.Grid.Cell.Size), na.rm = TRUE) # weighted mean of variance by area


    ### calculation of Moran's I - speed up by sf-package
    # checked agains ArcGIS Spatial Autocorrelation Morans I (Contiguity_edges_corners) - similar results
    # calculate Morans'I for shapefiles
    print("... Calculation of Morans'I")
    # moransI.shape <- readOGR(dsn = dirname(segments.scale.statistic), layer = tools::file_path_sans_ext(basename(segments.scale.statistic)), verbose = FALSE)
    # moransI.queen.nb <-poly2nb(moransI.shape) # create contiguity queen neughbours by Waller & Gotway
    moransI.shape_sf <- sf::st_read(dsn = paste0(getwd(), "/", dirname(segments.scale.statistic)),
                                    layer = tools::file_path_sans_ext(basename(segments.scale.statistic)),
                                    stringsAsFactors = FALSE, quiet = TRUE) # simple feature object, sf package

    moransI.shape <- as(moransI.shape_sf, "Spatial") # convert simple feature to spatialobjectdataframe
    poly2_nb_speed_up <- rgeos::gUnarySTRtreeQuery(moransI.shape) # speed up function for poly2nb
    moransI.queen.nb <-spdep::poly2nb(moransI.shape, foundInBox = poly2_nb_speed_up) # create contiguity queen neughbours by Waller & Gotway

    # get out NA's, set them to 0
    if(any(is.na(moransI.shape$Mean)))
    {
      moransI.shape$Mean[is.na(moransI.shape$Mean)] <- 0
    }


    ## relating to Espindola et al. 2006 all neighbours get weights of 1 - binary weights
    moransI.weights.binary <- spdep::nb2listw(moransI.queen.nb, style = "B", zero.policy = TRUE) # those with many neighbours are up-weighted compared to those with few
    moransI <- spdep::moran.test(x = moransI.shape$Mean, listw = moransI.weights.binary, zero.policy = TRUE, randomisation = TRUE, na.action = na.pass) # global Moran's I, if NA is there -> 0 weight

    process.time.run.extraction <- proc.time() - process.time.start.extraction
    print(paste0("------ Run of Extraction for Objective Function Parametern: " , process.time.run.extraction["elapsed"][[1]]/60, " Minutes ------"))



    # write to data frame
    if(Tool == "SAGA")
    {
      df.Objective.Function[nrow(df.Objective.Function)+1,] <- c(i, stat.intrasegment.variance, NA, moransI$estimate[[1]], NA, NA)
    }


    if(Tool == "GRASS")
    {
      if(Grass.Objective.Function.Method == "Threshold")
      {
        df.Objective.Function[nrow(df.Objective.Function)+1,] <- c(i, Grass.Segmentation.Minsize, stat.intrasegment.variance, NA, moransI$estimate[[1]], NA, NA)
      }

      if(Grass.Objective.Function.Method == "Minsize")
      {
        df.Objective.Function[nrow(df.Objective.Function)+1,] <- c(Grass.Segmentation.Threshold, i, stat.intrasegment.variance, NA, moransI$estimate[[1]], NA, NA)
      }

    } # end GRASS

    if(Tool == "GRASS Superpixels SLIC")
    {
      if(Grass.Objective.Function.Method == "Compactness")
      {
        df.Objective.Function[nrow(df.Objective.Function)+1,] <- c(Grass.SLIC.Superpixels, i, stat.intrasegment.variance, NA, moransI$estimate[[1]], NA, NA)
      }

      if(Grass.Objective.Function.Method == "Superpixels")
      {
        df.Objective.Function[nrow(df.Objective.Function)+1,] <- c(i, Grass.SLIC.Compactness, stat.intrasegment.variance, NA, moransI$estimate[[1]], NA, NA)
      }

    } # end GRASS Superpixels SLIC


    print("")
  } # end for


  # Final normalization for output -----------------------------------------------------------------------
  # Normalizing and Calculation of Objective Function by Camara et al. 2006 - Parameter Selection for
  # Region-Growing Image Segmentation Algorithms using Spatial Autocorrelation
  df.Objective.Function$Normalized.Intrasegment.Variance <- (max(df.Objective.Function$Intrasegment.Variance) - df.Objective.Function$Intrasegment.Variance) / diff(range(df.Objective.Function$Intrasegment.Variance))
  df.Objective.Function$Normalized.Morans.I <- (max(df.Objective.Function$Morans.I) - df.Objective.Function$Morans.I) / diff(range(df.Objective.Function$Morans.I))


  df.Objective.Function$Objective.Function <-  df.Objective.Function$Normalized.Intrasegment.Variance + df.Objective.Function$Normalized.Morans.I

  # write Objective.Function as csv if desired
  if(Objective.Function.save == TRUE)
  {


    if(is.null(Objective.Function.save.path))
    {
      if(Objective.Function.save.name == "")
      {
        Objective.Function.save.name <- "Objective.Function.csv"
      }

      save.path <- paste0(getwd(), "/", Objective.Function.save.path.default, "/" , Objective.Function.save.name)
    } else
    {
      save.path <- Objective.Function.save.path
    }


    write.csv(df.Objective.Function, save.path)

  }


  # get time of process
  process.time.run.scale <- proc.time() - process.time.start.scale
  print(paste0("------ Run of Scale Estimation: " , process.time.run.scale["elapsed"][[1]]/60, " Minutes ------"))

  # remove things out of memory
  gc()

  return(df.Objective.Function)

}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



