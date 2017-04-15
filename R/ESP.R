#' ESP - ESTIMATION OF SCALE PARAMETER
#'
#' The ESP function by Dragut et al. (2010) calculates optimal
#' segmentation scales by calculating, at first, the local variance for every segmented images
#' and then by comparing the results of different segmentation levels (ROC-LV).
#'
#' @param Scale.Input.Grid input grid for computing segmentation scale parameters
#' @param Scale.Input.Grid.Cell.Size cell size of input grid. Default: "1"
#' @param Scale.Statistic.Min.Size min size of area/polygon which is included in statistics (usefull for SAGA GIS segmentations). Default: "0"
#' @param ESP.save save estimations of function. Default: FALSE
#' @param ESP.save.path path where estimations are stored. Default: input path of \emph{segment.poly}
#' @param Scales containing scale parameters for loop-segmentation
#' @param Count amount of cells (Grid Statistics). Default:"1"
#' @param Min minimum value (Grid Statistics). Default:"0"
#' @param Max maximum value (Grid Statistics). Default:"0"
#' @param Range range of values (Grid Statistics). Default:"0"
#' @param Sum sum of values (Grid Statistics). Default:"0"
#' @param Mean mean of values (Grid Statistics). Default:"1"
#' @param Var variance of values (Grid Statistics). Default:"0"
#' @param Stddev standard deviation (Grid Statistics). Default:"1"
#' @param Quantile qunatile (Grid Statistics). Default:"0"
#' @param Grass.ESP.Method determining on which parameter the objective function (~loop) should be performed. Default: "Threshold"
#' @param other... see \code{\link{segmentation}}
#'
#'
#' @note
#' \itemize{
#'   \item DRAGUT, L., TIEDE, D. & S.R. LEVICK (2010): ESP: a tool to estimate
#'     scale parameter for multiresolution image segmentation of remotely
#'     sensed data. - International Journal of Geographical Information
#'     Science 24, 6, 859-871.
#' }
#'
#' @keywords ESP, USPO
#'
#'
#' @export
ESP <- function(Scale.Input.Grid, Scale.Input.Grid.Cell.Size = "1", Scale.Statistic.Min.Size  = "0", ESP.save = FALSE, ESP.save.path = NULL, Count = "1", Min = "0", Max = "0", Range = "0", Sum = "0", Mean = "0",
                Var = "0", Stddev = "1", Quantile = 0, Scales, Grass.ESP.Method = "Threshold", ...)
{


  # get start time of process
  process.time.start <- proc.time()

  # check input
  if (is.null(Segments.Grid) || is.null(Segments.Poly) || is.null(Input.Grid) || is.null(Scales) || is.null(Scale.Input.Grid) || (NoData == TRUE & is.null(Mask)))
  {
    stop("Some input parameters are missing!")
  }


  # create data table for ESP
  df.ESP <- data.frame("Scale Parameter" = double(), "Mean SD of LV" = double(), "ROC-LV" = double())
  ESP.save.path.default <- dirname(Segments.Poly)
  print("calculate optimal scale parameters")

  # Start loop through scales -----------------------------------------------------------------------
  for(i in Scales)
  {


    multiplier <- 1

    if(Seed.Method == "Fast Representativeness")
    {
      multiplier <- 100
    }


    print(paste0("Level of Generalisation|Threshold|Minsize: ", i))
    segments.poly <- paste0(file_path_sans_ext(Segments.Poly) , (i*multiplier),".", file_ext(Segments.Poly))
    # segments.scale.statistic <- paste0(file_path_sans_ext(segments.poly) , "_scaleStat.", file_ext(segments.poly))
    segments.scale.statistic <- segments.poly # update after changed dbf header


    # Segmentation -----------------------------------------------------------------------
    # perform segmentation
    if(Tool == "SAGA")
    {
      segmentation(Fast.Representativeness.LevelOfGeneralisation = i, Seed.Generation.Scale = as.character(i), ...)
    }
    else if(Tool == "GRASS")
    {
      if(Grass.ESP.Method == "Threshold")
      {
        segmentation(Grass.Segmentation.Threshold = as.character(i), ...)
      }

      if(Grass.ESP.Method == "Minsize")
      {
        segmentation(Grass.Segmentation.Minsize = i, ...)
      }

      if(Grass.ESP.Method == "Seeds")
      {
        segmentation(Fast.Representativeness.LevelOfGeneralisation = i, Seed.Generation.Scale = as.character(i), ...)

      }

    }


    # Calculation of grid statistics -----------------------------------------------------------------------
    # calculate grid statistics
    # rsaga.get.modules("shapes_grid", env = env)
    # rsaga.get.usage("shapes_grid", 2, env = env)
    # method: [0] standard | [1] shape wise, supports overlapping polygons
    print("Calculate Grid Statistics")
    rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = SOOC, param = list(
      GRIDS = Scale.Input.Grid, POLYGONS = segments.poly, METHOD = "0", NAMING = 1, RESULT = segments.scale.statistic,
      COUNT = Count, MIN = Min, MAX = Max, RANGE = Range, SUM = Sum, MEAN = Mean, VAR = Var, STDDEV = Stddev, QUANTILE = Quantile))

    # read statistic dbf
    print("Extraction for Estimate Scale Parameters")
    stat <- suppressMessages(read.dbf(paste0(file_path_sans_ext(segments.scale.statistic), ".dbf")))
    colnames(stat$dbf) <- c("segments_g", "ID", "NAME", "Cell_Count", "Stddev")
    write.dbf(stat, paste0(file_path_sans_ext(segments.scale.statistic), ".dbf")) # write dbf with better header

    ### calculation of weighted mean for that little region get less influcence
    stat$dbf <- stat$dbf[stat$dbf$Cell_Count > (as.numeric(Scale.Statistic.Min.Size) * as.numeric(Scale.Input.Grid.Cell.Size)),]
    stat.mean.SD <- weighted.mean(stat$dbf[[ncol(stat$dbf)]], stat$dbf[[ncol(stat$dbf)-1]] * as.numeric(Scale.Input.Grid.Cell.Size), na.rm = TRUE) # weighted mean of variance by area

    df.ESP[nrow(df.ESP)+1,] <- c(i, stat.mean.SD, NA)

    print("")
  } # end for


  # Calculation of ROC-LV -----------------------------------------------------------------------
  # calculate ROC-LV by Dragut, Tiede & Levick (2010), S. 863
  for(j in 1:nrow(df.ESP))
  {
    if((j+1) > nrow(df.ESP))
    {
      df.ESP$ROC.LV[j] <- 0
    }
    else
    {
      df.ESP$ROC.LV[j] <- (df.ESP$Mean.SD.of.LV[j] - df.ESP$Mean.SD.of.LV[j+1])/df.ESP$Mean.SD.of.LV[j+1]
    }

  }

  # write ESP as csv if desired
  if(ESP.save == TRUE)
  {


    if(is.null(ESP.save.path))
    {
      save.path <- paste0(getwd(), "/", ESP.save.path.default, "/ESP.csv")
    }
    else
    {
      save.path <- ESP.save.path
    }


    write.csv(df.ESP, save.path)

  }

  # get time of process
  process.time.run <- proc.time() - process.time.start
  print(paste0("------ Run of Scale Estimation: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))

  # remove things out of memory
  gc()

  return(df.ESP)

}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
