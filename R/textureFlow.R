#' Calculation of textural features in flow direction
#'
#' This function calculates textural features in flow direction according to STUMPF & KERLE (2011)
#'
#' @param grass.texture.flowDir full path to flow direction raster which must be performed by the D-Infinity method of TauDEM. Moreover, the angles must be in degree and already rotated (clockwise, NORTH 0°)
#' @param grass.texture.input GRASS GIS raster on which the textural features are computed (commonly a DTM). Therefore, the input must be already imported or created in GRASS GIS
#' @param saga.texture.shape full path of shapefile for the computation of grid statistics. Default: NULL
#' @param saga.texture.statistic full storing path of shapefile with object statistics. Default: NULL
#' @param grass.texture.window window size for the computation of textural features. Default: 3
#' @param grass.texture.distance distance between two samples. The distance must be smaller than the size of the moving window (>= 1). Default: 1
#' @param grass.texture.method the textural features that the user wants to compute. Default: "contrast,corr,var,sa,entr"
#' @param grass.texture.name name of texture for automatism. Number and names must fit to \emph{grass.texture.method}. Default: c("_Contr_", "_Corr_", "_VAR_", "_SA_", "_Entr_")
#' @param grass.categories.rules full stroing path of ruleset for reclassification of TauDEM flow direction. Default: tempdir() + rules_texture
#' @param grass.texture.save possibility to save texture in flow direction as .tif. Default: FALSE (tif's are removed at the end)
#' @param grass.texture.save.path (optional) output path of texture calculated in flow direction. Outputs are in .tif format. Default: tempdir()
#' @param grass.texture.statistics .csv with shapefile attributes including texture statistics
#' @param texture.output.name vector containing names for texture output, [1] for flow direction, [2] for perpendicular flow direction. Information on window size  is automatically added after "_". Default: c("textureFlowDir", "textureFlowDirPer")
#' @param show.output.on.console show output of GRASS and SAGA GIS modules on console. Default: FALSE
#' @param quiet no outputs in console. Default: TRUE
#'
#' @return
#' (optional) \linkS4class{data.table} with texture statistics of shapefile input
#'
#' @note
#' \itemize{
#'   \item the textural features MOC-1 and MOC-2 are not supported, yet
#'   \item calculation of textures is for all directions separately
#'   \item GRASS GIS session (mapset, location and raster for texture) must have been initialized before
#'   \item GRASS GIS texture orientations are counter-clockwise beginning from EAST, but output is calculated clockwise beginning from NORTH
#'   \item TO DO: ratio of texture in flow direction to texture in perpendicular flow direction
#'   \item STUMPF, A. & N. KERLE (2011): Object-oriented mapping of landslides using Random Forests. - Remote Sensing of Environment 115, 10, 2564-2577
#' }
#'
#' @keywords textural features, GRASS GIS, texture in flow direction
#'
#'
#' @export
textureFlow <- function(grass.texture.flowDir, grass.texture.input, saga.texture.shape = NULL, saga.texture.statistic = NULL, grass.texture.window = 3, quiet = TRUE, show.output.on.console = FALSE,
                        grass.texture.distance = 1, grass.texture.method = "contrast,corr,var,sa,entr", grass.texture.name = c("_Contr_", "_Corr_", "_Var_", "_SA_", "_Entr_"),
                        grass.categories.rules = "rules_texture", grass.texture.save = FALSE, grass.texture.save.path = tempdir(), texture.output.name = c("textureFlowDir", "textureFlowDirPer"))
{

  # get start time of process
  process.time.start <- proc.time()

  # calculate textures for all directions separately and for selected features, CARE: MOC-1 and MOC-2 did not work
  # print(parseGRASS("r.texture"))
  print("Calculate Textures")
  execGRASS("r.texture", flags = c("overwrite", "s", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    input = grass.texture.input, output = "texture", size = grass.texture.window,
    distance = grass.texture.distance, method = grass.texture.method))

  # load rotated flow direction of TauDEM into GRASS (must be calculated before!!!)
  execGRASS('r.in.gdal',  flags=c('o',"overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    input = grass.texture.flowDir,  output='flowDirInfRot'))

  ### categories flow direction into groups that coinced the GLCM directions
  # GRASS texture angle is counterclockwise beginning from East
  # 0 E-W | 45 NE - SW | 90 N-S | 135 NW-SE

  # print(parseGRASS("r.recode"))
  # N and S = 1 | 90, NE and SW = 2 | 45, E and W = 3 | 0, NW - SE = 4 | 135
  rules.categories <- c("0:22.5:1", "337.5:360:1",  "157.5:202.5:1",
                        "22.5:67.5:2", "202.5:247.5:2",
                        "67.5:112.5:3", "247.5:292.5:3",
                        "292.5:337.5:4", "112.5:157.5:4")
  # rules must be written
  if(grass.categories.rules == "rules_texture")
  {
    rule_path_name <- file.path(tempdir(), grass.categories.rules)
  } else
  {
    rule_path_name <- grass.categories.rules
  }

  write(rules.categories, file = rule_path_name)

  execGRASS("r.recode", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    input = "flowDirInfRot", output = "flowDirReClass", rules = rule_path_name))

  ### create GLCM by Flow Direction
  # print(parseGRASS("r.mapcalc"))
  print("Calculate Textures in Flow Direction (and perpendicular)")
  for(i in grass.texture.name)
  {
    # create expression for texture by flow direction
    expression.i <- paste0("if(flowDirReClass == 1, ",  paste0("texture", i, "90"),
                           ", if(flowDirReClass == 2, ", paste0("texture", i, "45"),
                           ", if(flowDirReClass == 3, ", paste0("texture", i, "0"),
                           ", if(flowDirReClass == 4, ", paste0("texture", i, "135", ", null()))))"))

    # create expression for texture by flow direction (perpendicular!)
    expression.i.per <- paste0("if(flowDirReClass == 1, ",  paste0("texture", i, "0"),
                               ", if(flowDirReClass == 2, ", paste0("texture", i, "135"),
                               ", if(flowDirReClass == 3, ", paste0("texture", i, "90"),
                               ", if(flowDirReClass == 4, ", paste0("texture", i, "45", ", null()))))"))


    # calculate glcm flow direction textures
    execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      expression = paste0(texture.output.name[1], "_", gsub("[_-]", '', i),  " = ", expression.i)))

    # calculate glcm perpendicular flow direction textures
    execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      expression = paste0(texture.output.name[2], "_", gsub("[_-]", '', i),  " = ", expression.i.per)))
  }


  ### loop through names and save every raster as sdat - SAGA format (optional)
  print("Save Textures in Flow Direction")

  # for ovoiding overwrites every raster get a unique name based on the GLCM calculation (window size and distance)
  j <- paste0(grass.texture.window, grass.texture.distance)
  for(i in grass.texture.name)
  {
    # print(parseGRASS("r.out.gdal"))
    execGRASS("r.out.gdal",  flags=c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input =  paste0(texture.output.name[1], "_", gsub("[_-]", '', i)),  output= paste0(grass.texture.save.path, "/", paste0(texture.output.name[1], j, "_", gsub("[_-]", '', i)), ".tif")))

    execGRASS("r.out.gdal",  flags=c("overwrite", "quiet"),  Sys_show.output.on.console = show.output.on.console, parameters = list(
      input =  paste0(texture.output.name[2], "_", gsub("[_-]", '', i)),  output= paste0(grass.texture.save.path, "/", paste0(texture.output.name[2], j, "_", gsub("[_-]", '', i)), ".tif")))
  }


  ### get object statistics - average
  # if shape is set, then start statistic for a shapefile
  if(!is.null(saga.texture.shape))
  {

    if(is.null(saga.texture.statistic))
    {
      saga.texture.statistic = saga.texture.shape
    }


    # print(parseGRASS("v.rast.stats"))
    print("Calculate Grid Statistics")

    counter <- 0
    vNames <- c() # vector containing relevant column names (for later)

    for(i in grass.texture.name)
    {
      # loop counter to know how many textures gonna be calculated

      # ... in flow direction
      # execGRASS("v.rast.stats",  flags=c("c", "quiet"), intern = FALSE, parameters = list(
      #   map =  grass.texture.shape,  raster = paste0("textureFlowDir_", gsub("[_-]", '', i)),
      #   column_prefix = paste0(gsub("[_-]", '', i), j), method = "average"))
      # vNames <- c(vNames, paste0(gsub("[_-]", '', i), j, "_average"))

      # ... perpendicular to flow direction
      # execGRASS("v.rast.stats",  flags=c("c", "quiet"), intern = FALSE, parameters = list(
      #   map =  grass.texture.shape,  raster = paste0("textureFlowDirPer_", gsub("[_-]", '', i)),
      #   column_prefix = paste0(gsub("[_-]", '', i), j, "P"), method = "average"))
      # vNames <- c(vNames, paste0(gsub("[_-]", '', i), j, "P_average"))

      ### calculate grid statistics
      # rsaga.get.modules("shapes_grid", env = env)
      # rsaga.get.usage("shapes_grid", 2, env = env)
      # method: [0] standard | [1] shape wise, supports overlapping polygons
      grid.i <- paste0(grass.texture.save.path, "/", paste0(texture.output.name[1], j, "_", gsub("[_-]", '', i)), ".tif")
      grid.i.p <- paste0(grass.texture.save.path, "/", paste0(texture.output.name[2], j, "_", gsub("[_-]", '', i)), ".tif")

      if(counter == 0)
      {
        saga.texture.i <- saga.texture.shape # first process will be with original shape
      } else
      {
        saga.texture.i <- saga.texture.statistic # next runs or processes will be with statistic shape
      }

      rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = show.output.on.console, param = list(
        GRIDS = grid.i, POLYGONS = saga.texture.i, METHOD = "0", NAMING = 1, RESULT = saga.texture.statistic,
        COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = "0"))

      vNames <- c(vNames, paste0(gsub("[_-]", '', i), j))

      rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = show.output.on.console, param = list(
        GRIDS = grid.i.p, POLYGONS = saga.texture.statistic, METHOD = "0", NAMING = 1, RESULT = saga.texture.statistic,
        COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = "0"))

      vNames <- c(vNames, paste0(gsub("[_-]", '', i), j, "P"))

      counter <- counter + 2
    }



    ### read dbfs
    print("Create Tables and Modify DBFs")

    saga.texture.shape.dbf <- suppressMessages(read.dbf(paste0(file_path_sans_ext(saga.texture.shape), ".dbf")))
    texture.stat <- suppressMessages(read.dbf(paste0(file_path_sans_ext(saga.texture.statistic), ".dbf")))
    n <- ncol(texture.stat$dbf)
    colnames(texture.stat$dbf)[(n-counter+1):n] <- vNames
    write.dbf(texture.stat, paste0(file_path_sans_ext(saga.texture.statistic), ".dbf")) # write dbf with better header

    # select data for output
    dt.texture.stat <- data.table(texture.stat$dbf[c("ID", vNames)])
    saga.texture.shape.dbf$dbf <- cbind(saga.texture.shape.dbf$dbf, texture.stat$dbf[vNames])
    write.dbf(saga.texture.shape.dbf, paste0(file_path_sans_ext(saga.texture.shape), ".dbf"))
  }

  # remove saved textures
  if(grass.texture.save == FALSE)
  {
    print("Delete temporary files")

    for(i in grass.texture.name)
    {
      fn <- paste0(grass.texture.save.path, "/", paste0(texture.output.name[2], j, "_", gsub("[_-]", '', i)), ".tiff")
      if(file.exists(fn)) file.remove(fn)
    }
  }


  # get time of process
  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) print(paste0("------ Run of Texture Flow Calculation: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))

  if(exists("dt.texture.stat"))
  {
    return(dt.texture.stat)
  }
}
