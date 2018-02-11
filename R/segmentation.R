#' Segmentation
#'
#' This function uses the open-source software of GRASS and SAGA GIS
#' to fulfill an imagery segmentation. It is possible to decide bet
#' ween SAGA and GRASS GIS segmentation.\cr \cr
#' In SAGA GIS, the tools SEED GENERATION or FAST REPRESENTATIVENESS
#' are used for computing seed points. Then, the SEEDED REGION GROWING al-
#' gorithm is used for the segmentation.\cr \cr
#' In GRASS GIS, the tool OBJECT-SEGMENTATION (i.segment) is used
#' for the segmentation. There is also the possibility to use
#' the new SLIC algorithm.\cr \cr
#' Moreover, there is the option to compute a generalisation of
#' segments by (multiple) majority filter (SAGA GIS) at the end.
#'
#' @param Tool GRASS, SAGA or GRASS Superixels SLIC. Definition of open-source software which will be used
#' @param Input.Grid vector containing grid(s) for segmentation. It is possible to add multiple gridsm, as well as different grids even in combination of SAGA and GRASS. By using SAGA and GRASS combination the following separation must be used: '<>' (SAGA before GRASS grids!)
#' @param Segments.Grid output path of raster with segments
#' @param Segments.Poly output path of polygon with segments
#' @param NoData input data contains NoData value. Default: FALSE
#' @param Mask mask raster to mask NoData from input. Default: NULL
#' @param show.output.on.console show output on console. Default: FALSE
#' @param env environment of RSAGA. Default: RSAGA::rsaga.env()
#' @param Generalisation.Flac performing (multiple) majority filter on segments. Default: FALSE
#' @param Generalization.Mode search mode by filtering: Default: "1" (circle, alternative: square)
#' @param Generalization.Threshold threshold for applying majority filters. Default: "0.0"
#' @param Seed.Method type of seed method for getting seeds. Default: "" (alternative: "Fast Representativeness", "Seed Generation",  "Superpixels SLIC")
#' @param Seed.Generation.Variance output raster with variance of seed generation. Default: temp
#' @param Seed.Generation.Points output of seed points as shapefile. Default: temp
#' @param Output.Seeds output of seed points as raster, used for segmentation. Default: temp
#' @param Seed.Generation.Type option of seed generation type. Default:"0" (minima of variance, alternative: maxima of variance)
#' @param Seed.Generation.Scale determining number of seed points in seed generation. Default: "10.0"
#' @param Saga.Output.Grid output of FAST REPRESENTATIVENESS in SAGA. Default: temp
#' @param Saga.Output.Lod output Lod of Representativeness Function in SAGA. Default: temp
#' @param Fast.Representativeness.LevelOfGeneralisation determining number of seed points. Default: "10.0"
#' @param Saga.Similarity output of similarity grid. Default: temp
#' @param Saga.Segments.Seeds.Table table of seeds information. Default: temp
#' @param Saga.Segmentation.Normalize normalisation during imagery segmentation. Default: "0"
#' @param Saga.Segmentation.Neighbourhood neighbourhood considered during imagery segmentation. Default: "0" (4, alternative: "1" for 8)
#' @param Saga.Segmentation.Method segmentation method during imagery segmentation. Default: "0" (feature space and position)
#' @param Saga.Segmentation.Sig.1 variance in feature space in imagery segmentation. Default: "1.0"
#' @param Saga.Segmentation.Sig.2 variance in position space in imagery segmentation. Default: "1.0"
#' @param Saga.Segmentation.Threshold similarity threshold for joining pixel to segments in imagery segmentation. Default: "0.0"
#' @param Saga.Segmentation.Refresh refresh image after imagery segmentation. Default: "0"
#' @param Saga.Segmentation.Leafsize parameter for speed optimation in imagery segmentation. Default: 256
#' @param Split split polygons to singlepart in vectorising grid classes. Default: "0"
#' @param AllVertices use all vertices by vectorising grid classes. Default: "FALSE"
#' @param Grass.Segmentation.Threshold similarity threshold for joining pixel to segments. Default: NULL, 0.0 is not allowed
#' @param Grass.Segmentation.Weighted option of weighing input grids in segmentation. Default: "FALSE"
#' @param Grass.Segmentation.Method type of GRASS Segmentation. Default: "region_growing"
#' @param Grass.Segmentation.Similarity distance measurement of similarity. Default: "euclidean"
#' @param Grass.Segmentation.Minsize minsize of segment. Default: 15
#' @param Grass.Segmentation.Memory memory to be used for segmentation. Default: 300
#' @param Grass.Segmentation.Iterations amount of allowed iterations. Default: 50
#' @param Grass.Segmentation.Seeds input of seeds raster. Enables bottom-up segmentation. Default: NULL
#' @param Segmentation.Boundary.Grid input of boundary raster. Enables top-down (or hierarchical) segmentation. Default: NULL. NULL values must be 0 (or any other not segment value!).
#' @param Grass.Segmentation.Goodness name for output goodness of fit estimate map. Default:"Grass.Segmentation.Goodness"
#' @param Grass.SLIC.Iter maximum number of iterations. Default: 10
#' @param Grass.SLIC.Superpixels approximate number of output super pixels. Default: 200
#' @param Grass.SLIC.Step distance (number of cells) between initial super pixel centers. A step size > 0 overrides the number of super pixels. Default: 0
#' @param Grass.SLIC.Compactness compactness. A larger value causes more compact superpixels. Default: 1.0
#' @param Grass.SLIC.Superpixels.MinSize minimum superpixel size. Default: 1
#' @param Grass.SLIC.Memory memory in MB. Default: 300
#' @param Grass.SLIC.Perturb Perturb initial super pixel centers. Percent of intitial superpixel radius. Default: 0, range: 0-100
#' @param burn.Boundary.into.Segments vector specifing if boundary grid is burned into segmentation (1) or seeds (2). Default: FALSE, maximum length: 2
#' @param estimateScaleParameter must be only be TRUE when scale parameter function is used. Default: FALSE
#' @param Mode.Filter.Flac re-assign objects of a specific size based on mode values to its surroundings (moving window). Default: FALSE
#' @param Mode.Filter.Size moving window size of mode-filter. Default: 3
#' @param Mode.Filter.Segment.MinSize objects smaller and equal to this size are selected for filtering. Default: 3
#' @param par.i run number. Default: ""
#'
#' @keywords segmentation, region growing, superpixels Simple Linear Iterative Clustering
#'
#'
#' @export
segmentation <- function(Tool, Segments.Grid, Segments.Poly, Input.Grid, Saga.Output.Grid = file.path(tempdir(), paste0("SagaRepresentativenessOutputGrid", par.i,  ".sgrd")), Saga.Output.Lod =  file.path(tempdir(), paste0("SagaRepresentativenessOutputLod", par.i, ".sgrd")),
                         Output.Seeds =  file.path(tempdir(), paste0("OutputSeed", par.i, ".sgrd")), Fast.Representativeness.LevelOfGeneralisation = "10.0", Saga.Similarity =  file.path(tempdir(),paste0("SagaSimilarity", par.i, ".sgrd")), Saga.Segments.Seeds.Table =  file.path(tempdir(), paste0("SagaSegmentsSeedsTable", par.i, ".mtab")),  Saga.Segmentation.Normalize = "0", Saga.Segmentation.Neighbourhood = "0", Saga.Segmentation.Method = "0",
                         Saga.Segmentation.Sig.1 = "1.0", Saga.Segmentation.Sig.2 = "1.0", Saga.Segmentation.Threshold = "0.0", Saga.Segmentation.Refresh = "0", Saga.Segmentation.Leafsize = 256, Split = "0",
                         Grass.Segmentation.Threshold = NULL, Grass.Segmentation.Weighted = FALSE, Grass.Segmentation.Method = "region_growing", Grass.Segmentation.Similarity = "euclidean", Grass.Segmentation.Minsize = 15, Grass.Segmentation.Memory = 300, Grass.Segmentation.Iterations = 50, Grass.Segmentation.Seeds = NULL, Grass.Segmentation.Neighbourhood = "0",
                         Segmentation.Boundary.Grid = NULL,  Grass.Segmentation.Goodness = paste0("Grass.Segmentation.Goodness", par.i), AllVertices = "FALSE", NoData = FALSE, Mask = NULL, NoData.Flag = -99999, show.output.on.console = FALSE, Seed.Method = "", Seed.Generation.Variance =  paste0(tempdir(), paste0("SeedGenerationVariance", par.i, ".sgrd")), Seed.Generation.Points =  file.path(tempdir(), paste0("SeedGenerationPoints", par.i, ".shp")),
                         Seed.Generation.Type = "0", Seed.Generation.Scale = "10.0", Generalisation.Flac = FALSE, Generalization.Mode = "1", Generalization.Radius = "1", Generalization.Threshold = "0.0", env = RSAGA::rsaga.env(),
                         Grass.SLIC.Iter = 10, Grass.SLIC.Superpixels = 200, Grass.SLIC.Step = 0, Grass.SLIC.Compactness = 1.0, Grass.SLIC.Superpixels.MinSize = 1, Grass.SLIC.Memory = 300, Grass.SLIC.Perturb = 0, burn.Boundary.into.Segments = FALSE,
                         estimateScaleParameter = FALSE, Mode.Filter.Flac = FALSE, Mode.Filter.Size = 7, Mode.Filter.Segment.MinSize = 3, par.i = "", ...)
{

 # browser()

  # get start time of process
  process.time.start <- proc.time()

  # check if input is missing
  if (is.null(Tool) || is.null(Segments.Grid) || is.null(Segments.Poly) || is.null(Input.Grid) || (NoData == TRUE & is.null(Mask)))
  {
    stop("Some input parameters are missing!")
  }


  # check input grids for separate SAGA and GRASS input if combined
  if(!is.na(match("<>", Input.Grid)))
  {
    pos <- match("<>", Input.Grid)
    Input.Grid.Saga <- Input.Grid[1:(pos-1)]
    Input.Grid.Saga <- paste(Input.Grid.Saga, collapse = ";")
    Input.Grid.GRASS <- Input.Grid[(pos+1):length(Input.Grid)]
  } else
  {
    # transform input grids to saga format, for grass no transformation
    Input.Grid.Saga <- paste(Input.Grid, collapse = ";")
    Input.Grid.GRASS <- Input.Grid
  }


  #  Seed Points------------------------------------------------------------------------
  if(Seed.Method == "Fast Representativeness")
  {
    # perform fast representativeness for getting seed points
    # rsaga.get.modules("statistics_grid", env = env)
    # rsaga.get.usage("statistics_grid", 0, env = env)
    print("perform fast representativeness for getting seed points")

    Input.Grid.Strings <- strsplit(Input.Grid.Saga, ";")

    if(length(Input.Grid.Strings[[1]]) == 1) # check if 2 grids are given
    {
      RSAGA::rsaga.geoprocessor(lib="statistics_grid", module = 0, env = env, show.output.on.console = show.output.on.console, param = list(
        INPUT = Input.Grid.Saga, RESULT = Saga.Output.Grid, RESULT_LOD = Saga.Output.Lod, SEEDS = Output.Seeds,
        LOD = Fast.Representativeness.LevelOfGeneralisation))
    }

    if(length(Input.Grid.Strings[[1]]) >= 2) # ... if yes, take first one for representativeness
    {
      print("... first grid of multiple is taken for seed points")
      Input.Grid1 <- Input.Grid.Strings[[1]][1]
      RSAGA::rsaga.geoprocessor(lib="statistics_grid", module = 0, env = env, show.output.on.console = show.output.on.console, param = list(
        INPUT = Input.Grid1, RESULT = Saga.Output.Grid, RESULT_LOD = Saga.Output.Lod, SEEDS = Output.Seeds,
        LOD = Fast.Representativeness.LevelOfGeneralisation))
    }

  } else if(Seed.Method == "Seed Generation")
  {
    # perform fast representativeness for getting seed points
    # rsaga.get.modules("imagery_segmentation", env = env)
    # rsaga.get.usage("imagery_segmentation", 2, env = env)
    # SEED_TYPE: [0] minima of variance |  [1] maxima of variance
    # METHOD: [0] band width smoothing | [1] band width search
    print("perform seed generation for getting seed points")

    RSAGA::rsaga.geoprocessor(lib="imagery_segmentation", module = 2, env = env, show.output.on.console = show.output.on.console, param = list(
      FEATURES = Input.Grid.Saga, VARIANCE = Seed.Generation.Variance, SEED_GRID =  Output.Seeds, SEED_POINTS = Seed.Generation.Points,
      SEED_TYPE = Seed.Generation.Type, METHOD = "0", NORMALIZE = "0", BAND_WIDTH = Seed.Generation.Scale))

  } else if(Seed.Method == "Superpixels SLIC")
  {
    print("GRASS Superpixels SLIC as Seed Method is selected")

  } else{

    print("No Seed Method was used!")
  }


  # check par.i
  # if(par.i != "")
  # {
  #   par.v <- gsub("[.]", "_", as.character(par.i))
  # } else {
  #   par.v <- par.i
  # }

  #  Segmentation ------------------------------------------------------------------------
  if(Tool != "SAGA" & Tool != "GRASS" & Tool != "GRASS Superpixels SLIC")
  {
    stop("Wrong tool software was selected. Please use 'SAGA' or 'GRASS' for segmentation procedure.")
  }

  #  Segmentation SAGA ------------------------------------------------------------------------
  if(Tool == "SAGA")
  {
    if(tools::file_ext(Segments.Grid) != "sgrd")
    {
      Segments.Grid.tmp <- file.path(tempdir(), paste0(tools::file_path_sans_ext(basename(Segments.Grid)), ".sgrd"))
    } else {
      Segments.Grid.tmp <- Segments.Grid
    }


    if(NoData == TRUE)
    {

      Mask.tmp.path <- Mask

      if(tools::file_ext(Mask) != "sgrd")
      {
        # RSAGA::rsaga.get.usage(lib="io_gdal", module = 1, env = env)
        # [4] GeoTIFF ... higher release: 2.3.1 format = "7"
        RSAGAUsage <- RSAGA::rsaga.get.usage(lib="io_gdal", module = 1, env = env, show = FALSE)
        formatTIFF <- gsub("\\D", "", grep('GeoTIFF', RSAGAUsage, value = TRUE))

        Mask.tmp.path <- file.path(tempdir(), "tmp_mask.sgrd")

        RSAGA::rsaga.geoprocessor(lib="io_gdal", module = 1, env = env, show.output.on.console = show.output.on.console, param = list(
          GRIDS = Mask, FILE = Mask.tmp.path, FORMAT = formatTIFF))
      }

      # rsaga.get.modules("grid_tools", env = env)
      # rsaga.get.usage("grid_tools", 24, env = env)
      print("... masking No Data Area in Saga Seeds")
      RSAGA::rsaga.geoprocessor(lib="grid_tools", module = 24, env = env, show.output.on.console = show.output.on.console, param = list(
        GRID = Output.Seeds, MASK = Mask.tmp.path, MASKED = Output.Seeds))
    }


    # perform SAGA seeded region growing
    # rsaga.get.modules("imagery_segmentation", env = env)
    # rsaga.get.usage("imagery_segmentation", 3, env = env)
    print("SAGA: Perform seeded region growing")

    RSAGA::rsaga.geoprocessor(lib="imagery_segmentation", module = 3, env = env, show.output.on.console = show.output.on.console, param = list(
      SEEDS = Output.Seeds, FEATURES= Input.Grid.Saga, SEGMENTS = Segments.Grid.tmp, SIMILARITY = Saga.Similarity, TABLE = Saga.Segments.Seeds.Table,
      NEIGHBOUR = Saga.Segmentation.Neighbourhood, METHOD = Saga.Segmentation.Method, SIG_1 = Saga.Segmentation.Sig.1, SIG_2 = Saga.Segmentation.Sig.2,
      THRESHOLD = Saga.Segmentation.Threshold, LEAFSIZE = Saga.Segmentation.Leafsize))
  # NORMALIZE = Saga.Segmentation.Normalize, REFRESH = Saga.Segmentation.Refresh

    if(tools::file_ext(Segments.Grid) != "sgrd")
    {
      # RSAGA::rsaga.get.usage(lib="io_gdal", module = 1, env = env)
      # [4] GeoTIFF ... higher release: 2.3.1 format = "7"
      RSAGAUsage <- RSAGA::rsaga.get.usage(lib="io_gdal", module = 1, env = env, show = FALSE)
      formatTIFF <- gsub("\\D", "", grep('GeoTIFF', RSAGAUsage, value = TRUE))

      RSAGA::rsaga.geoprocessor(lib="io_gdal", module = 1, env = env, show.output.on.console = show.output.on.console, param = list(
      GRIDS = Segments.Grid.tmp, FILE = Segments.Grid, FORMAT = formatTIFF))

    } else {
      Segments.Grid <- Segments.Grid.tmp
    }

  }



  #  Segmentations of GRASS ------------------------------------------------------------------------
  if(Tool == "GRASS" | Tool == "GRASS Superpixels SLIC")
  {
    # perform GRASS segmentation

    if(Tool == "GRASS" && (as.numeric(Grass.Segmentation.Threshold) <= 0 | as.numeric(Grass.Segmentation.Threshold) >=1))
    {
      print("Threshold value is wrong (0-1)! Default value of 0.25 is used.")
      Grass.Segmentation.Threshold <- "0.25"
    }


    # if a seed method is set, SAGA output.seed have to be transformed for GRASS input
    # UPDATE! SAGA Seed Output as GRASS Seed Input does not make sense!
    #         Therefore, the segmentation result will be used as seed Input!
    #         For the future as alternative: i.superpixels.slic
    #         https://grass.osgeo.org/grass72/manuals/addons/i.superpixels.slic.html
    if((Seed.Method  != "") & (Seed.Method != "Superpixels SLIC"))
    {

      # if(Seed.Method == "Seed Generation") # Fast Representativeness Seeds are already in the GRASS Format
      # {
      #   # reclassify seed output: all seeds prixels to value 1
      #   # rsaga.get.modules("grid_tools", env = env)
      #   # rsaga.get.usage("grid_tools", 15, env = env)
      #   # METHOD: [0] single | [1] range | [2] simple table | [3] user supplied table
      #   # RESULT_NODATA_CHOICE: [0] NoData value of input grid
      #   # SOPERATOR: [0] = | [1] < | [2] <= | [3] >= | [4] >
      #   rsaga.geoprocessor(lib="grid_tools", module = 15, env = env, show.output.on.console = show.output.on.console, param = list(
      #     INPUT = Output.Seeds, RESULT = Output.Seeds, METHOD = "0", OLD = "0.0", NEW = "1.0", SOPERATOR = "3",
      #     RESULT_NODATA_CHOICE = "0"))
      #
      #   # change data storage of reclassified seed output: from float to integer for GRASS
      #   # rsaga.get.modules("grid_tools", env = env)
      #   # rsaga.get.usage("grid_tools", 11, env = env)
      #   # TYPE: [1] unsigned 1 byte integer
      #   rsaga.geoprocessor(lib="grid_tools", module = 11, env = env, show.output.on.console = show.output.on.console, param = list(
      #     INPUT = Output.Seeds, OUTPUT = Output.Seeds, TYPE = "1", OFFSET = "0.0", SCALE = "1.0"))
      # }
      #
      # import the SAGA file in GRASS
      # print(parseGRASS("r.in.gdal"))
      # Grass.Segmentation.Seeds <- paste0(tools::file_path_sans_ext(Output.Seeds) , ".sdat")
      # execGRASS("r.in.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      #   input = Grass.Segmentation.Seeds, output = "Grass.Segmentation.Seeds"))


      # browser()

      if(tools::file_ext(Segments.Grid) != "sgrd")
      {
        Segments.Grid.tmp <- file.path(tempdir(), paste0(tools::file_path_sans_ext(basename(Segments.Grid)), ".sgrd"))
      } else {
        Segments.Grid.tmp <- Segments.Grid
      }

      if(NoData == TRUE)
      {
        # if(tools::file_ext(Mask) != "sgrd")
        # {
        #   # Mask.tmp <- Mask
        #   Mask.tmp <- raster::raster(Mask)
        #   Mask.tmp.path <- file.path(tempdir(), "tmp_mask.sgrd")
        #   raster::writeRaster(x = Mask.tmp, filename = Mask.tmp.path, NAflag = NoData.Flag, overwrite = TRUE)
        #   Mask <- Mask.tmp.path
        # }
        # # rsaga.get.modules("grid_tools", env = env)
        # # rsaga.get.usage("grid_tools", 24, env = env)
        # print("... masking No Data Area in Saga Seeds")
        # RSAGA::rsaga.geoprocessor(lib="grid_tools", module = 24, env = env, show.output.on.console = show.output.on.console, param = list(
        #   GRID = Output.Seeds, MASK = Mask, MASKED = Output.Seeds))

        Mask.tmp.path <- Mask

        if(tools::file_ext(Mask) != "sgrd")
        {
          # RSAGA::rsaga.get.usage(lib="io_gdal", module = 1, env = env)
          # [4] GeoTIFF ... higher release: 2.3.1 format = "7"
          RSAGAUsage <- RSAGA::rsaga.get.usage(lib="io_gdal", module = 1, env = env, show = FALSE)
          formatTIFF <- gsub("\\D", "", grep('GeoTIFF', RSAGAUsage, value = TRUE))

          Mask.tmp.path <- file.path(tempdir(), "tmp_mask.sgrd")

          RSAGA::rsaga.geoprocessor(lib="io_gdal", module = 1, env = env, show.output.on.console = show.output.on.console, param = list(
            GRIDS = Mask, FILE = Mask.tmp.path, FORMAT = formatTIFF))
        }

        # rsaga.get.modules("grid_tools", env = env)
        # rsaga.get.usage("grid_tools", 24, env = env)
        print("... masking No Data Area in Saga Seeds")
        RSAGA::rsaga.geoprocessor(lib="grid_tools", module = 24, env = env, show.output.on.console = show.output.on.console, param = list(
          GRID = Output.Seeds, MASK = Mask.tmp.path, MASKED = Output.Seeds))
      }

      # perform SAGA seeded region growing
      # rsaga.get.modules("imagery_segmentation", env = env)
      # rsaga.get.usage("imagery_segmentation", 3, env = env)
      print("SAGA: Perform seeded region growing for GRASS Seeds")
      RSAGA::rsaga.geoprocessor(lib="imagery_segmentation", module = 3, env = env, show.output.on.console = show.output.on.console, param = list(
        SEEDS = Output.Seeds, FEATURES= Input.Grid.Saga, SEGMENTS = Segments.Grid.tmp, SIMILARITY = Saga.Similarity, TABLE = Saga.Segments.Seeds.Table,
        NEIGHBOUR = Saga.Segmentation.Neighbourhood, METHOD = Saga.Segmentation.Method, SIG_1 = Saga.Segmentation.Sig.1, SIG_2 = Saga.Segmentation.Sig.2,
        THRESHOLD = Saga.Segmentation.Threshold, LEAFSIZE = Saga.Segmentation.Leafsize))
    # NORMALIZE = Saga.Segmentation.Normalize depricated?
    # REFRESH = Saga.Segmentation.Refresh ?BUG


      #  Burn Boundary into Seed ------------------------------------------------------------------------
      if(length(burn.Boundary.into.Segments) == 2 && burn.Boundary.into.Segments[2] == TRUE)
      {
        print("Burn Boundary into GRASS Seeds")

        Segments.Grid.r <- raster::raster(paste0(tools::file_path_sans_ext(Segments.Grid.tmp), ".sdat"))
        val <- raster::cellStats(x = Segments.Grid.r, stat = max)

        Segmentation.Boundary.Grid.r <- raster::raster(x = Segmentation.Boundary.Grid)
        Segments.Grid.r <- raster::overlay(x = Segmentation.Boundary.Grid.r, Segments.Grid.r,
                                        fun = function(x, y){return(ifelse(x > 0, x+val, y))})



        # Segments.Grid.rgdal <- paste0(tools::file_path_sans_ext(Segments.Grid.tmp), ".sdat")
        # Segments.Grid.rgdal.grid <- rgdal::readGDAL(fname = Segments.Grid.rgdal, silent = TRUE)
        # val <- max(Segments.Grid.rgdal.grid$band1, na.rm = TRUE)

        # formula for calculation
        # formula.expression <- paste0("ifelse(a > 0, (a+", val, "), b)")

        # boundary grid should have 0 for No Data values, or values that are not burned

        # "ifelse(a > 0, a, b)"
        # RSAGA::rsaga.geoprocessor(lib = "grid_calculus", env = env, module = 1, show.output.on.console = show.output.on.console, param=list(
        #   GRIDS =  paste(c(Segmentation.Boundary.Grid, Segments.Grid), collapse = ";"), RESULT = Segments.Grid, FORMULA = formula.expression,
        #   FNAME = "1", USE_NODATA = "1", TYPE = "5"))

        # rm(formula.expression, val, Segments.Grid.rgdal.grid, Segments.Grid.rgdal)
        rm(val, Segmentation.Boundary.Grid.r )
      }

      # change data storage of reclassified seed output: from float to integer for GRASS
      # rsaga.get.modules("grid_tools", env = env)
      # rsaga.get.usage("grid_tools", 11, env = env)
      # TYPE: [1] unsigned 1 byte integer, [5] unsigned 4 byte integer
      # RSAGA::rsaga.geoprocessor(lib="grid_tools", module = 11, env = env, show.output.on.console = show.output.on.console, param = list(
      #   INPUT = Segments.Grid.tmp, OUTPUT = Segments.Grid.tmp, TYPE = "5")) # OFFSET = "0.0", SCALE = "1.0"

      if(exists("Segments.Grid.r"))
      {
       rgdal::writeGDAL(dataset = as(Segments.Grid.r, "SpatialGridDataFrame"), fname = paste0(tools::file_path_sans_ext(Segments.Grid.tmp), ".tif"),
                         type = 'UInt32')

      } else {
        Segments.Grid.r <- raster::raster(Segments.Grid.tmp)

        rgdal::writeGDAL(dataset = as(Segments.Grid.r, "SpatialGridDataFrame"), fname = paste0(tools::file_path_sans_ext(Segments.Grid.tmp), ".tif"),
                         type = 'UInt32')
      }

      Grass.Segmentation.Seeds <- paste0(tools::file_path_sans_ext(Segments.Grid.tmp), ".tif")
      rgrass7::execGRASS("r.in.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        input = Grass.Segmentation.Seeds, output = paste0("Grass.Segmentation.Seeds", par.i)))

      rm(Segments.Grid.r, Segments.Grid.tmp)

    } # # # # # end if seed method

    ### test if multiple input files are set
    if(length(Input.Grid.GRASS) == 1) # 1 input grid
    {
      if(Grass.Segmentation.Neighbourhood == "0")
      {
        flag.Segmentation.GRASS <- c("overwrite", "quiet")
      }

      if(Grass.Segmentation.Neighbourhood == "1")
      {
        flag.Segmentation.GRASS <- c("overwrite", "quiet", "d")
      }

      if(tools::file_ext(Input.Grid.GRASS) == "sgrd")
      {
        inputGRASS <- paste0(tools::file_path_sans_ext(Input.Grid.GRASS) , ".sdat")
      } else{
        inputGRASS <- Input.Grid.GRASS
      }

      # print(parseGRASS("r.in.gdal"))
      rgrass7::execGRASS("r.in.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        input = inputGRASS, output = paste0("input.GRASS", par.i)))

      # print(parseGRASS("i.group"))
      rgrass7::execGRASS("i.group", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        group = paste0("GRASS.Segmentation.Group", par.i), paste0("input.GRASS", par.i))) # if this already exists it will be skipped

    } else # multiple input grids
    {

      if(Grass.Segmentation.Weighted == FALSE)
      {
        if(Grass.Segmentation.Neighbourhood == "0")
        {
          flag.Segmentation.GRASS <- c("overwrite", "quiet")
        }

        if(Grass.Segmentation.Neighbourhood == "1")
        {
          flag.Segmentation.GRASS <- c("overwrite", "quiet", "d")
        }
      } else
      {
        if(Grass.Segmentation.Neighbourhood == "0")
        {
          flag.Segmentation.GRASS <- c("w", "overwrite", "quiet")
        }

        if(Grass.Segmentation.Neighbourhood == "1")
        {
          flag.Segmentation.GRASS <- c("w", "overwrite", "quiet", "d")
        }
      }

      in.grids <- c()

      for(j in 1:length(Input.Grid.GRASS)) # read all grass input grids
      {
        # print(parseGRASS("r.in.gdal"))
        inputGRASS <- paste0(tools::file_path_sans_ext(Input.Grid.GRASS[j]) , ".sdat")
        in.grids[j] <- paste0("input.GRASS.", j)
        rgrass7::execGRASS("r.in.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
          input = inputGRASS, output = in.grids[j]))
      }

      # print(parseGRASS("i.group"))
      rgrass7::execGRASS("i.group", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        group = paste0("GRASS.Segmentation.Group", par.i), input = in.grids)) # if this already exists it will be skipped
    }


    if(Seed.Method == "Superpixels SLIC")
    {
      # https://grass.osgeo.org/grass72/manuals/addons/i.superpixels.slic.html
      # print(parseGRASS("i.superpixels.slic "))
      print("GRASS: Perform Superpixels SLIC for GRASS Seeds")

      # changed since UPDATE!
      # execGRASS("i.superpixels.slic", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      #   group = "GRASS.Segmentation.Group", output = "Grass.Segmentation.Seeds",
      #   iter = Grass.SLIC.Iter, k = Grass.SLIC.Superpixels, step = Grass.SLIC.Step, co = Grass.SLIC.Compactness, min = Grass.SLIC.Superpixels.MinSize))
      #
      rgrass7::execGRASS("i.superpixels.slic", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        group = paste0("GRASS.Segmentation.Group", par.i), output = paste0("Grass.Segmentation.Seeds", par.i), memory = Grass.SLIC.Memory, perturb = Grass.SLIC.Perturb,
        iterations = Grass.SLIC.Iter, num_pixels = Grass.SLIC.Superpixels, step = Grass.SLIC.Step, compactness = Grass.SLIC.Compactness, minsize = Grass.SLIC.Superpixels.MinSize))


      Grass.Segmentation.Seeds <- paste0("Grass.Segmentation.Seeds", par.i)

    }



    #  Segmentation GRASS ------------------------------------------------------------------------
    if(Tool == "GRASS")
    {
      print("GRASS: Perform Segmentation")
      # print(parseGRASS("i.segment"))
      if(!is.null(Grass.Segmentation.Seeds) & is.null(Segmentation.Boundary.Grid))
      {
        print("... with Seeds")
        rgrass7::execGRASS("i.segment", flags = flag.Segmentation.GRASS, Sys_show.output.on.console = show.output.on.console, parameters = list(
          group = paste0("GRASS.Segmentation.Group", par.i), output = paste0("output.segmentation.GRASS", par.i), threshold = as.numeric(Grass.Segmentation.Threshold),
          method = Grass.Segmentation.Method, similarity = Grass.Segmentation.Similarity, minsize = Grass.Segmentation.Minsize,
          memory = Grass.Segmentation.Memory, iterations = Grass.Segmentation.Iterations, seeds = paste0("Grass.Segmentation.Seeds", par.i)))
      }

      if(!is.null(Segmentation.Boundary.Grid))
      {
        if(tools::file_ext(Segmentation.Boundary.Grid) == "sgrd")
        {
          Grass.Segmentation.Boundary.Grid.GRASS <- paste0(tools::file_path_sans_ext(Segmentation.Boundary.Grid) , ".sdat")

        } else {
          Grass.Segmentation.Boundary.Grid.GRASS <-Segmentation.Boundary.Grid

        }

        # no data values must be 0
        rgrass7::execGRASS("r.in.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
          input = Grass.Segmentation.Boundary.Grid.GRASS, output = paste0("Grass.Segmentation.Boundary.Grid", par.i)))

        if(is.null(Grass.Segmentation.Seeds))
        {
          # browser()
          print("... with Boundary")
          rgrass7::execGRASS("i.segment", flags = flag.Segmentation.GRASS, Sys_show.output.on.console = show.output.on.console, parameters = list(
            group = paste0("GRASS.Segmentation.Group", par.i), output = paste0("output.segmentation.GRASS", par.i), threshold = as.numeric(Grass.Segmentation.Threshold),
            method = Grass.Segmentation.Method, similarity = Grass.Segmentation.Similarity, minsize = Grass.Segmentation.Minsize,
            memory = Grass.Segmentation.Memory, iterations = Grass.Segmentation.Iterations, bounds = paste0("Grass.Segmentation.Boundary.Grid", par.i)))
        } else
        {
          print("... with Boundary and Seeds")
          rgrass7::execGRASS("i.segment", flags = flag.Segmentation.GRASS, Sys_show.output.on.console = show.output.on.console, parameters = list(
            group = paste0("GRASS.Segmentation.Group", par.i), output = paste0("output.segmentation.GRASS", par.i), threshold = as.numeric(Grass.Segmentation.Threshold),
            method = Grass.Segmentation.Method, similarity = Grass.Segmentation.Similarity, minsize = Grass.Segmentation.Minsize,
            memory = Grass.Segmentation.Memory, iterations = Grass.Segmentation.Iterations, seeds = paste0("Grass.Segmentation.Seeds", par.i), bounds = paste0("Grass.Segmentation.Boundary.Grid", par.i)))
        }

      }

      if(is.null(Grass.Segmentation.Seeds) & is.null(Segmentation.Boundary.Grid))
      {
        rgrass7::execGRASS("i.segment", flags =flag.Segmentation.GRASS, Sys_show.output.on.console = show.output.on.console, parameters = list(
          group = paste0("GRASS.Segmentation.Group", par.i), output = paste0("output.segmentation.GRASS", par.i), threshold = as.numeric(Grass.Segmentation.Threshold),
          method = Grass.Segmentation.Method, similarity = Grass.Segmentation.Similarity, minsize = Grass.Segmentation.Minsize,
          memory = Grass.Segmentation.Memory, iterations = Grass.Segmentation.Iterations))

      }

    } # # # end tool grass segmentation


    #  Segmentation GRASS Superpixels SLIC ------------------------------------------------------------------------
    if(Tool == "GRASS Superpixels SLIC")
    {
      print("GRASS: Perform Superpixels SLIC")
      # execGRASS("i.superpixels.slic", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      #   group = "GRASS.Segmentation.Group", output = "output.segmentation.GRASS",
      #   iter = Grass.SLIC.Iter, k = Grass.SLIC.Superpixels, step = Grass.SLIC.Step, co = Grass.SLIC.Compactness, min = Grass.SLIC.Superpixels.MinSize))
      rgrass7::execGRASS("i.superpixels.slic", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        group = paste0("GRASS.Segmentation.Group", par.i), output = paste0("Grass.Segmentation.Seeds", par.i), memory = Grass.SLIC.Memory, perturb = Grass.SLIC.Perturb,
        iterations = Grass.SLIC.Iter, num_pixels = Grass.SLIC.Superpixels, step = Grass.SLIC.Step, compactness = Grass.SLIC.Compactness, minsize = Grass.SLIC.Superpixels.MinSize))
     } # # # end tool grass superpixels SLIC


    # print(parseGRASS("r.out.gdal"))
    if(tools::file_ext(Segments.Grid) == "sgrd")
    {
      rgrass7::execGRASS("r.out.gdal", flags = c("overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        input = paste0("output.segmentation.GRASS", par.i), output =  paste0(tools::file_path_sans_ext(Segments.Grid), ".sdat"), format = "SAGA"))
    } else {
      rgrass7::execGRASS("r.out.gdal", flags = c("overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        input = paste0("output.segmentation.GRASS", par.i), output =  Segments.Grid))
    }


    # remove segmentation items
    # print(parseGRASS("g.remove"))
    rgrass7::execGRASS("g.remove", flags = c("f", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      type = "group", name = paste0("GRASS.Segmentation.Group", par.i)))

    if(length(Input.Grid.GRASS) == 1)
    {

      rgrass7::execGRASS("g.remove", flags = c("f", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        type = "raster", name = paste0("input.GRASS", par.i)))
    } else
    {
      rgrass7::execGRASS("g.remove", flags = c("f", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        type = "raster", name = in.grids))
    }


    if(!is.null(Grass.Segmentation.Seeds))
    {
      rgrass7::execGRASS("g.remove", flags = c("f", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        type = "raster", name = paste0("Grass.Segmentation.Seeds", par.i)))
    }


  } # # # #  end of Segmentations of GRASS


  # # # # # check file extentsion of segmentation grid
  if(tools::file_ext(Segments.Grid) == "sgrd")
  {
    Segments.Grid.tmp <- raster::raster(paste0(tools::file_path_sans_ext(Segments.Grid), ".sdat"))
  } else {
    Segments.Grid.tmp <- raster::raster(Segments.Grid)
  }

    # browser()


  # Segments.Grid.tmp <- paste0(tools::file_path_sans_ext(Segments.Grid), ".sgrd")

  #  Burn Boundary into Segmentation ------------------------------------------------------------------------
  if(burn.Boundary.into.Segments[1] == TRUE)
  {
    print("Burn Boundary into Segments")

   # browser()

   if(tools::file_ext(Segmentation.Boundary.Grid) == "sgrd")
    {
     Segmentation.Boundary.Grid.tmp  <- raster::raster(paste0(tools::file_path_sans_ext(Segmentation.Boundary.Grid), ".sdat"))
    # Segmentation.Boundary.Grid.tmp <- Segmentation.Boundary.Grid
    } else {

      # Segmentation.Boundary.Grid.tmp <- paste0(tempdir(), "/", tools::file_path_sans_ext(basename(Segmentation.Boundary.Grid)), ".sdat")
#
#       RSAGA::rsaga.geoprocessor(lib = "io_gdal", module = 1, env = env, show.output.on.console = show.output.on.console, param = list(
#         GRIDS = Segmentation.Boundary.Grid, FILE =  paste0(tools::file_path_sans_ext(Segmentation.Boundary.Grid.tmp), ".sdat"), FORMAT =  "42"))

      Segmentation.Boundary.Grid.tmp <- raster::raster(Segmentation.Boundary.Grid)
      }

    # Segments.Grid.rgdal <- paste0(tools::file_path_sans_ext(Segments.Grid.tmp), ".sdat")
    # Segments.Grid.rgdal.grid <- rgdal::readGDAL(fname = Segments.Grid.rgdal, silent = TRUE)


    val <- Segments.Grid.tmp@data@max
    # val <- max(Segments.Grid.rgdal.grid$band1, na.rm = TRUE)

    # formula for calculation
    # formula.expression <- paste0("ifelse(a > 0, (a+", as.character(val, options(scipen = 999)), "), b)")

    # boundary grid should have 0 for No Data values, or values that are not burned
    # "ifelse(a > 0, a, b)"

    # Segmentation.Boundary.Grid.tmp <- paste0(tools::file_path_sans_ext(Segmentation.Boundary.Grid.tmp), ".sgrd")
    #
    # RSAGA::rsaga.geoprocessor(lib = "grid_calculus", env = env, module = 1, show.output.on.console = show.output.on.console, param=list(
    #   GRIDS =  paste(c(Segmentation.Boundary.Grid.tmp, Segments.Grid.tmp), collapse = ";"), RESULT = Segments.Grid.tmp, FORMULA = formula.expression,
    #   FNAME = "1", USE_NODATA = "1", TYPE = "5"))

    Segments.Grid.tmp <- raster::overlay(x = Segmentation.Boundary.Grid.tmp, y = Segments.Grid.tmp, fun = function(x, y){ifelse(x > 0, (x + val), y)}, forcefun = TRUE)

    rm(val, Segmentation.Boundary.Grid.tmp)
    # rm(formula.expression, val, Segments.Grid.rgdal.grid, Segments.Grid.rgdal)

  }

  #  Generalization ------------------------------------------------------------------------
  # perfrom generalization (majority filter)
  if(Generalisation.Flac == TRUE)
  {

    # browser()

    if(file_ext(Segments.Grid) != "sgrd")
    {
      Segments.Grid.r <- raster::raster(Segments.Grid)
      rgdal::writeGDAL(dataset = as(Segments.Grid.r, "SpatialGridDataFrame"),
                       fname = file.path(tempdir(), paste0(basename(tools::file_path_sans_ext(Segments.Grid)), ".sdat")),
                       drivername = "SAGA")


      Segments.Grid.tmp.path <- file.path(tempdir(), paste0(basename(tools::file_path_sans_ext(Segments.Grid)), ".sgrd"))

    } else {
      Segments.Grid.tmp.path <- Segments.Grid
    }

    # Segments.Grid.tmp.path <- paste0(tempdir(), "/", "Segments_Grid.tif")
    # raster::writeRaster(x = Segments.Grid.tmp, filename = Segments.Grid.tmp.path)
    #
    # RSAGA::rsaga.geoprocessor(lib = "io_gdal", module = 1, env = env, show.output.on.console = show.output.on.console, param = list(
    #          GRIDS = Segments.Grid.tmp.path, FILE =  paste0(tools::file_path_sans_ext(Segments.Grid.tmp.path), ".sdat"), FORMAT =  "42"))
    #
    # Segments.Grid.tmp.path  <- paste0(tools::file_path_sans_ext(Segments.Grid.tmp.path), ".sgrd")

    # rsaga.get.modules("grid_filter", env = env)
    # rsaga.get.usage("grid_filter", 6, env = env)
    # Mode: [0] Square | [1] Circle
    print("perform generalization")

    if(length(Generalization.Radius) == 1 && Generalization.Radius < 1)
    {
      print("Radius must be greater or equal to 1. Radius was set to 1.")
      Generalization.Radius = 1
    }

    if(length(Generalization.Radius) == 1) # 1 time use of majority filter
    {
      RSAGA::rsaga.geoprocessor(lib="grid_filter", module = 6, env = env, show.output.on.console = show.output.on.console, param = list(
        INPUT = Segments.Grid.tmp.path, RESULT = Segments.Grid.tmp.path, MODE = Generalization.Mode, RADIUS = Generalization.Radius,
        THRESHOLD = Generalization.Threshold))
    } else # use of majority filter
    {

      for(f in Generalization.Radius) # x time use of majority filter - first element is first filter procedure
      {
        RSAGA::rsaga.geoprocessor(lib="grid_filter", module = 6, env = env, show.output.on.console = show.output.on.console, param = list(
          INPUT = Segments.Grid.tmp.path, RESULT = Segments.Grid.tmp.path, MODE = Generalization.Mode, RADIUS = f,
          THRESHOLD = Generalization.Threshold))
      }

    }

    Segments.Grid.tmp <- raster::raster(paste0(tools::file_path_sans_ext(Segments.Grid.tmp.path), ".sdat"))

  }



  # Mode Filter -------------------------------------
    if(Mode.Filter.Flac == TRUE)
  {

    # browser()

    print("perform Mode Filter")
    Segments.Grid.tmp.path <- file.path(tempdir(), paste0("SegGridTmp", par.i, ".tif"))
    Segments.Grid.Mode.tmp.path  <- file.path(tempdir(), paste0("SegGridMode", par.i, ".tif"))

    # re-write segmentation raster
    rgrass7::writeRAST(x = as(Segments.Grid.tmp, "SpatialGridDataFrame"), vname = paste0("Segments.Grid.tmp", par.i), zcol = names(Segments.Grid.tmp), overwrite = TRUE,
                       flags = "quiet")

    # print(parseGRASS("r.to.vect"))
    # par.v
    #  "quiet"
    rgrass7::execGRASS("r.to.vect", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = paste0("Segments.Grid.tmp", par.i), output = paste0("Segments_Poly_Mode", par.i), type = "area"))


    # print(parseGRASS("v.to.rast"))
    # par.v
    #  "quiet"
    rgrass7::execGRASS("v.to.rast", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = paste0("Segments_Poly_Mode", par.i), output = paste0("Segments.Grid.tmp", par.i), type = "area", use = "cat"))


    # compute mode
    # print(parseGRASS("r.neighbors"))
    rgrass7::execGRASS("r.neighbors", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input =  paste0("Segments.Grid.tmp", par.i), output = paste0("Segments.Grid.Mode", par.i), method = "mode", size = Mode.Filter.Size))


    # write out data
    # print(parseGRASS("r.out.gdal"))
    rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = paste0("Segments.Grid.tmp", par.i), output = Segments.Grid.tmp.path))

    rgrass7::execGRASS("r.out.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input =  paste0("Segments.Grid.Mode", par.i), output = Segments.Grid.Mode.tmp.path))


    # read data
    Segments.Grid.tmp <- raster::raster(Segments.Grid.tmp.path)
    Segments.Grid.tmp.mode <- raster::raster(Segments.Grid.Mode.tmp.path)

    # browser()
    # calculate frequency table and select based on minsize
    Segments.Grid.tmp.freq <- as.data.table(raster::freq(x = Segments.Grid.tmp))
    Segments.Grid.tmp.remove <- Segments.Grid.tmp.freq[count <= Mode.Filter.Segment.MinSize,]$value

    # browser()
    if(length(Segments.Grid.tmp.remove) > 0)
    {
      # create and evaluate expression based on minsize
      # exp <- parse(text=paste0("x == ", Segments.Grid.tmp.remove, collapse = " | "))

      Segments.Grid.tmp.remove.r  <- raster::match(Segments.Grid.tmp, Segments.Grid.tmp.remove, nomatch = 0)
      # Segments.Grid.tmp.remove.r <- raster::calc(x = Segments.Grid.tmp, fun = function(x){return(ifelse(x > 0, 1, 0))}, forcefun = TRUE)


      # Segments.Grid.tmp.remove.r <- raster::calc(x = Segments.Grid.tmp, fun = function(x){return(ifelse(eval(exp), 1, 0))}, forcefun = TRUE)

      Segments.Grid.tmp <- raster::overlay(x = Segments.Grid.tmp, y = Segments.Grid.tmp.mode, z = Segments.Grid.tmp.remove.r,
                                fun =  function(x, y, z){return(ifelse(z > 0, y, x))})


      rm(Segments.Grid.tmp.mode, Segments.Grid.tmp.remove.r, Segments.Grid.tmp.freq)
    }

  }




  # mask no data area
  if(NoData == TRUE)
  {

    if(tools::file_ext(Mask) == "sgrd")
    {
      # Mask.tmp <- Mask
      Mask.tmp <- raster::raster(paste0(tools::file_path_sans_ext(Mask), ".sdat"))

    } else {
      # Mask.tmp <- paste0(tempdir(), "/", tools::file_path_sans_ext(basename(Mask)), ".sgrd")
      #
      # RSAGA::rsaga.geoprocessor(lib = "io_gdal", module = 1, env = env, show.output.on.console = show.output.on.console, param = list(
      #   GRIDS = Mask, FILE = paste0(tools::file_path_sans_ext(Mask.tmp), ".sdat"), FORMAT =  "42"))
      Mask.tmp <- raster::raster(Mask)
    }


    # rsaga.get.modules("grid_tools", env = env)
    # rsaga.get.usage("grid_tools", 24, env = env)
    # print("Masking No Data Area in Segments")
    # RSAGA::rsaga.geoprocessor(lib="grid_tools", module = 24, env = env, show.output.on.console = show.output.on.console, param = list(
    #   GRID = Segments.Grid.tmp, MASK = Mask.tmp, MASKED = Segments.Grid.tmp))
    Segments.Grid.tmp <- raster::overlay(Segments.Grid.tmp, Mask.tmp, fun = function(x, y){return(ifelse(is.na(y), NA, x))})
  }

  # browser()

  # write out data again
  if(tools::file_ext(Segments.Grid) == "sgrd")
  {

    Segments.Grid.tmp.path <- paste0(tempdir(), "/", "Segments_Grid.tif")
    raster::writeRaster(x = Segments.Grid.tmp, filename = Segments.Grid.tmp.path, overwrite=TRUE, NAflag = NoData.Flag)

    RSAGAUsage <- RSAGA::rsaga.get.usage(lib="io_gdal", module = 1, env = env,  show = FALSE)
    formatSAGA <- gsub("\\D", "", grep('SAGA GIS Binary', RSAGAUsage, value = TRUE))


    RSAGA::rsaga.geoprocessor(lib = "io_gdal", module = 1, env = env, show.output.on.console = show.output.on.console, param = list(
      GRIDS = Segments.Grid.tmp.path, FILE =  paste0(tools::file_path_sans_ext(Segments.Grid), ".sdat"), FORMAT =   formatSAGA))

  } else {
    # browser()
    if(file.exists(Segments.Grid) && !identical(Segments.Grid.tmp,  raster::raster(Segments.Grid)))
    {
         raster::writeRaster(x = Segments.Grid.tmp, filename = Segments.Grid, overwrite = TRUE, NAflag = NoData.Flag)
    }
  }

  #  browser()

  #  Vectorising Raster Objects ------------------------------------------------------------------------
  # vectorising grid classes
  # rsaga.get.modules("shapes_grid", env = env)
  # rsaga.get.usage("shapes_grid", 6, env = env)
  # CLASS_ALL = [1] all classes
  # ALLVERTICES: Keep Vertices on Straight Lines
  print("vectorising grid classes")
  # RSAGA::rsaga.geoprocessor(lib="shapes_grid", module = 6, env = env, show.output.on.console = show.output.on.console, param = list(
  #   GRID = Segments.Grid.tmp, POLYGONS = Segments.Poly, CLASS_ALL = "1", SPLIT = Split, ALLVERTICES = AllVertices))
  rgrass7::writeRAST(x = as(Segments.Grid.tmp, "SpatialGridDataFrame"), vname = paste0("Segments.Grid.tmp", par.i), zcol = names(Segments.Grid.tmp), overwrite = TRUE,
                     flags = "quiet")


  # print(parseGRASS("r.to.vect"))
  # par.v
  #  "quiet"
  rgrass7::execGRASS("r.to.vect", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
    input = paste0("Segments.Grid.tmp", par.i), output = paste0("Segments_Poly", par.i), type = "area"))

  if(estimateScaleParameter == FALSE)
  {
    # print(parseGRASS("v.out.ogr"))
    # par.v
    #  "quiet"
    rgrass7::execGRASS("v.out.ogr", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = paste0("Segments_Poly", par.i), output = Segments.Poly, type = "area", format = "ESRI_Shapefile"))

  }




  # clear temp data
  print("clear temp data")
  invisible(do.call(file.remove,list(list.files("Output/Temp", full.names=TRUE))))

  # get time of process
  process.time.run <- proc.time() - process.time.start
  print(paste0("------ Run of Segmentation: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))

}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
