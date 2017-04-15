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
#' @param Grass.Segmentation.Memor memory to be used for segmentation. Default: 300
#' @param Grass.Segmentation.Iterations amount of allowed iterations. Default: 50
#' @param Grass.Segmentation.Seeds input of seeds raster. Enables bottom-up segmentation. Default: NULL
#' @param Segmentation.Boundary.Grid input of boundary raster. Enables top-down (or hierarchical) segmentation. Default: NULL
#' @param Grass.Segmentation.Goodness name for output goodness of fit estimate map. Default:"Grass.Segmentation.Goodness"
#' @param Grass.SLIC.Iter maximum number of iterations. Default: 10
#' @param Grass.SLIC.Superpixels approximate number of output super pixels. Default: 200
#' @param Grass.SLIC.Step distance (number of cells) between initial super pixel centers. A step size > 0 overrides the number of super pixels. Default: 0
#' @param Grass.SLIC.Compactness compactness. A larger value causes more compact superpixels. Default: 1.0
#' @param Grass.SLIC.Superpixels.MinSize minimum superpixel size. Default: 1
#' @param Grass.SLIC.Memory memory in MB. Default: 300
#' @param Grass.SLIC.Perturb Perturb initial super pixel centers. Percent of intitial superpixel radius. Default: 0, range: 0-100
#' @param burn.Boundary.into.Segments vector specifing if boundary grid is burned into segmentation (1) or seeds (2). Default: FALSE, maximum length: 2
#'
#'
#' @keywords segmentation, region growing, superpixels Simple Linear Iterative Clustering
#'
#' @examples
#' segmentation()
#'
#' @export
segmentation <- function(Tool, Segments.Grid, Segments.Poly, Input.Grid, Saga.Output.Grid = paste0( segmentation.tmp.path, "SagaRepresentativenessOutputGrid.sgrd"), Saga.Output.Lod =  paste0( segmentation.tmp.path, "SagaRepresentativenessOutputLod.sgrd"),
                         Output.Seeds =  paste0(segmentation.tmp.path, "OutputSeed.sgrd"), Fast.Representativeness.LevelOfGeneralisation = "10.0", Saga.Similarity =  paste0( segmentation.tmp.path,"SagaSimilarity.sgrd"), Saga.Segments.Seeds.Table =  paste0( segmentation.tmp.path, "SagaSegmentsSeedsTable.mtab"),  Saga.Segmentation.Normalize = "0", Saga.Segmentation.Neighbourhood = "0", Saga.Segmentation.Method = "0",
                         Saga.Segmentation.Sig.1 = "1.0", Saga.Segmentation.Sig.2 = "1.0", Saga.Segmentation.Threshold = "0.0", Saga.Segmentation.Refresh = "0", Saga.Segmentation.Leafsize = 256, Split = "0",
                         Grass.Segmentation.Threshold = NULL, Grass.Segmentation.Weighted = FALSE, Grass.Segmentation.Method = "region_growing", Grass.Segmentation.Similarity = "euclidean", Grass.Segmentation.Minsize = 15, Grass.Segmentation.Memory = 300, Grass.Segmentation.Iterations = 50, Grass.Segmentation.Seeds = NULL, Grass.Segmentation.Neighbourhood = "0",
                         Segmentation.Boundary.Grid = NULL,  Grass.Segmentation.Goodness = "Grass.Segmentation.Goodness", AllVertices = "FALSE", NoData = FALSE, Mask = NULL, show.output.on.console = FALSE, Seed.Method = "", Seed.Generation.Variance =  paste0( segmentation.tmp.path, "SeedGenerationVariance.sgrd"), Seed.Generation.Points =  paste0(segmentation.tmp.path, "SeedGenerationPoints.shp"),
                         Seed.Generation.Type = "0", Seed.Generation.Scale = "10.0", Generalisation.Flac = FALSE, Generalization.Mode = "1", Generalization.Radius = "1", Generalization.Threshold = "0.0",
                         Grass.SLIC.Iter = 10, Grass.SLIC.Superpixels = 200, Grass.SLIC.Step = 0, Grass.SLIC.Compactness = 1.0, Grass.SLIC.Superpixels.MinSize = 1, Grass.SLIC.Memory = 300, Grass.SLIC.Perturb = 0, burn.Boundary.into.Segments = FALSE, ...)
{


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
      rsaga.geoprocessor(lib="statistics_grid", module = 0, env = env, show.output.on.console = show.output.on.console, param = list(
        INPUT = Input.Grid.Saga, RESULT = Saga.Output.Grid, RESULT_LOD = Saga.Output.Lod, SEEDS = Output.Seeds,
        LOD = Fast.Representativeness.LevelOfGeneralisation))
    }

    if(length(Input.Grid.Strings[[1]]) >= 2) # ... if yes, take first one for representativeness
    {
      print("... first grid of multiple is taken for seed points")
      Input.Grid1 <- Input.Grid.Strings[[1]][1]
      rsaga.geoprocessor(lib="statistics_grid", module = 0, env = env, show.output.on.console = show.output.on.console, param = list(
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

    rsaga.geoprocessor(lib="imagery_segmentation", module = 2, env = env, show.output.on.console = show.output.on.console, param = list(
      FEATURES = Input.Grid.Saga, VARIANCE = Seed.Generation.Variance, SEED_GRID =  Output.Seeds, SEED_POINTS = Seed.Generation.Points,
      SEED_TYPE = Seed.Generation.Type, METHOD = "0", NORMALIZE = "0", BAND_WIDTH = Seed.Generation.Scale))

  } else if(Seed.Method == "Superpixels SLIC")
  {
    print("GRASS Superpixels SLIC as Seed Method is selected")

  } else{

    print("No Seed Method was used!")
  }


  #  Segmentation ------------------------------------------------------------------------
  if(Tool != "SAGA" & Tool != "GRASS" & Tool != "GRASS Superpixels SLIC")
  {
    stop("Wrong tool software was selected. Please use 'SAGA' or 'GRASS' for segmentation procedure.")
  }

  #  Segmentation SAGA ------------------------------------------------------------------------
  if(Tool == "SAGA")
  {
    # perform SAGA seeded region growing
    # rsaga.get.modules("imagery_segmentation", env = env)
    # rsaga.get.usage("imagery_segmentation", 3, env = env)
    print("SAGA: Perform seeded region growing")

    rsaga.geoprocessor(lib="imagery_segmentation", module = 3, env = env, show.output.on.console = show.output.on.console, param = list(
      SEEDS = Output.Seeds, FEATURES= Input.Grid.Saga, SEGMENTS = Segments.Grid, SIMILARITY = Saga.Similarity, TABLE = Saga.Segments.Seeds.Table,
      NORMALIZE = Saga.Segmentation.Normalize, NEIGHBOUR = Saga.Segmentation.Neighbourhood, METHOD = Saga.Segmentation.Method, SIG_1 = Saga.Segmentation.Sig.1, SIG_2 = Saga.Segmentation.Sig.2,
      THRESHOLD = Saga.Segmentation.Threshold, REFRESH = Saga.Segmentation.Refresh, LEAFSIZE = Saga.Segmentation.Leafsize))
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
      # Grass.Segmentation.Seeds <- paste0(file_path_sans_ext(Output.Seeds) , ".sdat")
      # execGRASS("r.in.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      #   input = Grass.Segmentation.Seeds, output = "Grass.Segmentation.Seeds"))

      #
      # perform SAGA seeded region growing
      # rsaga.get.modules("imagery_segmentation", env = env)
      # rsaga.get.usage("imagery_segmentation", 3, env = env)
      print("SAGA: Perform seeded region growing for GRASS Seeds")
      rsaga.geoprocessor(lib="imagery_segmentation", module = 3, env = env, show.output.on.console = show.output.on.console, param = list(
        SEEDS = Output.Seeds, FEATURES= Input.Grid.Saga, SEGMENTS = Segments.Grid, SIMILARITY = Saga.Similarity, TABLE = Saga.Segments.Seeds.Table,
        NORMALIZE = Saga.Segmentation.Normalize, NEIGHBOUR = Saga.Segmentation.Neighbourhood, METHOD = Saga.Segmentation.Method, SIG_1 = Saga.Segmentation.Sig.1, SIG_2 = Saga.Segmentation.Sig.2,
        THRESHOLD = Saga.Segmentation.Threshold, REFRESH = Saga.Segmentation.Refresh, LEAFSIZE = Saga.Segmentation.Leafsize))


      #  Burn Boundary into Seed ------------------------------------------------------------------------
      if(length(burn.Boundary.into.Segments) == 2 && burn.Boundary.into.Segments[2] == TRUE)
      {
        print("Burn Boundary into GRASS Seeds")
        Segments.Grid.rgdal <- paste0(file_path_sans_ext(Segments.Grid), ".sdat")
        Segments.Grid.rgdal.grid <- rgdal::readGDAL(fname = Segments.Grid.rgdal, silent = TRUE)
        val <- max(Segments.Grid.rgdal.grid$band1, na.rm = TRUE)

        # formula for calculation
        formula.expression <- paste0("ifelse(a > 0, (a+", val, "), b)")

        # boundary grid should have 0 for No Data values, or values that are not burned
        # "ifelse(a > 0, a, b)"
        rsaga.geoprocessor(lib = "grid_calculus", env = env, module = 1, show.output.on.console = show.output.on.console, param=list(
          GRIDS =  paste(c(Segmentation.Boundary.Grid, Segments.Grid), collapse = ";"), RESULT = Segments.Grid, FORMULA = formula.expression,
          FNAME = "1", USE_NODATA = "1", TYPE = "5"))

        rm(formula.expression, val, Segments.Grid.rgdal.grid, Segments.Grid.rgdal)
      }



      # change data storage of reclassified seed output: from float to integer for GRASS
      # rsaga.get.modules("grid_tools", env = env)
      # rsaga.get.usage("grid_tools", 11, env = env)
      # TYPE: [1] unsigned 1 byte integer, [5] unsigned 4 byte integer
      rsaga.geoprocessor(lib="grid_tools", module = 11, env = env, show.output.on.console = show.output.on.console, param = list(
        INPUT = Segments.Grid, OUTPUT = Segments.Grid, TYPE = "5", OFFSET = "0.0", SCALE = "1.0"))

      Grass.Segmentation.Seeds <- paste0(file_path_sans_ext(Segments.Grid) , ".sdat")
      execGRASS("r.in.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        input = Grass.Segmentation.Seeds, output = "Grass.Segmentation.Seeds"))

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


      # print(parseGRASS("r.in.gdal"))
      inputGRASS <- paste0(file_path_sans_ext(Input.Grid.GRASS) , ".sdat")
      execGRASS("r.in.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        input = inputGRASS, output = "input.GRASS"))

      # print(parseGRASS("i.group"))
      execGRASS("i.group", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        group = "GRASS.Segmentation.Group", input = "input.GRASS")) # if this already exists it will be skipped

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
        inputGRASS <- paste0(file_path_sans_ext(Input.Grid.GRASS[j]) , ".sdat")
        in.grids[j] <- paste0("input.GRASS.", j)
        execGRASS("r.in.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
          input = inputGRASS, output = in.grids[j]))
      }

      # print(parseGRASS("i.group"))
      execGRASS("i.group", flags = c("quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        group = "GRASS.Segmentation.Group", input = in.grids)) # if this already exists it will be skipped
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
      execGRASS("i.superpixels.slic", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        group = "GRASS.Segmentation.Group", output = "Grass.Segmentation.Seeds", memory = Grass.SLIC.Memory, perturb = Grass.SLIC.Perturb,
        iterations = Grass.SLIC.Iter, num_pixels = Grass.SLIC.Superpixels, step = Grass.SLIC.Step, compactness = Grass.SLIC.Compactness, minsize = Grass.SLIC.Superpixels.MinSize))


      Grass.Segmentation.Seeds <- "Grass.Segmentation.Seeds"

    }

    #  Segmentation GRASS ------------------------------------------------------------------------
    if(Tool == "GRASS")
    {
      print("GRASS: Perform Segmentation")
      # print(parseGRASS("i.segment"))
      if(!is.null(Grass.Segmentation.Seeds) & is.null(Segmentation.Boundary.Grid))
      {
        print("... with Seeds")
        execGRASS("i.segment", flags = flag.Segmentation.GRASS, Sys_show.output.on.console = show.output.on.console, parameters = list(
          group = "GRASS.Segmentation.Group", output = "output.segmentation.GRASS", threshold = as.numeric(Grass.Segmentation.Threshold),
          method = Grass.Segmentation.Method, similarity = Grass.Segmentation.Similarity, minsize = Grass.Segmentation.Minsize,
          memory = Grass.Segmentation.Memory, iterations = Grass.Segmentation.Iterations, seeds = "Grass.Segmentation.Seeds"))
      }

      if(!is.null(Segmentation.Boundary.Grid))
      {
        Grass.Segmentation.Boundary.Grid.GRASS <- paste0(file_path_sans_ext(Segmentation.Boundary.Grid) , ".sdat")
        execGRASS("r.in.gdal", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
          input = Grass.Segmentation.Boundary.Grid.GRASS, output = "Grass.Segmentation.Boundary.Grid"))

        if(is.null(Grass.Segmentation.Seeds))
        {
          print("... with Boundary")
          execGRASS("i.segment", flags = flag.Segmentation.GRASS, Sys_show.output.on.console = show.output.on.console, parameters = list(
            group = "GRASS.Segmentation.Group", output = "output.segmentation.GRASS", threshold = as.numeric(Grass.Segmentation.Threshold),
            method = Grass.Segmentation.Method, similarity = Grass.Segmentation.Similarity, minsize = Grass.Segmentation.Minsize,
            memory = Grass.Segmentation.Memory, iterations = Grass.Segmentation.Iterations, bounds = "Grass.Segmentation.Boundary.Grid"))
        } else
        {
          print("... with Boundary and Seeds")
          execGRASS("i.segment", flags = flag.Segmentation.GRASS, Sys_show.output.on.console = show.output.on.console, parameters = list(
            group = "GRASS.Segmentation.Group", output = "output.segmentation.GRASS", threshold = as.numeric(Grass.Segmentation.Threshold),
            method = Grass.Segmentation.Method, similarity = Grass.Segmentation.Similarity, minsize = Grass.Segmentation.Minsize,
            memory = Grass.Segmentation.Memory, iterations = Grass.Segmentation.Iterations, seeds = "Grass.Segmentation.Seeds", bounds = "Grass.Segmentation.Boundary.Grid"))
        }

      }

      if(is.null(Grass.Segmentation.Seeds) & is.null(Segmentation.Boundary.Grid))
      {
        execGRASS("i.segment", flags =flag.Segmentation.GRASS, Sys_show.output.on.console = show.output.on.console, parameters = list(
          group = "GRASS.Segmentation.Group", output = "output.segmentation.GRASS", threshold = as.numeric(Grass.Segmentation.Threshold),
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
      execGRASS("i.superpixels.slic", flags = c("overwrite", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        group = "GRASS.Segmentation.Group", output = "Grass.Segmentation.Seeds", memory = Grass.SLIC.Memory, perturb = Grass.SLIC.Perturb,
        iterations = Grass.SLIC.Iter, num_pixels = Grass.SLIC.Superpixels, step = Grass.SLIC.Step, compactness = Grass.SLIC.Compactness, minsize = Grass.SLIC.Superpixels.MinSize))
     } # # # end tool grass superpixels SLIC


    # print(parseGRASS("r.out.gdal"))
    execGRASS("r.out.gdal", flags = c("overwrite"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      input = "output.segmentation.GRASS", output =  paste0(file_path_sans_ext(Segments.Grid), ".sdat"), format = "SAGA"))


    # remove segmentation items
    # print(parseGRASS("g.remove"))
    execGRASS("g.remove", flags = c("f", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
      type = "group", name = "GRASS.Segmentation.Group"))

    if(length(Input.Grid.GRASS) == 1)
    {

      execGRASS("g.remove", flags = c("f", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        type = "raster", name = "input.GRASS"))
    } else
    {
      execGRASS("g.remove", flags = c("f", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        type = "raster", name = in.grids))
    }


    if(!is.null(Grass.Segmentation.Seeds))
    {
      execGRASS("g.remove", flags = c("f", "quiet"), Sys_show.output.on.console = show.output.on.console, parameters = list(
        type = "raster", name = "Grass.Segmentation.Seeds"))
    }


  } # # # #  end of Segmentations of GRASS


  #  Burn Boundary into Segmentation ------------------------------------------------------------------------
  if(burn.Boundary.into.Segments[1] == TRUE)
  {
    print("Burn Boundary into Segments")

    Segments.Grid.rgdal <- paste0(file_path_sans_ext(Segments.Grid), ".sdat")
    Segments.Grid.rgdal.grid <- rgdal::readGDAL(fname = Segments.Grid.rgdal, silent = TRUE)
    val <- max(Segments.Grid.rgdal.grid$band1, na.rm = TRUE)

    # formula for calculation
    formula.expression <- paste0("ifelse(a > 0, (a+", as.character(val, options(scipen = 999)), "), b)")

    # boundary grid should have 0 for No Data values, or values that are not burned
    # "ifelse(a > 0, a, b)"
    rsaga.geoprocessor(lib = "grid_calculus", env = env, module = 1, show.output.on.console = show.output.on.console, param=list(
      GRIDS =  paste(c(Segmentation.Boundary.Grid, Segments.Grid), collapse = ";"), RESULT = Segments.Grid, FORMULA = formula.expression,
      FNAME = "1", USE_NODATA = "1", TYPE = "5"))

    rm(formula.expression, val, Segments.Grid.rgdal.grid, Segments.Grid.rgdal)

  }

  #  Generalization ------------------------------------------------------------------------
  # perfrom generalization (majority filter)
  if(Generalisation.Flac == TRUE)
  {
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
      rsaga.geoprocessor(lib="grid_filter", module = 6, env = env, show.output.on.console = show.output.on.console, param = list(
        INPUT = Segments.Grid, RESULT = Segments.Grid, MODE = Generalization.Mode, RADIUS = Generalization.Radius,
        THRESHOLD = Generalization.Threshold))
    } else # use of majority filter
    {

      for(f in Generalization.Radius) # x time use of majority filter - first element is first filter procedure
      {
        rsaga.geoprocessor(lib="grid_filter", module = 6, env = env, show.output.on.console = show.output.on.console, param = list(
          INPUT = Segments.Grid, RESULT = Segments.Grid, MODE = Generalization.Mode, RADIUS = f,
          THRESHOLD = Generalization.Threshold))
      }

    }

  }


  # mask no data area
  if(NoData == TRUE)
  {
    # rsaga.get.modules("grid_tools", env = env)
    # rsaga.get.usage("grid_tools", 24, env = env)
    print("Masking No Data Area in Segments")
    rsaga.geoprocessor(lib="grid_tools", module = 24, env = env, show.output.on.console = show.output.on.console, param = list(
      GRID = Segments.Grid, MASK = Mask, MASKED = Segments.Grid))

  }

  #  Vectorising Raster Objects ------------------------------------------------------------------------
  # vectorising grid classes
  # rsaga.get.modules("shapes_grid", env = env)
  # rsaga.get.usage("shapes_grid", 6, env = env)
  # CLASS_ALL = [1] all classes
  # ALLVERTICES: Keep Vertices on Straight Lines
  print("vectorising grid classes")
  rsaga.geoprocessor(lib="shapes_grid", module = 6, env = env, show.output.on.console = show.output.on.console, param = list(
    GRID = Segments.Grid, POLYGONS = Segments.Poly, CLASS_ALL = "1", SPLIT = Split, ALLVERTICES = AllVertices))


  # clear temp data
  print("clear temp data")
  invisible(do.call(file.remove,list(list.files("Output/Temp", full.names=TRUE))))

  # get time of process
  process.time.run <- proc.time() - process.time.start
  print(paste0("------ Run of Segmentation: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))

}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
