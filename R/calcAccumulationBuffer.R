#' Buffering of flow accumulations using thresholds
#'
#' This function buffers a flow accumulation input in dependence on given thesholds and distances
#' @param in.grid path of input grid (flow accumulation)
#' @param out.grid path of output grid
#' @param out.shp creation of an output shapefile. Default: NULL
#' @param return.sf return of simple feature. Default: FALSE
#' @param t vector containing thresholds for flow accumulation grid
#' @param buf.dist vector containing the buffer distances
#' @param working.path working path. Default: NULL
#' @param keep.intermediate keep intermediate data. Default: FALSE
#' @param noData value for no data input. Default: -99999
#' @param show.output.on.console show output on console. Default: FALSE
#' @param quiet no outputs in console. Default: TRUE
#'
#' @return
#' if \emph{out.shp} is set and \emph{return.sf} is set to TRUE, then a simple feature object is returned
#'
#' @note
#' length of \emph{t} and \emph{buf.dist} must be identical
#'
#' @keywords flow accumulation, buffer
#'
#'
#' @export
calcAccumulationBuffer <- function(in.grid, out.grid, out.shp = NULL, return.sf = FALSE, t, buf.dist, working.path = NULL, keep.intermediate = FALSE, noData = -99999, show.output.on.console = FALSE, quiet = TRUE)
{

  # get start time of process
  process.time.start <- proc.time()

  # check input
  if(identical(length(t), length(buf.dist)) == FALSE)
  {
    stop("Length of buffer and thresholds differ!")
  }

  # check for relative path
  if(!is.null(working.path))
  {
    working.path <- paste0(getwd(), "/")
  } else
  {
    working.path <- ""
  }


  # get default naming
  default_name <- basename(file_path_sans_ext(in.grid))

  # get name
  out.grid.i <- paste0(dirname(out.grid), "/", default_name, (t/1000), "k.sgrd")

  # start calculation of buffer
  if(quiet == FALSE) print("Calculation of Buffer")
  for(i in 1: length(t))
  {
    if(quiet == FALSE) print(paste0("... ", i, " of ", length(t)))

    # get formula
    formula.i <- paste0("ifelse(gt(a,", as.character(t[i], options(scipen = 999)), "), 1, (", noData, "))")

    # get cells higher then the threshold
    # rsaga.get.usage("grid_calculus", 1, env = env)
    rsaga.geoprocessor(lib = "grid_calculus", env = env, module = 1, show.output.on.console = show.output.on.console, param=list(
      GRIDS = in.grid, RESULT = out.grid.i[i], FORMULA = formula.i, FNAME = "1"))

    # start buffer
    rsaga.geoprocessor(lib="grid_tools", module = 8, env = env, show.output.on.console = show.output.on.console, param = list(
      FEATURES = out.grid.i[i], BUFFER = out.grid.i[i], DIST = buf.dist[i], BUFFERTYPE = "0"))

    # reclass values
    rsaga.geoprocessor(lib="grid_tools", module = 15, env = env, show.output.on.console = FALSE, param = list(
      INPUT = out.grid.i[i], RESULT = out.grid.i[i], METHOD = "0", OLD = "0.0", NEW = "1.0",
      SOPERATOR = "4", NODATAOPT = "1", NODATA = "0"))

    # start calculation of final buffer inside of the loop
    if(i == 1)
    {
     outgrid.tmp <- paste0(dirname(out.grid.i[i]), "/", default_name,"_temp.sgrd")

     rsaga.geoprocessor(lib = "grid_calculus", env = env, module = 1, show.output.on.console = show.output.on.console, param=list(
       GRIDS = out.grid.i[i], RESULT =  outgrid.tmp , FORMULA = "(a+0)", FNAME = "1"))

    } else
    {
      # adding buffer togehter

      if(i == length(t)) # last run
      {
        rsaga.geoprocessor(lib = "grid_calculus", env = env, module = 1, show.output.on.console = show.output.on.console, param=list(
          GRIDS = paste(c(outgrid.tmp, out.grid.i[i]),  collapse = ";"), RESULT = out.grid, FORMULA = "(a+b)", FNAME = "1"))


        # reset values of last run
        formula.last <- paste0("ifelse(gt(a,", 0, "), 1, (", noData, "))")

        rsaga.geoprocessor(lib = "grid_calculus", env = env, module = 1, show.output.on.console = show.output.on.console, param=list(
          GRIDS = out.grid, RESULT = out.grid, FORMULA =  formula.last, FNAME = "1"))


        if(keep.intermediate == FALSE)
        {
          files.2remove <- list.files(path = paste0(working.path, dirname(out.grid.i[i])), pattern = basename(file_path_sans_ext(out.grid.i[i])))
          files.2remove <- paste0(paste0(working.path, dirname(out.grid.i[i])), "/", files.2remove)

          files.2remove.1 <- list.files(path = paste0(working.path, dirname(out.grid.i[1])), pattern = basename(file_path_sans_ext(out.grid.i[1])))
          files.2remove.1 <- paste0(paste0(working.path, dirname(out.grid.i[1])), "/", files.2remove)

          # file.remove(files.2remove)
          unlink(c(files.2remove, files.2remove.1), force = TRUE)
        }

        # remove temp file
        files.2remove <- list.files(path = paste0(working.path, dirname(outgrid.tmp)), pattern = basename(file_path_sans_ext(outgrid.tmp)))
        files.2remove <- paste0(paste0(working.path, dirname(outgrid.tmp)), "/", files.2remove)

        # file.remove(files.2remove)
        unlink(files.2remove, force = TRUE, recursive = TRUE)

      } else
      {
        rsaga.geoprocessor(lib = "grid_calculus", env = env, module = 1, show.output.on.console = show.output.on.console, param=list(
          GRIDS = paste(c(outgrid.tmp, out.grid.i[i]),  collapse = ";"), RESULT = outgrid.tmp, FORMULA = "(a+b)", FNAME = "1"))

        # delete data if set
        if(keep.intermediate == FALSE)
        {
          files.2remove <- list.files(path = paste0(working.path, dirname(out.grid.i[i])), pattern = basename(file_path_sans_ext(out.grid.i[i])))
          files.2remove <- paste0(paste0(working.path, dirname(out.grid.i[i])), "/", files.2remove)

          # file.remove(files.2remove)
          unlink(files.2remove, force = TRUE, recursive = TRUE)
        }

      } # end of else of if(i == length(t))
    } # end of else of if(i == 1)


  } # end of for loop

  if(!is.null(out.shp))
  {
    if(quiet == FALSE) print("Convert grid to polygon")
    rsaga.geoprocessor(lib="shapes_grid", module = 6, env = env, show.output.on.console = FALSE, param = list(
      GRID =  out.grid, POLYGONS = out.shp))

    correctDBF(x = out.shp, end.n = 0, new.colnames = c("Stream", "ID", "NAME"))
  }

  if(return.sf == TRUE & !is.null(out.shp))
  {
    out.sf <- sf::st_read(dsn = paste(getwd(), dirname(out.shp), sep = "/"), layer = file_path_sans_ext(basename(out.shp)), quiet = TRUE, stringsAsFactors = FALSE)

    # get time of process
    process.time.run <- proc.time() - process.time.start
    if(quiet == FALSE) print(paste0("------ Run of calcAccumulationBuffer: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))

    return(out.sf)
  } else
  {
    # get time of process
    process.time.run <- proc.time() - process.time.start
    if(quiet == FALSE) print(paste0("------ Run of calcAccumulationBuffer: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))

  }
}

