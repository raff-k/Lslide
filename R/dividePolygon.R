#' Devide Polygon into N Equal Area Parts
#'
#' This function divides a polygon feature into n-party of approximatly the same area.
#' The division is performed along the y-axis (N-S). Therefore, polygons should be rotated
#' appropriate.
#'
#'
#' @param x object of class sf
#' @param angle Can be one single value, or a vector with the same length as the input object x. North is 0 degree, and East 90. Default: 0 degree for North-South. For East-West use 90 degree.
#' @param n number of parts
#' @param trans angle input is transformed, so that N-S cutting is appropriate. Default: TRUE.
#' @param btol tolerance value for lower and upper boundary in optim().Default: 0.0001
#' @param pgtol tolerance value to reduce the sum of residuals in optimizing function
#' @param maxit maximum iteration to convergence in finding optimal points. Default: 100
#' @param cores number of cores for parallel processing. Default: 1 (sequential)
#' @param quiet no outputs in console. Default: TRUE
#' @return
#' Geometry similar to sf input including divisions.
#'
#'
#' @note
#' \itemize{
#'   \item Algorithm does not find optimum if polygon split results in more objects than the given number
#'          n. This can happen when a polygon is very ragged.
#'   \item Algorithm cut polygon on the N-S-Axis (0-180 degree). Therefore, trans is used to adapt the rotation angle
#'         so that cutting is appropriate.
#'
#' }
#'
#'
#' @keywords rainfall tresholds, rainfall event, landslide, automatic appraoch
#'
#'
#' @export
#'
dividePolygon <- function(x, angle = 0, n, trans = TRUE, btol = 0.0001, pgtol = 1e-9, maxit = 100, cores = 1, quiet = TRUE)
{

  # get start time of process
  process.time.start <- proc.time()

  # # # # # # # #
  ## check input
  # # # # # # # #
  n.obj <- nrow(x) # number of features

  # check angle input
  if(length(angle) == 1)
  {
    angle <- rep(x = angle, n.obj) # create angle vector
  } else{
    if(length(angle) != n.obj){stop('Length of "angle" is not equal to "x" input features!')}
  }


  # transform angle input
  if(trans == TRUE)
  {
    # 0 to 90 degree
    angle <- ifelse(angle >= 0 & angle <= 90, 90 - angle, angle)

    # 90 to 180 degree
    angle <- ifelse(angle > 90 & angle < 180, 180 - (angle - 90), angle)

    # 180 to 270 degree
    angle <- ifelse(angle >= 180 & angle <= 270, 270 - angle, angle)

    # 270 to 360 degree
    angle <- ifelse(angle > 270 & angle <= 360, 180 - (angle - 270), angle)
  }

  # # # # # # # #
  ## create optimizer function
  # # # # # # # #

  # par - points of first lines guess
  # x.all - x extent - min and max values
  # n number of lines/parts
  # x.i sf object of run i
  opt.fun <- function(par, x.extent, n, x.i, splitByLines)
  {
    # browser()
    ## create lines
    # get points to optimize
    y.par <- par

    # split sf object by lines
    sfc.i <- splitByLines(x = x.extent, y = y.par, x.i = x.i)

    # get area of splitted polygons and squared residuals
    # sfc.i$A_RES <- ((sfc.i$A_FUNCT/n) - sfc.i$A_OPT)^2
    sfc.i$A_RES <- ((sfc.i$A_FUNCT/nrow(sfc.i)) - sfc.i$A_OPT)^2

    ## sum of squared residuals to optimize
    sum(sfc.i$A_RES)

  } #end of opt.fun

  # # # # # # # #
  ## rotate function (radiant)
  # # # # # # # #
  # x - radiant input
  rot <- function(x){matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2)}
  rot.inv <- function(x){matrix(c(cos(x), -sin(x), sin(x), cos(x)), 2, 2)}


  # # # # # # # #
  ## split sf object by lines
  # # # # # # # #
  # x - x coordiantes from extent
  # y - line coordiantes from guess
  # sf.i - sf-object to be split
  splitByLines <- function(x, y, x.i, final = FALSE)
  {
    # get line matrix
    m.xy <- expand.grid(x, y)
    m.lines <- lapply(X = seq(from = 1, to = n, by = 2), FUN = function(x, m.xy ){return(as.matrix(m.xy[c(x, x+1), ]))}, m.xy = m.xy)


    if(final)
    {
      # get lines
      lines <- sf::st_multilinestring(m.lines)
      x.i.dif <- list()
      x.i.run <- x.i

      for(i in 1:length(lines))
      {
        # get line string from multiline
        ls.i <- sf::st_geometry(sf::st_linestring(lines[[i]]))

        split.i <- sf::st_collection_extract(x = lwgeom::st_split(x = x.i.run, y = ls.i), type = "POLYGON")

        if((i+1) <= length(lines))
        {
          # get next line for intersection
          ls.i.next <- lines[[(i+1)]]

          # get interecting polygon that must be excluded from split
          inter <- unlist(sf::st_intersects(x = sf::st_linestring(ls.i.next), y = split.i))
          n.split <- 1:nrow(split.i)

          # write result of run
          x.i.dif[[i]] <- split.i[!n.split %in% inter,]

          # set polygon to be cut in next run
          x.i.run <- split.i[n.split %in% inter,]
         } else{

          # write result of run
            x.i.dif[[i]] <- split.i
        }
      } # end of for loop

      sfc.i <- do.call(rbind, x.i.dif) # combine final result of split

     } else {

      # get lines
      lines <- sf::st_sfc(sf::st_multilinestring(m.lines), crs = sf::st_crs(x.i))

      ## buffer line to get thin polygons
      lines.buf <- sf::st_buffer(x = lines, dist = 0.000001)  # create a very thin polygon

      ## splitting
      # calculate the differences of polygons and casting
      x.i.dif <- suppressWarnings(sf::st_difference(x = x.i, y = lines.buf))

      if(class(x.i.dif$geometry)[1] == "sfc_GEOMETRYCOLLECTION"){
        sfc.i <- sf::st_collection_extract(x = x.i.dif, "POLYGON")
      } else {
        sfc.i <- suppressWarnings(sf::st_cast(x = x.i.dif, to = "POLYGON"))
      }

    }

      ## calculate area
      sfc.i$A_OPT <- as.numeric(sf::st_area(x = sfc.i))

      return(sfc.i)
  } # end of function splitByLines




  # # # # # # # #
  ## divide features
  # # # # # # # #

  # calculate area
  x$A_FUNCT <- as.numeric(sf::st_area(x))
  x.cntrd <- sf::st_centroid(x)

  # get crs
  x.crs <- sf::st_crs(x)

  # init parallel
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl = cl, envir = environment(),
                          varlist = c('x', 'x.crs', 'angle', 'n', 'btol','pgtol', 'maxit', 'x.cntrd', 'rot', 'rot.inv', 'opt.fun', 'splitByLines'))

  parallel::clusterEvalQ(cl, library("sf", "lwgeom"))


  # x.div <- parallel::mclapply(mc.cores = mc.cores, X = 1:n.obj, fun = function(i, x, angle, n, btol, pgtol, maxit, x.cntrd, rot, opt.fun, splitByLines){
  x.div <- parallel::parLapply(cl = cl, X = 1:n.obj, fun = function(i, x.sf, crs, angle, n, btol, pgtol, maxit, x.cntrd, rot, rot.inv, opt.fun, splitByLines){
  # x.div <- lapply(X = 1:n.obj, FUN = function(i, x.sf, crs, angle, n, btol, pgtol, maxit, x.cntrd, rot, rot.inv, opt.fun, splitByLines){
      # browser()
      # get objects
        x.i <- x.sf[i,]

        # check if object is valid, and make it valuid
        if(suppressWarnings(!all(sf::st_is_valid(x.i)))){
          x.i <- lwgeom::st_make_valid(x.i)

          if(class(x.i$geometry)[1] == "sfc_GEOMETRYCOLLECTION")
          {
            x.i <- sf::st_collection_extract(x = x.i, type = "POLYGON")
          }

        }

        angle.i <- angle[i]
        x.cntrd.i <- sf::st_geometry(x.cntrd[i,])

        # rotate object
        x.i.rot <- (sf::st_geometry(x.i) - x.cntrd.i) * rot(angle.i*pi/180) + x.cntrd.i
        x.i.rot <- sf::st_sf(sf::st_set_geometry(x.i, NULL) , geometry = x.i.rot, crs = sf::st_crs(x.i))

        if(suppressWarnings(!all(sf::st_is_valid(x.i.rot)))){
          x.i.rot <- lwgeom::st_make_valid(x.i.rot)

          if(class(x.i.rot$geometry)[1] == "sfc_GEOMETRYCOLLECTION")
          {
            x.i.rot <- sf::st_collection_extract(x = x.i.rot, type = "POLYGON")
          }
        }

        # test
        # plot(sf::st_geometry(x.i))
        # plot(sf::st_geometry(x.i.rot), add = T)

        # get bounding box and start values
        bb.i <- sf::st_bbox(x.i.rot)

        y.start <- bb.i[2]
        y.end <- bb.i[4]
        y.par <- seq(from = y.start, to = bb.i[4], length.out = n + 1)[2:n]
        x.extent <- c(bb.i[1], bb.i[3])

        # # #
        # search for optimized splitting lines
        # # #
        y.par.opt <- optim(par = y.par, fn = opt.fun, method = "L-BFGS-B",
                          control = list(pgtol = pgtol, maxit = maxit),
                          lower = (y.start + btol), upper = (y.end - btol),
                          x.extent = x.extent, n = n, x.i = x.i.rot,
                          splitByLines = splitByLines)$par


        # # #
        # split polygon
        # # #
        x.i.rot.div <- splitByLines(x = x.extent, y = y.par.opt, x.i = x.i.rot, final = TRUE)
        # x.i.rot.div <- lwgeom::st_make_valid(sf::st_buffer(x.i.rot.div, dist = 0.000001))

        if(suppressWarnings(!all(sf::st_is_valid(x.i.rot.div))))
        {
          x.i.rot.div <- lwgeom::st_make_valid(x.i.rot.div)
        }

        # sf::st_intersects(x.i.rot.div)
        # plot(st_geometry(sf::st_intersection(x.i.rot.div)))

        # rotate back
        x.i.div <- (sf::st_geometry(x.i.rot.div) - x.cntrd.i) * rot.inv(angle.i*pi/180) + x.cntrd.i
        x.i.div <- sf::st_sf(sf::st_set_geometry(x.i.rot.div, NULL) , geometry = x.i.div, crs = sf::st_crs(x.i.rot.div))

        # set crs
        x.i.div <- sf::st_set_crs(x.i.div, crs)

        # set unique ID column
        x.i.div$ID_Split <- paste0(i, "_", c(1:nrow(x.i.div)))

        # check last time geometry
        if(suppressWarnings(!all(sf::st_is_valid(x.i.div))))
        {
          x.i.div <- lwgeom::st_make_valid(x.i.div)

          if(class(x.i.div$geometry)[1] == "sfc_GEOMETRYCOLLECTION")
          {
            x.i.div <- sf::st_collection_extract(x = x.i.div, type = "POLYGON")
          }
        }

        # test
        # plot(sf::st_geometry(x.i.div))
        # plot(sf::st_geometry(x.i.rot.div), add = T)

        return(x.i.div)

      }, x.sf = x, crs = x.crs, angle = angle, n = n, btol = btol, pgtol = pgtol, maxit = maxit, x.cntrd = x.cntrd, rot = rot, rot.inv = rot.inv, opt.fun = opt.fun, splitByLines = splitByLines)

  parallel::stopCluster(cl)

  # rbind simple features
  x.div <- do.call(rbind, x.div)
  # plot(sf::st_geometry(x.div))

  # get time of process
  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) cat("------ Run of dividePolygon: " , round(x = process.time.run["elapsed"][[1]]/60, digits = 3), " Minutes ------\n")


  return(x.div)
} # end of function dividePolygon



# setwd("d:/Users/qo23hel/Desktop/Test")
# x <- sf::st_read("test.shp")
# angle <- 90
# n = 3
# tol = 1
# cores = 2
# mc.cores = 2
# pgtol = 1e-9
# maxit = 100
# https://tereshenkov.wordpress.com/2017/09/10/dividing-a-polygon-into-a-given-number-of-equal-areas-with-arcpy/

