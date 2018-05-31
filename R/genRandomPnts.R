#' Generate random points based on distance constraints
#'
#' This function generate random points inside polygons based on distance constraints.
#' The function is based on spatstat::rSSI.
#'
#' @param x object of class sf
#' @param seed seed parameter for replication. Default: 123
#' @param dist inhibition distance of spatstat::rSSI. Default: NULL
#' @param n maximum number of points allowed. If n is finite, stop when the total number of points in the point pattern reaches n. If n is infinite (the default), stop only when it is apparently impossible to add any more points
#' @param maxit number of rejected proposals after which the algorithm should terminate. Default: 100
#' @param cores number of cores for parallel processing. Default: 1 (sequential)
#' @param quiet no outputs in console. Default: TRUE
#' @param ... additional parameters for spatstat::rSSI-function
#' @return
#' Randomly sampled points of class sf
#'
#'
#' @keywords random points, distance constraint
#'
#'
#' @export
#'
genRandomPnts <- function(x, seed = 123, dist = NULL, n = Inf, maxit = 100, quiet = TRUE, cores = 1, ...)
{

  # get start time of process
  process.time.start <- proc.time()

  # get crs
  crs <- sf::st_crs(x = x)

  # convert simple feature to spatial polygons
  x.sp <- x %>% as(., "Spatial") %>%  as(., "SpatialPolygons")


  # convert to owin object
  x.owin <- x.sp %>%
    slot(., "polygons") %>%
    lapply(X = ., FUN = function(x){sp::SpatialPolygons(list(x))}) %>%
    lapply(X = ., FUN = maptools::as.owin.SpatialPolygons) # spatstat::as.owin



  # generate random sampling with distant constraint (can be parallelized)
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl = cl, varlist = ls(), envir = environment())

  pts.ppp <- parallel::parLapply(cl = cl, X = 1:length(x.owin), fun = function(i, x.owin, r, n, quiet, seed, maxit, ...){
  # pts.ppp <- lapply(X = 1:length(x.owin), FUN = function(i, x.owin, r, n, quiet, seed, maxit, ...){


    # if(quiet == FALSE) cat("Run ", i, " of ", length(x.owin), "\n")
    set.seed(seed)

    pts.out <- spatstat::rSSI(r = r, n = n, giveup = maxit, win = x.owin[[i]], ...)

    return(list(pts.out, rep(i, pts.out$n)))

  }, quiet = quiet, x.owin = x.owin, r = dist, n = n, seed = seed, maxit = maxit, ...)

  parallel::stopCluster(cl)

  # back-conversion to simple feature
  pts.sf <- lapply(X = pts.ppp, FUN = '[[', 1) %>%
    lapply(X = ., FUN = function(x) sf::st_as_sfc(as(x, "SpatialPoints"))) %>%
    do.call(c, .)

  pts.sf <- sf::st_sf(ID = 1:length(pts.sf), geometry = pts.sf, crs = crs)


  # get intersected items
  pts.inter.x <- sf::st_intersects(x = pts.sf, y = x) %>% unlist

  if(length(pts.inter.x) != nrow(pts.sf))
  {
    warning("Some sample points are inside multiple polygons")
    if(length(unlist(lapply(X = pts.ppp, FUN = '[[', 2))) == nrow(pts.sf)){
    pts.sf$In <- unlist(lapply(X = pts.ppp, FUN = '[[', 2))
    }
  } else{
    pts.sf$In <- pts.inter.x
  }

  # get time of process
  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) cat(paste0("------ Run of genRandomPtsDist: " , round(x = process.time.run["elapsed"][[1]]/60, digits = 3), " Minutes ------\n"))

  return(pts.sf)
} # end of function genRandomPnts
