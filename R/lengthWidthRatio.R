#' Length-Width-Ratio
#'
#' This function calculates the length-width-ratio for SpatialPolygonsDataFrame objects.
#' The length-width vectors for the ratio are computed by a principal component analysis according to YI & MARSHALL (2000)
#'
#' @param spdf SpatialPolygonsDataFrame input
#' @param cores number of cores for parallel processing. Default: 1 (sequential)
#' @param quiet no outputs in console. Default: TRUE
#'
#' @note
#' \itemize{
#'   \item for the principal component analysis the \code{\link[stats]{princomp}} function is used
#'   \item SpatialPolygonsDataFrame input must have an \emph{ID} field with unique numbers
#'   \item YI, W.& S. MARSHALL (2000): Principal component analysis in application
#'         to object orientation. Geo-spatial Information Science, 3(3), 76-78.
#' }
#'
#'
#' @return
#' data.table containing the \emph{ID} field and a column (\emph{ratio}) with length-width-ratio
#'
#'
#'
#'
#' @keywords length-width-ratio, principal component analysis
#'
#'
#'
#'
#' @export
lengthWidthRatio <- function(spdf, cores = 1, quiet = TRUE)
{

  # get start time of process
  process.time.start <- proc.time()

  # read spdffile
  # spdf <- readspdfPoly(spdffile)

  # vector with Ratios
  # vRatio <- c()

  # init parallel
  switch(Sys.info()[[1]],
         Windows = type  <- "PSOCK",
         Linux = type  <- "FORK",
         Mac = type <- "FORK")

  cl <- parallel::makeCluster(cores, type = type)

  if(Sys.info()[[1]] == "Windows")
  {
    parallel::clusterExport(cl = cl, envir = environment(),
                            varlist = c('spdf'))

    parallel::clusterEvalQ(cl, library("sp", "rgeos"))
  }

  vRatio <- parallel::parSapply(cl = cl, X = 1:length(spdf@polygons), FUN = function(i, spdf){
  # vRatio  <- sapply(X = 1:length(spdf@polygons), FUN = function(i, spdf){
  # for(i in 1:length(spdf@polygons))
  # {

    ### perform principal component analysis
    pca.matrix <- matrix(ncol = 2)
    # retrieve polygons (nevessary if polygon consists of multiple parts)
    for(j in 1:length(spdf@polygons[[i]]@Polygons))
    {
      pca.matrix <- rbind(pca.matrix, spdf@polygons[[i]]@Polygons[[j]]@coords)
    }
    pca.matrix <- pca.matrix[-1:-1,]

    pca <- princomp(pca.matrix)

    eigenvector <- pca$loadings # eigenvectors
    eigenvalue <- pca$sdev^2 # retrieve eigenvalues

    # get length of vectors
    v1.x <- eigenvector[1] * eigenvalue[1] # last point of vector 1, x coordinate
    v1.y <- eigenvector[2] * eigenvalue[1] # last point of vector 1, y coordinate
    v2.x <- eigenvector[3] * eigenvalue[2] # last point of vector 2, x coordinate
    v2.y <- eigenvector[4] * eigenvalue[2] # last point of vector 2, y coordinate
    v1 <- sp::SpatialLines(list(sp::Lines(list(sp::Line(cbind(c(pca$center[1], pca$center[1] + v1.x), c(pca$center[2], pca$center[2] + v1.y)))), ID = "v1")))
    v2 <- sp::SpatialLines(list(sp::Lines(list(sp::Line(cbind(c(pca$center[1], pca$center[1] + v2.x), c(pca$center[2], pca$center[2] + v2.y)))), ID= "v2")))
    v1.length <- rgeos::gLength(v1)
    v2.length <- rgeos::gLength(v2)

    # calculate ratio: larger value is the numerator, s. ecognition reference book (2014, 285)
    if(v1.length >= v2.length)
    {
      ratio <- v1.length/v2.length
    }
    else
    {
      ratio <- v2.length/v1.length
    }

    # append ratio
    # vRatio <- c(vRatio, ratio)
    return(ratio)

  # } # end for
  }, spdf = spdf) # end lapply

  parallel::stopCluster(cl)


  dt.ratio <- data.frame(ID = spdf@data$ID)
  dt.ratio$ratio <- vRatio


  # get time of process
  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) cat(paste0("------ Run of LengthWidthRatio: " , round(x = process.time.run["elapsed"][[1]]/60, digits = 4), " Minutes ------ \n"))


  return(dt.ratio)
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
