#' Object orientation or main direction
#'
#' This function calculates the object orientation (or main direction) for SpatialPolygonsDataFrame objects
#' in xy-space. For this, a principal component analysis is performed according to YI & MARSHALL (2000)
#'
#' @param spdf \linkS4class{SpatialPolygonsDataFrame} input
#' @param angle angle taken for output. Default: "main" (alternative: "weighted main")
#' @param quiet no outputs in console. Default: TRUE
#'
#' @note
#' \itemize{
#'   \item for the principal component analysis the \code{\link[stats]{princomp}} function is used
#'   \item \linkS4class{SpatialPolygonsDataFrame} input must have an \emph{ID} field with unique numbers
#'   \item YI, W.& S. MARSHALL (2000): Principal component analysis in application
#'         to object orientation. Geo-spatial Information Science, 3(3), 76-78.
#' }
#'
#'
#' @return
#' data.table containing the \emph{ID} field, a column (\emph{angle}) containing the object orientations,
#'  and a column (\emph{angle_inv}) contaning the inverse object orientations. Both object orientations are necessairy
#'  due to the ambiguity of the orientation in xy-space without fix starting point
#'
#'
#'
#' @keywords object orientation, main direction, principal component analysis
#'
#' @examples
#' mainDirection()
#'
#' @export
mainDirection <- function(spdf, angle = "main", quiet = TRUE)
{

  # get start time of process
  process.time.start <- proc.time()


  # vector with angles
  vAngle <- c()

  for(i in 1:length(spdf@polygons))
  {

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

    # get orientation and length of vectors
    v1.x <- eigenvector[1] * eigenvalue[1] # last point of vector 1, x coordinate
    v1.y <- eigenvector[2] * eigenvalue[1] # last point of vector 1, y coordinate
    v2.x <- eigenvector[3] * eigenvalue[2] # last point of vector 2, x coordinate
    v2.y <- eigenvector[4] * eigenvalue[2] # last point of vector 2, y coordinate
    v1 <- SpatialLines(list(Lines(list(Line(cbind(c(pca$center[1], pca$center[1] + v1.x), c(pca$center[2], pca$center[2] + v1.y)))), ID = "v1")))
    v2 <- SpatialLines(list(Lines(list(Line(cbind(c(pca$center[1], pca$center[1] + v2.x), c(pca$center[2], pca$center[2] + v2.y)))), ID= "v2")))


    ### calculation of angles - atan2 gives values in radiant and must be converted to degree

    # calculate weighted main angle of object
    if(angle == "weighted main")
    {
      a1 <- as.numeric(atan2(v1.x, v1.y) * 180/pi)
      if(a1 < 0) a1 <- a1 + 180
      a2 <- as.numeric(atan2(v2.x, v2.y) * 180/pi)
      if(a2 < 0) a2 <- a2 + 180

      a <- weighted.mean(c(a1, a2), c(gLength(v1), gLength(v2)))
    }

    # calculate main angle of object
    if(angle == "main")
    {
      # check which vector is the longest for the main direction
      if(gLength(v1) >= gLength(v2))
      {
        vL.x <- v1.x
        vL.y <- v1.y

        vS.x <- v2.x
        vS.y <- v2.y
      } else
      {
        vL.x <- v2.x
        vL.y <- v2.y

        vS.x <- v1.x
        vS.y <- v1.y
      }

      a <- as.numeric(atan2(vL.x, vL.y) * 180/pi)
      if(a < 0) a <- a + 180
    }

    # append angle
    vAngle <- c(vAngle, a)
  } # end for

  dt.angle <- data.table(ID = spdf@data$ID)
  dt.angle$angle <- vAngle
  dt.angle$angle_inv <- vAngle + 180


  # get time of process
  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) print(paste0("------ Run of MainDirection: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))


  return(dt.angle)
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
