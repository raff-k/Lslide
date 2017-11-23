#' Conversion of raster layer to gray-scale image
#'
#' This function convert a \linkS4class{RasterLayer} input to a gray-scale image, so that it can be used within \code{\link[EBImage] functions.
#'
#' @param r \linkS4class{RasterLayer} to be converted either to an image or a matrix
#' @param NA.val.in replace value for no data (NA) in \emph{r}. Default: 0
#'
#' @return
#' \linkS4class{data.table} containing object-based textural features
#'
#'
#' @note
#' \itemize{
#'   \item \emph{r} is stretched to gray-scale (0-255) during process by linear normalization
#' }
#'
#'
#' @keywords EBImage, raster, gray-scale-image, matrix, conversion
#'
#'
#' @export
convR2GSM <- function(r, NA.val.in = 0)
{

  funReplace <- function(..., NA.val.in){ifelse(is.na(...), NA.val.in, ...)}
  r <-  raster::calc(r, function(x){funReplace(x, NA.val.in = NA.val.in)})


  # stretching to 0-255 range
  minV <- minValue(r)
  maxV <- maxValue(r)

  # the linear normalization
  # scale function
  funScale <- function(..., min, max){(... -min)*255/(max-min)+0}

  # scale grid
  r <- raster::calc(r, function(x){funScale(x, min = minV, max = maxV)})

  # dividing by 255 to fit haralick input
  r <- raster::calc(r, fun=function(x){return(x/255)})

  # convert to matrix
  r <- as.matrix(r)

  # return matrix
  return(r)
}
