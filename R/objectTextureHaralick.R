#' Caclulation of object-based textural features
#'
#' This function calculates textural features for objects by considering the values inside of the objects.
#' The calculation is based on the \code{\link[EBImage]{computeFeatures.haralick}} function
#'
#' @param r.texture \linkS4class{RasterLayer} for which the textural features are calculated
#' @param r.seg segments or objects in form of a \linkS4class{RasterLayer} (only 16 or 32 bit integer)
#' @param texture.haralick.nbins  size of gray level co-occurrence matrix. Default: 32
#' @param texture.haralick.scales definition of scales on which the texture is measured. Default: c(1,2)
#' @param NA.val.in replace value for no data (NA) in \emph{r.texture}. Default: 0
#' @param quiet no outputs in console. Default: TRUE
#'
#' @return
#' \linkS4class{data.table} containing object-based textural features
#'
#'
#' @note
#' \itemize{
#'   \item \emph{r.texture} is stretched to gray-scale (0-255) during process by linear normalization
#'   \item \emph{r.seg} must consist either of 16 bit or of 32 bit integer values
#'   \item the calculated Haralick features are rotation invariant (average over all four possible directions)
#'   \item a smaller texture.haralick.nbins can give more accurate estimates of the correlation, because the number of events per bin is higher.
#'   While a higher value will give more sensitivity
#'   \item see \code{\link[EBImage]{computeFeatures.haralick}}
#'   \item see \href{http://earlglynn.github.io/RNotes/package/EBImage/Haralick-Textural-Features.html}{here} for more details on Haralick textural features \emph{(last call: 13-04-2017)}
#' }
#'
#'
#' @keywords EBImage, computeFeatures.haralick, object-based textural features
#'
#' @examples
#' objectTextureHaralick()
#'
#' @export
objectTextureHaralick <- function(r.texture, r.seg, texture.haralick.nbins = 32, texture.haralick.scales = c(1,2), NA.val.in = 0, quiet = TRUE)
{

  # get start time of process
  process.time.start <- proc.time()

  # # # read data

  ## creation of gray scale image

  # replace NA
  # r.texture[is.na(r.texture)] <- NA.value
  funReplace <- function(..., NA.val.in){ifelse(is.na(...), NA.val.in, ...)}
  r.texture <-  raster::calc(r.texture, function(x){funReplace(x, NA.val.in = NA.val.in)})


  # stretching to 0-255 range
  minV <- minValue(r.texture)
  maxV <- maxValue(r.texture)

  # the linear normalization
  # scale function
  funScale <- function(..., min, max){(... -min)*255/(max-min)+0}

  # scale grid
  r.texture <- raster::calc(r.texture, function(x){funScale(x, min = minV, max = maxV)})

  # dividing by 255 to fit haralick input
  r.texture <- raster::calc(r.texture, fun=function(x){return(x/255)})

  # convert to matrix
  r.texture <- as.matrix(r.texture)


  ## creation of segments or objects input
  # r.seg[is.na(r.seg)] <- 0   # set no data to 0
  r.seg <- raster::calc(r.seg, fun=function(x){ifelse(is.na(x), 0, x)})

  # get unique IDs out of segments
  r.seg.pos <- unique(raster::values(r.seg))
  r.seg.pos <-  r.seg.pos[! r.seg.pos %in% 0] # remove 0

  # raster() image must be transformed to fit to readImage()
  r.seg <- as.matrix(r.seg) # get matrix out of data
  # r.seg<- t(r.seg) # same as aperm(x, c(2,1))

  # calculate texture
  texture <- computeFeatures.haralick(r.seg, r.texture, properties = FALSE, haralick.nbins = texture.haralick.nbins, haralick.scales = texture.haralick.scales)

  # create data.table for output
  dt.texture <- data.table(ID = c(1:nrow(texture)), stringsAsFactors = FALSE)
  dt.texture <- cbind(dt.texture, texture)
  colnames(dt.texture) <- gsub("\\.", "_", colnames(dt.texture))

  # kick out not segmentID columns
  dt.texture <- dt.texture[r.seg.pos]

  # get time of process
  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) print(paste0("------ Run of TextureObjectHaralick: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))


  # return object textures
  return(dt.texture)
}
