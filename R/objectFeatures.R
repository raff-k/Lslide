#' Caclulation of object features
#'
#' This function calculates object features based on the \code{\link[EBImage]{computeFeatures}} function
#'
#' @param r.seg segments or objects in form of a \linkS4class{RasterLayer} (only 16 or 32 bit integer)
#' @param quiet no outputs in console. Default: TRUE
#'
#' @return
#' \linkS4class{data.table} containing object features
#'
#'
#' @note
#' \itemize{
#'   \item as input only 16 or 32 bit integer values are allowed
#'   \item see \code{\link[EBImage]{computeFeatures}}
#' }
#'
#'
#' @keywords EBImage, computeFeatures, segmentation
#'
#' @examples
#' objectFeatures()
#'
#' @export
objectFeatures <- function(r.seg, quiet = TRUE)
{

  # get start time of process
  process.time.start <- proc.time()

  # read data
  # segmentsID <- raster(object.feature.segments) # raster containng the segments ID (16 or 32 bit integer)
  segmentsID <- as.matrix(r.seg) # get matrix out of data

  # calculate object features
  object.shape <- computeFeatures.shape(segmentsID, properties = FALSE)
  object.moment <- as.data.table(computeFeatures.moment(segmentsID, properties=FALSE))

  # change object angle - object angle scale is E: 90, S: 0, W: -90 and should be 90, 180, 270
  object.moment$m.theta <- object.moment$m.theta * 180/pi
  object.moment$m.angle <- with(object.moment, ifelse(m.theta > 0, 180 - m.theta, abs(m.theta)))


  # create data table for output
  dt.object.features <- data.table(ID = c(1:nrow(object.shape)), stringsAsFactors = FALSE)
  dt.object.features <- cbind(dt.object.features, object.shape, object.moment)


  # get time of process
  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) print(paste0("------ Run of ObjectFeatures: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))

  # return object textures
  return(dt.object.features)
}
