#' Calculation of object relations
#'
#' This function calculates statistics of \linkS4class{SpatialPolygonsDataFrame} object relations
#'
#' @param spdf \linkS4class{SpatialPolygonsDataFrame}, or full path of a shapefile
#' @param nb object neighborhood based on \code{\link[spdep]{poly2nb}}
#' @param var vector with number or name defining 1. field of evaluation (i.e. slope) and 2. field of weights (i.e. area). For example: c("Slope", "Area") or c(1,2)
#' @param quiet no outputs in console. Default: TRUE
#'
#' @return
#' \linkS4class{data.frame} containing object statistics
#'
#' @note
#' \itemize{
#'   \item relational object statistics according to ECOGNITION DEVELOPER (2014: 253-255):
#'   \itemize{
#'       \item \emph{m_wei} - mean value of selected features of an object and its neighbors  (of a selected class, see).
#'         For averaging, the feature values are weighted with the area of the objects.
#'       \item \emph{sd_nb} - standard deviation of selected features of an object and its neighbors
#'       It is also possible to calculate this feature with respect to a class, see
#'       \item \emph{m_dif} - mean difference between the feature value of an object and its neighbors (of a selected class, see)
#'       The feature values are weighted by the area of the respective objects
#'       \item \emph{m_dif_abs} - mean absolute difference
#'       \item \emph{m_dif_hv} - mean difference between the feature value of an object and the feature values of its neighbors (of a selected class, see),
#'        which have higher values than the object itself. The feature values are weighted by the area of the respective objects
#'       \item \emph{m_dif_lw} - mean difference between the feature value of an object and the feature values of its neighbors (of a selected class, see),
#'        which have lower values than the object itself. The feature values are weighted by the area of the respective objects
#'       \item \emph{ratio_nb} - proportion between the feature value of an object and the mean feature value of its neighbors (of a selected class, see)
#'         For averaging the feature values are weighted with the area of the corresponding objects
#'       \item \emph{sum_nb} - sum of the feature values of the neighbors (of a selected class, see )
#'       \item \emph{num_nb} - number of neighbors of (of a selected class, see )
#'       \item \emph{min_nb} - minimum value of the feature values of an object and its neighbors (of a selected class, see)
#'       \item \emph{max_nb} - maximum value of the feature values of an object and its neighbors (of a selected class, see)
#'   }
#'   \item ECOGNITION DEVELOPER (2014) Reference Book. Trimble Documentation, München, Germany
#' }
#'
#'
#' @keywords relation of objects, spdep, object-oriented image analysis
#'
#'
#' @export
#'
relationalObjectFunction <- function(spdf, nb, var, quiet = TRUE)
{
  if(length(nb) != nrow(spdf))
  {
    stop("mismatch between amount of attributes and neighbors!")
  }

  # get start time of process
  process.time.start <- proc.time()

  # create empty table
  result.table <- data.frame(ID = integer(),
                             Value = double(),
                             Weight = double(),
                             mean.w = double(),
                             sd.nb  = double(),
                             mean.dif = double(),
                             mean.dif.abs = double(),
                             mean.dif.to.hv = double(),
                             mean.dif.to.lw = double(),
                             ratio.nb = double(),
                             sum.nb = double(),
                             num.nb = double(),
                             min.nb = double(),
                             max.nb = double(),
                             stringsAsFactors=FALSE)


  # loop through every object, start with first in data frame of spdf, then second, and so on
  for(i in 1:nrow(spdf))
  {
    # espdftract neighbors
    nb.i <- i
    nb.i <- c(nb.i, nb[[i]]) # neighbors including selected object
    nb.n <- c(nb[[i]]) # neighbors espdfcluding selected object

    # get data frame of neighbors
    nb.df.i <- as.data.frame(spdf[nb.i,]) # including selected object
    nb.df.n <- as.data.frame(spdf[nb.n,]) # espdfcluding selected object
    value.i <- spdf[[var[1]]][i] # value of i
    weight.i <- spdf[[var[2]]][i] # weight of i

    # polygons with NA as value or without any neighbor will be ignored
    if(is.na(value.i) | nb.n[1] == 0)
    {
      # polygon value is NA
      if(is.na(value.i))
      {
        results.stat <- c(i, rep(NA, 13))
      }

      # polygon has no neighbor
      if(nb.n[1] == 0)
      {
        results.stat <- c(i, rep(NA, 10), 0, NA, NA)
      }

      # put results in table
      result.table <- rbind(result.table, results.stat)

      # jump loop
      nespdft
    }

    ####### calculate statistics, adapted by eCognition (2014) reference book, p. 253 ff.
    mean.w <- weighted.mean(nb.df.i[[var[1]]], nb.df.i[[var[2]]], na.rm = TRUE) # weighted mean
    sd.nb <- sd(nb.df.i[[var[1]]], na.rm = TRUE) # standard deviation
    if(is.na(sd.nb)){sd.nb <- 0}

    ### important parameter: postive means higher | negative means smaller
    mean.dif <-  weighted.mean((value.i - nb.df.n[[var[1]]]), nb.df.n[[var[2]]], na.rm = TRUE) # Calculates the mean difference between the feature value of an image object and its neighbors of a selected class. Note that the feature values are weighted by area
    mean.dif.abs <- abs(mean.dif) # absolute value
    mean.dif.to.hv <- weighted.mean(value.i - nb.df.n[[var[1]]][nb.df.n[[var[1]]] > value.i], nb.df.n[[var[2]]][nb.df.n[[var[1]]] > value.i], na.rm = TRUE) # the mean difference between the feature value of an image object and the feature values of its neighbors, which have higher values than the object itself, weighted by area
    if(is.na(mean.dif.to.hv)){ mean.dif.to.hv <- 0}
    mean.dif.to.lw <- weighted.mean(value.i - nb.df.n[[var[1]]][nb.df.n[[var[1]]] < value.i], nb.df.n[[var[2]]][nb.df.n[[var[1]]] < value.i], na.rm = TRUE) # the mean difference between the feature value of an image object and the feature values of its neighbors, which have lower values than the object itself, weighted by area
    if(is.na(mean.dif.to.lw)){ mean.dif.to.lw <- 0}
    ###

    ratio.nb <- value.i / weighted.mean(nb.df.n[[var[1]]], nb.df.n[[var[2]]], na.rm = TRUE)  # ratio between the feature value of an image object and the weighted mean feature value of its neighbors
    sum.nb <- sum(nb.df.n[[var[1]]], na.rm = TRUE) # Calculates the sum of the feature values of the neighbors
    if(nb.n[1] == 0){num.nb <- 0}else {num.nb <- length(nb.n)} # Calculates the number of neighbors
    min.nb <- min(nb.df.i[[var[1]]], na.rm = TRUE) # Returns the minimum value of the feature values of an image object and its neighbors
    max.nb <- max(nb.df.i[[var[1]]], na.rm = TRUE) # Returns the maximum value of the feature values of an image object and its neighbors

    # results to vector
    results.stat <- c(i, value.i, weight.i, mean.w, sd.nb, mean.dif, mean.dif.abs, mean.dif.to.hv, mean.dif.to.lw,
                      ratio.nb, sum.nb, num.nb, min.nb, max.nb)

    # put results in table
    result.table <- rbind(result.table, results.stat)
  }

  # rename table columns
  names(result.table) <- c("ID", "Value", "Weight","m_wei", "sd_nb", "m_dif", "m_dif_abs", "m_dif_hv",
                           "m_dif_lw", "ratio_nb",  "sum_nb", "num_nb", "min_nb", "max_nb")

  # get time of process
  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) print(paste0("------ Run of RelationalFunction: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))

  # return output
  return(result.table)
}
