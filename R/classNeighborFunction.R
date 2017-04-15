#' Calculation of relations from a class to its neighbor objects
#'
#' This function calculates statistics for every object of a specific class to its neighbor objects.
#' It is computationally identical to the \code{\link{relationalObjectFunction}}.
#' However, the function additionally enables the possibility to perform statistics for neighbors in a specific direction,
#' likewise the \code{\link{relationalClassFunction}}, and to select class objects by expressions.
#'
#'
#' @param spdf \linkS4class{SpatialPolygonsDataFrame}, or full path of a shapefile
#' @param nb object neighborhood based on \code{\link[spdep]{poly2nb}}. Otherwise, neighborhood is computed by spdep::poly2nb(..., Queens = TRUE)
#' @param var vector with number or name defining 1. field of evaluation (i.e. slope) and 2. field of weights (i.e. area). For example: c("Slope", "Area") or c(1,2)
#' @param class.var defination of class. Can be string or integer, or an expression in form of a list. For examples: class = as.integer(1), class = list("> 2.5", "expression"), class = list(c("> 2", "&", "< 7"), "expression")
#' @param class.as.neighbors class object is also neighbor object. Default: FALSE
#' @param only.neighbors object statistics are computed only for neighbors. Default: TRUE (FALSE means computation for the entiry data)
#' @param calc.nb.flow include statistics for neighbors in flow direction. Default: FALSE
#' @param col.flow column specifying direction. Default: NULL
#' @param bb input of bounding box geometry. Default: NULL
#' @param bb.intersection method for intersection with bounding box. Default: "rgeos" (alternative: "sp")
#' @param return.bb.flow return geoemtry of counding box. Default: FALSE
#' @param other... see \code{\link{getBoundingBox}}
#' @param quiet no outputs in console. Default: TRUE
#'
#' @return
#' \linkS4class{data.frame} containing object statistics, and optional the created geometry of \code{\link{getBoundingBox}}
#'
#' @note
#' \itemize{
#'   \item class statistics to neighbor objects acccording to ECOGNITION DEVELOPER (2014: 253-255, 375-381):
#'   \itemize{
#'       \item see \code{\link{relationalObjectFunction}}
#'       \item \emph{~_nb_FlAl} - statistics for all objects located in flow direction
#'       \item \emph{~_nb_FlMa} - statistics for direct neighbors located in flow direction
#'       }
#'   \item ECOGNITION DEVELOPER (2014) Reference Book. Trimble Documentation, München, Germany
#'   \item see \code{\link{getBoundingBox}}
#'   \item \linkS4class{SpatialPolygonsDataFrame} objects with NA values are ignored
#'   \item \linkS4class{SpatialPolygonsDataFrame} objects without neighbors are ignored
#'   \item \linkS4class{SpatialPolygonsDataFrame} input must have a valid \emph{Class} and \emph{ID} field.
#'      \emph{ID} field must be contained of unique numbers
#' }
#'
#'
#' @keywords relation of class objects to their neighbor objects, spdep, object-oriented image analysis
#'
#'
#' @export
#'
classNeighborFunction <- function(spdf, nb, class.var, var, calc.nb.flow = FALSE, col.flow = NULL, bb.intersection = "rgeos", class.as.neighbors = FALSE, bb = NULL,
                                  return.bb.flow = FALSE, scale.factor.flow = c(1, 1), centroid.flow = TRUE, set.centroid.flow =  "inverse", k.centroid.flow = 2, k.flow = 2, scale.side.flow = "small", quiet = TRUE)
{
  # get start time of process
  process.time.start <- proc.time()

  if(calc.nb.flow == TRUE)
  {
    if(is.null(col.flow) & is.null(bb))
    {
      stop("column containing flow direction values must be specified")
    }

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
                               mean.w.flow.a.i = double(),
                               sd.nb.flow.a.i = double(),
                               mean.w.flow.a.n = double(),
                               sd.nb.flow.a.n = double(),
                               mean.w.flow.m.i = double(),
                               sd.nb.flow.m.i = double(),
                               mean.w.flow.m.n = double(),
                               sd.nb.flow.m.n = double(),
                               mean.dif.a = double(),
                               mean.dif.abs.a = double(),
                               mean.dif.to.hv.a = double(),
                               mean.dif.to.lw.a = double(),
                               mean.dif.m = double(),
                               mean.dif.abs.m = double(),
                               mean.dif.to.hv.m = double(),
                               mean.dif.to.lw.m = double(),
                               ratio.nb.a = double(),
                               sum.nb.a = double(),
                               num.nb.a = double(),
                               min.nb.a = double(),
                               max.nb.a = double(),
                               ratio.nb.m = double(),
                               sum.nb.m = double(),
                               num.nb.m = double(),
                               min.nb.m = double(),
                               max.nb.m = double(),
                               stringsAsFactors=FALSE)


  } else
  {
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
  }

  # read shapefile (maptools required)
  if(!class(spdf) == "SpatialPolygonsDataFrame")
  {
    sf <- sf::st_read(dsn = dirname(spdf), layer = basename(file_path_sans_ext(spdf)), quiet = TRUE)
    x <- as(sf, "Spatial")

    rm(sf, spdf)
  } else
  {
    x <- spdf
    rm(spdf)
  }

  ### check class attributes
  class.type <- class(x@data$Class) # retrieve type of class column in data frame

  # check input of class variable
  if(length(class.var) == 1) # numbers and factors are allowed
  {
    # check if type is the same
    if(class.type == class(class.var))
    {
      # get IDs from related class objects
      ID.class <- subset(x, x@data$Class == class.var)@data$ID
      # ID.class <- x[x@data$Class == class.var,]@data$ID # does not work with NAs
    } else # if the type is different. stop!
    {
      stop("Mismatch in class types. Check type of class input variable and Class column!")
    }
  } else # if the input is an expression
  {
    if(length(class.var) == 2 & class.var[2] == "expression") # check expression
    {
      class.tmp <- class.var[1]
      var.tmp <- var[1]

      # create expression
      if(class(var[1]) == "character")
      {
        class.var <- paste0(as.character(var[1]), " ", class.tmp)
        var[1] <- paste0('"', var[1], '"') # parse string that it fits in expression
      }

      if(length(class.tmp[[1]]) == 1)
      {
        # build expression
        expression <- paste0('x@data[x@data[', as.character(var[1]), "] ", class.tmp, ",]")

        # set class for table
        class.var <- paste0(as.character(var.tmp), " ", class.tmp)
      } else if(length(class.tmp[[1]]) == 3)
      {
        # split expression
        exp_first <- class.tmp[[1]][1]
        exp_operator <-class.tmp[[1]][2]
        exp_last <- class.tmp[[1]][3]

        # build expression
        expression <- paste0('x@data[x@data[', as.character(var[1]), "] ", exp_first, " ",
                             exp_operator, ' x@data[', as.character(var[1]), "] ", " ", exp_last, ",]")

        # set class for table
        class.var <- paste0(as.character(var.tmp), " ", exp_first, " ", exp_operator, " ", exp_last)
      } else
      {
        stop("Wrong expression input.")
      }

      # get IDs from exression result
      ID.class <- eval(parse(text = expression))$ID

      if(length(ID.class) == length(x@data$ID) | length(ID.class) == 0)
      {
        stop("Expression is invalid. Either no results or all ID's selected.")
      }

      # reset var
      var[1] <- var.tmp

    } else # if the type is different. stop!
    {
      stop("Class input invalid!")
    }
  }



  # check and clean geometry (cleangeo required)
  # x <- createSPComment(x)
  # x <- checkGeometry(x)

  # check neighborhood and create it if neighborhood does not exists
  if(missing(nb))
  {
    nb_speed_up <- rgeos::gUnarySTRtreeQuery(x) # speed up function for poly2nb
    nb <- spdep::poly2nb(x, queen = TRUE, foundInBox = nb_speed_up) # neighborhood based on queen continuity
  }

  # extract neighbors
  ID.class.pos <- match(ID.class, x@data$ID)


  # Calculate Neighbor in Flow Direction ----------------------------------------------
  if(calc.nb.flow == TRUE)
  {

    if(is.null(bb))
    {
      if(quiet == FALSE) print(paste0("Calculate Flow Bounding Box"))

      # check invalid input for flow column
      flow.v <- x[ID.class.pos,]@data[[col.flow]]
      ID.class.pos.excl <- c()

      if(any(is.na(flow.v)) | any(is.null(flow.v)) | any(flow.v == -9999) |  any(flow.v == "-9999"))
      {
        ID.class.pos.excl <- which(is.na(flow.v) | flow.v == "-9999" | is.null(flow.v) | flow.v == -9999)

        # get geometry of class
        x.class <- x[ID.class.pos[-ID.class.pos.excl],]
      } else
      {
        # get geometry of class
        x.class <- x[ID.class.pos,]
      }


      # create bounding boxes in flow direction of class
      # source("../R/module_BoundingBox.R")
      bb <- getBoundingBox(shape = x.class, col.name = col.flow, scale.factor = scale.factor.flow, k.centroid = k.centroid.flow,
                           k = k.flow, centroid = centroid.flow, set.centroid = set.centroid.flow, scale.side = scale.side.flow, quiet = TRUE)
    } else
    {
      bb <- bb
      ID.class.pos.excl <- c()
    }

    if(bb.intersection == "sp")
    {
      intersection.sp <- sp::over(bb, x, returnList = TRUE)
      nb.flow <- lapply(intersection.sp, function(x) as.numeric(row.names(x[1]))) # extract neighbor items as indices of row name
      rm(intersection.sp) # can become very large
    }

    if(bb.intersection == "rgeos")
    {
      nb.flow <- rgeos::gIntersects(bb, x, byid = TRUE, returnDense = FALSE)

      # intersection.rgeos <- rgeos::gBinarySTRtreeQuery(bb, x)
      # intersection.rgeos[unlist(lapply(intersection.rgeos , is.null))] <- 0
      # intersection.rgeos <- as.logical(unlist(intersection.rgeos))
      # nb.flow <- x[intersection.rgeos,]@data$ID
    }

    # add empty neighbor for invalid flow input
    if(length(ID.class.pos.excl) != 0)
    {
      names(nb.flow) <- ID.class.pos[-ID.class.pos.excl]

      # adding empty neighbor hoods to the list
      for(l in 1:length(ID.class.pos.excl))
      {
        nb.flow[length(nb.flow)+1] <- 0
        names(nb.flow)[length(nb.flow)] <- ID.class.pos.excl[l]
      }

      # order list to fit ID.class
      nb.flow <- nb.flow[order(as.numeric(names(nb.flow)))]
    }

    # length(nb.ID.class) == length(nb.flow) should be TRUE
    nb.ID.class.i <- nb[ID.class.pos]

    # counter for following loop
    counter.flow <- 1

  } # end of calculation of neighbors in flow direction



  if(quiet == FALSE) print(paste0("Calculate Statistics ..."))
  if(quiet == FALSE) print(paste0("... number of iterations: " , length(ID.class.pos)))


  #  start looping through class ------------------------------------------------------------------------
  for(i in ID.class.pos)
  {
    ID_i <- x@data[i,]$ID

    # extract neighbors
    nb.i <- c(i, nb[[i]]) # neighbors including selected object
    nb.n <- c(nb[[i]]) # neighbors excluding selected object


    if(class.as.neighbors == FALSE)
    {
      nb.i <- nb.i[!nb.i %in% ID.class.pos] # remove class objects
      nb.i <- c(i, nb.i)

      nb.n <- nb.n[!nb.n %in% ID.class.pos]
    }


    # if flow direction is set
    if(calc.nb.flow == TRUE)
    {
      nb.flow.all.i <- unique(c(i, nb.flow[[counter.flow]]))  # neighbors including selected object
      nb.flow.all.n <- nb.flow.all.i[nb.flow.all.i != i] # neighbors excluding selected object

      if(class.as.neighbors == FALSE)
      {
        nb.flow.all.i <- nb.flow.all.i[!nb.flow.all.i %in% ID.class.pos] # remove class objects
        nb.flow.all.i <- c(i, nb.flow.all.i)

        nb.flow.all.n  <- nb.flow.all.n[!nb.flow.all.n %in% ID.class.pos]
      }

      if(length(nb.flow.all.n) == 0)
      {
        nb.flow.all.n <- 0
      }

      nb.flow.match.i <- Reduce(intersect, list(unique(c(nb.flow[[counter.flow]], i)), nb.i)) # neighbors including selected object
      nb.flow.match.n <- Reduce(intersect, list(nb.flow[[counter.flow]], nb.n)) # neighbors excluding selected object

      if(class.as.neighbors == FALSE)
      {
        nb.flow.match.i <- nb.flow.match.i[!nb.flow.match.i %in% ID.class.pos] # remove class objects
        nb.flow.match.i <- c(i, nb.flow.match.i)

        nb.flow.match.n  <- nb.flow.match.n[!nb.flow.match.n %in% ID.class.pos]
      }


      if(length(nb.flow.match.n) == 0)
      {
        nb.flow.match.n <- 0
      }

      counter.flow <- counter.flow + 1

    } # end if flow direction

    # get data frame of neighbors
    nb.df.i <- as.data.frame(x[nb.i,]) # including selected object
    nb.df.n <- as.data.frame(x[nb.n,]) # excluding selected object

    # flow direction
    if(calc.nb.flow == TRUE)
    {
      nb.df.flow.a.i <- as.data.frame(x[nb.flow.all.i,]) # including selected object
      nb.df.flow.a.n <- as.data.frame(x[nb.flow.all.n,]) # excluding selected object

      nb.df.flow.m.i <- as.data.frame(x[nb.flow.match.i,]) # including selected object
      nb.df.flow.m.n <- as.data.frame(x[nb.flow.match.n,]) # excluding selected object
    }

    value.i <- x[[var[1]]][i] # value of i
    weight.i <- x[[var[2]]][i] # weight of i

    # polygons with NA as value or without any neighbor will be ignored
    if(is.na(value.i) || is.null(nb.n)) # || is.na(nb.n[1]) || nb.n[1] == 0)
    {
      # polygon value is NA
      if(is.na(value.i))
      {
        if(calc.nb.flow == TRUE)
        {
          results.stat <- c(i, rep(NA, 39))
        } else
        {
          results.stat <- c(i, rep(NA, 13))
        }
      }

      # polygon has no neighbor
      if(is.na(nb.n) || is.na(nb.n[1]) || nb.n[1] == 0)
      {
        if(calc.nb.flow == TRUE)
        {
          results.stat <- c(ID_i, rep(NA, 10), 0, rep(NA, 28))
        } else
        {
          results.stat <- c(ID_i, rep(NA, 10), 0, NA, NA)
        }
      }

      # put results in table
      result.table <- rbind(result.table, results.stat)

      # jump loop
      next
    } # end of check NA

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


    # flow features -----------------------------------------------------------
    ####### calculate statistics, adapted by eCognition (2014) reference book, p. 253 ff.
    if(calc.nb.flow == TRUE)
    {
      if(nb.flow.match.n[1] == 0 | nb.flow.all.n[1] == 0)
      {
        # set NAs
        results.stat <- c(ID_i, value.i, weight.i, mean.w, sd.nb, mean.dif, mean.dif.abs, mean.dif.to.hv, mean.dif.to.lw,
                          ratio.nb, sum.nb, num.nb, min.nb, max.nb, rep(NA, 26))

        # put results in table
        result.table <- rbind(result.table, results.stat)

        # jump loop
        next
      }
      # weighted mean and std for all intersected polygons in flow direction - i inclusive class, n exclusive class
      mean.w.flow.a.i <- weighted.mean(nb.df.flow.a.i[[var[1]]], nb.df.flow.a.i[[var[2]]], na.rm = TRUE) # weighted mean
      sd.nb.flow.a.i <- sd(nb.df.flow.a.i[[var[1]]], na.rm = TRUE) # standard deviation
      if(is.na(sd.nb.flow.a.i)){sd.nb.flow.a.i <- 0}

      mean.w.flow.a.n <- weighted.mean(nb.df.flow.a.n[[var[1]]], nb.df.flow.a.n[[var[2]]], na.rm = TRUE) # weighted mean
      sd.nb.flow.a.n <- sd(nb.df.flow.a.n[[var[1]]], na.rm = TRUE) # standard deviation
      if(is.na(sd.nb.flow.a.n)){sd.nb.flow.a.n <- 0}

      # weighted mean and std for all neighbors in flow direction - i inclusive class, n exclusive class
      mean.w.flow.m.i <- weighted.mean(nb.df.flow.m.i[[var[1]]], nb.df.flow.m.i[[var[2]]], na.rm = TRUE) # weighted mean
      sd.nb.flow.m.i <- sd(nb.df.flow.m.i[[var[1]]], na.rm = TRUE) # standard deviation
      if(is.na(sd.nb.flow.m.i)){sd.nb.flow.m.i <- 0}

      mean.w.flow.m.n <- weighted.mean(nb.df.flow.m.n[[var[1]]], nb.df.flow.m.n[[var[2]]], na.rm = TRUE) # weighted mean
      sd.nb.flow.m.n <- sd(nb.df.flow.m.n[[var[1]]], na.rm = TRUE) # standard deviation
      if(is.na(sd.nb.flow.m.n)){sd.nb.flow.m.n <- 0}

      ### important parameter: postive means higher | negative means smaller
      # all flow intersections
      mean.dif.a <-  weighted.mean((value.i - nb.df.flow.a.n[[var[1]]]), nb.df.flow.a.n[[var[2]]], na.rm = TRUE) # Calculates the mean difference between the feature value of an image object and its neighbors of a selected class. Note that the feature values are weighted by area
      mean.dif.abs.a <- abs(mean.dif.a) # absolute value
      mean.dif.to.hv.a <- weighted.mean(value.i - nb.df.flow.a.n[[var[1]]][nb.df.flow.a.n[[var[1]]] > value.i], nb.df.flow.a.n[[var[2]]][nb.df.flow.a.n[[var[1]]] > value.i], na.rm = TRUE) # the mean difference between the feature value of an image object and the feature values of its neighbors, which have higher values than the object itself, weighted by area
      if(is.na(mean.dif.to.hv.a)){ mean.dif.to.hv.a <- 0}
      mean.dif.to.lw.a <- weighted.mean(value.i - nb.df.flow.a.n[[var[1]]][nb.df.flow.a.n[[var[1]]] < value.i], nb.df.flow.a.n[[var[2]]][nb.df.flow.a.n[[var[1]]] < value.i], na.rm = TRUE) # the mean difference between the feature value of an image object and the feature values of its neighbors, which have lower values than the object itself, weighted by area
      if(is.na(mean.dif.to.lw.a)){ mean.dif.to.lw.a <- 0}


      # flow neighbors
      mean.dif.m <-  weighted.mean((value.i - nb.df.flow.m.n[[var[1]]]), nb.df.flow.m.n[[var[2]]], na.rm = TRUE) # Calculates the mean difference between the feature value of an image object and its neighbors of a selected class. Note that the feature values are weighted by area
      mean.dif.abs.m <- abs(mean.dif.m) # absolute value
      mean.dif.to.hv.m <- weighted.mean(value.i - nb.df.flow.m.n[[var[1]]][nb.df.flow.m.n[[var[1]]] > value.i], nb.df.flow.m.n[[var[2]]][nb.df.flow.m.n[[var[1]]] > value.i], na.rm = TRUE) # the mean difference between the feature value of an image object and the feature values of its neighbors, which have higher values than the object itself, weighted by area
      if(is.na(mean.dif.to.hv.m)){ mean.dif.to.hv.m <- 0}
      mean.dif.to.lw.m <- weighted.mean(value.i - nb.df.flow.m.n[[var[1]]][nb.df.flow.m.n[[var[1]]] < value.i], nb.df.flow.m.n[[var[2]]][nb.df.flow.m.n[[var[1]]] < value.i], na.rm = TRUE) # the mean difference between the feature value of an image object and the feature values of its neighbors, which have lower values than the object itself, weighted by area
      if(is.na(mean.dif.to.lw.m)){ mean.dif.to.lw.m <- 0}


      # # # other parameter # # #

      # all flow intersections
      ratio.nb.a <- value.i / weighted.mean(nb.df.flow.a.n[[var[1]]], nb.df.flow.a.n[[var[2]]], na.rm = TRUE)  # ratio between the feature value of an image object and the weighted mean feature value of its neighbors
      sum.nb.a <- sum(nb.df.flow.a.n[[var[1]]], na.rm = TRUE) # Calculates the sum of the feature values of the neighbors
      if(nb.flow.all.n[1] == 0){num.nb.a <- 0}else {num.nb.a <- length(nb.flow.all.n)} # Calculates the number of neighbors
      min.nb.a <- min(nb.df.flow.a.i[[var[1]]], na.rm = TRUE) # Returns the minimum value of the feature values of an image object and its neighbors
      max.nb.a <- max(nb.df.flow.a.i[[var[1]]], na.rm = TRUE) # Returns the maximum value of the feature values of an image object and its neighbors


      # flow neighbors
      ratio.nb.m <- value.i / weighted.mean(nb.df.flow.m.n[[var[1]]], nb.df.flow.m.n[[var[2]]], na.rm = TRUE)  # ratio between the feature value of an image object and the weighted mean feature value of its neighbors
      sum.nb.m <- sum(nb.df.flow.m.n[[var[1]]], na.rm = TRUE) # Calculates the sum of the feature values of the neighbors
      if(nb.flow.match.n[1] == 0){num.nb.m <- 0}else {num.nb.m <- length(nb.flow.match.n)} # Calculates the number of neighbors
      min.nb.m <- min(nb.df.flow.m.i[[var[1]]], na.rm = TRUE) # Returns the minimum value of the feature values of an image object and its neighbors
      max.nb.m <- max(nb.df.flow.m.i[[var[1]]], na.rm = TRUE) # Returns the maximum value of the feature values of an image object and its neighbors

      # results to vector
      results.stat <- c(ID_i, value.i, weight.i, mean.w, sd.nb, mean.dif, mean.dif.abs, mean.dif.to.hv, mean.dif.to.lw,
                        ratio.nb, sum.nb, num.nb, min.nb, max.nb,
                        mean.w.flow.a.i, sd.nb.flow.a.i, mean.w.flow.a.n, sd.nb.flow.a.n, mean.w.flow.m.i, sd.nb.flow.m.i, mean.w.flow.m.n, sd.nb.flow.m.n,
                        mean.dif.a, mean.dif.abs.a, mean.dif.to.hv.a, mean.dif.to.lw.a , mean.dif.m, mean.dif.abs.m, mean.dif.to.hv.m, mean.dif.to.lw.m,
                        ratio.nb.a, sum.nb.a, num.nb.a, min.nb.a, max.nb.a, ratio.nb.m, sum.nb.m, num.nb.m, min.nb.m, max.nb.m)
    } else
    {
      # results to vector
      results.stat <- c(ID_i, value.i, weight.i, mean.w, sd.nb, mean.dif, mean.dif.abs, mean.dif.to.hv, mean.dif.to.lw,
                        ratio.nb, sum.nb, num.nb, min.nb, max.nb)
    }



    # put results in table
    result.table <- rbind(result.table, results.stat)

  } # end loop ID class

  if(calc.nb.flow == TRUE)
  {
    # rename table columns
    names(result.table) <- c("ID", "Value", "Weight","m_wei", "sd_nb", "m_dif", "m_dif_abs", "m_dif_hv",
                             "m_dif_lw", "ratio_nb",  "sum_nb", "num_nb", "min_nb", "max_nb",
                             "m_w_FlAlI", "sd_nb_FlAlI", "m_w_FlAlN", "sd_nb_FlAlN", "m_w_FlMaI", "sd_nb_FlMaI", "m_w_FlMaN", "sd_nb_FlMaN",
                             "m_dif_FlAl", "m_dif_abs_FlAl", "m_dif_hv_FlAl", "m_dif_lw_FlAl" , "m_dif_FlMa", "m_dif_abs_FlMa", "m_dif_hv_FlMa", "m_dif_lw_FlMA",
                             "rat_nb_FlAl", "sum_nb_FlAl", "num_nb_FlAl", "min_nb_FlAl", "max_nb_FlAl", "rat_nb_FlMa", "sum_nb_FlMa", "num_nb_FlMa", "min_nb_FlMa", "max_nb_FlMa")
  } else
  {
    # rename table columns
    names(result.table) <- c("ID", "Value", "Weight","m_wei", "sd_nb", "m_dif", "m_dif_abs", "m_dif_hv",
                             "m_dif_lw", "ratio_nb",  "sum_nb", "num_nb", "min_nb", "max_nb")
  }

  # join output tables
  # return.table <- merge(x@data, result.table, all = TRUE, by = "ID")
  return.table <- result.table

  # return output
  if(calc.nb.flow == TRUE)
  {
    if(return.bb.flow == TRUE)
    {
      # return.bb.flow.spdf <- Reduce(rgeos::gUnion, bb.list)

      # get time of process
      process.time.run <- proc.time() - process.time.start
      if(quiet == FALSE) print(paste0("------ Run of ClassNeighborFunction: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))

      return(list(return.table, bb))

    } else
    {
      # get time of process
      process.time.run <- proc.time() - process.time.start
      if(quiet == FALSE) print(paste0("------ Run of ClassNeighborFunction: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))

      return(return.table)
    }

  } else
  {
    # get time of process
    process.time.run <- proc.time() - process.time.start
    if(quiet == FALSE) print(paste0("------ Run of ClassNeighborFunction: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))

    return(return.table)
  }
}
