#' Calculation of object relations to a class
#'
#' This function calculates statistics from neighbors of a given class to the given class (input is \linkS4class{SpatialPolygonsDataFrame}).
#' Moreover, it is possible to calculate statistics for neighbors in a specific direction of a class object by setting
#' the parameter for the \emph{getBoundingBox} in the function call.
#' In addition, if the neighbors are attributed to another second class, also a class-to-class relationship can be performed.
#' Furthermore, attributing the first and the second class must not been necessarily hard coded, instead it is also possible to use expressions for the class-to-class relationship.
#' However, here it can happen that some objects are assigned to both classes. As default, these objects are excluded from the statistic, but they can be included by setting the specific argument \emph{(class.as.neighbors)}.
#' Moreover, the function is also capable to compute statistics on angle inputs by setting the function argument of \emph{var.is.angle} to TRUE
#'
#' @param spdf \linkS4class{SpatialPolygonsDataFrame}, or full path of a shapefile
#' @param nb object neighborhood based on \code{\link[spdep]{poly2nb}}. Otherwise, neighborhood is computed by spdep::poly2nb(..., Queens = TRUE)
#' @param var vector with number or name defining 1. field of evaluation (i.e. slope) and 2. field of weights (i.e. area). For example: c("Slope", "Area") or c(1,2)
#' @param var.is.angle variable for statistc represents angles. Other mathematical computation. Default: FALSE
#' @param class.var defination of class. Can be string or integer, or an expression in form of a list. For examples: class = as.integer(1), class = list("> 2.5", "expression"), class = list(c("> 2", "&", "< 7"), "expression")
#' @param class.var2 defination of a second class, see \emph{class.var}. Enables class-to-class relationship. Default: NULL
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
#'   \item object statistics to class according to ECOGNITION DEVELOPER (2014: 253-255, 375-381):
#'   \itemize{
#'       \item see \code{\link{relationalObjectFunction}}
#'       \item \emph{ang_m_cl} - mean angle of an object in relation to its surrounding class objects
#'       \item \emph{ang_d_cl} - mean difference angle of an object in relation to its surrounding class objects
#'       \item \emph{ang_d_cl_a} - absolute mean difference angle
#'       \item \emph{bor_cl_a} - absolute border of an object shared with neighboring objects of a defined class
#'       \item \emph{bor_cl_rel} - relative border of an object shared with neighboring objects of a defined class
#'       \item \emph{dist_cl} - distance of a not-neighbor object of a defined class to a defined class object
#'       \item \emph{flow_nb_cl} - object is located in flow direction (TRUE/FALSE)
#'   }
#'   \item ECOGNITION DEVELOPER (2014) Reference Book. Trimble Documentation, München, Germany
#'   \item see \code{\link{getBoundingBox}}
#'   \item \linkS4class{SpatialPolygonsDataFrame} objects with NA values are ignored
#'   \item \linkS4class{SpatialPolygonsDataFrame} objects without neighbors are ignored
#'   \item \linkS4class{SpatialPolygonsDataFrame} input must have a valid \emph{Class} and \emph{ID} field.
#'      \emph{ID} field must be contained of unique numbers
#' }
#'
#'
#' @keywords relation of objects to a class, spdep, object-oriented image analysis
#'
#' @examples
#'relationalClassFunction()
#'
#' @export
#'
relationalClassFunction <- function(spdf, nb, class.var, class.var2 = NULL, var, var2 = NULL, var.is.angle = FALSE, class.as.neighbors = FALSE, only.neighbors = TRUE, centroid.for.distance = TRUE,
                                    calc.nb.flow = FALSE, col.flow = NULL, bb.intersection = "rgeos", return.bb.flow = FALSE, scale.factor.flow = c(1, 1), bb = NULL,
                                    centroid.flow = FALSE, set.centroid.flow = "inverse", k.centroid.flow = 2, k.flow = 2, scale.side.flow = "small", quiet = TRUE)
{
  # get start time of process
  process.time.start <- proc.time()

  # read shapefile (sf required)
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

  # Get Class 1 ---------------------------------
  ### check class.var attributes
  class.type <- class(x@data$Class) # retrieve type of class column in data frame

  # check input of class.var variable
  if(length(class.var) == 1) # numbers and factors are allowed
  {
    # check if type is the same
    if(class.type == class(class.var))
    {
      # get IDs from related class objects
      ID.class <- subset(x, x@data$Class == class.var)@data$ID
      # ID.class <- x[x@data$Class == class.var,]@data$ID
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

        # set class.var for table
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

        # set class.var for table
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


  # If set: Get Class 2 ---------------------------------
  if(!is.null(class.var2))
  {
    # check if var is set, if not var2 is var
    if(is.null(var2))
    {
      var2 <- var
    }

    ### check class.var2 attributes
    class.type <- class(x@data$Class) # retrieve type of class column in data frame

    # check input of class.var2 variable
    if(length(class.var2) == 1) # numbers and factors are allowed
    {
      # check if type is the same
      if(class.type == class(class.var2))
      {
        if(quiet == FALSE) print(paste0("From Class to Class Relationship is performed"))
        # get IDs from related class objects
        ID.class2 <- subset(x, x@data$Class == class.var2)@data$ID
        From.To.Class <- TRUE

        # ID.class <- x[x@data$Class == class.var2,]@data$ID
      } else # if the type is different. stop!
      {
        stop("Mismatch in class types. Check type of class input variable and Class column!")
      }
    } else # if the input is an expression
    {
      if(length(class.var2) == 2 & class.var2[2] == "expression") # check expression
      {
        class.tmp <- class.var2[1]
        var.tmp <- var2[1]

        # create expression
        if(class(var2[1]) == "character")
        {
          class.var2 <- paste0(as.character(var2[1]), " ", class.tmp)
          var2[1] <- paste0('"', var2[1], '"') # parse string that it fits in expression
        }

        if(length(class.tmp[[1]]) == 1)
        {
          # build expression
          expression <- paste0('x@data[x@data[', as.character(var2[1]), "] ", class.tmp, ",]")

          # set class.var2 for table
          class.var2 <- paste0(as.character(var.tmp), " ", class.tmp)
        } else if(length(class.tmp[[1]]) == 3)
        {
          # split expression
          exp_first <- class.tmp[[1]][1]
          exp_operator <-class.tmp[[1]][2]
          exp_last <- class.tmp[[1]][3]

          # build expression
          expression <- paste0('x@data[x@data[', as.character(var2[1]), "] ", exp_first, " ",
                               exp_operator, ' x@data[', as.character(var2[1]), "] ", " ", exp_last, ",]")

          # set class.var2 for table
          class.var2 <- paste0(as.character(var.tmp), " ", exp_first, " ", exp_operator, " ", exp_last)
        } else
        {
          stop("Wrong expression input.")
        }

        if(quiet == FALSE) print(paste0("From Class to Class Relationship is performed"))
        # get IDs from exression result
        ID.class2 <- eval(parse(text = expression))$ID
        From.To.Class <- TRUE

        if(length(ID.class2) == length(x@data$ID) | length(ID.class2) == 0)
        {
          stop("Expression is invalid. Either no results or all ID's selected.")
        }

        # reset var
        var2[1] <- var.tmp

      } else # if the type is different. stop!
      {
        stop("Class input invalid!")
      }
    }
  }


  # Get Geometries ---------------------------------
  # extract geoemtries by class ---- NAs not permitted in row index
  # x.extr <- x[x@data$ID == ID.class,] # extract polygons
  x.extr <- x[which(x@data$ID %in% ID.class),]
  # x.extr.union <- gUnaryUnion(x[x@data$ID[ID.class],]) # union or merge polygons to one object
  # x.extr.union <- rgeos::gUnaryUnion(x[which(x@data$ID %in% ID.class),])
  x.extr.union <- rgeos::gUnaryUnion(x.extr)


  # modify x.extra.union
  x.extr.union <- gSimplify(x.extr.union, tol = 0.00001)
  x.extr.union <- gBuffer(x.extr.union, byid=TRUE, width = 0)

  if(gIsValid(x.extr.union) == FALSE)
  {
    # check and clean geometry (cleangeo required)
    x.extr.union<- createSPComment(x.extr.union)
    x.extr.union <- checkGeometry(x.extr.union)
  }

  # create spatial lines for intersection later
  x.extr.union.l <- as(x.extr.union, "SpatialLines")

  if(gIsValid(x.extr.union.l) == FALSE)
  {
    # check and clean geometry (cleangeo required)
    x.extr.union.l <- createSPComment(x.extr.union.l)
    x.extr.union.l <- checkGeometry(x.extr.union.l)
  }

  # gLength(x.extr.union.l) == gLength(x.ext.b)  # should be TRUE
  x.ext.b <- rgeos::gBoundary(x.extr.union) # get boundaries of unioned object


  # Set Neighborhood --------------------------------
  if(missing(nb)) # spdep required
  {
    nb_speed_up <- rgeos::gUnarySTRtreeQuery(x) # speed up function for poly2nb
    nb <- spdep::poly2nb(x, queen = TRUE, foundInBox = nb_speed_up) # neighborhood based on queen continuity
  }

  # extract neighbors from ID.class
  ID.class.pos <- match(ID.class, x@data$ID) # get position of IDs in data frame for neighborhood
  nb.class <- unlist(nb[ID.class.pos]) # extract all neighbors from the class objects to a vector

  if(class.as.neighbors == TRUE)
  {
    nb.class <- nb.class[!nb.class %in% 0] # remove objects with zero neighbors
  }

  if(class.as.neighbors == FALSE)
  {
    nb.remove <- c(0, ID.class.pos) # remove the class object itselfes and elements with 0 neighbors
    nb.class <- nb.class[!nb.class %in% nb.remove] # remove objects with zero neighbors
  }

  nb.class.uni <- unique(nb.class) # remove duplicats in neighbor objects

  # stop executing if there is no neighbor detected
  if(!length(nb.class.uni))
  {
    stop("Selected class or expression has no neighbor relationship.")
  }


  # if From Class to Class is set
  if(exists("From.To.Class") && From.To.Class == TRUE)
  {
    if(only.neighbors == FALSE)
    {
      if(quiet == FALSE) warning("By from class to class calculation is the analysis for the whole data set disabled (only.neighbors = TRUE)")

      only.neighbors <- TRUE # analysis for whole data frame is disabled
    }

    ID.class2.pos <- match(ID.class2, x@data$ID) # get position of IDs of class2 in data frame for neighborhood

    # extract neighbors from ID.class that are of class2
    nb.class2.uni <- Reduce(intersect, list(ID.class2.pos, nb.class.uni))

    # IMPORTANT: nb.class
    nb.class2 <- Reduce(intersect, list(ID.class2.pos, nb.class))
  }

  # recheck
  # stop executing if there is no neighbor detected
  if(!length(nb.class.uni))
  {
    stop("Selected class or expression has no neighbor relationship.")
  }



  # Calculate Neighbor in Flow Direction ----------------------------------------------
  if(calc.nb.flow == TRUE)
  {
    if(is.null(bb))
    {
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

    # extract all neighbors from the class objects to a vector
    nb.flow.unlist <- unique(unlist(nb.flow))

    # remove class IDs and 0
    nb.flow.unlist <- nb.flow.unlist[!nb.flow.unlist %in% c(0, ID.class.pos)]

  } # end of calculation of neighbors in flow direction





  # get data frame of neighbors
  # class.df.nb <- nb.df.incl.class <- as.data.frame(x[c(ID.class.pos, nb.class.uni),]) # including class objects
  # nb.df.excl.class <- as.data.frame(x[nb.class.uni,]) # excluding class objects, only neighbors

  # create data frame from class
  class.df <- as.data.frame(x@data[ID.class.pos,]) # data frame including only class objects


  # start operations --------------------------------------
  ### relational class function created in adpation of ecognition (2014) reference book, p. 375 ff.


  ### Number of
  # Number of neighbors of neighbors belonging to the selected class in a certain distance (in pixels) around the image object.
  nb.class.uni.o <- nb.class.uni[order(nb.class.uni)]
  nb.num.nb <- nb[nb.class.uni.o]
  nb.num <- unlist(lapply(nb.num.nb , function(x){length(x[])}))

  nb.class.output <- data.frame(Index = nb.class.uni.o, num_nb = nb.num, stringsAsFactors=FALSE)
  nb.class.output$ID <- x@data$ID[nb.class.output$Index] # neighbor ID


  if(exists("From.To.Class") && From.To.Class == TRUE)
  {
    # create table
    nb.class2.uni.o <- nb.class2.uni[order(nb.class2.uni)]
    nb.num.nb2 <- nb[nb.class2.uni.o]
    nb.num2 <- unlist(lapply(nb.num.nb2 , function(x){length(x[])}))

    nb.class2.output <- data.frame(Index = nb.class2.uni.o, num_nb = nb.num2, stringsAsFactors=FALSE)
    nb.class2.output$ID <- x@data$ID[nb.class2.output$Index] # neighbor ID


    # inner join
    nb.class.output <- merge(nb.class2.output, nb.class.output, by = "Index")

    # modify and rename table
    nb.class.output[c(2,3)] <- NULL
    names(nb.class.output)[2:3] <- c("num_nb", "ID")

    # overwrite nb.class
    nb.class <- nb.class2
  }



  ### Calculate statistics and border relations for neighbor objects of selected class
  if(quiet == FALSE) print(paste0("Calculate Statistics for class neighbors"))
  # Border to
  # The absolute border of an image object shared with neighboring objects of a defined classification.
  # extract neighbors
  x.extr.nb <- x[nb.class.output$Index,]

  # modify extr.nb
  x.extr.nb <- gSimplify(x.extr.nb, tol = 0.00001)
  x.extr.nb <- gBuffer(x.extr.nb, byid = TRUE, width = 0)

  # check geometry
  if(gIsValid(x.extr.nb) == FALSE)
  {
    # check and clean geometry (cleangeo required)
    x.extr.nb <- createSPComment(x.extr.nb)
    x.extr.nb <- checkGeometry(x.extr.nb, quiet = quiet)
  }

  # get back data frame
  row.names(x.extr.nb) <- as.character(nb.class.output$Index)
  x.extr.nb.spdf <- SpatialPolygonsDataFrame(x.extr.nb, x[nb.class.output$Index,]@data)

  # convert to SpatialLines
  x.extr.nb.l <- as(x.extr.nb, "SpatialLines")

  # check geometry
  if(gIsValid(x.extr.nb.l) == FALSE)
  {
    # check and clean geometry (cleangeo required)
    x.extr.nb.l <- checkGeometry(x.extr.nb.l, quiet = quiet)
  }

  # convert to SpatialLinesDataFrame
  row.names(x.extr.nb.l) <- row.names(x[nb.class.output$Index,]@data)
  x.extr.nb.sldf <- SpatialLinesDataFrame(x.extr.nb.l, x[nb.class.output$Index,]@data)

  ### get border length
  # gIntersection
  borTo.intersect <- rgeos::gIntersection(x.extr.nb.sldf, x.extr.union.l, byid = TRUE,
                                          id = row.names(x.extr.nb.sldf), drop_lower_td = TRUE)

  # get absolut border length
  borTo.abs <- rgeos::gLength(borTo.intersect, byid = TRUE) # id is index of id in the original data frame |  in meter

  # check of invalid numbers
  if(any(is.na(borTo.abs)) | any(is.null(borTo.abs)))
  {
    borTo.abs <- ifelse(is.na(borTo.abs), 0, borTo.abs)
    borTo.abs <- ifelse(is.null(borTo.abs), 0, borTo.abs)
  }

  # check of missing numbers, can happens because of Queens neighborhood! (sharing only an edge point)
  if(length(borTo.abs) != length(nb.class.output$Index))
  {
    nb.missing <- setdiff(as.character(nb.class.output$Index), names(borTo.abs))
    nb.add <- rep(0, length(nb.missing))
    names(nb.add) <- nb.missing

    # add missing values to border vales
    borTo.abs <- c(borTo.abs, nb.add)
    borTo.abs <- borTo.abs[order(as.numeric(names(borTo.abs)))]
  }

  # http://stackoverflow.com/questions/20351624/rgeos-gintersection-in-loop-takes-too-long-to-clip-path-networ
  # borTo.abs.sldf <- SpatialLinesDataFrame(borTo.intersect , x[nb.class.output$Index,]@data)
  # rgdal::writeOGR(borTo.abs.sldf, dsn = paste(getwd(), dirname(seg1.sel), sep = "/"), overwrite_layer = TRUE,
  #                 layer = "line", driver = "ESRI Shapefile")



  ### Relative Border to
  # The feature Rel. Border To determines the relative border length an object shares with
  # the image border. It describes the ratio of the shared border length of an image object
  # (with a neighboring image object assigned to a defined class) to the total border length. If
  # the relative border of an image object to image objects of a certain class is 100, the image
  # object is totally embedded in them.
  b <- rgeos::gLength(x.extr.nb.sldf, byid = TRUE)

  borTo.rel <- borTo.abs/b*100 # as ratio 0-1

  if(quiet == FALSE) print(paste0("... finished border calculation"))

  # start loop through neighbors -----------------------

  # initialize empty vectors for result
  # borTo <- c() # absolute border
  # borToRel <- c() # relative border
  w.mean.nb.class <- c() # weighted mean of class
  sd.nb.class <- c() # standard deviation
  mean.dif.nb.class <- c() # Mean difference
  mean.dif.nb.class.abs <- c() # absolute value of Mean difference
  mean.dif.to.hv.nb.class <-c() # # Mean difference to higher values
  mean.dif.to.lw.nb.class  <- c() # Mean difference to lower values
  ratio.nb.class <- c() # ratio
  sum.nb.class <- c() # sum
  num.nb.class <- c() # number
  min.nb.class <- c() # min value
  max.nb.class <- c() # max value


  if(var.is.angle == TRUE)
  {
    angle.mean <- c() # mean angle of all
    angle.mean.class <- c() # mean angle of class
    angle.dif <- c() # difference between class angle and mean angle of class
    angle.dif.abs <- c() # absolute difference between class angle and mean angle of class
  }


  # if flow direction is set
  if(calc.nb.flow == TRUE)
  {
    flow.class <- c() # max value
  }

  # i is position in data frame of x
  for(i in nb.class.output$Index)
  {

    ID_i <- nb.class.output$ID[nb.class.output$Index == i]

    value.nb.i <- x[[var[1]]][i] # value of i
    weight.nb.i <- x[[var[2]]][i] # weight value of i

    # if feature has no data it is ignored
    if(is.null(value.nb.i) || is.na(value.nb.i))
    {
      # fill vector with NAs
      #borTo.abs[which(nb.class.output$Index == i)] <- NA # absolute border
      #borTo.rel[which(nb.class.output$Index == i)] <- NA # relative border
      w.mean.nb.class <- c(w.mean.nb.class, NA) # weighted mean of class
      sd.nb.class <- c(sd.nb.class, NA) # standard deviation
      mean.dif.nb.class <- c(mean.dif.nb.class, NA) # Mean difference
      mean.dif.nb.class.abs <- c(mean.dif.nb.class.abs, NA) # absolute value of Mean difference
      mean.dif.to.hv.nb.class <-c(mean.dif.to.hv.nb.class, NA) # # Mean difference to higher values
      mean.dif.to.lw.nb.class  <- c(mean.dif.to.lw.nb.class, NA) # Mean difference to lower values
      ratio.nb.class <- c(ratio.nb.class , NA) # ratio
      sum.nb.class <- c(sum.nb.class, NA) # sum
      num.nb.class <- c(num.nb.class, NA) # number
      min.nb.class <- c(min.nb.class, NA) # min value
      max.nb.class <- c(max.nb.class, NA) # max value

      # if flow direction is set
      if(calc.nb.flow == TRUE)
      {
        flow.class <- c(flow.class, NA)
      }

      if(var.is.angle == TRUE)
      {
        angle.mean <- c(angle.mean, NA) # mean angle of all
        angle.mean.class <- c(angle.mean.class, NA) # mean angle of class
        angle.dif <- c(angle.dif, NA) # difference between class angle and mean angle of class
        angle.dif.abs <- c(angle.dif.abs, NA) # absolute difference between class angle and mean angle of class
      }


      # jump loop
      next
    }

    # check class.df.nb
    ####### calculate statistics, adapted by eCognition (2014) reference book, p. 253 ff.
    # get neighbors of neighbor objects of selected class
    nb.i <- nb[[i]] # get all neighbors of object
    nb.i <- nb.i[nb.i %in% ID.class.pos] # remove all neighbor objects of i that are NOT class objects - class objects remain!
    value.nb.i.v <- c(value.nb.i, x@data[[var[1]]][nb.i]) # vector with feature values and its class objects
    value.i.class <- c(x@data[[var[1]]][nb.i]) # vector with class objects values of feature i
    weight.nb.i.v <- c(weight.nb.i, x@data[[var[2]]][nb.i]) # vector with feature weight and its class objects weights
    weight.i.class <- c(x@data[[var[2]]][nb.i]) # vector with class objects weights of feature i

    if(var.is.angle == FALSE)
    {
      # statistic calculation
      w.mean.nb.class.temp <- weighted.mean(value.nb.i.v, weight.nb.i.v, na.rm = TRUE) # weighted mean
      if(is.na(w.mean.nb.class.temp)){w.mean.nb.class.temp <- 0}
      w.mean.nb.class <- c(w.mean.nb.class, w.mean.nb.class.temp)

      sd.nb.class.temp <- sd(value.nb.i.v, na.rm = TRUE) # standard deviation of the feature and the class objects
      if(is.na(sd.nb.class.temp)){sd.nb.class.temp <- 0}
      sd.nb.class <- c(sd.nb.class, sd.nb.class.temp)

      ### important parameter: postive means higher | negative means smaller
      mean.dif.nb.class.r <-  weighted.mean((value.nb.i - value.i.class), weight.i.class, na.rm = TRUE) # Calculates the mean difference between the feature value of an image object and its neighbors of a selected class. Note that the feature values are weighted by area
      mean.dif.nb.class <- c(mean.dif.nb.class, mean.dif.nb.class.r)
      mean.dif.nb.class.abs <- c(mean.dif.nb.class.abs, abs(mean.dif.nb.class.r)) # absolute value
      mean.dif.to.hv.nb.class.r <- weighted.mean(value.nb.i - value.i.class[value.i.class > value.nb.i], weight.i.class[value.i.class > value.nb.i], na.rm = TRUE) # the mean difference between the feature value of an image object and the feature values of its neighbors of a specific class, which have higher values than the object itself, weighted by area
      if(is.na(mean.dif.to.hv.nb.class.r)){ mean.dif.to.hv.nb.class <- c(mean.dif.to.hv.nb.class, 0)}else{ mean.dif.to.hv.nb.class <- c(mean.dif.to.hv.nb.class, mean.dif.to.hv.nb.class.r)}
      mean.dif.to.lw.nb.class.r <- weighted.mean(value.nb.i - value.i.class[value.i.class < value.nb.i], weight.i.class[value.i.class < value.nb.i], na.rm = TRUE) # the mean difference between the feature value of an image object and the feature values of its neighbors of a specific class, which have lower values than the object itself, weighted by area
      if(is.na(mean.dif.to.lw.nb.class.r)){ mean.dif.to.lw.nb.class <- c(mean.dif.to.lw.nb.class, 0)}else{ mean.dif.to.lw.nb.class <- c(mean.dif.to.lw.nb.class, mean.dif.to.lw.nb.class.r)}
      ###

      ratio.nb.class <- c(ratio.nb.class, value.nb.i / weighted.mean(value.i.class, weight.i.class, na.rm = TRUE))  # ratio between the feature value of an image object and the weighted mean feature value of its neighbors of a specific class
      sum.nb.class <- c(sum.nb.class, sum(value.i.class, na.rm = TRUE)) # Calculates the sum of the feature values of the neighbors of a specific class
      num.nb.class <- c(num.nb.class, length(nb.i)) # Calculates the number of neighbors of a specific class
      min.nb.class <- c(min.nb.class, min(value.nb.i.v, na.rm = TRUE)) # Returns the minimum value of the feature values of an image object and its neighbors of a specific class
      max.nb.class <- c(max.nb.class, max(value.nb.i.v, na.rm = TRUE)) # Returns the maximum value of the feature values of an image object and its neighbors of a specific class

      # if flow direction is set
      if(calc.nb.flow == TRUE)
      {
        if(is.element(i, nb.flow.unlist) == TRUE)
        {
          flow.class <- c(flow.class, 1)
        } else
        {
          flow.class <- c(flow.class, 0)
        }
      }

    } # end of if var.is.angle == FALSE

    if(var.is.angle == TRUE)
    {
      # # #  calculate mean angle
      # The negative on angle deals with the fact that we are changing from counterclockwise to clockwise.
      # The +90 deals with the offset of ninety degrees. And lastly we need to mod by 360 to keep our angle in the desired range
      # Angles in R are counterclockwise from E:0|360, N 90, W: 180, S: 270
      value.nb.i.tf <- (-value.nb.i + 90) %% 360
      value.nb.i.tf <- value.nb.i.tf *pi/180

      # kick out no data
      value.i.class <- value.i.class[!is.na(value.i.class) & value.i.class != -9999]


      # start calculation of mean angle
      if(length(value.i.class) == 0)
      {
        angle.mean <- c(angle.mean, NA) # mean angle of all
        angle.mean.class <- c(angle.mean.class, NA) # mean angle of class
        angle.dif <- c(angle.dif, NA) # difference between class angle and mean angle of class
        angle.dif.abs <- c(angle.dif.abs, NA) # absolute difference between class angle and mean angle of class
      } else
      {
        # transform to R specific angles
        value.i.class <- (-value.i.class + 90) %% 360
        value.i.class <- value.i.class *pi/180

        angle.m <- ((atan2(mean(sin(c(value.nb.i.tf, value.i.class))), mean(cos(c(value.nb.i.tf, value.i.class)))) * (-180)/pi) + 90) %% 360
        angle.class <- ((atan2(mean(sin(value.i.class)), mean(cos(value.i.class))) * (-180)/pi) + 90) %% 360
        angle.d <- value.nb.i - angle.class
        angle.d.a <- abs(angle.d)

        angle.mean <- c(angle.mean, angle.m) # mean angle of all
        angle.mean.class <- c(angle.mean.class, angle.class) # mean angle of class
        angle.dif <- c(angle.dif, angle.d) # difference between class angle and mean angle of class
        angle.dif.abs <- c(angle.dif.abs, angle.d.a)
      }



      # fill vector with NAs
      w.mean.nb.class <- c(w.mean.nb.class, NA) # weighted mean of class
      mean.dif.nb.class <- c(mean.dif.nb.class, NA)
      sd.nb.class <- c(sd.nb.class, NA) # standard deviation
      mean.dif.nb.class.abs <- c(mean.dif.nb.class.abs, NA) # absolute value of Mean difference
      mean.dif.to.hv.nb.class <-c(mean.dif.to.hv.nb.class, NA) # # Mean difference to higher values
      mean.dif.to.lw.nb.class  <- c(mean.dif.to.lw.nb.class, NA) # Mean difference to lower values
      ratio.nb.class <- c(ratio.nb.class , NA) # ratio
      sum.nb.class <- c(sum.nb.class, NA) # sum
      num.nb.class <- c(num.nb.class, NA) # number
      min.nb.class <- c(min.nb.class, NA) # min value
      max.nb.class <- c(max.nb.class, NA) # max value

      # if flow direction is set
      if(calc.nb.flow == TRUE)
      {
        flow.class <- c(flow.class, NA)
      }
    }


  } # end loop


  # put results to table
  nb.class.output$bor_cl_a <- borTo.abs
  nb.class.output$bor_cl_rel <- borTo.rel
  nb.class.output$sd_nb_cl <- sd.nb.class

  if(var.is.angle == TRUE)
  {
    nb.class.output$ang_m <- angle.mean
    nb.class.output$ang_m_cl <- angle.mean.class
    nb.class.output$ang_d_cl <- angle.dif
    nb.class.output$ang_d_cl_a <- angle.dif.abs
  }
  nb.class.output$m_w <-  w.mean.nb.class
  nb.class.output$m_d_nb_cl <- mean.dif.nb.class
  nb.class.output$m_d_nb_cl_a <- mean.dif.nb.class.abs
  nb.class.output$m_d_hv_nb_cl <- mean.dif.to.hv.nb.class
  nb.class.output$m_d_lw_nb_cl <- mean.dif.to.lw.nb.class
  nb.class.output$rat_nb_cl <- ratio.nb.class
  nb.class.output$sum_nb_cl <- sum.nb.class
  nb.class.output$num_nb_cl <- num.nb.class
  nb.class.output$min_nb_cl <- min.nb.class
  nb.class.output$max_nb_cl <- max.nb.class

  if(calc.nb.flow == TRUE)
  {
    nb.class.output$flow_nb_cl <- flow.class
  }

  if(quiet == FALSE) print(paste0("... finished neighbor to class calculation"))


  if(only.neighbors == FALSE & var.is.angle == FALSE)
  {
    #### Calculate statistics and distances for all objects to class objects
    # create empty table
    result.table <- data.frame(ID = integer(),
                               sd.class  = double(),
                               mean.dif.class = double(),
                               mean.dif.abs.class = double(),
                               mean.dif.to.hv.class = double(),
                               mean.dif.to.lw.class = double(),
                               stringsAsFactors=FALSE)

    results.stat <- c() # empty vector for loop results

    # loop through all objects and get shorest distance to class objects
    if(quiet == FALSE) print(paste0("Calculate Statistics for all objects"))

    ### Distance to
    # The distance (in pixels) of the image object concerned to the closest image ob-
    # ject assigned to a defined class.
    # modify extr.nb
    x.simply <- gSimplify(x, tol = 0.00001)
    x.simply <- gBuffer(x.simply, byid = TRUE, width = 0)

    # check geometry
    if(gIsValid(x.simply) == FALSE)
    {
      # check and clean geometry (cleangeo required)
      x.simply <- createSPComment(x.simply)
      x.simply <- checkGeometry(x.simply, quiet = quiet)
    }

    if(centroid.for.distance == TRUE)
    {
      x.simply.centroid <- rgeos::gCentroid(x.simply, byid = TRUE)

      # calc distances
      dist.ToClass <- rgeos::gDistance(x.simply.centroid, x.extr.union, byid = TRUE) # in meter, objects itselfes and neighbors have the value 0
    }

    if(centroid.for.distance == FALSE)
    {
      # calc distances
      dist.ToClass <- rgeos::gDistance(x.simply, x.extr.union, byid = TRUE) # in meter, objects itselfes and neighbors have the value 0
    }

    dist.ToClass.v <- as.vector(dist.ToClass)

    if(length(dist.ToClass.v) != length(x))
    {
      print("Something went wrong in calculation of distances. Length differ to original data set!")
    }

    if(quiet == FALSE) print(paste0("... finished distance calculation"))

    # start looping through objects
    # i is first position in data frame of x, then second and so on
    for(i in 1:length(x@data$ID))
    {
      ID_i <- x@data[i,]$ID

      value.i <- x[[var[1]]][i] # value of i of first object

      # check if ID == ID.class OR value.i is NA then go forward | NAs to dist
      # remember: ID.class.pos is Index not the actual ID
      if(is.element(ID_i, x@data$ID[ID.class.pos]) | is.na(value.i))
      {
        # put NAs into vector
        results.stat <- c(ID_i, rep(NA, 5))

        # put results in table
        result.table <- rbind(result.table, results.stat)

        # jump loop
        next
      }

      ### Standard deviation
      sd.class.i <- c(value.i, class.df[[var[1]]]) # put class values and actual i value into a vector
      sd.class <- sd(sd.class.i , na.rm = TRUE) # calculate standard deviation
      if(is.na(sd.class)){sd.class <- 0}

      ### Mean Difference to
      # The mean difference of the layer L mean value of the image object concerned to the layer L
      # mean value of all image objects assigned to a defined class
      mean.dif.class <-  weighted.mean((value.i -  class.df[[var[1]]]),  class.df[[var[2]]], na.rm = TRUE) # Calculates the mean difference between the feature value of an image object and its neighbors of a selected class. Note that the feature values are weighted by area
      mean.dif.abs.class <- abs(mean.dif.class) # absolute value
      mean.dif.to.hv.class <- weighted.mean(value.i - class.df[[var[1]]][class.df[[var[1]]] > value.i], class.df[[var[2]]][class.df[[var[1]]] > value.i], na.rm = TRUE) # the mean difference between the feature value of an image object and the feature values of a selected class, which have higher values than the object itself, weighted by area
      if(is.na(mean.dif.to.hv.class)){ mean.dif.to.hv.class <- 0}
      mean.dif.to.lw.class <- weighted.mean(value.i - class.df[[var[1]]][class.df[[var[1]]] < value.i], class.df[[var[2]]][class.df[[var[1]]] < value.i], na.rm = TRUE) # the mean difference between the feature value of an image object and the feature values of a selected class, which have lower values than the object itself, weighted by area
      if(is.na(mean.dif.to.lw.class)){ mean.dif.to.lw.class <- 0}

      results.stat <- c(ID_i, sd.class, mean.dif.class, mean.dif.abs.class, mean.dif.to.hv.class, mean.dif.to.lw.class) # loop results into vector

      # put results into table
      result.table <- rbind(result.table, results.stat)

    } # end loop

    if(quiet == FALSE) print(paste0("... finished all objects to class calculation"))

    result.table$dist <-  dist.ToClass.v

    # rename table columns
    names(result.table) <- c("ID", "sd_cl", "m_dif_cl", "m_dif_abs_cl",
                             "m_dif_hv_cl", "m_dif_lw_cl", "dist_cl")
  } # of check only.neighbors

  # create class data frame
  return.table <- as.data.frame(x@data$ID)
  names(return.table) <- "ID"
  return.table$Cl <- x@data$Class
  return.table$ClRel <- NA
  return.table$ClRel[nb.class.uni] <- class.var # real neighbors
  return.table$ClRel[ID.class.pos] <- "Selection" # selected polygons
  return.table$ClRel[is.na(return.table$ClRel)] <- "Other" # not relevant polygons

  if(exists("From.To.Class") && From.To.Class == TRUE)
  {
    return.table$ClTo <- NA
    return.table$ClTo[ID.class2.pos] <- "Selection" # selected polygons
    return.table$ClTo[nb.class.uni] <- class.var2 # real neighbors
    return.table$ClTo[is.na(return.table$ClTo)] <- "Other" # not relevant polygons
  }

  return.table$Val <- x@data[[var[1]]]
  # colnames(return.table)[4]
  return.table$Wei <- x@data[[var[2]]]

  # join output tables
  if(only.neighbors == FALSE)
  {
    return.table <- merge(return.table, result.table, all = TRUE, by = "ID")
  }

  return.table <- merge(return.table, nb.class.output, all = TRUE, by = "ID")


  # get time of process
  process.time.run <- proc.time() - process.time.start
  print(paste0("------ Run of RelationalClassFunction: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))

  # return output
  if(calc.nb.flow == TRUE & return.bb.flow == TRUE)
  {
    return(list(return.table, bb))
  } else
  {
    return(return.table)
  }

}

