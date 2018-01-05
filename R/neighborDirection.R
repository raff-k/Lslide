#' Calculation of the direction of every object to its neighbors
#'
#' This function neighborDirection calculate the direction in degree to its surrounding neighbors or a selection
#'
#' @param spdf \linkS4class{SpatialPolygonsDataFrame}, or full path of a shapefile
#' @param spdf.bb \linkS4class{SpatialPolygonsDataFrame} from which the bounding boxes are created
#' @param col.name column name for angles (in degree, N = 0, E = 90, S = 180, W = 270)
#' @param modus "bb" for using a boundig box or "nb" for using object neighborhood \emph{nb}. Modus \emph{bb} returns data.frame based on spdf. Modus \emph{nb} returns list based on nb.
#' @param tol tolerance for select neighbor of specific range. Default: 360
#' @param sp.PoS \linkS4class{SpatialPoints} to which the direction is calculated. As standard they are calculated from spdf input. Default: NULL
#' @param nb object neighborhood based on \code{\link[spdep]{poly2nb}}, for modus "nb. Default: Null
#' @param bb bounding box for modus "bb". Default: Null
#' @param quiet no outputs in console. Default: TRUE
#'
#' @return
#' In the case of modus \emph{nb} a list is returned. If modus \emph{bb} is set, a data.frame based on the input spdf is returned
#'
#' @note
#' \itemize{
#'   \item \linkS4class{SpatialPolygonsDataFrame} input must have an \emph{ID} field with unique numbers
#' }
#'
#'
#' @keywords neighbor direction, object-oriented image analysis
#'
#'
#' @export
neighborDirection <- function(spdf, col.name, modus = "nb", tol = 360, spdf.bb = NULL, sp.PoS = NULL, nb = NULL, bb = NULL, quiet = TRUE, ...)
{

  # get start time of process
  process.time.start <- proc.time()

  if(is.null(sp.PoS ))
  {
    if(quiet == FALSE){cat("Point on surface...\n")}
    sp.PoS <- rgeos::gPointOnSurface(spgeom = spdf, byid = TRUE)
    spdf.PoS <- SpatialPointsDataFrame(sp.PoS, data.frame(ID = row.names(spdf)))
  }

  if(modus == "bb" & is.null(bb))
  {
    if(is.null(spdf.bb)){stop("Input for spdf.bb is missing!\n")}

    if(quiet == FALSE){cat("Creation of bounding box...\n")}
    bb <- getBoundingBox(shape = spdf.bb, col.name = col.name, ...)
  } else {
    spdf.bb <- bb
  }

  if(modus == "nb" & is.null(nb))
  {
    if(quiet == FALSE){cat("Creation of neighborhood...\n")}
    nb_speed_up <- rgeos::gUnarySTRtreeQuery(spdf) # speed up function for poly2nb
    nb <- spdep::poly2nb(spdf, queen = TRUE, foundInBox = nb_speed_up) # neighborhood based on queen continuity
  }


  # initialize variables for loop
  if(modus == "bb"){obj.inTol.T <- c()} # variable containing information about which neighbors are in tolerance to angle

  obj.nb.angle <- list()

  # get object that intersect with class object
  if(modus == "bb")
    {
    if(quiet == FALSE){cat("Intersect geometries...\n")}
    obj.inter <- rgeos::gIntersects(spgeom1 = spdf.bb, spgeom2 = spdf, byid = TRUE, returnDense = FALSE)


    # index.class <- as.numeric(row.names(spdf[which(spdf@data$ID == spdf.bb@data$ID),]))
    # xy.class <- sp::coordinates(sp.PoS[index.class,])
    index.class <- match(spdf.bb@data$ID, spdf@data$ID)
    xy.class <- sp::coordinates(sp.PoS[index.class,])

    xy.class.angle <- spdf[index.class,]@data[[col.name]]

    # intersected objects and remove classes
    obj.inter.unlCl <- unique(unlist(obj.inter))
    obj.inter.unlCl <- setdiff(obj.inter.unlCl, index.class)

    # get subset of points intersecting the class bounding box
    obj.inter.pts <- sp.PoS[obj.inter.unlCl,]

    # extract cooridates of points
    xy.obj.inter <- sp::coordinates(obj.inter.pts)

    # get amount of polygons in shape
    run <- length(obj.inter)
  }

  if(modus == "nb")
  {
    # align variable names
    obj.inter <- nb
    index.class <- as.numeric(row.names(spdf))
    xy.class.angle <- spdf@data[[col.name]]
    xy.class <- sp::coordinates(sp.PoS)

    # extract cooridates of points
    xy.obj.inter <- sp::coordinates(spdf.PoS)

    # get amount of polygons in shape
    run <- length(spdf)
  }

  # start loop ------
  if(quiet == FALSE){cat("Start loop...\n")}

  for(i in 1:run)
  {
    # print(i)
    if(quiet == FALSE)
    {
      if(i == round(run * 0.25)) {cat("... 25 % done\n")}
      if(i == round(run * 0.5)) {cat("... 50 % done\n")}
      if(i == round(run * 0.75)) {cat("... 75 % done\n")}
    }

    # get object that intersect with class object
    # obj.inter <- unique(unlist(rgeos::gIntersects(spgeom1 = spdf.bb[i,], spgeom2 = spdf, byid = TRUE, returnDense = FALSE)))

    # remove class object itself
    # index.class <- as.numeric(row.names(spdf[which(spdf@data$ID == spdf.bb[i,]@data$ID),]))
    # obj.inter.sub <- setdiff(obj.inter[[i]], index.class)
    obj.inter.i <- obj.inter[[i]]
    obj.inter.i <- setdiff(obj.inter.i, index.class[i])

    if(!is.null(obj.inter.i) & xy.class.angle[i] != -9999 & !is.na(xy.class.angle[i]))
    {
      # get cooridante from actual i-class
      # xy.class <- sp::coordinates(sp.PoS[index.class,])
      # xy.class.angle <- spdf[index.class,]@data[col.name][1,1]
      xy.class.angle.i <- xy.class.angle[i]

      # get subset of points intersecting the class bounding box
      # obj.inter.pts <- sp.PoS[obj.inter,]

      # extract cooridates of points
      # xy.obj.inter <- sp::coordinates(obj.inter.pts)
      xy.obj.inter.i <- xy.obj.inter[match(obj.inter.i, as.numeric(row.names(xy.obj.inter))),]

      if(class(xy.obj.inter.i ) == "matrix")
      {
        # substract class cooridante from matrix
        xy.obj.inter.i[,1] <- xy.obj.inter.i[,1] - xy.class[i,1] # x coordiante
        xy.obj.inter.i[,2] <- xy.obj.inter.i[,2] - xy.class[i,2] # y coordiante

        xy.obj.angle.i <- (atan2(x = xy.obj.inter.i[, 1], y = xy.obj.inter.i[, 2]) * ((-180)/pi) + 90) %% 360


      } else{
        xy.obj.inter.i[1] <- xy.obj.inter.i[1] - xy.class[i,1] # x coordiante
        xy.obj.inter.i[2] <- xy.obj.inter.i[2] - xy.class[i,2] # y coordiante

        xy.obj.angle.i <- (atan2(x = xy.obj.inter.i[1], y = xy.obj.inter.i[2]) * ((-180)/pi) + 90) %% 360

      }




      # # #  calculate mean angle
      # The negative on angle deals with the fact that we are changing from counterclockwise to clockwise.
      # The +90 deals with the offset of ninety degrees. And lastly we need to mod by 360 to keep our angle in the desired range
      # Angles in R are counterclockwise from E:0|360, N 90, W: 180, S: 270
      # get direction
      # xy.obj.angle.i <- (atan2(x = xy.obj.inter.i[, 1], y = xy.obj.inter.i[, 2]) * ((-180)/pi) + 90) %% 360
      # browser()
      # compare with angle direction and selection
      anglediff <- (xy.class.angle.i - xy.obj.angle.i + 180 + 360) %% 360 - 180

      names(anglediff) <- obj.inter.i

      xy.obj.inTol <- sapply(X = anglediff, FUN = function(x, tol)
                        {
                          ifelse(x <= tol && x >= -tol, 1, 0)
                      }, tol = tol)

      if(length(which(xy.obj.inTol == 1)) >= 1)
      {
        # add results to variables
        if(modus == "bb")
        {
          obj.inTol.T.sub <- as.numeric(unique(names(xy.obj.inTol[which(xy.obj.inTol == 1)])))
          # if(length(obj.inTol.T.sub) == 0){obj.inTol.T.sub <- NA}

          obj.inTol.T <- c(unique(obj.inTol.T), obj.inTol.T.sub) # indices!
        }

        xy.obj.angle.sub <- xy.obj.angle.i[which(xy.obj.inTol == 1)]
        obj.nb.angle[[i]] <- list(Object = index.class[i], NeighborDirection = xy.obj.angle.sub)
      } else {
        obj.nb.angle[[i]] <- list(Object = index.class[i], NeighborDirection = NA)
      }


    } else {
      obj.nb.angle[[i]] <- list(Object = index.class[i], NeighborDirection = NA)
      }# end of if

  } # end of loop


  # get time of process
  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE){cat("------ Run of ClassNeighborFunction: " , process.time.run["elapsed"][[1]]/60, " Minutes ------")}

  if(modus == "bb")
  {
    df <- data.frame(Index = row.names(spdf), ID = spdf@data$ID, NbInDir = 0, NbDir = NA)
    df[obj.inTol.T, ]$NbInDir <- 1
    df[unlist(lapply(obj.nb.angle, `[[`, 1)), ]$NbDir <- lapply(obj.nb.angle, `[[`, 2)

    return(df)
  }

  if(modus == "nb")
  {
    return(obj.nb.angle)
  }


} # end of function


# # # # # # # #
# EXAMPLE ----------------------


# setwd("D:/Users/rAVer/Desktop/LSlide New")
#
# source("module_BoundingBox.R")
#
# # save.image("NBDirection.RData")
# load("NBDirection.RData")
#
# pacman::p_load("sf", "sp", "rgdal", "rgeos", "spdep", "maptools")
#
# # spdf <- rgdal::readOGR(dsn = getwd(), layer = "L3_seg_sel_out")
# spdf <- sf::st_read(dsn = paste0(getwd(), "/Input"), layer = "L3_seg_sel_out")
# # seg1.sf <- sf::st_transform(seg1.sf, sp::CRS(paste0("+init=", tolower(epsg.code)))@projargs)
#
# # transform to SP
# spdf <- as(spdf, 'Spatial')
# scarp <- subset(spdf, spdf@data$Class == 1) # 640
#
# # spdf.bb <- getBoundingBox(shape = scarp, col.name = "Flow")
#
#
#
# df <- NeighborDirection(spdf = spdf, spdf.bb = scarp, tol = 45, col.name = "Flow", modus = "bb", quiet = FALSE)
# nb.dir <- NeighborDirection(spdf = spdf, col.name = "Flow", modus = "nb", quiet = FALSE)
#
#
#
#
# # writeSpatialShape(x = bb, fn = paste(getwd(), "/Output/bb.shp", sep = "/"))
# rgdal::writeOGR(obj = spdf.PoS, dsn = paste0(getwd(), "/Output"), layer = "spdf_PoS", driver = "ESRI Shapefile")
# rgdal::writeOGR(obj = bb, dsn = paste0(getwd(), "/Output"), layer = "bb", driver = "ESRI Shapefile")
