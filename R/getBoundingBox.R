#' Rotated and scaled bounding box
#'
#' This function creates a bounding box for every SpatialPolygonsDataFrame object.
#' Moreover, the bounding box is rotated and scaled in a direction given by the user. \cr \cr
#'
#' The calculation is performed by trigonometric functions. \cr
#' R takes angles in the following anti-clockwise sense: \cr EAST 0°, NORTH 90°, WEST 180°, SOUTH 270°. \cr \cr
#' However, the angles in the output have the following clockwise sense: \cr
#' NORTH 0°, EAST 90°, SOUTH 180°, WEST 270°.
#
#'
#' @param spdf SpatialPolygonsDataFrame. SPDF must contains a valid \emph{ID} field and a field containing direction angles
#' @param col.name name or number of column containing angles in degree
#' @param scale.factor scale factor for bounding box. This can also be vector with two elements: 1. for first side (default: small), and 2. second side (default: longside). Default: 2
#' @param k factor controlling the distance between first two coordinates related to the side and centroid. Default: 2
#' @param k.centroid.flow same as \emph{k}, used when \emph{centroid} is set to FALSE
#' @param projection projection of \emph{spdf} in CRS-format. If \emph{spdf} is projected, projection is read from \emph{spdf}. Default: NA
#' @param centroid centroid is taken for the bounding box. Default: FALSE
#' @param set.centroid set intial starting point of bounding box with respect to the centroid. Default: "inverse" (alternative: "direction")
#' @param scale.side side which is scaled in direction. Default: "long" (alternative: small)
#' @param quiet show output on console. Default: FALSE
#'
#' @keywords bounding box, scaled, rotated, direction
#'
#' @note
#' \itemize{
#'  \item SpatialPolygonsDataFrame must have a valid \emph{ID} field and a column specifying
#'         the angles in which the bounding box is rotated
#'  \item Angles will be transformed to clockwise sense: NORTH 0°, EAST 90°, SOUTH 180°, WEST 270°.
#'  }
#'
#' @examples
#' getBoundingBox()
#'
#'
#' @export
getBoundingBox <- function(spdf, col.name, scale.factor = 2, k = 2, projection = NA, centroid = FALSE, set.centroid = "inverse", k.centroid = 2, scale.side = "long", quiet = TRUE)
{

  # get start time of process
  process.time.start <- proc.time()

  # get angle out of data frame and remove NA's
  angle <- spdf@data[, col.name]
  angle <- angle[!is.na(angle)]

  # get amount of polygons in spdf
  len <- length(spdf@polygons)

  # check if both length are equal
  if(length(angle) != len)
  {
    stop("Length of angles and polygons differ!")
  }

  # get projection of spatial objects
  if(is.projected(spdf))
  {
    projection <- proj4string(spdf)
  }

  # check scale factors
  if(length(scale.factor) == 1)
  {
    f.scale.factor <- scale.factor
    s.scale.factor <- scale.factor

  } else if(length(scale.factor) == 2)
  {
    f.scale.factor <- scale.factor[1]
    s.scale.factor <- scale.factor[2]

  } else
  {
    stop("Invalid input of scale factor!")
  }

  # empty list for SpatialPolygons
  sp_list <- list()

  ### start creation of bounding boxes rotated and scaled in flow direction
  # start loop through polygons ---------------------------------------------
  for(i in 1:len)
  {

    if(is.na(angle[i]) | angle[i] == -9999 |  angle[i] == "-9999")
    {
        stop("Invalid angle input")
    }

    # The negative on angle deals with the fact that we are changing from counterclockwise to clockwise.
    # The +90 deals with the offset of ninety degrees. And lastly we need to mod by 360 to keep our angle in the desired range
    # Angles in R are counterclockwise from E:0|360, N 90, W: 180, S: 270
    angle.rev <- (-angle[i] + 90) %% 360


    ### create bounding box
    # get information out of bounding box ---------------------------------------------
    bbox <- sp::bbox(spdf[i,])

    # retrieve points and length out of bounding box
    xMin <- bbox[1,1]
    yMin <- bbox[2,1]
    xMax <- bbox[1,2]
    yMax <- bbox[2,2]

    yLength <- yMax - yMin
    xLength <- xMax - xMin

    if(yLength >= xLength)
    {
      l.side <- yLength
      s.side <- xLength
    } else
    {
      l.side <- xLength
      s.side <- yLength
    }

    # create coordinates
    x <- c(xMin, xMax, xMax, xMin)
    y <- c(yMin, yMin, yMax, yMax)

    # l1 = SpatialLines(list(Lines(list(Line(cbind(x[1:2], y[1:2]))), ID = "l1")))
    # l2 = SpatialLines(list(Lines(list(Line(cbind(x[2:3], y[2:3]))), ID = "l2")))
    # rgeos::gLength(l1)
    # rgeos::gLength(l2)

    ### create SpatialPolygon
    bbox.sp <- SpatialPolygons(list(Polygons(list(Polygon(cbind(x,y))), paste("p", i))))
    proj4string(bbox.sp) <- projection

    ### get center point of bbox
    center.sp <- sp::coordinates(bbox.sp) # center of bbox
    # center.sp.r <- rgeos::gCentroid(bbox.sp) # center by rgeos
    # center.sp <- sp::coordinates(spdf[i,]) # center of polygon

    ### rotate and scale
    # start of rotation and scale ---------------------------------------------
    # retrieve inverse angles to set the center more to the polygon (against flow direction)

    if(centroid == FALSE)
    {
      if(set.centroid ==  "inverse")
      {
        angle.inv <- (angle.rev - 180 + 360) %% 360
        pC.x <- center.sp[[1]] + s.side/k.centroid * cos(angle.inv *pi/180) # from center in rotational direction
        pC.y <- center.sp[[2]] + s.side/k.centroid * sin(angle.inv *pi/180) # from center in rotational direction
        center <- c(pC.x, pC.y)

      } else if(set.centroid == "direction")
      {
        pC.x <- center.sp[[1]] + s.side/k.centroid * cos(angle.rev *pi/180) # from center in rotational direction
        pC.y <- center.sp[[2]] + s.side/k.centroid * sin(angle.rev *pi/180) # from center in rotational direction
        center <- c(pC.x, pC.y)

      } else
      {
        stop("Wrong Input fpr set.centroid ...")
      }
    }


    if(centroid == TRUE)
    {
      center <- c(center.sp[[1]], center.sp[[2]])
    }

    # retrieve perpendicular angles for small side (+-90 deg) --- keep them in 0 - 360 deg by +360 and mod 360
    angle.per.1 = (angle.rev - 90 + 360) %% 360
    angle.per.2 = (angle.rev + 90 + 360) %% 360

    if(scale.side == "long")
    {
      # get small side coordinates first
      p1.x <- center[[1]] + s.side/k * f.scale.factor * cos(angle.per.1*pi/180)
      p1.y <- center[[2]] + s.side/k * f.scale.factor * sin(angle.per.1*pi/180)

      p2.x  <- center[[1]] + s.side/k * f.scale.factor * cos(angle.per.2*pi/180)
      p2.y  <- center[[2]] + s.side/k * f.scale.factor * sin(angle.per.2*pi/180)

      # get long side now
      p3.x <- p1.x + l.side * s.scale.factor * cos(angle.rev*pi/180)
      p3.y <- p1.y + l.side * s.scale.factor * sin(angle.rev*pi/180)

      p4.x <- p2.x + l.side * s.scale.factor * cos(angle.rev*pi/180)
      p4.y <- p2.y + l.side * s.scale.factor * sin(angle.rev*pi/180)
    }


    if(scale.side == "small")
    {
      # get long side coordinates first
      p1.x <- center[[1]] + l.side/k * f.scale.factor * cos(angle.per.1*pi/180)
      p1.y <- center[[2]] + l.side/k * f.scale.factor * sin(angle.per.1*pi/180)

      p2.x  <- center[[1]] + l.side/k * f.scale.factor * cos(angle.per.2*pi/180)
      p2.y  <- center[[2]] + l.side/k * f.scale.factor * sin(angle.per.2*pi/180)

      # get small side now
      p3.x <- p1.x + s.side * s.scale.factor * cos(angle.rev*pi/180)
      p3.y <- p1.y + s.side * s.scale.factor * sin(angle.rev*pi/180)

      p4.x <- p2.x + s.side * s.scale.factor * cos(angle.rev*pi/180)
      p4.y <- p2.y + s.side * s.scale.factor * sin(angle.rev*pi/180)
    }


    ### create SpatialPolygon from created points | x: xMin, xMax, xMax, xMin   y: yMin, yMin, yMax, yMax
    bbox.rot.sp <- SpatialPolygons(list(Polygons(list(Polygon(cbind(c(p1.x, p3.x, p4.x, p2.x), c(p1.y, p3.y, p4.y, p2.y)))), paste("p", i))))
    proj4string(bbox.rot.sp) <- projection

    # put spatial polygon into list
    sp_list[[i]] <- bbox.rot.sp
  } # end loop

  # create SpatialPolygonsDataFrame
  # create SpatialPolygonsDataFrame from created SpatialPolygons ---------------------------------------------
  sp_list_joined <- SpatialPolygons(lapply(sp_list, function(x){x@polygons[[1]]}))
  spdf <- SpatialPolygonsDataFrame(Sr = sp_list_joined, data = spdf@data["ID"], match.ID = FALSE)
  proj4string(spdf) <- projection

  # get running time of process
  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) print(paste0("------ Run of getBoundingBox: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))


  return(spdf)
} # end function
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #







# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# HOW TO USE FUNCTIONS AND OTHERS----------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# set working directory
# # setwd("E:/Masterarbeit/R/PCA Example")
# setwd("E:/Masterarbeit/Output/Representative_Area")
#
# # # load spdffile
# # spdf <- "scarp.shp"
# spdf <- "Landslide_Example.shp"
# spdf_sf <- st_read(dsn = getwd(), layer = file_path_sans_ext(basename(spdf)), stringsAsFactors = FALSE, quiet = TRUE)
# spdf_sp <- as(spdf_sf, "Spatial") # transform simple feature to spatial
#
# # # get bounding boxes
# bbox <- getBoundingBox(spdf_sp, col.name = 8)
#
# # plot SPDFs
# plot(spdf_sp)
# plot(bbox, col = "red", add = T, pch = 1)
#
# # write spdffile
# sf <- st_as_sf(bbox) # transform spatial to simple feature
# st_write(sf, "scarpBBox.shp", quiet = TRUE)
