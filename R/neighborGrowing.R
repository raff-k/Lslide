#' Growing of neighbors
#'
#' This function grows neighbor \linkS4class{SpatialPolygonsDataFrame} objects together
#'
#' @param spdf \linkS4class{SpatialPolygonsDataFrame}, or full path of a shapefile
#' @param nb object neighborhood based on \code{\link[spdep]{poly2nb}}
#' @param ID.start vector containing \emph{ID}'s where the growing starts. Default: NULL
#' @param ID.class \emph{ID}'s of a specific class in order to allow growing only for class members. Default: NULL
#' @param return.input \emph{spdf} is returned with neighbor information. Default: TRUE
#' @param return.gUnaryUnionNeighbors return united neighbors based on \code{\link[rgeos]{gUnaryUnion}}. Default: TRUE
#' @param quiet no outputs in console. Default: TRUE
#'
#' @return
#' list containing input \linkS4class{SpatialPolygonsDataFrame} with information on growing neighbors, geometry with united neighbors, or both (default)
#'
#' @note
#' \itemize{
#'   \item \linkS4class{SpatialPolygonsDataFrame} input must have an \emph{ID} field with unique numbers
#' }
#'
#'
#' @keywords neighbor growing, spdep, object-oriented image analysis
#'
#' @examples
#' neighborGrowing()
#'
#' @export
neighborGrowing <- function(spdf, nb, ID.start = NULL, ID.class = NULL, return.input = TRUE, return.gUnaryUnionNeighbors = TRUE, quiet = TRUE)
{
  # get start time of process
  process.time.start <- proc.time()

  # create recursive function
  getNeighbors <- function(nb, nb.stillToCheck = c(), nb.check, grow, all.found)
  {
    if (all.found == TRUE)
    {
      # return neighbor indices that grow together
      return(grow)

    } else
    {
      # get neighbor indices that grow together
      grow <- c(grow, nb.check)

      # get new neighborhood
      nb.stillToCheck <- unique(c(nb.stillToCheck, nb[[nb.check]]))

      # get neighbors that still have to be checked
      nb.stillToCheck <- setdiff(nb.stillToCheck, grow)

      # check break condition
      if(length(nb.stillToCheck) == 0 || nb.stillToCheck == 0)
      {
        return(getNeighbors(nb, NULL, NULL, grow, TRUE))
      } else
      {
        return(getNeighbors(nb, nb.stillToCheck, nb.stillToCheck[1], grow, all.found))
      }
    }
  } # end of getNeighbors



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


  if(missing(nb)) # spdep required
  {
    nb_speed_up <- rgeos::gUnarySTRtreeQuery(x) # speed up function for poly2nb
    nb <- spdep::poly2nb(x, queen = TRUE, foundInBox = nb_speed_up) # neighborhood based on queen continuity
  }

  # create new column
  x@data$NeighborGrowing <- NA

  # set expression limit for recursive function
  options(expressions = length(nb) + 1)

  # modify neighborhood depending on input
  if(!is.null(ID.class) & is.null(ID.start))
  {
    # get position of ID in data frame
    ID.class.pos <- match(ID.class, x@data$ID)

    # delete neighbors that are not of the given class (set to 0)
    nb <- lapply(nb, function(x){ if(length(x[x %in% ID.class.pos]) > 0) {x[x %in% ID.class.pos]} else {0} }) # remove neighbors that are not class

    # get neighborhoood for class polygons
    nb.mod <- nb[ID.class.pos] # only class polygons are considered
    names(nb.mod) <- ID.class.pos

  } else if(!is.null(ID.start) & is.null(ID.class))
  {
    ID.start.pos <- match(ID.start, x@data$ID)
    nb.mod <- nb[ID.start.pos] # only start polygons are considered
    names(nb.mod) <- ID.start.pos

  } else if(!is.null(ID.start) & !is.null(ID.class))
  {
    ID.class.pos <- match(ID.class, x@data$ID)
    ID.start.pos <- match(ID.start, x@data$ID)

    # delete neighbors that are not of the given class (set to 0)
    nb <- lapply(nb, function(x){ if(length(x[x %in% ID.class.pos]) > 0) {x[x %in% ID.class.pos]} else {0} }) # remove neighbors that are not class

    nb.mod <- nb[ID.start.pos] # only start polygons are considered
    names(nb.mod) <- ID.start.pos
  } else
  {
    nb.mod <- nb
    names(nb.mod) <- c(1:length(nb))
  }

  # vector containing all checked neighbors
  checked.nb <- c()

  counter.elements <- 1

  # start loop through neighborhoods --------------------------------------------------------------
  for(i in 1:(length(nb.mod)))
  {

    # if(i == 13) stop()

    if(as.numeric(names(nb.mod)[i]) %in% checked.nb)
    {
      # go to next loop
      next
    }


    # initial neighborhood
    nb_i <- nb.mod[[i]]


    # vector contaning elements that grow together
    grow <- c(as.numeric(names(nb.mod)[i]))

    if(nb_i[1] != 0)
    {
      if(length(nb_i) >=2)
      {
        grow <- getNeighbors(nb = nb, nb.stillToCheck =  nb_i[-1], nb.check = nb_i[1], grow = grow, all.found = FALSE)
      } else
      {
        grow <- getNeighbors(nb = nb, nb.check = nb_i[1], grow = grow, all.found = FALSE)
      }
    }

    # set number for growing neighbors that belong together in column
    x@data[grow,]$NeighborGrowing <- rep(counter.elements, length(grow))

    # increment counter
    counter.elements <- counter.elements + 1

    # set checked neighbors to vector
    checked.nb <- unique(c(checked.nb, grow))

  } # end of loop through neighbors
  # # #

  if(is.null(ID.start) & is.null(ID.class))
  {
    # check results
    if(any(is.na(x@data$NeighborGrowing)) == TRUE)
    {
      print("Take care! Something went wrong. Data contains NA values!")
    }
  }

  # build union of neighbors
  if(return.gUnaryUnionNeighbors == TRUE)
  {
    # dissolve polygons by NeighborGrowing ID
    dis.sp <- rgeos::gUnaryUnion(spgeom = x, id = x@data$NeighborGrowing)

    # create data table
    df.dis <- data.frame(unique(x@data$NeighborGrowing))

    if(any(is.na(df.dis)))
    {
      df.dis <- na.omit(df.dis)
    }

    colnames(df.dis) <- "NG"
    row.names(df.dis) <- row.names(dis.sp)

    # create SpatialPolygonsDataFrame
    dis.spdf <- SpatialPolygonsDataFrame(dis.sp, df.dis)
  }


  # get time of process
  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) print(paste0("------ Run of NeighborGrowing: " , process.time.run["elapsed"][[1]]/60, " Minutes ------"))


  # return results depeding of boolean variables set in the function
  if(return.input == TRUE & return.gUnaryUnionNeighbors == TRUE)
  {
    return(list(x, dis.spdf))
  }

  if(return.input == TRUE & return.gUnaryUnionNeighbors == FALSE)
  {
    return(x)
  }

  if(return.input == FALSE & return.gUnaryUnionNeighbors == TRUE)
  {
    return(dis.spdf)
  }
}
