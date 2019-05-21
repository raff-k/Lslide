#' Geometry check and repair of SpatialPolygonsDataFrame
#'
#' This function checks the geometry of a \linkS4class{SpatialPolygonsDataFrame} input. Invalid geometries are repaired by the \code{\link[cleangeo]{clgeo_Clean}} function.
#'
#' @param spdf \linkS4class{SpatialPolygonsDataFrame} input
#' @param bug.correct.rgeos correct rgeos bug, see \href{https://gis.stackexchange.com/questions/163445/r-solution-for-topologyexception-input-geom-1-is-invalid-self-intersection-er}{TopologyException} \emph{(last call: 13-04-2017)}. Default: TRUE
#' @param bug.simpl.tol tolerance in rgeos::gSimplify. Default: 0.00001
#' @param bug.buf.width buffer width in rgeos::gBuffer. Default: 0
#' @param quiet no outputs in console. Default: TRUE
#'
#' @return
#' repaired \linkS4class{SpatialPolygonsDataFrame}
#'
#'
#' @keywords invalid geometry, rgeos, cleangeo
#'
#'
#' @export
checkGeometry <- function(spdf, bug.correct.rgeos = TRUE, bug.simpl.tol = 0.00001, bug.buf.width = 0, quiet = TRUE)
{
  # get report
  report <- cleangeo::clgeo_CollectionReport(spdf)

  if(bug.correct.rgeos == TRUE)
  {
    df.spdf <- spdf@data
    spdf <-  rgeos::gSimplify(spdf , tol = bug.simpl.tol)
    spdf  <- rgeos::gBuffer(spdf, byid = TRUE, width = bug.buf.width)

    row.names(spdf) <- row.names(df.spdf )

    spdf <- SpatialPolygonsDataFrame(spdf , df.spdf)
  }


  if(is.element(FALSE, report$valid)) # if a geometrie is not correct, clean it!
  {
    if(quiet == FALSE) print(paste0("Geometry is corrupt... start cleaning..."))
    spdf <- cleangeo::clgeo_Clean(spdf) # clean geometry

    # check again
    report <- cleangeo::clgeo_CollectionReport(spdf)

    # control
    if(is.element(FALSE, report$valid))
    {
      if(quiet == FALSE) print(paste0("Geometry is still corrupt. Take care!"))
    } else
    {
      # double check (rgeos required)
      report.rgeos <- sapply(slot(spdf, "polygons"), function(spdf){gIsValid(SpatialPolygons(Srl = list(spdf)))})

      if(identical(report.rgeos, report$valid) == FALSE)
      {
        if(quiet == FALSE) print(paste0("Geometry sucessfully corrected."))
      } else
      {
        if(quiet == FALSE) print(paste0("Geometry is still corrupt. Take care!"))
      }

    }
  } # end if

  # return geometry
  return(spdf)
}
