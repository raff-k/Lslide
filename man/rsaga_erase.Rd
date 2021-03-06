% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sf_utils.R
\name{rsaga_erase}
\alias{rsaga_erase}
\title{Erase one geometry from another using Saga GIS}
\usage{
rsaga_erase(x, y, method = "1", split = "0", attributes = "1",
  env.rsaga = RSAGA::rsaga.env(), check.geom = TRUE, quiet = TRUE)
}
\arguments{
\item{x}{object of class sf. First element: Should be either line or polygon}

\item{y}{object of class sf. Second element: Always polygon.}

\item{method}{method of erase. Either "1": , or "2": Line-Polygon Intersection. Default: "1"}

\item{split}{Set to "1", if multi-part polygons should be splitted to single-part polygons. Default: "0"}

\item{attributes}{attributes inherited to intersection result. [0] polygon, [1] line, [2] line and polygon. Default: "1"}

\item{env.rsaga}{SAGA GIS environemnt. Default: RSAGA::rsaga.env()}

\item{check.geom}{If set to TRUE then geometry is checked with sf::st_is_valid(). If there are invalid geometries, geometries are repaired using lwgeom::st_make_valid(). Default: TRUE}
}
\value{
Geometry of class sfc
}
\description{
This function erase one geometry from another. The projection must be identical.
}
\keyword{erase}
\keyword{feature,}
\keyword{simple}
