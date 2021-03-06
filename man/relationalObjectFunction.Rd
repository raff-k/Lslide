% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/relationalObjectFunction.R
\name{relationalObjectFunction}
\alias{relationalObjectFunction}
\title{Calculation of object relations}
\usage{
relationalObjectFunction(spdf, nb, var, quiet = TRUE)
}
\arguments{
\item{spdf}{\linkS4class{SpatialPolygonsDataFrame}, or full path of a shapefile}

\item{nb}{object neighborhood based on \code{\link[spdep]{poly2nb}}}

\item{var}{vector with number or name defining 1. field of evaluation (i.e. slope) and 2. field of weights (i.e. area). For example: c("Slope", "Area") or c(1,2)}

\item{quiet}{no outputs in console. Default: TRUE}
}
\value{
\linkS4class{data.frame} containing object statistics
}
\description{
This function calculates statistics of \linkS4class{SpatialPolygonsDataFrame} object relations
}
\note{
\itemize{
  \item relational object statistics according to ECOGNITION DEVELOPER (2014: 253-255):
  \itemize{
      \item \emph{m_wei} - mean value of selected features of an object and its neighbors  (of a selected class, see).
        For averaging, the feature values are weighted with the area of the objects.
      \item \emph{sd_nb} - standard deviation of selected features of an object and its neighbors
      It is also possible to calculate this feature with respect to a class, see
      \item \emph{m_dif} - mean difference between the feature value of an object and its neighbors (of a selected class, see)
      The feature values are weighted by the area of the respective objects
      \item \emph{m_dif_abs} - mean absolute difference
      \item \emph{m_dif_hv} - mean difference between the feature value of an object and the feature values of its neighbors (of a selected class, see),
       which have higher values than the object itself. The feature values are weighted by the area of the respective objects
      \item \emph{m_dif_lw} - mean difference between the feature value of an object and the feature values of its neighbors (of a selected class, see),
       which have lower values than the object itself. The feature values are weighted by the area of the respective objects
      \item \emph{ratio_nb} - proportion between the feature value of an object and the mean feature value of its neighbors (of a selected class, see)
        For averaging the feature values are weighted with the area of the corresponding objects
      \item \emph{sum_nb} - sum of the feature values of the neighbors (of a selected class, see )
      \item \emph{num_nb} - number of neighbors of (of a selected class, see )
      \item \emph{min_nb} - minimum value of the feature values of an object and its neighbors (of a selected class, see)
      \item \emph{max_nb} - maximum value of the feature values of an object and its neighbors (of a selected class, see)
  }
  \item ECOGNITION DEVELOPER (2014) Reference Book. Trimble Documentation, München, Germany
}
}
\keyword{analysis}
\keyword{image}
\keyword{object-oriented}
\keyword{objects,}
\keyword{of}
\keyword{relation}
\keyword{spdep,}
