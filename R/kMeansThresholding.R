#' Thresholds by k-Means clustering
#'
#' This function calculates thresholds by using k-Means clustering. The number of clusters is automatically derived by the Schwarz Bayesian criterion
#'
#' @param df input data frame
#' @param n.clust number of initial cluster (Default: NULL)
#' @param G see \code{\link[mclust]{Mclust}}. Default: 25
#' @param size.data percentage of input data used in anaylsis. Default: NULL
#' @param size.sample sample in initialization of \code{\link[mclust]{Mclust}} in order to speed up the process. Default: 5000, to use full data set to NULL
#' @param seed set seed. Default: 123
#' @param iter.max number of iterlations in \code{\link[mclust]{Mclust}}. Default: 5000
#'
#' @return vector containing the cluster centers
#'
#' @note
#' \itemize{
#'   \item see \href{https://cran.r-project.org/web/packages/mclust/mclust.pdf}{mclust} \emph{(last call: 13-04-2017)}
#'   \item see \href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/kmeans.html}{kmeans} \emph{(last call: 13-04-2017)}
#' }
#'
#' @keywords k-Means thresholding, k-Means clustering, Schwarz Bayesian criterio
#'
#' @examples
#' kMeansThresholding()
#'
#' @export
kMeansThresholding <- function(df, n.clust = NULL, G = 25, size.data = NULL, size.sample = 5000, seed = 123, iter.max = 5000)
{

  if(any(is.na(df)) == TRUE)
  {
    if(class(df) == "data.frame" || class(df)[1] == "data.table")
    {
      df <- df[complete.cases(df),]
    } else
    {
      df <- df[complete.cases(df)]
    }

  }

  m <- as.matrix(df)


  if(!is.null(size.data))
  {
    if(!findInterval(size.data, c(1, 99) ) == 1)
    {
      stop("Wrong input for size.data. Only numbers in the range of 1 - 99 are allowed!")
    }

    set.seed(seed)
    s <- sample(1:nrow(m), size = (nrow(m) * size.data/100))
    m <-  as.matrix(m[s, ])
  }


  if(is.null(n.clust))
  {
    if(!is.null(size.sample) && nrow(m) >= size.sample)
    {
      set.seed(seed)
      #        User      System verstrichen
      #       277.97        1.23      281.60
      n.clust <- mclust::Mclust(m, G = 1:G, initialization = list(subset = sample(1:nrow(m), size = size.sample)))
    }

    if(is.null(size.sample) || nrow(m) < size.sample)
    {
      n.clust <- mclust::Mclust(m, G = 1:G,  na.action=na.exclude)
    }

    n.clust.best <- dim(n.clust$z)[2]
  } else
  {
    n.clust.best <- n.clust
  }

  kmeans.cluster  <- kmeans(x = m, centers = n.clust.best, iter.max = iter.max)

  if( kmeans.cluster$ifault == 4)
  {
    kmeans.cluster <- kmeans(x = m, centers = n.clust.best, algorithm = "MacQueen", iter.max = iter.max)
  }


  return(kmeans.cluster$centers)
}
