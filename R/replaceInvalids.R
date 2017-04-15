#' Replace invalid fields in a data.frame or data.table
#'
#' This function replaces invalid fields (NA, NAN, and NULL) in a data.frame and data.table. These fields can lead
#' to erroneous attribute tables in desktop GIS software.
#'
#' @param x input data.frame or data.table
#' @param replace.value value that is used instead of NA, NAN, or NULL. Default: -9999
#'
#' @note
#' \itemize{
#'   \item see \href{http://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table}{Fastest way to replace NAs in a large data.table} \emph{(last call: 13-04-2017)}
#' }
#'
#' @keywords data.frame, data.table, NA, NAN, NULL, replace
#'
#' @examples
#' replaceInvalids()
#'
#' @export
replaceInvalids <- function(x, replace.value = -9999)
{
  # NA and NAN
  for(j in seq_len(ncol(x)))
    set(x, which(is.na(x[[j]])), j, replace.value)

  # NULL
  for(j in seq_len(ncol(x)))
    set(x, which(is.null(x[[j]])), j, replace.value)
}
