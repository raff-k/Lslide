#' Correct .dbf of file
#'
#' This function enables the possibility to overwrite the header of an existing .dbf-file.
#' Some applications, such as SAGA GIS, store shapefiles in formats that are not conform with the
#' ESRI-shapefile format. By this function, "corrupt" headers can be replaced.
#'
#' @param x full path to .dbf-file, but without extension
#' @param end.n start number of adjusting. Default: number of columns
#' @param adjust.n number of columns that are renamed. Default: 0
#' @param new.colnames vector with new names
#'
#' @note
#' length of vector with new names must be identical to number of adjusted columns
#'
#' @keywords dbf, correct dbf
#'
#'
#' @export
correctDBF <- function(x, end.n = length(colnames(d$dbf)), adjust.n = 0, new.colnames)
{
  d <- suppressMessages(shapefiles::read.dbf(paste0(file_path_sans_ext(x), ".dbf")))
  colnames(d$dbf) <- c(colnames(d$dbf)[0:(end.n-adjust.n)], new.colnames)
  shapefiles::write.dbf(d, paste0(file_path_sans_ext(x), ".dbf")) # write dbf with better header
  rm(d)
}
