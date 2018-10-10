#' Erase one geometry from another
#'
#' This function erase one geometry from another. The projection must be identical.
#'
#' @param x object of class sf
#' @param y object of class sf
#' @return
#' Geometry of class sfc
#'
#'
#' @keywords simple feature, erase
#'
#'
#' @export
#'
st_erase = function(x, y) st_difference(x, st_union(st_combine(y))) # erase y from x




#' Dissolve geometry
#'
#' This function dissovle a geometry based on a field. Optionnally, field statistics can be computed.
#'
#' @param x object of class sf
#' @param by field for dissolving. Default: NULL
#' @param ... Optional: field statistcs
#' @return
#' dissolved geometry of class sf
#'
#'
#' @keywords simple feature, dissolve
#'
#'
#' @export
#'
st_dissolve = function(x, by = NULL, ...) x %>% dplyr::group_by(.dots = by) %>% dplyr::summarise(...)



#' Row bind simple features
#'
#' This function combines simple feature with different columns using dplr::bind_rows().
#'
#' @param x list of objects of class sf
#' @param geom vector containing the geometries of x. Order must be the same as in list.
#' @return
#' sf object with combined data.frames
#'
#'
#' @keywords simple feature, rbind
#'
#'
#' @export
#'
st_rbind = function(x, geom) x %>% lapply(X = ., FUN = function(x) sf::st_set_geometry(x = x, value = NULL)) %>%
  dplyr::bind_rows() %>% sf::st_set_geometry(x = ., value = geom)



#' Calculate the perimeter of a polygon
#'
#' This function calculates the perimeter of a polygon.
#'
#' @param x object of class sf
#' @return
#' Vector with perimeter values
#'
#'
#' @keywords simple feature, perimeter
#'
#'
#' @export
#'
st_perimeter = function(x) x %>% st_cast("MULTILINESTRING") %>% sf::st_length(.)



#' Calculate the shape indices of a polygon
#'
#' This function calculates the shape index of a polygon.
#'
#' @param x object of class sf
#' @return
#' Vector with shape index and interior edge ratio values
#'
#'
#' @keywords simple feature, shape index, interior edge ratio
#'
#'
#' @export
#'
st_shape_indices = function(x){
  perimeter <- x %>% st_perimeter(.)
  area <- x %>% sf::st_area(.) %>% units::drop_units(.)

  shape_index <- perimeter / (2 * sqrt(pi * area))
  interior_edge_ratio <- perimeter / area

  return(list(shape_index = shape_index, interior_edge_ratio = interior_edge_ratio))
}

