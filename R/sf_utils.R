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
st_perimeter = function(x) x %>% st_cast("MULTILINESTRING") %>% sf::st_length(.) %>% units::drop_units(.)



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




#' Calculate the effective mesh size
#'
#' This function calculates the effective mesh size.
#'
#' @param geom.frag polygon of class sf representing the fragmentation geometry
#' @param geom.boundary polygon of class sf representing subregions, e.g. administrative boundaries
#' @param total.area method how area is accumulated. Either 0 or 1. 0 (default) means total area of subregion polygon. 1 means total area as sum of meshs
#' @param method analysis method, relevant if geom.boundary is present. 0: all methods, 1: cutting-out (CUT) procedure, 2: cross-boundary connections (CBC) procedure (default)
#' @param do.preProcessing If TRUE (default), the input of geom.frag is, first, dissolved to single part feature, and second, splitted to multi-parts. By this step it is assured, that polygon connected to each other are summarized
#' @param quiet If set to FALSE, actual state is printed to console. Default: TRUE.
#' @return
#' Vector with shape index and interior edge ratio values
#'
#'
#' @keywords simple feature, mesh, effective mesh size
#'
#'
#' @export
#'
sf_mesh = function(geom.frag, geom.boundary = NULL, total.area = 0, method = 2, do.preProcessing = TRUE, quiet = TRUE)
{
  # get start time of process
  process.time.start <- proc.time()

  ## check input
  if(missing(geom.frag)){ stop('Input of "geom.frag" is missing!')}
  if(!missing(geom.frag) && !("sf" %in% class(geom.frag))){ stop('Input of "geom.frag" is not of class "sf"!')}
  if(!is.null(geom.boundary) && !("sf" %in% class(geom.boundary))){ stop('Input of "geom.boundary" is not of class "sf"!')}
  if(is.null(geom.boundary)){
    warning('Input of "geom.boundary" is null. Method of "total.area" is set to 1. Parameter "method" is masked.')
    total.area <- 1
    method <- -1
  }


  ## check validity of geometries
  if(all(sf::st_is_valid(geom.frag))){ stop('Input of "geom.frag" contains not valid geometries. Please try lwgeom::st_make_valid().')}
  if(all(sf::st_is_valid(geom.boundary))){ stop('Input of "geom.boundary" contains not valid geometries. Please try lwgeom::st_make_valid().')}


  if(do.preProcessing)
  {
    cat("... union geometries to a single geometry with resolved boundaries \n")
    geom.frag <- geom.frag %>% sf::st_union(.)

    cat("... split multi-parts to single-parts polygon \n")
    geom.frag <- geom.frag %>% st_cast(., "POLYGON") %>% sf::st_sf(ID_FRAG = 1:length(.), geometry = .)
  } else {
    geom.frag$ID_FRAG <- 1:nrow(geom.frag) ## add unique IDs
  }

  if(!is.null(geom.boundary)){ geom.boundary$ID_BOUNDS <- 1:nrow(geom.boundary) } # ## add unique IDs

  ## add Area in m_sq
  geom.frag$A_FRAG <- sf::st_area(geom.frag) %>% units::drop_units(.)
  if(!is.null(geom.boundary)){ geom.boundary$A_BOUNDS <- sf::st_area(geom.boundary) %>% units::drop_units(.) }


  ## subset data
  geom.frag <- geom.frag[, c("ID_FRAG", "A_FRAG", "geometry")]
  if(!is.null(geom.boundary)){geom.boundary <- geom.boundary[, c("ID_BOUNDS", "A_BOUNDS", "geometry")]}

  ## defining functions
  # intern_getArea <- function(x, y)
  # {
  #   # browser()
  if(!is.null(geom.boundary))
  {
    ## get intersection
    inter <- sf::st_intersection(x = geom.boundary, y = geom.frag)
    inter <- inter %>% st_cast(., "MULTIPOLYGON") %>% sf::st_cast(., "POLYGON") # cast is necessairy to split multi-polygons
    inter$A_FRAG_INTER <- sf::st_area(inter) %>% units::drop_units(.) # overwrite area of fragments

    ## calculation of mesh indices
    df.inter <- sf::st_set_geometry(x = inter, value = NULL) %>% data.table::as.data.table(.)

    df.inter.multi <-  df.inter[, list(A_BOUNDS = unique(A_BOUNDS), # splitted multi-parts are joined together
                                    A_FRAG = unique(A_FRAG),
                                    A_FRAG_INTER = sum(A_FRAG_INTER, na.rm = TRUE)), by = list(ID_BOUNDS, ID_FRAG)]



    # cutting-out (CUT) procedure
    # mesh.CUT <- df.inter[, list(Fg = unique(A_BOUNDS),
    #                             Fi = sum(A_FRAG_INTER^2, na.rm = TRUE),
    #                             Fi_count = length(!is.na(A_FRAG_INTER)),
    #                             Fi_sum = sum(A_FRAG_INTER, na.rm = TRUE),
    #                             Fi_min = min(A_FRAG_INTER, na.rm = TRUE),
    #                             Fi_max = max(A_FRAG_INTER, na.rm = TRUE),
    #                             Fi_CBC1 = sum(A_FRAG_INTER*A_FRAG, na.rm = TRUE),
    #                             # FiErg_sum = sum(A_FRAG, na.rm = TRUE),
    #                             Fi_CBC2 = sum(2*A_FRAG_INTER*A_FRAG-A_FRAG_INTER^2, na.rm = TRUE)), by = ID_BOUNDS] %>%
    #             dplyr::mutate(., mEff_CUT = Fi/Fg) %>%
    #             dplyr::mutate(., mEff_CBC1 = Fi_CBC1/Fg) %>%
    #             dplyr::mutate(., mEff_CBC2 = Fi_CBC2/Fg)

    mesh.CUT <- df.inter[, list(Fg = unique(A_BOUNDS),
                                Fi = sum(A_FRAG_INTER^2, na.rm = TRUE),
                                Fi_count = length(!is.na(A_FRAG_INTER)),
                                Fi_sum = sum(A_FRAG_INTER, na.rm = TRUE),
                                Fi_min = min(A_FRAG_INTER, na.rm = TRUE),
                                Fi_max = max(A_FRAG_INTER, na.rm = TRUE)), by = ID_BOUNDS] %>%
                dplyr::mutate(., mEff_CUT = Fi/Fg)
    View(mesh.CUT)

    #### here to GO ON!
    mesh.CUT2 <- df.inter.multi[, list(Fg = unique(A_BOUNDS),
                                Fi = sum(A_FRAG_INTER^2, na.rm = TRUE),
                                Fi_count = length(!is.na(A_FRAG_INTER)),
                                Fi_sum = sum(A_FRAG_INTER, na.rm = TRUE),
                                Fi_min = min(A_FRAG_INTER, na.rm = TRUE),
                                Fi_max = max(A_FRAG_INTER, na.rm = TRUE)), by = ID_BOUNDS] %>%
      dplyr::mutate(., mEff_CUT = Fi/Fg)
    View(mesh.CUT2)
    View(mesh.CBC)
    # cross-boundary connections (CBC) procedure

  } else {


  }
    # out <- lapply(X = ID, FUN = function(i, inter){return(inter[which(inter$ID_BOUNDS == i), ])}, inter = inter)
    # names(out) <- ID
   #  return(inter)
   # } # end of cutting-out (CUT) procedure

 #  out.m1 <- intern_getArea(x = geom.boundary, y = geom.frag)


  # 2: cross-boundary connections (CBC) procedure
  # intern_m2 <- function(x, y)
  # {
  #   browser()
  #   inter <- sf::st_intersects(x = x, y = y)
  #   out <- lapply(X = 1:length(inter), FUN = function(i, y, inter)
  #   {
  #     inter.i <- inter[[i]]
  #     if(length(inter.i) > 0)
  #     {
  #       tmp <- y[inter.i, ]
  #       tmp$ID_BOUNDS <- i
  #       return(tmp)
  #     } else {
  #       return(NULL)
  #     }
  #   }, y = y, inter = inter) %>% do.call(what = rbind, args = .) # end of lapply function
  #
  #   # names(out) <- ID
  #   return(out)
  # } # end of cross-boundary connections (CBC) procedure
  #
  # intern_m2(x = geom.boundary, y = geom.frag)


  # # # CALCULATE MESH SIZE INDICES

  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) cat("------ Run of sf_mesh: " , round(process.time.run["elapsed"][[1]]/60, digits = 4), " Minutes \n")
  return(df.result)
}
