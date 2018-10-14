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
st_erase = function(x, y) sf::st_difference(x, sf::st_union(sf::st_combine(y))) # erase y from x



#' Plot sf geometry
#'
#' This function simply plots the geometry of a sf object.
#'
#' @param x object of class sf
#' @return
#' plot
#'
#'
#' @keywords simple feature, plot
#'
#'
#' @export
#'
st_plot = function(x) plot(sf::st_geometry(x))




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
st_perimeter = function(x) x %>% st_cast("MULTILINESTRING") %>% sf::st_length(.) %>% as.numeric(.) # units::drop_units(.)



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
  area <- x %>% sf::st_area(.) %>% as.numeric(.) # units::drop_units(.)

  shape_index <- perimeter / (2 * sqrt(pi * area))
  interior_edge_ratio <- perimeter / area

  return(list(shape_index = shape_index, interior_edge_ratio = interior_edge_ratio))
}





#' Calculate integration index
#'
#' This function calculates an index defining the intergration of new land use area into conisting land use area
#' of the same land use class conisdering two different time steps.
#'
#' @param geom.old object of class sf representing a land use class of a previous time step
#' @param geom.new object of class sf representing a land use class of a following time step
#' @param geom.boundary polygon of class sf representing subregions, e.g. administrative boundaries
#' @param tol tolerance value for overlapping area1 m square
#' @param return.geom If set to TRUE, intermediate geometries are returned as well. Default: FALSE
#' @param quiet show output on console. Default: FALSE
#' @return
#' Vector with integration index (completly integrated: 2/3 < R < 1; good integrated: 1/3 < R < 2/3; low integrated: 0 < R < 1/3; not integrated: 0)
#'
#'
#' @keywords simple feature, shape index, interior edge ratio
#'
#'
#' @export
#'
st_integration_index = function(geom.old, geom.new, geom.boundary = NULL, tol = 0.1, ignr.overlap = FALSE, return.geom = FALSE, quiet = FALSE){

  # get start time of process
  process.time.start <- proc.time()

  ## check input
  if(missing(geom.old) || missing(geom.new)){ stop('Input is missing!')}

  ## check validity of geometries
  if(!all(sf::st_is_valid(geom.old))){ stop('Input of "geom.old" contains not valid geometries. Please try lwgeom::st_make_valid().')}
  if(!all(sf::st_is_valid(geom.new))){ stop('Input of "geom.new" contains not valid geometries. Please try lwgeom::st_make_valid().')}
  if(!is.null(geom.boundary) && !all(sf::st_is_valid(geom.boundary))){ stop('Input of "geom.boundary" contains not valid geometries. Please try lwgeom::st_make_valid().')}

  if(!quiet) cat("... union input geometries \n")
  geom.old <- sf::st_union(x = geom.old)
  geom.new <- sf::st_union(x = geom.new)

  # # # # START CALCULATION OF INTEGRATION INDEX
  ## check for overlapping polygon

  ## common area
  if(!quiet) cat("... intersection of input geometries \n")
  inter <- suppressWarnings(sf::st_intersection(x = geom.old, y = geom.new) %>% sf::st_collection_extract(x = ., type = c("POLYGON")))

  ## new area
  if(!quiet) cat('... erase intersection from "geom.new" (this can take a while!) \n')
  erase <- suppressWarnings(Lslide::st_erase(x = geom.new, y = inter) %>%
                            sf::st_collection_extract(x = ., type = c("POLYGON")) %>%
                            sf::st_cast(x = ., to = "POLYGON") %>%
                            .[which(x = as.numeric(sf::st_area(.)) >= tol)])

  if(!quiet) cat('... conversion to lines \n')
  line.erase <-  sf::st_cast(x = erase, to = "MULTILINESTRING")
  line.inter <- sf::st_cast(x = inter, to = "MULTILINESTRING")

  ## common border
  if(!quiet) cat('... find border lines by intersection \n')
  unique.border <- suppressWarnings(sf::st_intersection(x = line.inter, y = line.erase) %>%
                    sf::st_collection_extract(x = ., type = c("LINESTRING")))

  ## check boundary constraints
  if(!is.null(geom.boundary))
  {
    geom.boundary$ID_BOUNDS <- 1:nrow(geom.boundary)
    geom.boundary <- geom.boundary[, c("ID_BOUNDS", "geometry")]

    if(!quiet) cat('... intersection with boundaries \n')
    unique.border <- suppressWarnings(sf::st_intersection(x = geom.boundary, y = unique.border) %>%
                            st_cast(., "MULTILINESTRING") %>% sf::st_cast(., "LINESTRING")) # cast is necessairy to split multi-object

    dt.unique.border <- unique.border %>%
                        sf::st_set_geometry(x = ., value = NULL) %>%
                        data.table::as.data.table(.)


    erase <- suppressWarnings(sf::st_intersection(x = geom.boundary, y = erase) %>%
                              st_cast(., "MULTIPOLYGON") %>% sf::st_cast(., "POLYGON"))
    dt.erase  <- erase  %>%
                sf::st_set_geometry(x = ., value = NULL) %>%
                data.table::as.data.table(.)


    if(!quiet) cat('... get statistics and calculate index \n')
    dt.unique.border$L <- sf::st_length(x = unique.border) %>% as.numeric() # units::drop_units(x = .)
    dt.erase$P <- Lslide::st_perimeter(x = erase)

    dt.result.UB <- dt.unique.border[,list(L = sum(L, na.rm = TRUE)), by = ID_BOUNDS]
    dt.result.E <- dt.erase[,list(P = sum(P, na.rm = TRUE)), by = ID_BOUNDS]

    dt.result <- merge(x = dt.result.E, y = dt.result.UB, by = "ID_BOUNDS") %>%
                  dplyr::mutate(.data = ., R = L/P)

  } else {
    if(!quiet) cat('... get statistics and calculate index \n')

        dt.result <- data.table::data.table(P = sum(Lslide::st_perimeter(x = erase), na.rm = TRUE),
                            L = sum((sf::st_length(x = unique.border) %>% as.numeric(.)), na.rm = TRUE)) %>% # units::drop_units(x = .))) %>%
                 dplyr::mutate(.data = ., R = L/P)

  } # end of if-else boundary check


  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) cat("------ Run of st_integration_index: " , round(process.time.run["elapsed"][[1]]/60, digits = 4), " Minutes \n")

  if(return.geom)
  {
    return(list(integration_index = dt.result, unique_border = unique.border, new_area = erase))
  } else {
    return(dt.result)
  }
} # end of function st_integration_index





#' Calculate the effective mesh size
#'
#' This function calculates the effective mesh size.
#'
#' @param geom.frag polygon of class sf representing the fragmentation geometry
#' @param geom.boundary polygon of class sf representing subregions, e.g. administrative boundaries
#' @param total.area Numeric value representing size for area. Only to use if geom.boundary is not present. Value must match with the conversion constant c (default hectare). Default: NULL
#' @param conv constant to convert original square m output. Default: 10000 to convert to hectare. If set to 1, than meter square is the result.
#' @param do.preProcessing If TRUE (default), the input of geom.frag is, first, dissolved to single part feature, and second, splitted to multi-parts. By this step it is assured, that polygon connected to each other are summarized
#' @param return.geom If set to TRUE, intermediate geometries are returned as well. Default: FALSE
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
st_mesh = function(geom.frag, geom.boundary = NULL, total.area = NULL, conv = 10000, do.preProcessing = TRUE, return.geom = FALSE, quiet = TRUE)
{
  # get start time of process
  process.time.start <- proc.time()

  ## check input
  if(missing(geom.frag)){ stop('Input of "geom.frag" is missing!')}
  if(!missing(geom.frag) && !("sf" %in% class(geom.frag))){ stop('Input of "geom.frag" is not of class "sf"!')}
  if(!is.null(geom.boundary) && !("sf" %in% class(geom.boundary))){ stop('Input of "geom.boundary" is not of class "sf"!')}
  if(is.null(geom.boundary) && (is.null(total.area) || !is.numeric(total.area))){
   stop('If input of "geom.boundary" is null, then input of "total.area" must be a numeric value representing the size of the study area.')
  }


  ## check validity of geometries
  if(!all(sf::st_is_valid(geom.frag))){ stop('Input of "geom.frag" contains not valid geometries. Please try lwgeom::st_make_valid().')}
  if(!is.null(geom.boundary) && !all(sf::st_is_valid(geom.boundary))){ stop('Input of "geom.boundary" contains not valid geometries. Please try lwgeom::st_make_valid().')}


  if(do.preProcessing)
  {
    if(!quiet) cat("... union geometries to a single geometry with resolved boundaries \n")
    geom.frag <- geom.frag %>% sf::st_union(.)

    if(!quiet) cat("... split multi-parts to single-parts polygon \n")
    geom.frag <- geom.frag %>% st_cast(., "POLYGON") %>% sf::st_sf(ID_FRAG = 1:length(.), geometry = .)
  } else {
    geom.frag$ID_FRAG <- 1:nrow(geom.frag) ## add unique IDs
  }

  if(!is.null(geom.boundary)){ geom.boundary$ID_BOUNDS <- 1:nrow(geom.boundary) } # ## add unique IDs

  ## add Area in m_sq
  geom.frag$A_FRAG <- sf::st_area(geom.frag) %>% as.numeric() %>% `/` (conv) # units::drop_units(.)
  if(!is.null(geom.boundary)){ geom.boundary$A_BOUNDS <- sf::st_area(geom.boundary) %>% as.numeric() %>% `/` (conv)}


  ## subset data
  geom.frag <- geom.frag[, c("ID_FRAG", "A_FRAG", "geometry")]
  if(!is.null(geom.boundary)){geom.boundary <- geom.boundary[, c("ID_BOUNDS", "A_BOUNDS", "geometry")]}


  # # # CALCULATE MESH SIZE INDICES
  ## start calculation
  if(!is.null(geom.boundary))
  {
    ## get intersection
    if(!quiet) cat("... intersection to boundary \n")
    inter <- suppressWarnings(sf::st_intersection(x = geom.boundary, y = geom.frag))
    inter <- suppressWarnings(inter %>% st_cast(., "MULTIPOLYGON") %>% sf::st_cast(., "POLYGON")) # cast is necessairy to split multi-polygons
    inter$A_FRAG_INTER <- sf::st_area(inter) %>% as.numeric() %>% `/` (conv) # overwrite area of fragments, # units::drop_units(.)

    ## calculation of mesh indices
    df.inter <- sf::st_set_geometry(x = inter, value = NULL) %>% data.table::as.data.table(.)

    df.inter.multi <-  df.inter[, list(A_BOUNDS = unique(A_BOUNDS), # splitted multi-parts are joined together
                                    A_FRAG = unique(A_FRAG),
                                    A_FRAG_INTER = sum(A_FRAG_INTER, na.rm = TRUE)), by = list(ID_BOUNDS, ID_FRAG)]

    if(!quiet) cat("... get statistics \n")
    ## cutting-out (CUT) procedure
    mesh.CUT <- df.inter[, list(Fg = unique(A_BOUNDS),
                                Fi = sum(A_FRAG_INTER^2, na.rm = TRUE),
                                Fi_count = length(!is.na(A_FRAG_INTER)),
                                Fi_sum = sum(A_FRAG_INTER, na.rm = TRUE),
                                Fi_min = min(A_FRAG_INTER, na.rm = TRUE),
                                Fi_max = max(A_FRAG_INTER, na.rm = TRUE),
                                Fi_ave = mean(A_FRAG_INTER, na.rm = TRUE),
                                Ci = sum((A_FRAG_INTER/A_BOUNDS)^2, na.rm = TRUE)), by = ID_BOUNDS] %>%
                      dplyr::mutate(., Di = 1-Ci)  %>%
                      dplyr::mutate(., Si = 1/Ci)  %>%
                      dplyr::mutate(., mEff_CUT = Fi/Fg)

    # View(mesh.CUT)

    ## cross-boundary connections (CBC) procedure
    mesh.CBC <- df.inter.multi[, list(Fg = unique(A_BOUNDS),
                                FiErg_count = length(!is.na(A_FRAG)),
                                FiErg_sum = sum(A_FRAG, na.rm = TRUE),
                                FiErg_min = min(A_FRAG, na.rm = TRUE),
                                FiErg_max = max(A_FRAG, na.rm = TRUE),
                                FiErg_ave = mean(A_FRAG, na.rm = TRUE),
                                CiErg = sum((A_FRAG/A_BOUNDS)^2, na.rm = TRUE),
                                Fi_CBC1 = sum(A_FRAG_INTER*A_FRAG, na.rm = TRUE),
                                Fi_CBC2 = sum(2*A_FRAG_INTER*A_FRAG-A_FRAG_INTER^2, na.rm = TRUE)), by = ID_BOUNDS] %>%
                dplyr::mutate(., DiErg = 1-CiErg)  %>%
                dplyr::mutate(., SiErg = 1/CiErg)  %>%
                dplyr::mutate(., mEff_CBC1 = Fi_CBC1/Fg) %>%
                dplyr::mutate(., mEff_CBC2 = Fi_CBC2/Fg)

    # View(mesh.CBC)

    df.result <- merge(x = mesh.CUT, y = mesh.CBC[, -2], by = "ID_BOUNDS")

  } else {

    if(!quiet) cat("... get statistics \n")
    df.result <- sf::st_set_geometry(x = geom.frag, value = NULL) %>% data.table::as.data.table(.)
    df.result <- df.result[, list(Fg = total.area,
                                  Fi = sum(A_FRAG^2, na.rm = TRUE),
                                  Fi_count = length(!is.na(A_FRAG)),
                                  Fi_sum = sum(A_FRAG, na.rm = TRUE),
                                  Fi_min = min(A_FRAG, na.rm = TRUE),
                                  Fi_max = max(A_FRAG, na.rm = TRUE),
                                  Fi_ave = mean(A_FRAG, na.rm = TRUE),
                                  Ci = sum((A_FRAG/total.area)^2, na.rm = TRUE)),] %>%
              dplyr::mutate(., Di = 1-Ci)  %>%
              dplyr::mutate(., Si = 1/Ci)  %>%
              dplyr::mutate(., mEff_CUT = Fi/total.area)

  }

  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) cat("------ Run of st_mesh: " , round(process.time.run["elapsed"][[1]]/60, digits = 4), " Minutes \n")

  if(return.geom)
  {
    if(is.null(geom.boundary))
    {
      return(list(mesh = df.result, geom.frag = geom.frag))
    } else {
      return(list(mesh = df.result, geom.frag = geom.frag, geom.inter = inter))
    }
  } else {
    return(df.result)
  }
}
