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
st_erase = function(x, y, precision = 0, do.subset = TRUE)
{
  if(do.subset)
  {
    inter <- sf::st_intersects(x = y, y = x)  %>% unlist(.) %>% unique(.)
    x.remain <- x[-inter,]
    x <- x[inter,]
  }

  if(precision != 0)
  {
    x <- x %>% sf::st_set_precision(x = ., precision = precision) %>% lwgeom::st_make_valid(.)
    y <- y %>% sf::st_combine(.) %>% sf::st_union(.) %>%
               sf::st_set_precision(x = ., precision = precision) %>% lwgeom::st_make_valid(.)
  }

  out <- sf::st_difference(x = x, y = y) # erase y from x

  if(do.subset)
  {
    out <- rbind(out, x.remain)
  }

  return(out)
}



#' Erase one geometry from another using Saga GIS
#'
#' This function erase one geometry from another. The projection must be identical.
#'
#' @param x object of class sf. First element: Should be either line or polygon
#' @param y object of class sf. Second element: Always polygon.
#' @param method method of erase. Either "1": , or "2": Line-Polygon Intersection. Default: "1"
#' @param split Set to "1", if multi-part polygons should be splitted to single-part polygons. Default: "0"
#' @param attributes attributes inherited to intersection result. [0] polygon, [1] line, [2] line and polygon. Default: "1"
#' @param env.rsaga SAGA GIS environemnt. Default: RSAGA::rsaga.env()
#' @param check.geom If set to TRUE then geometry is checked with sf::st_is_valid(). If there are invalid geometries, geometries are repaired using lwgeom::st_make_valid(). Default: TRUE
#' @return
#' Geometry of class sfc
#'
#'
#' @keywords simple feature, erase
#'
#'
#' @export
#'
rsaga_erase = function(x, y, method = "1", split = "0", attributes = "1", env.rsaga = RSAGA::rsaga.env(), check.geom = TRUE, quiet = TRUE)
{
  path.x <- file.path(tempdir(), "tmp_x.shp")
  path.y <- file.path(tempdir(), "tmp_y.shp")
  path.result <- file.path(tempdir(), "tmp_result.shp")

  sf::st_write(obj = x, dsn = path.x, delete_layer = TRUE, quiet = quiet)
  sf::st_write(obj = y, dsn = path.y, delete_layer = TRUE, quiet = quiet)

  if(method == "1")
  {
    # RSAGA::rsaga.get.usage(lib = "shapes_polygons", module = 15, env = env.rsaga)
    RSAGA::rsaga.geoprocessor(lib = "shapes_polygons", module = 15, env = env.rsaga, show.output.on.console = !quiet, param = list(
      A = path.x, B = path.y, RESULT = path.result, SPLIT = split))
  }

  if(method == "2")
  {
    # RSAGA::rsaga.get.usage(lib = "shapes_lines", module = 3, env = env.rsaga)
    # ATTRIBUTES: [1] line
    RSAGA::rsaga.geoprocessor(lib = "shapes_lines", module = 3, env = env.rsaga, show.output.on.console = !quiet, param = list(
      LINES = path.x, POLYGONS = path.y, ATTRIBUTES = attributes, DIFFERENCE = path.result))
  }

  ## read data
  out <- sf::st_read(dsn = path.result, quiet = quiet)

  ## check validity
  if(check.geom && !all(sf::st_is_valid(out)))
  {
    warning('Some invalid geometries by "rsaga_erase". Try to correct geomeries using lwgeom::st_make_valid()!')
    out <- lwgeom::st_make_valid(x = out)

    if(method == "1")
    {
      out <- suppressWarnings(out %>% sf::st_collection_extract(x = ., type = c("POLYGON")))
    }

    if(method == "2")
    {
      out <- suppressWarnings(out %>%  sf::st_collection_extract(x = ., type = c("LINESTRING")))
    }
  }
  return(out)
}




#' Return geometry from bounding box
#'
#' This function calculate the geometry based on a bounding box.
#'
#' @param x object of class sf
#' @param extent vector containing numeric values, in the following order: xmin, xmax, ymax, ymin
#' @return
#' Geometry of class sfc
#'
#'
#' @keywords simple feature, geometry of bounding box
#'
#'
#' @export
#'
st_bbox_geom = function(x, extent = NULL)
{
  if(!is.null(x) & !is.null(extent))
  {
    stop("Input conflict: either x or extent!")
  }

  if(is.null(extent))
  {
    out <- sf::st_bbox(x) %>% sf::st_as_sfc(.)
  } else {
    names(extent) <- c("xmin", "xmax", "ymax", "ymin")
    out <- sf::st_bbox(extent) %>% sf::st_as_sfc(.)
  }
  return(out)
} # end of function st_bbox_geom


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
st_perimeter = function(x, ...) suppressWarnings(x %>% sf::st_cast(x = ., to = "MULTILINESTRING", ...) %>% sf::st_length(.) %>% as.numeric(.)) # units::drop_units(.)



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
#' @param tol tolerance value for overlapping area m square
#' @param precision precision for process. See sf::st_set_precision(). Default: 0
#' @param env.rsaga environment of SAGA GIS. If st_erase fails then SAGA GIS erase is used. Default: NULL, but in function call if not set: RSAGA::rsaga.env()
#' @param use.saga use SAGA GIS for erase process. Default: FALSE
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
st_integration_index = function(geom.old, geom.new, geom.boundary = NULL, tol = 0.1, precision = 0,
                                ignr.overlap = FALSE, env.rsaga = NULL, use.saga = FALSE, return.geom = FALSE, quiet = FALSE){

  # get start time of process
  process.time.start <- proc.time()

  ## check input
  if(missing(geom.old) || missing(geom.new)){ stop('Input is missing!')}

  ## check and set precision
  if(!quiet) cat("... check and set st_precision \n")
  if(sf::st_precision(geom.old) != precision)
  {
    geom.old <- geom.old %>% sf::st_set_precision(x = ., precision = precision) %>%
      lwgeom::st_make_valid(.) %>%
      sf::st_collection_extract(x = ., type = "POLYGON")
  }

  if(sf::st_precision(geom.new) != precision)
  {
    geom.new <- geom.new %>% sf::st_set_precision(x = ., precision = precision) %>%
      lwgeom::st_make_valid(.) %>%
      sf::st_collection_extract(x = ., type = "POLYGON")
  }

  ## check validity of geometries
  if(!all(sf::st_is_valid(geom.old))){ stop('Input of "geom.old" contains not valid geometries. Please try lwgeom::st_make_valid().')}
  if(!all(sf::st_is_valid(geom.new))){ stop('Input of "geom.new" contains not valid geometries. Please try lwgeom::st_make_valid().')}
  if(!is.null(geom.boundary) && !all(sf::st_is_valid(geom.boundary))){ stop('Input of "geom.boundary" contains not valid geometries. Please try lwgeom::st_make_valid().')}

  # if(!quiet) cat("... union input geometries \n")
  # geom.old <- sf::st_union(x = geom.old) %>% sf::st_cast(., "POLYGON") %>% sf::st_set_precision(x = ., precision = precision)
  # geom.new <- sf::st_union(x = geom.new) %>% sf::st_cast(., "POLYGON") %>% sf::st_set_precision(x = ., precision = precision)

  # # # # START CALCULATION OF INTEGRATION INDEX
  ## check for overlapping polygon

  ## common area
  if(!quiet) cat("... intersection of input geometries \n")
  inter <- suppressWarnings(sf::st_intersection(x = geom.old, y = geom.new) %>%
                              sf::st_collection_extract(x = ., type = c("POLYGON")))#  %>%
                              # sf::st_set_precision(x = ., precision = precision)

  ## new area
  if(!quiet) cat('... erase intersection from "geom.new" (this can take a while!) \n')

  if(use.saga)
  {
    if(is.null(env.rsaga))
    {
      env.rsaga <-  RSAGA::rsaga.env()
    }

    erase <- Lslide::rsaga_erase(x = geom.new, y = inter, method = "1",
                                 split = "1", env.rsaga = env.rsaga) %>%
            .[which(x = as.numeric(sf::st_area(.)) >= tol),]

  } else {
      erase <- tryCatch({
                      suppressWarnings(Lslide::st_erase(x = geom.new, y = inter, precision = precision) %>%
                            sf::st_collection_extract(x = ., type = c("POLYGON")) %>%
                            sf::st_cast(x = ., to = "POLYGON") %>%
                            .[which(x = as.numeric(sf::st_area(.)) >= tol),])
                    }, error = function(e){
                      warning(paste('SAGA GIS is used due to error in st_erase():', e))
                              if(is.null(env.rsaga))
                              {
                                env.rsaga <-  RSAGA::rsaga.env()
                              }
                       tmp.erase <- Lslide::rsaga_erase(x = geom.new, y = inter, method = "1",
                                                        split = "1", env.rsaga = env.rsaga) %>%
                                    .[which(x = as.numeric(sf::st_area(.)) >= tol),]
                      return(tmp.erase)
                      })
    } # end of use.saga


  if(!quiet) cat('... conversion to lines \n')
  line.erase <-  sf::st_cast(x = erase, to = "MULTILINESTRING") %>% sf::st_set_precision(x = ., precision = precision)
  line.inter <- sf::st_cast(x = inter, to = "MULTILINESTRING") %>% sf::st_set_precision(x = ., precision = precision)

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
                            sf::st_collection_extract(x = ., type = c("LINESTRING")) %>%
                            st_cast(., "MULTILINESTRING") %>% sf::st_cast(., "LINESTRING")) # cast is necessairy to split multi-object

    dt.unique.border <- unique.border %>%
                        sf::st_set_geometry(x = ., value = NULL) %>%
                        data.table::as.data.table(.)


    erase <- suppressWarnings(sf::st_intersection(x = geom.boundary, y = erase) %>%
                              sf::st_collection_extract(x = ., type = c("POLYGON")) %>%
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
  geom.frag$A_FRAG <- sf::st_area(geom.frag) %>% as.numeric() %>% "/" (conv) # units::drop_units(.)
  if(!is.null(geom.boundary)){ geom.boundary$A_BOUNDS <- sf::st_area(geom.boundary) %>% as.numeric() %>% "/" (conv)}


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
    inter$A_FRAG_INTER <- sf::st_area(inter) %>% as.numeric() %>% "/" (conv) # overwrite area of fragments, # units::drop_units(.)

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
} # end of function st_mesh






#' Calculate urban sprawl
#'
#' This function calculates the urban sprawl of urban area.
#'
#' @param geom.urban polygon of class sf representing the fragmentation geometry
#' @param geom.boundary polygon of class sf representing subregions, e.g. administrative boundaries
#' @param dist vector containing distance between lines in x and y direction. Default: c(100, 100) [m]
#' @param trans transformation function {x-1+1/(x+trans.k)} with x as free line and trans.k a constant
#' @param trans.k constant in km for transformation function trans. Default: 1
#' @param tol tolerance value for intersection with erased lines. Buffering procedure is used.Default: 0.1 [m]
#' @param precision precision for process. See sf::st_set_precision(). Default: 0
#' @param extent Numeric value representing extent for area. Format of vector: c(xmin, xmax, ymax, ymin) . Default: NULL
#' @param force.extent If TRUE extent is used instead of geom.boundary (if both are present). Default: FALSE
#' @param do.preProcessing If TRUE (default), the input of geom.frag is, first, dissolved to single part feature, and second, splitted to multi-parts. By this step it is assured, that polygon connected to each other are summarized
#' @param return.geom If set to TRUE, intermediate geometries are returned as well. Default: FALSE
#' @param env.rsaga environment of SAGA GIS. If st_erase fails then SAGA GIS erase is used. Default: NULL, but in function call if not set: RSAGA::rsaga.env()
#' @param use.saga use SAGA GIS for erase process. Default: FALSE
#' @param quiet If set to FALSE, actual state is printed to console. Default: TRUE.
#' @return
#'  strong urban sprawl: 40-50%, less urban sprawl: 80-90%
#'
#'
#' @keywords simple feature, urban sprawl
#'
#'
#' @export
#'
st_urban_sprawl = function(geom.urban, geom.boundary = NULL, dist = c(100, 100), trans = function(x, trans.k){x-1+1/(x+trans.k)}, trans.k = 1, tol = 0.1, extent = NULL, force.extent = FALSE,
                           precision = 0, do.preProcessing = TRUE, env.rsaga = NULL, use.saga = FALSE, return.geom = FALSE, quiet = TRUE)
{
  # get start time of process
  process.time.start <- proc.time()

  ## check input
  if(missing(geom.urban)){ stop('Input of "geom.urban" is missing!')}
  if(!missing(geom.urban) && !("sf" %in% class(geom.urban))){ stop('Input of "geom.urban" is not of class "sf"!')}
  if(!is.null(geom.boundary) && !("sf" %in% class(geom.boundary))){ stop('Input of "geom.boundary" is not of class "sf"!')}
  if(is.null(geom.boundary) && is.null(extent)){
    warning('If input of "geom.boundary" and "extent" is null. Fish net is created using the extent of "geom.urban"!')
  }
  if(is.null(geom.boundary) & is.null(extent) & force.extent){
    stop('If "force.extent" is TRUE, than extent should be set!"!')
  }
  if(!is.null(geom.boundary) & !is.null(extent)){
    warning('Fish net is created using extent of "geom.boundary". "Extent" is skipped. For "extent" use "force.extent"!')
  }

  if(sf::st_precision(geom.urban) != precision)
  {
    geom.urban <- geom.urban %>% sf::st_set_precision(x = ., precision = precision) %>%
                                 lwgeom::st_make_valid(.) %>%
                                 sf::st_collection_extract(x = ., type = "POLYGON")
  }


  ## check validity of geometries
  if(!all(sf::st_is_valid(geom.urban))){ stop('Input of "geom.urban" contains not valid geometries. Please try lwgeom::st_make_valid().')}
  if(!is.null(geom.boundary) && !all(sf::st_is_valid(geom.boundary))){ stop('Input of "geom.boundary" contains not valid geometries. Please try lwgeom::st_make_valid().')}

  ## add unique ID and subset data
  if(!is.null(geom.boundary)){ geom.boundary$ID_BOUNDS <- 1:nrow(geom.boundary) } # ## add unique IDs
  if(!is.null(geom.boundary)){geom.boundary <- geom.boundary[, c("ID_BOUNDS", "geometry")]}


  ## ... create boundary for fishnet
  if(!is.null(geom.boundary) & !force.extent){bbox.fishnet <- Lslide::st_bbox_geom(x = geom.boundary)
  } else if(!is.null(geom.boundary) & force.extent){bbox.fishnet <- Lslide::st_bbox_geom(extent = extent)
  } else if(is.null(geom.boundary) & !is.null(extent)){bbox.fishnet <- Lslide::st_bbox_geom(extent = extent)
  } else {bbox.fishnet <- Lslide::st_bbox_geom(x = geom.urban) }


  ## create and subset fish net
  if(!quiet) cat("... create fishnet \n")
  fishnet <- Lslide::st_make_grid_lines(x = bbox.fishnet, cellsize = dist)

  if(!is.null(geom.boundary) & !force.extent){
    if(!quiet) cat("... subset fishnet to area of interest \n")
    fishnet <- suppressWarnings(sf::st_intersection(x = fishnet, y = geom.boundary)) %>%
                sf::st_collection_extract(x = ., type = "LINESTRING", warn = FALSE)
  } else {
    fishnet <- suppressWarnings(sf::st_intersection(x = fishnet, y = bbox.fishnet)) %>%
      sf::st_collection_extract(x = ., type = "LINESTRING", warn = FALSE)
  }

  # fishnet <- fishnet[, c("ID", "geometry")]


  ## do post-processing of urban geometry
  if(do.preProcessing)
  {
    if(!quiet) cat("... union geometries to a single geometry with resolved boundaries \n")
    geom.urban <- geom.urban %>% sf::st_union(.)

    if(!quiet) cat("... split multi-parts to single-parts polygon \n")
    geom.urban <- geom.urban %>% st_cast(., "POLYGON") %>% sf::st_sf(ID_URBAN = 1:length(.), geometry = .)
  } else {
    geom.urban$ID_URBAN <- 1:nrow(geom.urban) ## add unique IDs
  }

  ## check intersection
  inter.check <- sf::st_intersects(x = geom.urban, y = fishnet)
  inter.empty <- which(sapply(inter.check, function(x) length(x) == 0))
  if(length(inter.empty) > 0)
  {
    inter.A <- geom.urban[inter.empty,] %>% sf::st_area(.) %>% as.numeric(.)
    inter.num <- length(inter.empty)

    warning('Some urban polygons are not intersected by fishnet: ', inter.num,
            ' | area min: ', min(inter.A, na.rm = TRUE), ' - max: ', max(inter.A, na.rm = TRUE), ' [m_sqr] \n')
  }


  ## erase urban area from fishnet using SAGA GIS
  if(!quiet) cat("... erase urban area from fishnet \n")

  if(use.saga)
  {
    if(is.null(env.rsaga))
    {
      env.rsaga <-  RSAGA::rsaga.env()
    }
    erase <-  Lslide::rsaga_erase(x = fishnet, y = geom.urban, method = "2", env.rsaga = env.rsaga)
  } else {
    erase <- suppressWarnings(Lslide::st_erase(x = fishnet, y = geom.urban, precision = precision) %>%
                                sf::st_collection_extract(x = ., type = "LINESTRING"))
  }

  ## split to single part
  erase.single <- erase %>% sf::st_cast(x = ., to = 'LINESTRING', warn = FALSE)


  ## selection of lines
  if(!quiet) cat("... selection of lines \n")
  lines.urban <- sf::st_intersects(x = geom.urban %>%
                                     sf::st_buffer(x = ., dist = tol),
                                   y = erase.single) %>%
                unlist(.) %>% unique(.)


  if(!is.null(geom.boundary) & !force.extent)
  {
    lines.boundary <- sf::st_intersects(x = geom.boundary %>%
                                          sf::st_boundary(x = .) %>%
                                          sf::st_buffer(x = ., dist = tol),
                                        y = erase.single) %>%
                      unlist(.) %>% unique(.)
  } else {
    lines.boundary <- sf::st_intersects(x = bbox.fishnet %>%
                                          sf::st_boundary(x = .) %>%
                                          sf::st_buffer(x = ., dist = tol),
                                        y = erase.single) %>%
                      unlist(.) %>% unique(.)
  }


  # lines.urban  %>% length(.)
  # lines.boundary %>% length(.)

  check.missing <- setdiff(c(1:nrow(erase.single)),
                     (c(lines.urban, lines.boundary)  %>% unique(.)))

  if(length(check.missing) > 0)
  {
    warning("Some lines are neither intersected by boundary nor by urban area: ", length(check.missing), " in total \n")
  }

  erase.single$Type <- NA

  # TPYE 1: Lines between urban area (red lines)
  type.urban <- setdiff(x = lines.urban, y = lines.boundary)
  erase.single$Type[type.urban] <- 1

  # TPYE 2: Lines not touching urban area (blue lines)
  type.boundary <- setdiff(x = lines.boundary, y = lines.urban)
  erase.single$Type[type.boundary] <- 2

  # TYPE 3: Lines between boundary and urban area (green lines)
  type.bound.urb <- intersect(x = lines.urban, y = lines.boundary)
  erase.single$Type[type.bound.urb] <- 3

  # summary(erase.single$Type)
  # c(type.bound.urb, type.boundary, type.urban) %>% unique(.) %>% length(.)

  # browser()

  ## group single lines back to multi-lines
  if((!is.null(geom.boundary) & !force.extent))
  {
    erase.final <- st_dissolve(x = erase.single[which(erase.single$Type == 2 | erase.single$Type == 3),],
                               by = list("Type", "ID_FNET", "ID_BOUNDS"))

    erase.final <- rbind(erase.final, erase.single[which(erase.single$Type == 1),])

  } else {
    erase.final <- st_dissolve(x = erase.single[which(erase.single$Type == 2 | erase.single$Type == 3),],
                               by = list("Type", "ID_FNET"))

    erase.final <- rbind(erase.final, erase.single[which(erase.single$Type == 1),])
  }

  erase.final <- erase.final %>% sf::st_collection_extract(x = ., type = "LINESTRING")

  erase.final$L <- sf::st_length(x = erase.final) %>% as.numeric(.) %>% "/" (1000) # convert from meter to km!
  erase.final$L_trans <- ifelse(erase.final$Type == 2, erase.final$L, trans(x = erase.final$L, trans.k = trans.k))


  ## summarizing final data
  if(!quiet) cat("... get statistics: urban sprawl \n")

  if((!is.null(geom.boundary) & !force.extent))
  {
    df.result <- erase.final %>% sf::st_set_geometry(x = ., value = NULL) %>%
                                 data.table::as.data.table(.) %>%
                                  .[,list(L = sum(L, na.rm = TRUE),
                                          L_trans = sum(L_trans, na.rm = TRUE)), by = ID_BOUNDS] %>%
                                  dplyr::mutate(.data = ., FFE = L_trans/L*100)
  } else {
    df.result <- erase.final %>% sf::st_set_geometry(x = ., value = NULL) %>%
                                 data.table::as.data.table(.) %>%
                                 .[,list(L = sum(L, na.rm = TRUE),
                                            L_trans = sum(L_trans, na.rm = TRUE)),] %>%
                                 dplyr::mutate(.data = ., FFE = L_trans/L*100)
  }


  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) cat("------ Run of urban_sprawl: " , round(process.time.run["elapsed"][[1]]/60, digits = 4), " Minutes \n")

  if(return.geom)
  {
    return(list(FFE = df.result, geom = erase.final))
  } else {
    return(df.result)
  }

} # end of function st_urban_sprawl





#' Make multiline fishnet
#'
#' This function creates a fishnet, consisting of multilines in x and y direction.
#'
#' @param x object of class sf
#' @return
#' Geometry of class sfc
#'
#'
#' @keywords simple feature, fishnet, multilines
#'
#'
#' @export
#'
st_make_grid_lines = function(x, cellsize)
{
  # ... create fishnet using corners (lower left corner is starting point)
  fishnet <- sf::st_make_grid(x = x, cellsize = cellsize, what = "corners")

  # ... get sf::st_make_grid function parameters (automatically calculated in function call)
  offset <- sf::st_bbox(obj =  x)[1:2]
  nx <- ceiling((sf::st_bbox(obj = x)[3] - offset[1])/cellsize[1]) + 1 # add 1 to fit overall columns
  ny <- ceiling((sf::st_bbox(obj = x)[4] - offset[2])/cellsize[2]) + 1 # add 1 to fit overall rows

  # ... create lines
  linesX <- lapply(1:ny, function(i, nx, fishnet){
    if(i == 1)
    {
      index <- 1:nx
    } else {
      index <- ((i-1)*nx+1):(i*nx)
    }

    line <- fishnet[index] %>% st_combine(.) %>%
      sf::st_multilinestring(x = ., dim = "XY") %>%
      st_sfc(.) %>%
      sf::st_sf(ID_FNET = paste0("X", i), geometry = .)

    return(line)
  }, nx = nx, fishnet = fishnet) %>% do.call(what = rbind, args = .)

  seqY <- seq(from = 1, to = length(fishnet), by = nx)
  linesY <- lapply(0:(nx-1), FUN = function(i, seqY, fishnet){

    if(i == 0){
      index <- seqY
    } else {
      index <- seqY+i
    }

    line <- fishnet[index] %>% st_combine(.) %>%
      sf::st_multilinestring(x = ., dim = "XY") %>%
      st_sfc(.) %>%
      sf::st_sf(ID_FNET = paste0("Y", i), geometry = .)

    return(line)
  }, seqY = seqY, fishnet) %>% do.call(what = rbind, args = .)



  fishnet.multi <- rbind(linesX, linesY)
  sf::st_crs(fishnet.multi) <- sf::st_crs(fishnet)

  if(!all(sf::st_is_valid(fishnet.multi)))
  {
    warning("Fishnet contains invalid geometries!")
  }

  return(fishnet.multi)
} # end of function st_make_grid_lines
