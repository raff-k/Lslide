#' Get rain event of a landslide
#'
#' This function calcuates different precipitation characteristics for a specific time-series:
#' total precipitation, number of rainfall events, weighted mean intensitiy of rainfall events (normalized by MAP, RD or RDN),
#' cumulative critical event rainfall (normalized by MAP, RD or RDN), maximum rainfall during critical rainfall event,
#' duration of critical rainfall event, critical rainfall intensitiy (normalized by MAP, RD or RDN), rainfall at day of failure (start date),
#' rainfall intensity at day of failure (start date), maximum rainfall at day of failure (start date).
#'
#' @param x vector containing precipitation
#' @param type type of precipitation data, i.e. "daily" or "hourly". Default: "daily"
#' @param rainThresh rainfall below this threshold is considered as dry. This threshold should be adapted according to type. Default: 0.2
#' @param rainOffLength length of dry period to separate rainfall events. This parameter should be adapted according to type. Default: 2 (days), for "hourly" data 48-96 h can be found in literature
#' @param RD average number of rainy days in a year, proxy for locate climate conditions. Default: NULL
#' @param MAP mean annual precipitation, the long-term yearly average precipitation, see CRU - climate research units for number. Default: NULL
#' @param RDN a climatic index that provides better description (or proxy) than the MAP for the occurence of extreme storm events (Guzzetti et al. 2006: 247). Default: MAP/RD
#' @return
#' vector containing rain characeristics (see description).
#'
#'
#' @note
#' \itemize{
#'   \item Guzzetti, F., Peruccacci, S., Rossi, M., & Stark, C. P. (2007). Rainfall thresholds for the initiation of landslides in central and southern Europe. Meteorology and atmospheric physics, 98(3-4), 239-267.
#'   \item Rossi, M., Luciani, S., Valigi, D., Kirschbaum, D., Brunetti, M. T., Peruccacci, S., & Guzzetti, F. (2017). Statistical approaches for the definition of landslide rainfall thresholds and their uncertainty using rain gauge and satellite data. Geomorphology, 285, 16-27.
#' }
#'
#'
#' @keywords rainfall tresholds, rainfall event, landslide, rainfall characeristics
#'
#'
#' @export

getRainEventData <- function(x, type = "daily", rainThresh = 0.2, rainOffLength = 2, RD = NULL, MAP = NULL, RDN = MAP/RD)
{
  if(type != "daily" & type != "hourly")
  {
    stop('"type"must be either "daily" or "hourly"')
  }

  if(rainThresh < 0)
  {
    stop('"rainThresh" must be a positive value')
  }

  if(rainOffLength < 0)
  {
    stop('"rainOffLength" must be a positive value')
  }

  if((!is.null(RD) && RD == 0) | (!is.null(MAP) && MAP == 0))
  {
    stop('"RD" or "MAP" are nor allowed to be "0"')
  }


  # reverse x, meaning that precipitation before event is at first position
  x <- rev(x)

  # x.tmp <- x
  # # convert data to daily basis
  # if(type == "daily") # TO DO: if using date, it could be aggregated by dates(?)
  # {
  #   #  get number of days
  #   d.num <- length(x)/24
  #
  #   # create group for splitting
  #   split.goup <- rep(x = 1:d.num, each = 24)
  #
  #   # split data into days
  #   x.split <- split(x = x, f = split.goup)
  #
  #   # aggregate values
  #   x.tmp <- sapply(X = x.split, FUN = sum, na.rm = TRUE)
  # }

  ## find position of rainfall under threshold
  x.pos.rainThresh <- which(x <= rainThresh)

  ## set elements under threshold to 0 mm (dry!)
  # x.tmp[x.pos.rainThresh] <- 0

  ## get consecutive positions under rainTreshold
  x.pos.rainThresh.U <- split(seq_along(along.with = x.pos.rainThresh), cumsum(c(0, diff(x.pos.rainThresh) > 1))) # https://stackoverflow.com/questions/36041576/how-do-i-find-ranges-of-successive-numbers-in-a-vector-in-r
  x.pos.rainThresh.ULen <- unlist(x.pos.rainThresh.U[which(lengths(x.pos.rainThresh.U) >= rainOffLength)])
  x.pos.rainThresh.index <- x.pos.rainThresh[x.pos.rainThresh.ULen]

  ## get rain events
  x.pos.rainEvent <- setdiff(seq_along(along.with = x), x.pos.rainThresh.index)
  x.pos.rainEvent.Len <- split(seq_along(along.with = x.pos.rainEvent), cumsum(c(0, diff(x.pos.rainEvent) > 1)))
  x.pos.rainEvent.index <- lapply(x.pos.rainEvent.Len, function(x, y){y[x]}, y = x.pos.rainEvent) # get original indices


  ### get precipitation characteristics

  ## general
  precip.tot <- sum(x, na.rm = TRUE)

  # browser()

  ## rain events
  # total sum of precipitation
  RE.sum <- sapply(X = x.pos.rainEvent.index, FUN = function(X, precip) {
    sum(precip[X], na.rm = TRUE)
  }, precip = x)

  if(!is.null(MAP)){RE.sum.MAP <- RE.sum/MAP}
  if(!is.null(RD)){RE.sum.RD <- RE.sum/RD}
  if(!is.null(MAP) & !is.null(RD)){RE.sum.RDN <- RE.sum/RDN}

  # maximum of precipitation of each rain event
  RE.max <- sapply(X = x.pos.rainEvent.index, FUN = function(X, precip) {
    max(precip[X], na.rm = TRUE)
  }, precip = x)

  # duration of rain events
  RE.dur <- sapply(X = x.pos.rainEvent.index, FUN = length)

  # number of rain events
  rainEvent.nb <- length(x.pos.rainEvent.index)

  # weighted mean intensity of rain events
  RE.wID <- weighted.mean(x = RE.sum/RE.dur, w = RE.dur, na.rm = TRUE)
  if(!is.null(MAP)){RE.wID.MAP <- RE.wID/MAP}
  if(!is.null(RD)){RE.wID.RD <- RE.wID/RD}
  if(!is.null(MAP) & !is.null(RD)){RE.wID.RDN <- RE.wID/RDN}


  ## critical rainfall event
  # critical event rainfall
  cRE.sum <- RE.sum[[1]]
  if(!is.null(MAP)){cRE.sum.MAP <- RE.sum.MAP[[1]]}
  if(!is.null(RD)){cRE.sum.RD <- RE.sum.RD[[1]]}
  if(!is.null(MAP) & !is.null(RD)){cRE.sum.RDN <- RE.sum.RDN[[1]]}

  # critical maximal event rainfall
  cRE.max <- RE.max[[1]]

  # critical event rainfall duration
  cRE.dur <- RE.dur[[1]]

  # critical event rainfall intensity
  cRE.ID <- cRE.sum/cRE.dur
  if(!is.null(MAP)){cRE.ID.MAP <- cRE.sum.MAP/cRE.dur}
  if(!is.null(RD)){cRE.ID.RD <- cRE.sum.RD/cRE.dur}
  if(!is.null(MAP) & !is.null(RD)){cRE.ID.RDN <- cRE.sum.RDN/cRE.dur}


  ## rainfall at the date of failure
  if(type == "daily")
  {
    cRD <- x[1]
    cRD.ID <- x[1]
    cRD.max <- x[1]
  }

  if(type == "hourly")
  {
    cRD <- sum(x[c(1:24)], na.rm = TRUE)
    cRD.ID <- cRD/length(which(x[c(1:24)] > 0))
    cRD.max <- max(x[c(1:24)], na.rm = TRUE)
  }


  ## return data
  if(!is.null(MAP) & is.null(RD)){
    res <- c(precip.tot, rainEvent.nb, RE.wID, RE.wID.MAP, cRE.sum, cRE.sum.MAP, cRE.max, cRE.dur, cRE.ID, cRE.ID.MAP, cRD, cRD.ID, cRD.max)
    names(res) <- c("P_total", "RE_number", "RE_weightIntensity", "RE_weightIntensity_MAP", "cRE_sum", "cRE_sum_MAP", "cRE_max", "cRE_duration", "cRE_Intensity", "cRE_Intensity_MAP", "cRD_sum", "cRD_ID", "cRD_max")
  } else if(!is.null(RD) & is.null(MAP)){
    res <- c(precip.tot, rainEvent.nb, RE.wID, RE.wID.RD, cRE.sum, cRE.sum.RD, cRE.max, cRE.dur, cRE.ID, cRE.ID.RD, cRD, cRD.ID, cRD.max)
    names(res) <- c("P_total", "RE_number", "RE_weightIntensity", "RE_weightIntensity_RD", "cRE_sum", "cRE_sum_RD", "cRE_max", "cRE_duration", "cRE_Intensity", "cRE_Intensity_RD", "cRD_sum", "cRD_ID", "cRD_max")
  } else if(!is.null(RD) & !is.null(MAP)){
    res <- c(precip.tot, rainEvent.nb, RE.wID, RE.wID.MAP, RE.wID.RD, cRE.sum, cRE.sum.MAP, cRE.sum.RD, cRE.max, cRE.dur, cRE.ID, cRE.ID.MAP, cRE.ID.RD, cRD, cRD.ID, cRD.max)
    names(res) <- c("P_total", "RE_number", "RE_weightIntensity", "RE_weightIntensity_MAP", "RE_weightIntensity_RD",
                    "cRE_sum", "cRE_sum_MAP", "cRE_sum_RD", "cRE_max", "cRE_duration", "cRE_Intensity", "cRE_Intensity_MAP", "cRE_Intensity_RD", "cRD_sum", "cRD_ID", "cRD_max")

  } else {
    res <- c(precip.tot, rainEvent.nb, RE.wID, cRE.sum, cRE.max, cRE.dur, cRE.ID, cRD, cRD.ID, cRD.max)
    names(res) <- c("P_total", "RE_number", "RE_weightIntensity", "cRE_sum", "cRE_max", "cRE_duration", "cRE_Intensity", "cRD_sum", "cRD_ID", "cRD_max")

  }

  return(res)
} # end of function
