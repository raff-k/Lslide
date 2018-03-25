#' Get rain event of a landslide
#'
#' This function calcuates different precipitation characteristics for a specific time-series:
#' total precipitation, number of rainfall events, weighted mean intensitiy of rainfall events (normalized by MAP, RD or RDN),
#' cumulative critical event rainfall (normalized by MAP, RD or RDN), maximum rainfall during critical rainfall event,
#' duration of critical rainfall event, critical rainfall intensitiy (normalized by MAP, RD or RDN), rainfall at day of failure (start date),
#' rainfall intensity at day of failure (start date), maximum rainfall at day of failure (start date).
#'
#' @param x vector containing precipitation
#' @param dates vector containing dates. The length of dates must be similar to the length of x. Default: NULL
#' @param timesteps time period or length of observation. The rev(x) and rev(dates) are subsetted to this length according to date.of.failure if set. I.e. 24 for hourly or 1 for daily data. Default: NULL
#' @param date.of.failure date of failure. If set data is subsetted to this date. Must be of class "POSIXct" or "POSIXt". Default: NULL
#' @param sub.RainEvent examine potential sub-rain-events of critical rainfall event. Default: TRUE
#' @param all.RainEvent if TRUE, all rain events in data are extracted. By setting this option, no critival rain event metrics are computed. Default: FALSE
#' @param cumu.RainFall vector containing time intervals for cumulative rainfall. I.e. c(24, 48, 96) for 1, 2 and 4 days aggregation. Default: NULL
#' @param return.DataFrame only the rain events are returned as a data.frame. Default: FALSE
#' @param S1.rainThresh isolated rainfall measurements below this thresholds are removed in the first step. Default: 0.2
#' @param S3.rainThresh exclusion of irrelevant rainfall sub-events under and equal to this threshold (third step). Default: 1 [mm]
#' @param S1.rainOffLength dry periods between isolated rain events in the first step. Default: c(3, 6) (hours). When dates is NULL, then the smallest values is used for separation.
#' @param S2.rainOffLength dry periods between rainfall sub-events in the second step. Default: c(6, 12) (hours). When dates is NULL, then the smallest values is used for separation.
#' @param S4.rainOffLength dry periods between rainfall sub-events in the second step. Default: c(6, 12) (hours). When dates is NULL, then the smallest values is used for separation.
#' @param RD average number of rainy days in a year, proxy for locate climate conditions. Default: NULL
#' @param MAP mean annual precipitation, the long-term yearly average precipitation, see CRU - climate research units for number. Default: NULL
#' @param RDN a climatic index that provides better description (or proxy) than the MAP for the occurence of extreme storm events (Guzzetti et al. 2006: 247). Default: MAP/RD
#' @param index.month.warm.season month indices of the warm season. First element is start, and second element represents the end (all including). Only relevant when dates are set. Default: c(4, 10) (including April, including Ocotober)
#'
#' @return
#' vector containing rainfall metrics (see description). If return.DataFrame is TRUE a data.frame is returned containing similar
#' rain metrics for all rain events.
#'
#'
#' @note
#' \itemize{
#'   \item thresholds are oriented at hourly data
#'   \item timesteps with precipitation equal 0 are included (see Melillo et al. 2015: 314)
#'   \item Guzzetti, F., Peruccacci, S., Rossi, M., & Stark, C. P. (2007). Rainfall thresholds for the initiation of landslides in central and southern Europe. Meteorology and atmospheric physics, 98(3-4), 239-267.
#'   \item Rossi, M., Luciani, S., Valigi, D., Kirschbaum, D., Brunetti, M. T., Peruccacci, S., & Guzzetti, F. (2017). Statistical approaches for the definition of landslide rainfall thresholds and their uncertainty using rain gauge and satellite data. Geomorphology, 285, 16-27.
#'   \item Melillo, M., Brunetti, M. T., Peruccacci, S., Gariano, S. L., & Guzzetti, F. (2015). An algorithm for the objective reconstruction of rainfall events responsible for landslides. Landslides, 12(2), 311-320.
#' }
#'
#'
#' @keywords rainfall tresholds, rainfall event, landslide, rainfall metrics
#'
#'
#' @export
getRainEventData <- function(x, dates = NULL, timesteps = NULL, date.of.failure = NULL, sub.RainEvent = TRUE, all.RainEvent = FALSE, cumu.RainFall = NULL, return.DataFrame = FALSE,
                             S1.rainThresh = 0.2, S3.rainThresh = 1, S1.rainOffLength = c(3, 6), S2.rainOffLength = c(6, 12), S4.rainOffLength = c(48, 96),
                             RD = NULL, MAP = NULL, RDN = MAP/RD, index.month.warm.season = c(4, 10))
{

  # # # # # # # # # CHECK POTENTIAL ERRORS # # # # # # # # #


  if(all.RainEvent)
  {
    if(!is.null(timesteps) || !is.null(date.of.failure) || is.null(cumu.RainFall))
    {
      warnings('The function parameters "timesteps", "date.of.failure", and "cumu.RainFall" are set to NULL')

      timesteps <- NULL
      date.of.failure <- NULL
      cumu.RainFall <- NULL
    }
  } #  end of if all.RainEvent


  if(!is.null(dates) && length(dates) != length(x))
  {
    stop('Length of "x" and "dates" must be identical')
  }

  if(length(index.month.warm.season) != 2 | max(index.month.warm.season) > 12 | min(index.month.warm.season) < 1)
  {
    Stop('Function argument "index.month.warm.season" is wrongly defined')
  }

  # if(type != "daily" & type != "hourly")
  # {
  #   stop('"type"must be either "daily" or "hourly"')
  # }

  if(S1.rainThresh < 0 | S3.rainThresh < 0)
  {
    stop('"rainThresh" must be a positive value')
  }

  if(length(S1.rainOffLength) > 2)
  {
    stop('"rainOffLength" should contain maximum 2 elements')
  }

  if(length(which(S1.rainOffLength < 0)) > 0 | length(which(S2.rainOffLength < 0)) > 0 | length(which(S4.rainOffLength < 0)) > 0)
  {
    stop('"rainOffLength" thresholds must contain positive numbers')
  }


  if((!is.null(RD) && RD == 0) | (!is.null(MAP) && MAP == 0))
  {
    stop('"RD" or "MAP" are nor allowed to be "0"')
  }

  if(!is.null(date.of.failure) && class(date.of.failure)[1] != "POSIXct" & class(date.of.failure)[1] != "POSIXt")
  {
    stop('Date of failure must be of class POSIXct or POSIXt')
  }


  if(!is.null(timesteps) && (length(timesteps) > 1 | timesteps[1] < 0 | length(x) < timesteps[1]))
  {
    stop('Something wrong in "timesteps". Only a single positive integer number is accepted, which is equal or smaller to "x"')
  }

  if(!is.null(cumu.RainFall) && (length(which(cumu.RainFall <= 0)) > 0 | max(cumu.RainFall) > length(x)))
  {
    stop('The calculation of cumulative rainfall should have at least have 1 time step intervall indicating by positive integer numbers.
         In addition, the cumulative rainfall should not exceeded the length of the precipitation vector')
  }






  # # # # # # # # # START ALGORITHM # # # # # # # # #


  # reverse x, meaning that precipitation before event is at first position
  x <- rev(x)
  x.input <- x

  # reverse dates
  dates <- rev(dates)



  ## subset data to date of failure
  if(!is.null(date.of.failure))
  {
    if(!is.null(dates))
    {
      # get subset indices
      DoF.sub <- which(dates <= date.of.failure)

      # subset dates
      dates <- dates[DoF.sub]

      # subset precipitation
      x <- x[DoF.sub]
    } else{
      warning('date of failure is set, but no dates avaiable. Therefore, subsetting is skipped')
    }
  }


  ## subset data to timesteps
  if(!is.null(timesteps))
  {
    if(timesteps > length(x))
    {
      stop("After subsetting data to the date of failure, there are more timesteps than data")
    }

    # subset precipitation
    x <- x[1:timesteps]

    # subset dates
    if(!is.null(dates)){dates <- dates[1:timesteps]}
  }

  if(!is.null(cumu.RainFall) && any(timesteps < cumu.RainFall))
  {
    stop('The time intervall for cumulative rain fall exceeded timesteps')
  }


  ## extract months of dates and extent dates
  ## align thresholds
  if(!is.null(dates))
  {
    # get months of date
    dates <- list(dates, as.numeric(format(dates, "%m")))
  } else {
    S1.rainOffLength <- min(S1.rainOffLength) # threshold of first step
    S2.rainOffLength <- min(S2.rainOffLength) # threshold of second step
    S4.rainOffLength <- min(S4.rainOffLength) # threshold of fourth step
  }



  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  ## Step 1: Detection and exclusion of isolated rainfall measurements --------------------------
  # ... find isolated rainfall (rainThresh must be fix to 0)
  S1.isolRF <- findRainFallPosition(x = x, dates = dates, rainThresh = c(0, 0), rainOffLength = S1.rainOffLength,
                                    op.rainThresh = ">", op.rainOffLength = "<=",
                                    index.month.warm.season = index.month.warm.season)

  # ... set irrelevant rainfall depening on threshold to 0
  x[S1.isolRF[which(x[S1.isolRF] <= S1.rainThresh)]] <- 0 # THRESHOLD MUST BE DEFINED!!!!!!




  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  ## Step 2: Identification of rainfall sub-events ----------------------------------------------
  S2.DryPeriods <- findRainFallPosition(x = x, dates = dates, rainThresh = c(0,0), rainOffLength = S2.rainOffLength,
                                        op.rainThresh = "==", op.rainOffLength = ">=",
                                        index.month.warm.season = index.month.warm.season)
  S2.RainEvents <- findRainEvent(x = x, x.pos.dryPeriods = S2.DryPeriods)


  # ... sum of precipitation of sub-rain-events
  S2.RE.sum <- sapply(X = S2.RainEvents, FUN = function(X, precip) {
    sum(precip[X], na.rm = TRUE)
  }, precip = x)

  # ... duration of sub-rain-events
  # S2.RE.dur <- sapply(X = S2.RainEvents, FUN = length) # DO WE NEED THIS?




  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  ## Step 3: Exclusion if irrelevant rainfall sub-events ----------------------------------------
  S3.RainEvents <- S2.RainEvents[which(S2.RE.sum > S3.rainThresh)]
  S3.RE.exclusion <- unlist(S2.RainEvents[which(S2.RE.sum <= S3.rainThresh)])

  # ... set rainfall to 0
  x[S3.RE.exclusion] <- 0




  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  ## Step 4/5: Identification of rainfall events --------------------------------------------------
  S4.DryPeriods <- findRainFallPosition(x = x, dates = dates, rainThresh = c(0,0), rainOffLength = S4.rainOffLength,
                                        op.rainThresh = "==", op.rainOffLength = ">=",
                                        index.month.warm.season = index.month.warm.season)
  S4.RainEvents <- findRainEvent(x = x, x.pos.dryPeriods = S4.DryPeriods)

  ## ... get rainfall metrics
  # general
  precip.tot <- sum(x.input, na.rm = TRUE)
  names(precip.tot) <- "sum_total"

  res <- precip.tot


  # cumulative rainfall
  if(!is.null(cumu.RainFall))
  {
    precip.cum <- sapply(X = cumu.RainFall, function(X, precip){sum(precip[1:X], na.rm = TRUE)}, precip = x)
    names(precip.cum) <- paste0("sum_cumu_", cumu.RainFall)

    res <- c(res, precip.cum)
  }

  if(all.RainEvent)
  {
    # all rain event data
    cERM <- calcEventRainfallMetrics(x = x, dates = dates, list.RainEvents = S4.RainEvents, modus = "sub", RD = RD, MAP = MAP, RDN = RDN)

    # ... gsub s in names
    names(cERM) <- gsub(pattern = "^(s)", replacement = "", x = names(cERM))

  } else {
    # critical rain event
    cERM <- calcEventRainfallMetrics(x = x, dates = dates, list.RainEvents = S4.RainEvents, modus = "main", RD = RD, MAP = MAP, RDN = RDN)

  }

  # merge
  res <- c(res, cERM)





  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  ## Step 6: Rainfall measurements for events with landslides --------------------------------------------------

  if(sub.RainEvent)
  {

    if(all.RainEvent)
    {
      # ... find all sub events
      cERM.S6.SubEvents <- lapply(X = 1:length(S4.RainEvents), FUN = function(i, S4RE, S3RE, precip, dates, RD, MAP, RDN)
        {

          # browser()
          # ... get max and min index
          # ... ... find nearest min index
          S3RE.minIndices <- sapply(X = S3RE, FUN = min)
          S4RE.minIndex.pos <- unname(which.min(abs(S3RE.minIndices - min(S4RE[[i]]))))


          S4RE.minIndex <- min(S3RE[[S4RE.minIndex.pos]])
          S4RE.maxIndex <- max(S4RE[[i]])

          # ... check for sub-rainfall events
          S6SRE.check.min <- sapply(X = S3RE, FUN = function(x, minIndex){min(x) >= minIndex}, minIndex = S4RE.minIndex)
          S6SRE.check.max <- sapply(X = S3RE, FUN = function(x, maxIndex){max(x) <= maxIndex}, maxIndex = S4RE.maxIndex)
          S6SE <- S3RE[which(S6SRE.check.min & S6SRE.check.max)]

          # ... calculate rainfall metrics
          cERM.sub <- calcEventRainfallMetrics(x = precip, dates = dates, list.RainEvents = S6SE, modus = "sub", RD = RD, MAP = MAP, RDN = RDN)
          names(cERM.sub) <- paste0(names(cERM.sub), "_", i)

          return(cERM.sub)


      }, S3RE = S3.RainEvents, S4RE = S4.RainEvents, precip = x, dates = dates, RD = RD, MAP = MAP, RDN = RDN)

      res <- c(res, unlist(cERM.S6.SubEvents))

      # ... flatten list
      # S6.SubEvents.fl <- unlist(S6.SubEvents, recursive = FALSE)
      # S6.SubEvents.Names <- gsub(pattern = "\\." , replacement = "_" , names(S6.SubEvents.fl))

      # ... calculate sub-event rainfall metrics
      # cERM.sub <- calcEventRainfallMetrics(x = x, dates = dates, list.RainEvents =  S6.SubEvents.fl, modus = "sub", RD = RD, MAP = MAP, RDN = RDN)

      # ... rename items
      # cERM.sub.NumNames <- paste0(paste0(c(1:length(S4.RainEvents)), "_"), unlist(sapply(X = 1:length(S6.SubEvents.LenName), FUN = function(x, n){rep(x, n[x])}, n = S6.SubEvents.LenName)))

    } else {

        ## ... find sub-events of critical rainfall event
        # get largest index of first rain event
        S4.RainEvents.maxIndex <- max(S4.RainEvents[[1]])

        # ... check for sub-rainfall events
        S6.SubEvents.check <- sapply(X = S3.RainEvents, FUN = function(x, maxIndex){max(x) <= maxIndex}, maxIndex = S4.RainEvents.maxIndex)
        S6.SubEvents <- S3.RainEvents[which(S6.SubEvents.check)]

        # if(length(S6.SubEvents) > 1) # > 1, because 1 would be similar to the critical rainfall event
        # {
          cERM.sub <- calcEventRainfallMetrics(x = x, dates = dates, list.RainEvents = S6.SubEvents, modus = "sub", RD = RD, MAP = MAP, RDN = RDN)
          res <- c(res, cERM.sub)
        # }
    } # end of if(all.RainEvent)
  } # end if sub.RainEvent




  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  ## Return data --------------------------------------------------------------------------------

  if(return.DataFrame)
  {
    # browser()

    # if(sub.RainEvent && length(S6.SubEvents) > 1)
    if(sub.RainEvent)
    {
      if(all.RainEvent)
      {
        # ... remove of first two items: "sum_total", "RE_number", "RE_weightIntensity"
        # remove items
        res <- res[!names(res) %in% c("sum_total", "RE_number", "RE_weightIntensity",
                                      "RE_weightIntensity_MAP", "RE_weightIntensity_RD", "RE_weightIntensity_RDN")]


        # ... get and remove number of sub-events
        extr.sRENum.pos <- grep(pattern = "sRE_number_.", x = names(res)) # ... get positions
        extr.sRENum <- res[extr.sRENum.pos] # ... get number of sub events
        res <- res[-extr.sRENum.pos] # remove this variable from results

        extr.sREWeight.pos <- grep(pattern = "sRE_weightIntensity_.", x = names(res)) # ... get positions
        extr.sREWeight <- res[extr.sREWeight.pos] # ... get number of sub events
        res <- res[-extr.sREWeight.pos] # remove this variable from results

        # ... ... create data.frame for weighted Intensity
        sREWeight.ColName <- unique(gsub(pattern = "_[[:digit:]]", replacement = "", x = names(extr.sREWeight)))
        sREWeight.df <- data.frame(matrix(extr.sREWeight, ncol = length(sREWeight.ColName)))
        colnames(sREWeight.df) <- sREWeight.ColName


        # ... get start indices of sub-events
        startSub <- grep(pattern = "sRE_sum_1.", x = names(res))[1]


        # ... create main data frame
        # ... ... names
        res.df.main.ColName <- names(res)[1:(startSub-1)] # get names
        res.df.main.ColName <- unique(gsub(pattern = "_[[:digit:]]", replacement = "", x = res.df.main.ColName)) # adapt colnames

        # ... ... data.frame
        res.df.main <- data.frame(matrix(res[1:(startSub-1)], ncol = ((startSub-1)/length(S4.RainEvents))))
        colnames(res.df.main) <- res.df.main.ColName # re-name columns
        res.df.main$sRE_number <- extr.sRENum # add number of sub-events
        res.df.main <- cbind(res.df.main, sREWeight.df) # weightesIntensitiy data.frame
        res.df.main$ID <- c(1:nrow(res.df.main))*100 # id for ordering later


        # ... create sub data frame
        # ... ... names
        res.df.sub.ColName <- names(res)[startSub:length(res)] # get names
        res.df.sub.ColName <- unique(gsub(pattern = "^s|_[[:digit:]]", replacement = "", x = res.df.sub.ColName)) # adapt colnames

        # ... ... data.frame
        res.df.sub <- data.frame(matrix(res[startSub:length(res)], ncol = ((startSub-1)/length(S4.RainEvents))))
        colnames(res.df.sub) <- res.df.sub.ColName # re-name columns
        res.df.sub$ID <- unlist(sapply(X = 1:length(extr.sRENum), FUN = function(i, x){i*100 + seq(1:extr.sRENum[i])}, x = extr.sRENum)) # get ID variables

        # ... row-bind data frames
        res.df <- dplyr::bind_rows(res.df.main, res.df.sub)

        # re-order data.frame based on ID
        res.df <- res.df[with(res.df, order(ID)), ]
        row.names(res.df)[which(!is.na(res.df$sRE_number))] <- paste0("RE_", res.df[which(!is.na(res.df$sRE_number)), ]$ID)
        row.names(res.df)[which(is.na(res.df$sRE_number))] <- paste0("sRE_", res.df[which(is.na(res.df$sRE_number)), ]$ID)

        # ... finally remove ID
        res.df$ID <- NULL


      } else {

          # get main event items
          items <- c("sRE_number", "sRE_weightIntensity", "sRE_weightIntensity_MAP", "sRE_weightIntensity_RD", "sRE_weightIntensity_RDN",
                     "sum_total", names(res)[grep(pattern = "sum_cumu", names(res))], "RE_total",
                     "RE_number", "RE_weightIntensity", "RE_weightIntensity_MAP", "RE_weightIntensity_RD", "RE_weightIntensity_RDN")

          # remove items
          res.df <- res[!names(res) %in% items]

          # rename
          names(res.df) <- substring(names(res.df), 2)

          ## create data frame
          startSub <- which(names(res.df) == "RE_sum_1")

          res.colNames <- names(res.df)[1: (startSub-1)]

          # ... create first main then sub data.frame
          res.df.main <- data.frame(matrix(res.df[1:(startSub-1)], ncol = length(c(1:(startSub-1)))))
          res.df.sub <- data.frame(matrix(res.df[startSub:length(res.df)], ncol = length(c(1:(startSub-1)))))

          # ... combine both data.frames
          if(length(S6.SubEvents) > 1) # if there is only 1 sub-event, the event is not returned!!!!
          {
            res.df <- rbind(res.df.main, res.df.sub)
          } else {
            res.df <- res.df.main
          }


          # ... set names
          colnames(res.df) <-  res.colNames

          # ... set items to data frame
          items.res <- res[which(names(res) %in% items)]
          items.res <- items.res[items[which(items %in% names(res))]] # just re-ordering

          # ... create empty data.frame and add values
          items.df <- data.frame(matrix(NA, nrow = nrow(res.df), ncol = length(items.res)))
          items.df[1, ] <- items.res
          colnames(items.df) <- names(items.res)

          # ... add new columns to result
          res.df <- cbind(res.df, items.df)

          if(length(S6.SubEvents) > 1) # if there is only 1 sub-event, the event is not returned!!!!
          {
            rownames(res.df) <- c("cRE", paste0("sRE", c(1:length(S6.SubEvents))))
          } else {
            rownames(res.df) <- "cRE"
          }
      } # end of if all.RainEvent
    } else {

      if(all.RainEvent)
      {
        # ... remove of first two items: "sum_total", "RE_number"
        # remove items
        res <- res[!names(res) %in% c("sum_total", "RE_number", "RE_weightIntensity",
                                      "RE_weightIntensity_MAP", "RE_weightIntensity_RD", "RE_weightIntensity_RDN")]

        # get colnames
        res.df.ColName <- unique(gsub(pattern = "^s|_[[:digit:]]", replacement = "", x = names(res))) # adapt colnames

        # ... ... data.frame
        res.df <- data.frame(matrix(res, ncol = length(res.df.ColName)))
        colnames(res.df) <- res.df.ColName # re-name columns
        rownames(res.df) <- paste0("RE", seq(100, 100*nrow(res.df), 100))

      } else {

        # re-order vector
        items.fst <- which(names(res) == "cRE_sum")
        res <- res[c(items.fst:length(res), 1:(items.fst-1))]

        res.df <- data.frame(matrix(data = res, ncol = length(res)))
        names(res.df) <- gsub(pattern = "^(c)", replacement = "", x = names(res))
      } # end of if a..RainEvent

    } # end of if-else: sub.RainEvent & length(S6.SubEvents) > 1

    # modified data.frame output
    return(res.df)
  } else {


    # standard vector output
    return(res)
  }


} # end of function getRainEventData()
















#' calcEventRainfallMetrics
#'
#' This function calcuates different precipitation metrics.
#'
#' @param x vector containing precipitation
#' @param list.RainEvents list containing indices of rain events
#' @param modus modus of calculation of rain metrics. "sub" of sub-rainfall events or "main" for critical rainfall event
#' @param ... for RD, MAP or RDN
#' @return
#' vector containing rainfall metrics
#'
#' @export
calcEventRainfallMetrics <- function(x, dates, list.RainEvents, modus, RD = RD, MAP = MAP, RDN = RDN)
{

  ## cumulative rainfall
  # total sum of precipitation of each rain event
  RE.sum <- sapply(X = list.RainEvents, FUN = function(X, precip) {
    sum(precip[X], na.rm = TRUE)
  }, precip = x)


  # maximum of precipitation of each rain event
  RE.max <- sapply(X = list.RainEvents, FUN = function(X, precip) {
    max(precip[X], na.rm = TRUE)
  }, precip = x)


  ## duration of rain events
  RE.dur <- sapply(X = list.RainEvents, FUN = length)


  ## number of rain events
  rainEvent.nb <- length(list.RainEvents)


  ### ... MAIN rainfall metrics
  if(modus == "main")
  {

    ## naming of variables
    names(rainEvent.nb) <- "RE_number"

    ## overall sum of rain events
    RE.tot <- sum(RE.sum)
    names(RE.tot) <- "RE_total"

    # normalizations
    if(!is.null(MAP)){RE.sum.MAP <- RE.sum/MAP}
    if(!is.null(RD)){RE.sum.RD <- RE.sum/RD}
    if(!is.null(MAP) & !is.null(RD)){RE.sum.RDN <- RE.sum/RDN}


    ## weighted mean intensity of rain events
    RE.wID <- weighted.mean(x = RE.sum/RE.dur, w = RE.dur, na.rm = TRUE)
    names(RE.wID) <- "RE_weightIntensity"

    # normalizations
    if(!is.null(MAP)){RE.wID.MAP <- RE.wID/MAP; names(RE.wID.MAP) <- "RE_weightIntensity_MAP"}
    if(!is.null(RD)){RE.wID.RD <- RE.wID/RD; names(RE.wID.RD) <- "RE_weightIntensity_RD"}
    if(!is.null(MAP) & !is.null(RD)){RE.wID.RDN <- RE.wID/RDN; names(RE.wID.RDN) <- "RE_weightIntensity_RDN"}


    ## critical rainfall event
    # critical event rainfall
    cRE.sum <- RE.sum[[1]]
    names(cRE.sum) <- "cRE_sum"

    # normalizations
    if(!is.null(MAP)){cRE.sum.MAP <- RE.sum.MAP[[1]]; names(cRE.sum.MAP) <- "cRE_sum_MAP"}
    if(!is.null(RD)){cRE.sum.RD <- RE.sum.RD[[1]]; names(cRE.sum.RD) <- "cRE_sum_RD"}
    if(!is.null(MAP) & !is.null(RD)){cRE.sum.RDN <- RE.sum.RDN[[1]]; names(cRE.sum.RDN) <- "cRE_sum_RDN"}


    ## critical maximal event rainfall
    cRE.max <- RE.max[[1]]
    names(cRE.max) <- "cRE_max"


    ## critical event rainfall duration
    cRE.dur <- RE.dur[[1]]
    names(cRE.dur) <- "cRE_duration"


    ## critical event rainfall intensity
    cRE.ID <- cRE.sum/cRE.dur
    names(cRE.ID) <- "cRE_Intensity"

    # normalizations
    if(!is.null(MAP)){cRE.ID.MAP <- cRE.sum.MAP/cRE.dur; names(cRE.ID.MAP) <- "cRE_Intensity_MAP"}
    if(!is.null(RD)){cRE.ID.RD <- cRE.sum.RD/cRE.dur; names(cRE.ID.RD) <- "cRE_Intensity_RD"}
    if(!is.null(MAP) & !is.null(RD)){cRE.ID.RDN <- cRE.sum.RDN/cRE.dur; names(cRE.ID.RDN) <- "cRE_Intensity_RDN"}


    ## get range of rain event
    # cRE.range <- paste0(range(list.RainEvents[[1]],  na.rm = TRUE), collapse = ":")
    cRE.range <- range(list.RainEvents[[1]],  na.rm = TRUE)
    cRE.range.start <- min(cRE.range)
    cRE.range.end <- max(cRE.range)

    names(cRE.range.start) <- "cRE_range_start"
    names(cRE.range.end) <- "cRE_range_end"

    ## get dates out of range of rain event
    if(!is.null(dates))
    {
      cRE.date.start <- as.numeric(gsub("-|[[:space:]]", "", format(dates[[1]][cRE.range.start], "%Y-%m-%d-%H")))
      names(cRE.date.start) <- "cRE_date_start"

      # browser()
      cRE.date.end <- as.numeric(gsub("-|[[:space:]]", "", format(dates[[1]][cRE.range.end], "%Y-%m-%d-%H")))
      names(cRE.date.end) <- "cRE_date_end"
    }
  } # end of main rainfall metrics



  ### ... SUB rainfall metrics
  if(modus == "sub")
  {
    ## naming of variables
    names(rainEvent.nb) <- "sRE_number"
    names(RE.sum) <-  paste0("sRE_sum_", c(1:rainEvent.nb))
    names(RE.max) <-  paste0("sRE_max_", c(1:rainEvent.nb))
    names(RE.dur) <-  paste0("sRE_dur_", c(1:rainEvent.nb))


    ## rainfall intensity
    RE.ID <- RE.sum/RE.dur
    names(RE.ID) <- paste0("sRE_Intensity_", c(1:rainEvent.nb))

    # normalizations c(rainEvent.nb, RE.sum, RE.max, RE.dur, RE.ID)
    if(!is.null(MAP))
    {
      RE.sum.MAP <- RE.sum/MAP
      names(RE.sum.MAP) <- paste0("sRE_sum_MAP_", c(1:rainEvent.nb))

      RE.ID.MAP <- RE.sum.MAP/RE.dur
      names(RE.ID.MAP) <- paste0("sRE_Intensity_MAP_", c(1:rainEvent.nb))
    }

    if(!is.null(RD))
    {
      RE.sum.RD <- RE.sum/RD
      names(RE.sum.RD) <- paste0("sRE_sum_RN_", c(1:rainEvent.nb))

      RE.ID.RD <- RE.sum.RD/RE.dur
      names(RE.ID.RD) <- paste0("sRE_Intensity_RD_", c(1:rainEvent.nb))
    }

    if(!is.null(MAP) & !is.null(RD))
    {
      RE.sum.RDN <- RE.sum/RDN
      names(RE.sum.RDN) <- paste0("sRE_sum_RND_", c(1:rainEvent.nb))

      RE.ID.RDN <- RE.sum.RDN/RE.dur
      names(RE.ID.RDN) <- paste0("sRE_Intensity_RDN_", c(1:rainEvent.nb))
    }


    ## weighted mean intensity of sub rain events
    RE.wID <- weighted.mean(x = RE.sum/RE.dur, w = RE.dur, na.rm = TRUE)
    names(RE.wID) <- "sRE_weightIntensity"

    # normalizations
    if(!is.null(MAP)){RE.wID.MAP <- RE.wID/MAP; names(RE.wID.MAP) <- "sRE_weightIntensity_MAP"}
    if(!is.null(RD)){RE.wID.RD <- RE.wID/RD; names(RE.wID.RD) <- "sRE_weightIntensity_RD"}
    if(!is.null(MAP) & !is.null(RD)){RE.wID.RDN <- RE.wID/RDN; names(RE.wID.RDN) <- "sRE_weightIntensity_RDN"}



    ## get range of rain events
    # RE.range <- sapply(X = list.RainEvents, FUN = function(X) { paste0(range(x,  na.rm = TRUE), collapse = ":")})
    RE.range.start <- sapply(X = list.RainEvents, FUN = min, na.rm = TRUE)
    names(RE.range.start ) <- paste0("sRE_range_start_", c(1:rainEvent.nb))

    RE.range.end <- sapply(X = list.RainEvents, FUN = max, na.rm = TRUE)
    names(RE.range.end) <- paste0("sRE_range_end_", c(1:rainEvent.nb))


    ## get dates out of range of rain event
    if(!is.null(dates))
    {
      RE.date.start <- as.numeric(gsub("-|[[:space:]]", "", format(dates[[1]][RE.range.start], "%Y-%m-%d-%H")))
      # RE.date.start <- as.numeric(gsub("-|[[:space:]]", "", substring(dates[[1]][RE.range.start], 1, 13)))
      names(RE.date.start) <- paste0("sRE_date_start_", c(1:rainEvent.nb))

      RE.date.end <- as.numeric(gsub("-|[[:space:]]", "", format(dates[[1]][RE.range.end], "%Y-%m-%d-%H")))
      # RE.date.end <- as.numeric(gsub("-|[[:space:]]", "", substring(dates[[1]][RE.range.end], 1, 13)))
      names(RE.date.end) <- paste0("sRE_date_end_", c(1:rainEvent.nb))
    }

  } # end of sub rainfall metrics



  ### ... return data
  if(!is.null(MAP) & is.null(RD)){

    if(modus == "main")
    {
      res <- c(RE.tot, rainEvent.nb, RE.wID, RE.wID.MAP, cRE.sum, cRE.sum.MAP, cRE.max, cRE.dur, cRE.ID, cRE.ID.MAP, cRE.range.start, cRE.range.end)
      # names(res) <- c("RE_total", "RE_number", "RE_weightIntensity", "RE_weightIntensity_MAP", "cRE_sum", "cRE_sum_MAP", "cRE_max", "cRE_duration", "cRE_Intensity", "cRE_Intensity_MAP", "cRE_range_start", "cRE_range_end")
      if(!is.null(dates)){res <- c(res, cRE.date.start, cRE.date.end)}
    }

    if(modus == "sub")
    {
      res <- c(rainEvent.nb, RE.sum, RE.sum.MAP, RE.max, RE.dur, RE.ID, RE.ID.MAP, RE.wID, RE.wID.MAP, RE.range.start, RE.range.end)
      # names(res) <- c("sRE_number", paste0("sRE_sum_", c(1:rainEvent.nb)), paste0("sRE_sum_MAP_", c(1:rainEvent.nb)),
      #                 paste0("sRE_max_", c(1:rainEvent.nb)), paste0("sRE_dur_", c(1:rainEvent.nb)),
      #                 paste0("sRE_Intensity_", c(1:rainEvent.nb)), paste0("sRE_Intensity_MAP_", c(1:rainEvent.nb)),
      #                 paste0("sRE_range_start_", c(1:rainEvent.nb)), paste0("sRE_range_end_", c(1:rainEvent.nb)))
      if(!is.null(dates)){res <- c(res, RE.date.start, RE.date.end)}
      #if(length(list.RainEvents) <= 1){names(res) <- gsub(pattern = "_1", replacement = "", names(res))}
    }

  } else if(!is.null(RD) & is.null(MAP)){

    if(modus == "main")
    {
      res <- c(RE.tot, rainEvent.nb, RE.wID, RE.wID.RD, cRE.sum, cRE.sum.RD, cRE.max, cRE.dur, cRE.ID, cRE.ID.RD, cRE.range.start, cRE.range.end)
      # names(res) <- c("RE_total", "RE_number", "RE_weightIntensity", "RE_weightIntensity_RD", "cRE_sum", "cRE_sum_RD", "cRE_max", "cRE_duration", "cRE_Intensity", "cRE_Intensity_RD", "cRE_range_start", "cRE_range_end")
      if(!is.null(dates)){res <- c(res, cRE.date.start, cRE.date.end)}
    }

    if(modus == "sub")
    {
      res <- c(rainEvent.nb, RE.sum, RE.sum.RD, RE.max, RE.dur, RE.ID, RE.ID.RD, RE.wID,  RE.wID.RD, RE.range.start, RE.range.end)
      # names(res) <- c("sRE_number", paste0("sRE_sum_", c(1:rainEvent.nb)), paste0("sRE_sum_RD_", c(1:rainEvent.nb)),
      #                 paste0("sRE_max_", c(1:rainEvent.nb)), paste0("sRE_dur_", c(1:rainEvent.nb)),
      #                 paste0("sRE_Intensity_", c(1:rainEvent.nb)), paste0("sRE_Intensity_RD_", c(1:rainEvent.nb)),
      #                 paste0("sRE_range_start_", c(1:rainEvent.nb)), paste0("sRE_range_end_", c(1:rainEvent.nb)))
      if(!is.null(dates)){res <- c(res, RE.date.start, RE.date.end)}
      #if(length(list.RainEvents) <= 1){names(res) <- gsub(pattern = "_1", replacement = "", names(res))}
    }

  } else if(!is.null(RD) & !is.null(MAP)){

    if(modus == "main")
    {
      res <- c(RE.tot, rainEvent.nb, RE.wID, RE.wID.MAP, RE.wID.RD, RE.wID.RDN, cRE.sum, cRE.sum.MAP, cRE.sum.RD, cRE.sum.RDN, cRE.max, cRE.dur, cRE.ID, cRE.ID.MAP, cRE.ID.RD, cRE.ID.RDN, cRE.range.start, cRE.range.end)
      # names(res) <- c("RE_total", "RE_number", "RE_weightIntensity", "RE_weightIntensity_MAP", "RE_weightIntensity_RD", "RE_weightIntensity_RDN",
      #                 "cRE_sum", "cRE_sum_MAP", "cRE_sum_RD", "cRE_sum_RDN", "cRE_max", "cRE_duration", "cRE_Intensity", "cRE_Intensity_MAP", "cRE_Intensity_RD", "cRE_Intensity_RDN", "cRE_range_start", "cRE_range_end")
      if(!is.null(dates)){res <- c(res, cRE.date.start, cRE.date.end)}
    }

    if(modus == "sub")
    {
      res <- c(rainEvent.nb, RE.sum, RE.sum.MAP, RE.sum.RD, RE.sum.RDN, RE.max, RE.dur, RE.ID, RE.ID.MAP, RE.ID.RD, RE.ID.RDN, RE.wID, RE.wID.MAP,  RE.wID.RD, RE.wID.RDN, RE.range.start, RE.range.end)
      # names(res) <- c("sRE_number", paste0("sRE_sum_", c(1:rainEvent.nb)), paste0("sRE_sum_MAP_", c(1:rainEvent.nb)), paste0("sRE_sum_RD_", c(1:rainEvent.nb)), paste0("sRE_sum_RDN_", c(1:rainEvent.nb)),
      #                 paste0("sRE_max_", c(1:rainEvent.nb)), paste0("sRE_dur_", c(1:rainEvent.nb)),
      #                 paste0("sRE_Intensity_", c(1:rainEvent.nb)), paste0("sRE_Intensity_MAP_", c(1:rainEvent.nb)), paste0("sRE_Intensity_RD_", c(1:rainEvent.nb)), paste0("sRE_Intensity_RDN_", c(1:rainEvent.nb)),
      #                 paste0("sRE_range_start_", c(1:rainEvent.nb)), paste0("sRE_range_end_", c(1:rainEvent.nb)))
      if(!is.null(dates)){res <- c(res, RE.date.start, RE.date.end)}
      #if(length(list.RainEvents) <= 1){names(res) <- gsub(pattern = "_1", replacement = "", names(res))}
    }

  } else {

    if(modus == "main")
    {
      res <- c(RE.tot, rainEvent.nb, RE.wID, cRE.sum, cRE.max, cRE.dur, cRE.ID, cRE.range.start, cRE.range.end)
      # names(res) <- c("RE_total", "RE_number", "RE_weightIntensity", "cRE_sum", "cRE_max", "cRE_duration", "cRE_Intensity", "cRE_range_start", "cRE_range_end")
      if(!is.null(dates)){res <- c(res, cRE.date.start, cRE.date.end)}
    }

    if(modus == "sub")
    {
      res <- c(rainEvent.nb, RE.sum, RE.max, RE.dur, RE.ID, RE.wID, RE.range.start, RE.range.end)
      # names(res) <- c("sRE_number", paste0("sRE_sum_", c(1:rainEvent.nb)), paste0("sRE_max_", c(1:rainEvent.nb)),
      #                 paste0("sRE_dur_", c(1:rainEvent.nb)), paste0("sRE_Intensity_", c(1:rainEvent.nb)),
      #                 paste0("sRE_range_start_", c(1:rainEvent.nb)), paste0("sRE_range_end_", c(1:rainEvent.nb)))
      if(!is.null(dates)){res <- c(res, RE.date.start, RE.date.end)}
      #if(length(list.RainEvents) <= 1){names(res) <- gsub(pattern = "_1", replacement = "", names(res))}
    }

  } # end of if - else if - else statement

    return(res)

} # end of calcEventRainfallMetrics












#' findRainFallPosition
#'
#' This function return indices of a specific selection.
#'
#' @param x vector containing precipitation
#' @param dates vector containing months of dates. The length of dates must be similar to the length of x
#' @param rainThresh list containing indices of rain events
#' @param rainOffLength modus of calculation of rain metrics. "sub" of sub-rainfall events or "main" for critical rainfall event
#' @param op.rainThresh operator for rain-threshold: x OP rainThresh
#' @param op.rainOffLength operator for rainOff-threshold: lengths of event OP rainOffLength-Threshold
#' @return
#' vector containing indices of x corresponding to specific selection.
#'
#' @export
findRainFallPosition <- function(x, dates, rainThresh, rainOffLength, op.rainThresh, op.rainOffLength,
                                 index.month.warm.season)
{
  # browser()

  ## find position of rainfall under/over/equal to threshold | x OP rainThresh
  if(!is.null(dates))
  {
    ## find position of dates correpsonding to season
    pos.warm <- which(dates[[2]] >= index.month.warm.season[1] & dates[[2]] <= index.month.warm.season[2])
    pos.cold <- which(dates[[2]] < index.month.warm.season[1] | dates[[2]] > index.month.warm.season[2])

    ## thresholding x based on season and season threshold
    x.pos.rainThresh.warm <- which(do.call(op.rainThresh, list(x[pos.warm], rainThresh[1])))
    x.pos.rainThresh.cold <- which(do.call(op.rainThresh, list(x[pos.cold], rainThresh[2])))

    # x.pos.rainThresh <- sort(c(pos.warm[x.pos.rainThresh.warm], pos.cold[x.pos.rainThresh.cold]))

  } else {
    x.pos.rainThresh <- which(do.call(op.rainThresh, list(x, rainThresh[1])))
  }

  if(!is.null(dates))
  {
    ## thresholding x based on season and season threshold
    # warm
    if(length(x.pos.rainThresh.warm) > 0)
    {
      x.pos.rainThresh.C.warm <- split(seq_along(along.with = x.pos.rainThresh.warm), cumsum(c(0, diff(x.pos.rainThresh.warm) > 1)))
      x.pos.rainThresh.CLen.warm <- unlist(x.pos.rainThresh.C.warm[which(do.call(op.rainOffLength, list(lengths(x.pos.rainThresh.C.warm), rainOffLength[1])))])

      x.pos.rainThresh.index.warm <- pos.warm[x.pos.rainThresh.warm[x.pos.rainThresh.CLen.warm]]
    } else {
      x.pos.rainThresh.index.warm <- NULL
    }


    # cold
    if(length(x.pos.rainThresh.cold) > 0)
    {
      x.pos.rainThresh.C.cold <- split(seq_along(along.with = x.pos.rainThresh.cold), cumsum(c(0, diff(x.pos.rainThresh.cold) > 1)))
      x.pos.rainThresh.CLen.cold <- unlist(x.pos.rainThresh.C.cold[which(do.call(op.rainOffLength, list(lengths(x.pos.rainThresh.C.cold), rainOffLength[2])))])

      x.pos.rainThresh.index.cold <- pos.cold[x.pos.rainThresh.cold[x.pos.rainThresh.CLen.cold]]
    } else {
      x.pos.rainThresh.index.cold <- NULL
    }

    x.pos.rainThresh.index <- sort(c(x.pos.rainThresh.index.warm, x.pos.rainThresh.index.cold))
  } else {

    ## get consecutive positions under rainTreshold | lengths of event OP rainOffLength-Threshold
    # ## get consecutive positions under rainTreshold | lengths of event OP rainOffLength-Threshold
    x.pos.rainThresh.C <- split(seq_along(along.with = x.pos.rainThresh), cumsum(c(0, diff(x.pos.rainThresh) > 1)))
    x.pos.rainThresh.CLen <- unlist(x.pos.rainThresh.C[which(do.call(op.rainOffLength, list(lengths(x.pos.rainThresh.C), rainOffLength)))])

    x.pos.rainThresh.index <- x.pos.rainThresh[x.pos.rainThresh.CLen]
  }


  return(x.pos.rainThresh.index)

} # end of function findIsolatedRainEvent








#' findRainEvent
#'
#' This function return a list of rain events.
#'
#' @param x vector containing precipitation
#' @param x.pos.dryPeriods indices of x containing dry periods. Result of findRainFallPosition.
#' @return
#' list containing rain events as specific indices of x.
#'
#' @export
findRainEvent <- function(x = x, x.pos.dryPeriods)
{
  # browser()

  # get rainy days from dry period
  x.pos.rainEvent <- setdiff(seq_along(along.with = x), x.pos.dryPeriods)

  # create consecutive number of rainy days
  x.pos.rainEvent.Len <- split(seq_along(along.with = x.pos.rainEvent), cumsum(c(0, diff(x.pos.rainEvent) > 1)))

  # get original index for rainy days
  x.pos.rainEvent.index <- lapply(x.pos.rainEvent.Len, function(x, y){y[x]}, y = x.pos.rainEvent) # get original indices

  return(x.pos.rainEvent.index)
}
