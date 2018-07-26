#' Find optimal hyper-parameters for function highPassThresholding
#'
#' This functions optimize the hyper-parameter of the function highPassThresholding using standard
#' classification error measurements. As comparison a valid inventory in form of a binary classification
#' must be given.
#'
#'
#' @param x raster object or input path of file
#' @param inventory binarized raster inventory of same extent and resolution as x
#' @param range.scale.factor range of different scale factors
#' @param range.threshold range of different thresholds
#' @param ... see highPassThresholding()
#' @param cores number of cores for parallel processing. Default: 1 (sequential)
#' @param NoData = -99999
#' @param path.runfile full path for run file with .txt as ending. Default: NULL
#' @param env.rsaga
#' @param quiet no outputs in console. Default: TRUE
#' @return
#'  data.frame with clssification measurements ().
#'
#' @note
#' \itemize{
#'   \item SAGA GIS must be installed.
#'
#' }
#'
#' @keywords high-pass filter, SAGA GIS
#'
#'
#' @export
#'
optHiPassThresh <- function(x, inventory, range.scale.factor, range.threshold, path.runfile = NULL, env.rsaga, path.save = tempdir(), NoData = -99999, ..., cores = 1, quiet = TRUE)
{
  ## get start time of process
  process.time.start <- proc.time()

  if(!raster::compareRaster(x, inventory, stopiffalse = FALSE))
  {
    stop('Input raster "x" and "inventory" are different either in extent, number of rows and columns, projection, resolution, and/or origin')
  }

  if(!is.null(path.runfile) && tools::file_ext(path.runfile) != "txt")
  {
    stop('Format of "path.runfile" is not correct. Should be end with .txt!')
  } else if(!is.null(path.runfile)) {
    cat(paste0("run", "\t", "total", "\t", "scale",  "\t", "threshold", "\t", "TN", "\t", "FN", "\t", "FP", "\t", "TP", "\t", "TPR", "\t", "TNR", "\t", "FNR", "\t", "FPR", "\t", "acc", "\t", "rndm_acc", "\t", "f1score", "\t", "quality", "\t", "Kappa", "\n"), file = path.runfile, append = TRUE)
  }

  if(class(x) == "RasterLayer")
  {
    path.x <- file.path(tempdir(), "input.sgrd")
    raster::writeRaster(x = x, filename = paste0(tools::file_path_sans_ext(path.x), ".sdat"), overwrite = TRUE, NAflag = NoData)
  } else {
    path.x <- x
  }


  ## init run function
  intern.runFun <- function(i, combi, input, inventory, path.runfile, env.rsaga = env.rsaga, path.save, NoData, ...){
    # result <- parallel::parApply(cl = cl, X = combi, MARGIN = 1, FUN = function(x, input, inventory, path.runfile, ...){

    # browser()
    combi.i <- combi[i,]

    # init variables
    scale.i <- combi.i[[1]]
    threshold.i <- combi.i[[2]]

    tryCatch({

      # get high-pass image
      hipass <- Lslide::hiPassThresh(x = input, scale.factor = scale.i, threshold = threshold.i, do.use.temp.HPF = TRUE, env.rsaga = env.rsaga, path.save = path.save, NoData = NoData, ...)
      hipass.val <- unique(raster::values(hipass))

      if(length(hipass.val) == 1 && is.na(hipass.val)) {stop("Raster contains only NAs")}

      # resetNull data
      hipass <- raster::calc(x = hipass, fun = function(x){ifelse(is.na(x), 0, 2)})

      # overlay data
      result <- raster::overlay(x = hipass, y = inventory, fun = sum, na.rm = TRUE)

      # get statistics

      # ... create data frames
      stat <- table(raster::values(result))
      df.scale <- data.frame(scale = scale.i, threshold = threshold.i)
      df.stat <- data.frame(matrix(ncol = length(stat), nrow = 1))
      colnames(df.stat) <- names(stat)
      df.stat[1, ] <- stat
      df.stat <- cbind(df.scale, df.stat)

      # ... # 0: TRUE NEGATIVES (TN), 3: TRUE POSITIVES (TP), 1: FALSE NEGATIVES, 2: FALSE POSITIVES
      pos.TN <- grep(pattern = "0", x = colnames(df.stat))
      pos.FN <- grep(pattern = "1", x = colnames(df.stat))
      pos.FP <- grep(pattern = "2", x = colnames(df.stat))
      pos.TP <- grep(pattern = "3", x = colnames(df.stat))

      if(length(pos.TN) > 0) {
        colnames(df.stat)[pos.TN] <- "TN"
      } else {
        df.stat$TN <- 0
      }

      if(length(pos.FN) > 0) {
        colnames(df.stat)[pos.FN] <- "FN"
      } else {
        df.stat$FN <- 0
      }

      if(length(pos.FP) > 0) {
        colnames(df.stat)[pos.FP] <- "FP"
      } else {
        df.stat$FP <- 0
      }

      if(length(pos.TP) > 0) {
        colnames(df.stat)[pos.TP] <- "TP"
      } else {
        df.stat$TP <- 0
      }

      # ... confusion matrix error measurements
      df.stat$TPR <- tryCatch({df.stat$TP/(df.stat$TP + df.stat$FN)}, error = function(e){return(NA)})  # sensitivity, recall, hit rate, or true positive rate (TPR)
      df.stat$TNR <- tryCatch({df.stat$TN/(df.stat$TN + df.stat$FP)}, error = function(e){return(NA)})  # specificity or true negative rate (TNR)
      df.stat$FNR <- tryCatch({df.stat$FN/(df.stat$FN + df.stat$TP)}, error = function(e){return(NA)})  # miss rate or false negative rate (FNR)
      df.stat$FPR <- tryCatch({df.stat$FP/(df.stat$FP + df.stat$TN)}, error = function(e){return(NA)})  # miss rate or false negative rate (FPR)
      df.stat$acc <- tryCatch({(df.stat$TP + df.stat$TN)/(df.stat$TP + df.stat$TN + df.stat$FP + df.stat$FN)}, error = function(e){return(NA)})  # accuracy (ACC)

      # ... hard-coded data.frame numbers here!
      df.stat$rndm_acc <- tryCatch({((as.numeric(df.stat$TN + df.stat$FP) * as.numeric(df.stat$TN + df.stat$FN)) + (as.numeric(df.stat$FN + df.stat$TP) * as.numeric(df.stat$FP + df.stat$TP)))/(as.numeric(sum(df.stat[1, 3:6])) * as.numeric(sum(df.stat[1, 3:6])))}, error = function(e){return(NA)})
      df.stat$f1score <- tryCatch({(2*df.stat$TP)/((2*df.stat$TP) + df.stat$FP + df.stat$FN)}, error = function(e){return(NA)})  # F1 score is the harmonic mean of precision and sensitivity

      df.stat$quality <- tryCatch({(df.stat$TP)/(df.stat$TP + df.stat$FP + df.stat$FN)}, error = function(e){return(NA)})  # Tarolli et al. 2012:77, Heipke et al. (1997)


      # J. Richard Landis and Gary G. Koch - The Measurement of Observer Agreement for Categorical Data, Biometrics, Vol. 33, No. 1 (Mar., 1977), pp. 159-174.
      # http://standardwisdom.com/softwarejournal/2011/12/confusion-matrix-another-single-value-metric-kappa-statistic/
      df.stat$Kappa <- tryCatch({(df.stat$acc - df.stat$rndm_acc)/(1 - df.stat$rndm_acc)}, error = function(e){return(NA)})


      if(!is.null(path.runfile))
      {
        cat(paste0(i, "\t", nrow(combi), "\t", scale.i, "\t", threshold.i, "\t", df.stat$TN, "\t", df.stat$FN, "\t", df.stat$FP, "\t", df.stat$TP,
                   "\t", df.stat$TPR, "\t", df.stat$TNR, "\t", df.stat$FNR, "\t", df.stat$FPR, "\t", df.stat$acc, "\t", df.stat$rndm_acc, "\t", df.stat$f1score, "\t", df.stat$quality, "\t", df.stat$Kappa, "\n"), file = path.runfile, append=TRUE)
      }

      return(df.stat)

    }, error = function(e){
      if(!is.null(path.runfile))
      {
        cat(paste0(i, "\t", nrow(combi), "\t", scale.i, "\t", threshold.i, "\t", NA, "\t", NA, "\t", NA, "\t", NA, "\t", NA, "\t", NA,
                 "\t", NA, "\t", NA, "\t", NA, "\t", NA, "\t", NA, "\t", NA, "\t", NA, "\n"), file = path.runfile, append=TRUE)
      }

      df.scale <- data.frame(scale = scale.i, threshold = threshold.i)
      df.stat <- data.frame(matrix(ncol = 13, nrow = 1))
      colnames(df.stat) <- c("TN",  "FN",  "FP",  "TP",  "TPR",  "TNR",  "FNR",  "FPR",  "acc",  "rndm_acc",  "f1score",  "quality",  "Kappa")
      df.stat[1, ] <- NA
      df.stat <- cbind(df.scale, df.stat)

      return(df.stat)
    })
 } # end of runFun


  ## init high-pass function
  intern.hiPassFun <- function(path.x, path.save, env.rsaga, scale.factor, show.output.on.console = FALSE)
  {
    scale.txt <- gsub(pattern = "\\.", replacement = "", x = as.character(scale.factor))
    path.hipass <- file.path(path.save, paste0("hipass_", scale.txt, ".sgrd"))

    # RSAGA::rsaga.get.usage(lib = "grid_filter", module = 11, env = env.rsaga)
      RSAGA::rsaga.geoprocessor(lib = "grid_filter", module = 11, env = env.rsaga, show.output.on.console = show.output.on.console, param = list(
        GRID = path.x, HIPASS = path.hipass, SCALE = scale.factor))
  } # end of intern.hiPassFun




  ## get all combinations
  combi <- expand.grid(threshold = range.threshold, scale = range.scale.factor)[, c(2, 1)] # scale | threshold

  # sequential processing
  if(cores == 1){
    if(quiet == FALSE) cat("... start finding optimal hyper-parameters (sequentiel) \n")
    if(quiet == FALSE) cat("... ... generation of high-pass filtered images \n")
    noOutput <- lapply(X = range.scale.factor, FUN = intern.hiPassFun, path.x = path.x,  env.rsaga = env.rsaga, path.save = path.save)

    if(quiet == FALSE) cat("... ... generation of thesholding images \n")
    result <- lapply(X = 1:nrow(combi), FUN = intern.runFun, combi = combi, input = path.x, inventory = inventory, path.runfile = path.runfile, env.rsaga = env.rsaga, path.save = path.save, NoData = NoData, ...)
  }


  # parallel processing
  if(cores > 1){

    if(quiet == FALSE) cat("... start finding optimal hyper-parameters (parallel) \n")

    ## init parallel
    switch(Sys.info()[[1]],
           Windows = type  <- "PSOCK",
           Linux = type  <- "FORK",
           Mac = type <- "FORK")

    cl <- parallel::makeCluster(cores, type = type)

    if(Sys.info()[[1]] == "Windows")
    {
      parallel::clusterExport(cl = cl, varlist = names(environment()), envir = environment())
    }

    if(quiet == FALSE) cat("... ... generation of high-pass filtered images \n")
    noOutput <- parallel::parLapply(cl = cl, X = range.scale.factor, fun = intern.hiPassFun, path.x = path.x,  env.rsaga = env.rsaga, path.save = path.save)

    if(quiet == FALSE) cat("... ... generation of thesholding images \n")
    result <- parallel::parLapply(cl = cl, X = 1:nrow(combi), fun = intern.runFun, combi = combi, input = path.x, inventory = inventory, path.runfile = path.runfile, env.rsaga = env.rsaga, path.save = path.save, NoData = NoData, ...) # end of (par)apply

    parallel::stopCluster(cl)
  }

  # ... row bind lists
  result <- dplyr::bind_rows(result)


  process.time.run <- proc.time() - process.time.start
  if(quiet == FALSE) cat("------ Run of optHiPassThresh: " , round(x = process.time.run["elapsed"][[1]]/60, digits = 3), " Minutes ------\n")

  return(result)

} # end of function optHiPassThresh
