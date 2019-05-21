#' @title Get rain event of a landslide
#'
#' @description This function calcuates different precipitation characteristics for a specific time-series:
#' total precipitation, number of rainfall events, weighted mean intensitiy of rainfall events (normalized by MAP, RD or RDN),
#' cumulative critical event rainfall (normalized by MAP, RD or RDN), maximum rainfall during critical rainfall event,
#' duration of critical rainfall event, critical rainfall intensitiy (normalized by MAP, RD or RDN), rainfall at day of failure (start date),
#' rainfall intensity at day of failure (start date), maximum rainfall at day of failure (start date).
#'
#' @param Re vector containing the rain event variable, e.g. cumulated event rainfall (in mm) or intensity (mm/h)
#' @param D vector containing the duration of the rainfall events
#' @param method method to compute threshold. Either "LS" for least square, or "NLS" for non-linear least squares Method. Default: "nls"
#' @param prob.threshold exceedance probability level. Default: 0.05 (5 [percent] )
#' @param log10.transform log-transformation of input vectors Re and/or D. Default: TRUE
#' @param bootstrapping If TRUE bootstrapping is performed. Default: TRUE
#' @param R the number of bootstrap replicates, see boot::boot() for more information. Default: 1000
#' @param seed replicable bootstrapping. Default: 123
#' @param use.integralError for estimating x of prob.threshold, the entire function is integrated first to estimate the bias: (1 - INTEGRAL)/2. Default. TRUE
#' @param ... more options passed to boot:boot() function, i.e. parallel for paralell processing
#'
#' @return vector containing rainfall metrics (see description). If return.DataFrame is TRUE a data.frame is returned containing similar rain metrics for all rain events.
#'
#'
#' @note
#' \itemize{
#'   \item Brunetti, M. T., Peruccacci, S., Rossi, M., Luciani, S., Valigi, D., & Guzzetti, F. (2010). Rainfall thresholds for the possible occurrence of landslides in Italy. Natural Hazards and Earth System Sciences, 10(3), 447.
#'   \item Peruccacci, S., Brunetti, M. T., Luciani, S., Vennari, C., & Guzzetti, F. (2012). Lithological and seasonal control on rainfall thresholds for the possible initiation of landslides in central Italy. Geomorphology, 139, 79-90.
#'   \item Rossi, M., Luciani, S., Valigi, D., Kirschbaum, D., Brunetti, M. T., Peruccacci, S., & Guzzetti, F. (2017). Statistical approaches for the definition of landslide rainfall thresholds and their uncertainty using rain gauge and satellite data. Geomorphology, 285, 16-27.
#'   \item Guzzetti, F., Peruccacci, S., Rossi, M., & Stark, C. P. (2007). Rainfall thresholds for the initiation of landslides in central and southern Europe. Meteorology and atmospheric physics, 98(3-4), 239-267.
#' }
#'
#' @keywords rainfall tresholds, rainfall event, landslide, automatic appraoch
#'
#'
#' @export
getRainThreshFreqM <- function(Re, D, method = "NLS", prob.threshold = 0.05, log10.transform = FALSE,
                               bootstrapping = TRUE, R = 1000, seed = 123, use.integralError = TRUE, ...){

  ## log transformation of input vectors
    if(log10.transform)
  {
    Re <- log(x = Re, base = 10) # cumultive precipitation of rain event
    D <- log(x = D, base = 10) # duration of rain event
  }



  if(method == "LS"){
    ## Least Squares Method (LS)
    # ReD.LS <- stats::lsfit(y = Re, x = D) # delivers same result as lm()
    if(bootstrapping){
      # ... create bootstrap function
      boot.LS <- function(formula, data, indices){
        # subset data
        data <- data[indices,]

        # boot model
        LM.boot <- stats::lm(formula, data)

        return(coef(LM.boot))
      }

      # set seed for replication
      set.seed(seed)

      # bootstrapping with R replications
      ReD.LS.boot <- boot::boot(data = data.frame(Re = Re, D = D), statistic = boot.LS,
                      R = R, formula = Re ~ D)

      # boot::boot.ci(ReD.LS.boot, conf= c(.05, 0.5, .95))

      # ... according to Rossi et al. (2017: 20)
      # compute quantiles: 5th for min, 50 as median for best fit, 95th for max
      ReD.LS.boot.conf <- apply(X = ReD.LS.boot$t, MARGIN = 2, FUN = quantile, probs = c(.05, .50, .95), na.rm=TRUE)
      # ReD.LS.boot.stat <- apply(X = ReD.LS.boot$t, MARGIN = 2, FUN = median, na.rm=TRUE)

      # ... intercept alpha from model
      # T50 <- ReD.LS.boot.stat[1]
      T5 <- ReD.LS.boot.conf[1, 1]
      T50 <- ReD.LS.boot.conf[2, 1]
      T95 <- ReD.LS.boot.conf[3, 1]


      # ... slope gamma from model
      # gamma <- ReD.LS.boot.stat[2]
      gamma <- c(ReD.LS.boot.conf[2, 2], ReD.LS.boot.conf[1, 2], ReD.LS.boot.conf[3, 2])
      names(gamma) <- c("gamma_Tx_median", "gamma_Tx_min", "gamma_Tx_max")

      # Getting residuals of fit
      # ... defining linear function
      linFunct <- function(x, intercept, gamma){return((intercept + gamma * x))}

      ## ... fitting median model and get residuals and error estimates
      # ReD.Model.fit <- linFunct(x = D, intercept = ReD.LS.boot.stat[1], gamma = ReD.LS.boot.stat[2])
      ReD.Model.fit <- linFunct(x = D, intercept = ReD.LS.boot.conf[2, 1], gamma = ReD.LS.boot.conf[2, 2])

      ReD.Model.Res <- Re - ReD.Model.fit

      # ... computing error estimates using http://pages.mtu.edu/~fmorriso/cm3215/UncertaintySlopeInterceptOfLeastSquaresFit.pdf
      # n <- length(D)
      # SSE <- sum(ReD.Model.Res^2, na.rm = TRUE) # Error sum of squares
      # SST <- sum((Re - mean(Re))^2, na.rm = TRUE) # Total sum of squares
      # SSR <- SST - SSE # Regression sum of squares
      # Sxx <- sum((D - mean(D))^2, na.rm = TRUE)
      # Sxy <- sum((D - mean(D)) * (Re - mean(Re)), na.rm = TRUE)
      # SDyx <- sqrt(SSE/(n-2)) # Standard Deviation of y(x)
      # SDm <- sqrt(SDyx^2/Sxx) # Standard Deviation of Slope
      # SDb <- sqrt(SDyx^2 * (1/n + ((mean(D, na.rm = TRUE)^2)/Sxx))) # Standard Deviation of Intercept

      # ... get 5th and 95th percentile of gamma (slope) and intercept
      # ... gamma (slope)
      # ReD.Model.gamma.quant <- ReD.LS.boot.stat[2] + qt(c((.05/2), .95+(0.05/2)), df = (n-2)) * SDm

      # ... intercept
      # ReD.Model.interc.quant <- ReD.LS.boot.stat[1] + qt(c((.05/2), .95+(0.05/2)), df = (n-2)) * SDb
      # T5 <- ReD.Model.interc.quant[1]
      # T95 <- ReD.Model.interc.quant[2]

      ## ... fitting min and max model and get residuals
      # min
      ReD.Model.fit.min <- linFunct(x = D, intercept = ReD.LS.boot.conf[1, 1], gamma = ReD.LS.boot.conf[1, 2])
      # ReD.Model.fit.min <- linFunct(x = D, intercept = ReD.Model.interc.quant[1], gamma = ReD.Model.gamma.quant[1])
      ReD.Model.Res.min <- Re - ReD.Model.fit.min

      # max
      ReD.Model.fit.max <- linFunct(x = D, intercept = ReD.LS.boot.conf[3, 1], gamma = ReD.LS.boot.conf[3, 2])
      # ReD.Model.fit.max <- linFunct(x = D, intercept = ReD.Model.interc.quant[2], gamma = ReD.Model.gamma.quant[2])
      ReD.Model.Res.max <- Re - ReD.Model.fit.max

      # plot(Re ~ D)
      # abline(lm(Re ~ D))
      # abline(a = ReD.LS.boot.stat[1], b = ReD.LS.boot.stat[2], col = "red")
      # abline(a = ReD.LS.boot.stat[1] - ReD.Model.interc.quant[1], b = ReD.LS.boot.stat[2], col = "orange")

    } else {
      ReD.Model <- stats::lm(Re ~ D)
      ReD.Model.fit <- ReD.Model$fitted.values
      ReD.Model.Res <- ReD.Model$residuals
      # ReD.LS <- stats::lsfit(y = Re, x = D)

      # intercept alpha from model
      T50 <- coef(ReD.Model)[[1]]

      # slope gamma from model
      gamma <- coef(ReD.Model)[[2]]
      names(gamma) <- "gamma"

    }# end of if bootstrapping

  } else if(method == "NLS") { # end of if method == "LS"

      ## Nonlinear Least Squares Method (NLS)

      # ... optimizer function
      opt.NLS <- function(par, x, y)
      {
        t <- par[1]
        alpha <- par[2]
        gamma <- par[3]
        y.fit <- t + alpha * (x^gamma)
        sum((y - y.fit)^2)
      }

      # ... optimize parameter
      par.opt.NLS <- optim(x = D, y = Re, par = c(mean(D), sd(D), 1),
                           fn = opt.NLS, method = "Nelder-Mead", #, "BFGS",
                           control = list(pgtol = 1e-9, maxit = 10000))$par

      if(bootstrapping)
      {

        # ... create bootstrap function
        boot.NLS <- function(formula, data, start, control, indices)
        {
          # subset data
          data <- data[indices,]

          # boot model
          NLS.boot <- minpack.lm::nlsLM(formula = formula, data = data,
                                      start = start, control = control)

          return(coef(NLS.boot))
        } # end of boot.NLS

        # set seed for replication
        set.seed(seed)

        # bootstrapping with R replications
        ReD.NLS.boot <- boot::boot(data = data.frame(y = Re, x = D), statistic = boot.NLS,
                                 R = R, formula = y ~ t + a * (x^gamma), control = list(maxiter = 500),
                                 start = list(t = par.opt.NLS[1], a = par.opt.NLS[2], gamma = par.opt.NLS[3]))


        # ... according to Rossi et al. (2017: 20)
        # compute quantiles: 5th for min, 50 as median for best fit, 95th for max
        ReD.NLS.boot.conf <- apply(X = ReD.NLS.boot$t, MARGIN = 2, FUN = quantile, probs = c(.05, .50, .95), na.rm=TRUE)
        T5 <- ReD.NLS.boot.conf[1, 1]
        T50 <- ReD.NLS.boot.conf[2, 1]
        T95 <- ReD.NLS.boot.conf[3, 1]

        gamma <- c(ReD.NLS.boot.conf[1, 2], ReD.NLS.boot.conf[2, 2], ReD.NLS.boot.conf[3, 2])
        names(gamma) <- c("gamma_Tx_min", "gamma_Tx_median", "gamma_Tx_max")


        # Getting residuals of fit
        # ... defining linear function
        powerLawFunct <- function(x, t, a, gamma){return((t + a * (x^gamma)))}

        ## ... fitting median model and get residuals
        ReD.Model.fit <- powerLawFunct(x = D, t = ReD.NLS.boot.conf[2, 1], a = ReD.NLS.boot.conf[2, 2], gamma = ReD.NLS.boot.conf[2, 3])
        ReD.Model.Res <- Re - ReD.Model.fit

        ## ... fitting min and max model and get residuals
        # min
        ReD.Model.fit.min <- powerLawFunct(x = D, t = ReD.NLS.boot.conf[1, 1], a = ReD.NLS.boot.conf[1, 2], gamma = ReD.NLS.boot.conf[1, 3])
        ReD.Model.Res.min <- Re - ReD.Model.fit.min

        # max
        ReD.Model.fit.max <- powerLawFunct(x = D, t = ReD.NLS.boot.conf[3, 1], a = ReD.NLS.boot.conf[3, 2], gamma = ReD.NLS.boot.conf[3, 3])
        ReD.Model.Res.max <- Re - ReD.Model.fit.max

      } else{

        # get model
        ReD.NLS <- nls(formula = y ~ t + a * (x^gamma), data = data.frame(x = D, y = Re),
                                    start = list(t = par.opt.NLS[1], a = par.opt.NLS[2], gamma = par.opt.NLS[3]),
                                    control = list(maxiter = 500))


        ReD.NLS <- minpack.lm::nlsLM(formula = y ~ t + a * (x^gamma), data = data.frame(x = D, y = Re),
                                    start = list(t = par.opt.NLS[1], a = par.opt.NLS[2], gamma = par.opt.NLS[3]),
                                    control = nls.lm.control(maxiter = 500))

        ReD.NLS.fit <- predict(ReD.NLS)
        # ReD.Model.fit <- powerLawFunct(x = D, t = coef(ReD.NLS)[1], a = coef(ReD.NLS)[2], gamma = coef(ReD.NLS)[3])

      }
  } else { # end of if method == "NLS"
    stop('Selected method is not supported. Please use "NLS" or "LS"!')
  } # end of if method

  # PLOT NLS
  # plot(y = Re, x = D)
  # a<-coef(ReD.NLS)[1]
  # b<-coef(ReD.NLS)[2]
  # k<-coef(ReD.NLS)[3]
  # x <- seq(from = 0, to = 5, by = 0.001)
  # lines(x = x, a+b*x^k,col='red')
  # ggplot2::ggplot(data.frame(D = D, ReD.NLS.fit = ReD.NLS.fit),
  #                 ggplot2::aes(D,ReD.NLS.fit)) + ggplot2::geom_point() + ggplot2::geom_smooth()




  # IS THAT NECESAIRY?
  ## standardize residuals
  # ReD.Model.SDRes <- ReD.Model$residuals
  # ReD.Model.SDRes <- scale(ReD.Model$residuals)[, 1]


  ## check normality of residuals (limit: large sample sizes easily produce significant results from small deviations from normality)
  # if(shapiro.test(ReD.Model.Res)$p.value <= 0.005){warning("Residuals are highly significantly not normal distributed")}
  # if(bootstrapping && shapiro.test(ReD.Model.Res.min)$p.value <= 0.005){warning("Residuals are highly significantly not normal distributed")}
  # if(bootstrapping &&shapiro.test(ReD.Model.Res.max)$p.value <= 0.005){warning("Residuals are highly significantly not normal distributed")}



  ## Probability Densitiy Function (PDF) with Kernel Density Estimation (KDE) with a Gaussian function
  ReD.PDF <- stats::density(x = ReD.Model.Res, kernel = "gaussian")
  ReD.PDF.df <- data.frame(x = ReD.PDF$x, y = ReD.PDF$y, type = "PDF")


  if(bootstrapping)
  {
    # min
    ReD.PDF.min <- stats::density(x = ReD.Model.Res.min, kernel = "gaussian")
    ReD.PDF.df.min <- data.frame(x = ReD.PDF.min$x, y = ReD.PDF.min$y, type = "PDF")
    # plot(ReD.PDF.min)

    # max
    ReD.PDF.max <- stats::density(x = ReD.Model.Res.max, kernel = "gaussian")
    ReD.PDF.df.max <- data.frame(x = ReD.PDF.max$x, y = ReD.PDF.max$y, type = "PDF")
    # plot(ReD.PDF.max)
  }



  ## Modelling PDF using a Gaussian Function
  GaussianFunction <- function(ReD.PDF)
  {
    # ... using Nonlinear regression model (NLS)
    # # https://stats.stackexchange.com/questions/220109/fit-a-gaussian-to-data-with-r-with-optim-and-nls

    # ... optimizer function
    opt.Gauss <- function(par, x, y)
    {
      m <- par[1]
      sd <- par[2]
      k <- par[3]
      rhat <- k * exp(-0.5 * ((x - m)/sd)^2)
      sum((y - rhat)^2)
    }

    # ... optimize parameter
    par.opt.Gauss <- stats::optim(x = ReD.PDF$x, y = ReD.PDF$y, par = c(mean(ReD.PDF$x), sd(ReD.PDF$x), 1), fn = opt.Gauss,
                           method = "BFGS")


    # ... comute NLS using a Gaussian Function and fit values
    ReD.PDF.NLS <- stats::nls(y ~ k*exp(-1/2*(x-mu)^2/sigma^2),
                      start = c(mu = par.opt.Gauss$par[1], sigma = par.opt.Gauss$par[2], k = par.opt.Gauss$par[3]),
                      data = data.frame(x = ReD.PDF$x, y = ReD.PDF$y), control = list(maxiter = 10000, reltol=1e-9))


    # ... return model and fitted values
    return(list(stats::predict(ReD.PDF.NLS), ReD.PDF.NLS))

  } # end of Gaussian Function


  ReD.PDF.result <- GaussianFunction(ReD.PDF = ReD.PDF)
  ReD.PDF.fit <- ReD.PDF.result[[1]]
  ReD.PDF.NLS <- ReD.PDF.result[[2]]

  if(bootstrapping)
  {
    ReD.PDF.result.min <- GaussianFunction(ReD.PDF = ReD.PDF.min)
    ReD.PDF.fit.min <- ReD.PDF.result.min[[1]]
    ReD.PDF.NLS.min <- ReD.PDF.result.min[[2]]

    ReD.PDF.result.max <- GaussianFunction(ReD.PDF = ReD.PDF.max)
    ReD.PDF.fit.max <- ReD.PDF.result.max[[1]]
    ReD.PDF.NLS.max <- ReD.PDF.result.max[[2]]
  }

  # plotting density curves
  # ReD.PDF.fit.df <- data.frame(x = ReD.PDF$x, y = ReD.PDF.fit, type = "Gaussian fit")
  # ReD.PDF.df <- rbind(ReD.PDF.df, ReD.PDF.fit.df)
  #
  # ggExample <- ggplot2::ggplot(data = ReD.PDF.df, ggplot2::aes(x = x, y = y, color = type, linetype = type)) +
  #   ggplot2::geom_line(size = 0.9) +
  #   ggplot2::scale_color_manual(values = c("royalblue4", "black")) +
  #   ggplot2::scale_linetype_manual(values = c(2, 1))
  #
  # ggExample


  ## Extracting the x procent proceeding probability
  # https://stackoverflow.com/questions/44313022/how-to-compute-confidence-interval-of-the-fitted-value-via-nls
  # https://stackoverflow.com/questions/37455512/predict-x-values-from-simple-fitting-and-annoting-it-in-the-plot?rq=1
  # https://stackoverflow.com/questions/43322568/predict-x-value-from-y-value-with-a-fitted-model

  # integrate.xy is not 1 for the whole function. Therehore, the error is estimated first and divided by 2 for both sides
  if(use.integralError)
  {
    integral.error <-  (1 - sfsmisc::integrate.xy(x = ReD.PDF$x, fx = ReD.PDF.fit, a = min(ReD.PDF$x), b = max(ReD.PDF$x)))/2
  } else{
    integral.error <- 0
  }


  integrFunctToThresh <- function(ReD.PDF, ReD.PDF.fit, thresh, lower, upper, par, integral.error)
  {
    # ... optimizer function for integration
    opt.Integral <- function(par, x, y, thresh)
    {
      b <- par[1]
      abs(sfsmisc::integrate.xy(x = x, fx = y, b = b) - thresh)
    }


    # ... find x parameter to area threshold of integral
    par.opt <- optim(thresh = (thresh - integral.error), x = ReD.PDF$x, y = ReD.PDF.fit,
                          lower = lower, upper = upper,
                          par = par, # start value
                          fn = opt.Integral, method = "L-BFGS-B",
                          control = list(pgtol = 1e-9, maxit = 10000, ndeps = 1e-9))$par

    # sfsmisc::integrate.xy(x = ReD.PDF$x, fx = ReD.PDF.fit, b = par.opt.Mean)
    # predict(object = ReD.PDF.NLS, newdata = data.frame(x = par.opt.Mean))

    # return correspond to x which fulfills threshold
    return(par.opt)
  }

  fit.Gauss.50 <- integrFunctToThresh(ReD.PDF = ReD.PDF, ReD.PDF.fit = ReD.PDF.fit, thresh = 0.5, integral.error = integral.error,
                                      lower = (mean(ReD.PDF$x) - sd(ReD.PDF$x)), upper = (mean(ReD.PDF$x) + sd(ReD.PDF$x)), par = mean(ReD.PDF$x))

  fit.Gauss.x <- integrFunctToThresh(ReD.PDF = ReD.PDF, ReD.PDF.fit = ReD.PDF.fit, thresh = prob.threshold, integral.error = integral.error,
                                     lower = (min(ReD.PDF$x) + 0.0001), upper = mean(ReD.PDF$x), par = mean(c(min(ReD.PDF$x), mean(ReD.PDF$x)))) # if b becomes the min value, then in sfsmisc::integrate.xy a = b fails!


  if(bootstrapping)
  {
    # min
    fit.Gauss.50.min <- integrFunctToThresh(ReD.PDF = ReD.PDF.min, ReD.PDF.fit = ReD.PDF.fit.min, thresh = 0.5, integral.error = integral.error,
                                        lower = (mean(ReD.PDF.min$x) - sd(ReD.PDF.min$x)), upper = (mean(ReD.PDF.min$x) + sd(ReD.PDF.min$x)), par = mean(ReD.PDF.min$x))

    fit.Gauss.x.min <- integrFunctToThresh(ReD.PDF = ReD.PDF.min, ReD.PDF.fit = ReD.PDF.fit.min, thresh = prob.threshold, integral.error = integral.error,
                                       lower = (min(ReD.PDF.min$x) + 0.0001), upper = mean(ReD.PDF.min$x), par = mean(c(min(ReD.PDF.min$x), mean(ReD.PDF.min$x))))

    # max
    fit.Gauss.50.max <- integrFunctToThresh(ReD.PDF = ReD.PDF.max, ReD.PDF.fit = ReD.PDF.fit.max, thresh = 0.5, integral.error = integral.error,
                                            lower = (mean(ReD.PDF.max$x) - sd(ReD.PDF.max$x)), upper = (mean(ReD.PDF.max$x) + sd(ReD.PDF.max$x)), par = mean(ReD.PDF.max$x))

    fit.Gauss.x.max <- integrFunctToThresh(ReD.PDF = ReD.PDF.max, ReD.PDF.fit = ReD.PDF.fit.max, thresh = prob.threshold, integral.error = integral.error,
                                         lower = (min(ReD.PDF.max$x) + 0.0001), upper = mean(ReD.PDF.max$x), par = mean(c(min(ReD.PDF.max$x), mean(ReD.PDF.max$x))))
  }



  # par.opt.Mean <- optim(thresh = (0.5 - integral.error), x = ReD.PDF$x, y = ReD.PDF.fit,
  #                       lower = (mean(ReD.PDF$x) - sd(ReD.PDF$x)),
  #                       upper = (mean(ReD.PDF$x) + sd(ReD.PDF$x)),
  #                       par = mean(ReD.PDF$x), # start value
  #                       fn = opt.Integral, method = "L-BFGS-B",
  #                       control = list(pgtol = 1e-9, maxit = 10000, ndeps = 1e-9))$par



  # ... ... (2) integration for probablilty level threshold
  # fit.Gauss.x  <- optim(thresh = (prob.threshold - integral.error), x = ReD.PDF$x, y = ReD.PDF.fit,
  #                           lower = (min(ReD.PDF$x) + 0.0000001), #  if b becomes the min value, then in sfsmisc::integrate.xy a = b fails!
  #                           upper = mean(ReD.PDF$x),
  #                           par = mean(min(ReD.PDF$x), (mean(ReD.PDF$x) - sd(ReD.PDF$x))), #  start value
  #                           fn = opt.Integral, method = "L-BFGS-B",
  #                           control = list(pgtol = 1e-9, maxit = 10000, ndeps = 1e-9))$par


  # sfsmisc::integrate.xy(x = ReD.PDF$x, fx = ReD.PDF.fit, a = min(ReD.PDF$x), b = par.opt.probLevel)
  # ggExample + ggplot2::geom_segment(inherit.aes = FALSE, color = "red", size = 0.8,
  #                                   ggplot2::aes(x = fit.Gauss.x , xend = fit.Gauss.x,
  #                                                y = 0, yend = predict(object = ReD.PDF.NLS, newdata = data.frame(x = fit.Gauss.x))))

  ## get Intercept of probablilty level
  # Tx is the curve parallel to the best-fit line T50 (slope = gamma), with intercept alpha.x= alpha.50−alpha.opt
  Tx <- T50 - (fit.Gauss.50 - fit.Gauss.x)
  # alpha <- exp(Tx) # re-transform log alpha to alpha ?

  if(bootstrapping)
  {
    Tx.min <- T5 - (fit.Gauss.50.min - fit.Gauss.x.min)
    Tx.max <- T95 - (fit.Gauss.50.max - fit.Gauss.x.max)
  }


  ### LS
  ## create data frame for ggplot2
  # df.fit.LS <- data.frame(x = D, y = Re, y_fit = ReD.Model.fit,
  #                  y_fit_min = ReD.Model.fit.min, y_fit_max = ReD.Model.fit.max,
  #                  y_fit_Lx = linFunct(x = D, intercept = Tx, gamma = gamma),
  #                  y_fit_Lx_min = linFunct(x = D, intercept = Tx.min, gamma = ReD.Model.gamma.quant[1]),
  #                  y_fit_Lx_max = linFunct(x = D, intercept = Tx.max, gamma = ReD.Model.gamma.quant[2]))
  #
  #
  # ggplot2::ggplot(data = df.fit.LS) +
  #   ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y)) + # add Re-D point pairs
  #   # ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y_fit)) + # add fit of regression line (median fit)
  #   ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y_fit_Lx), col = "green") + # add fit of regression line (median fit)
  #   ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y_fit_Lx_min), col = "blue") + # add fit of regression line (median fit)
  #   ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y_fit_Lx_max), col = "red") + # add fit of regression line (median fit)
  #   ggplot2::xlim(0, 6) + ggplot2::ylim(-1, 6)
  #
  #
  # # ... plot best fit with 5-95th quantile
  # ggplot2::ggplot(data = df.fit.LS) +
  #   ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y)) + # add Re-D point pairs
  #   ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y_fit_min), col = "blue") + # add fit of regression line (median fit)
  #   ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y_fit), col = "black") +
  #   ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y_fit_max), col = "red") +
  #   ggplot2::xlim(0, 6) + ggplot2::ylim(-1, 6)
  #
  #
  # ### NLS
  # ## create data frame for ggplot2
  # df.fit.NLS <- data.frame(x = D, y = Re, y_fit = ReD.Model.fit,
  #                         y_fit_min = ReD.Model.fit.min, y_fit_max = ReD.Model.fit.max,
  #                         y_fit_Lx = powerLawFunct(x = D, t = Tx, a = ReD.NLS.boot.conf[2, 2], gamma = ReD.NLS.boot.conf[2, 3]),
  #                         y_fit_Lx_min = powerLawFunct(x = D, t = Tx.min, a = ReD.NLS.boot.conf[1, 2], gamma = ReD.NLS.boot.conf[1, 3]),
  #                         y_fit_Lx_max = powerLawFunct(x = D, t = Tx.max, a = ReD.NLS.boot.conf[3, 2], gamma = ReD.NLS.boot.conf[3, 3]))
  #
  #
  # ggplot2::ggplot(data = df.fit.NLS) +
  #   ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y)) + # add Re-D point pairs
  #   ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_min), col = "orange", method = "loess", size = 0.7) +
  #   ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_max), col = "yellow", method = "loess", size = 0.7) +
  #   ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit), col = "black", method = "loess", size = 0.7)+ # add fit of regression line (median fit)
  #   ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx), col = "green", method = "loess", size = 0.7) + # add fit of regression line (median fit)
  #   ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx_min), col = "blue", method = "loess", size = 0.7) + # add fit of regression line (median fit)
  #   ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx_max), col = "red", method = "loess", size = 0.7) + # add fit of regression line (median fit)
  #   ggplot2::xlim(0, 6) + ggplot2::ylim(-1, 6)
  #
  #
  # ggplot2::ggplot(data = df.fit.NLS) +
  #   ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y)) + # add Re-D point pairs
  #   ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit), col = "black", method = "loess", size = 0.7) +
  #   ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_min), col = "orange", method = "loess", size = 0.7) +
  #   ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_max), col = "yellow", method = "loess", size = 0.7) +
  #   ggplot2::xlim(0, 6) + ggplot2::ylim(-1, 6)
  #
  #
  # ggplot2::ggplot(data = df.fit.NLS) +
  #   ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y)) + # add Re-D point pairs
  #   ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx), col = "green", method = "loess", size = 0.7) + # add fit of regression line (median fit)
  #   ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx_min), col = "blue", method = "loess", size = 0.7) + # add fit of regression line (median fit)
  #   ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx_max), col = "red", method = "loess", size = 0.7) + # add fit of regression line (median fit)
  #   ggplot2::xlim(0, 6) + ggplot2::ylim(-1, 6)
  #
  #
  # ## combine both LS and NSL
  # ggplot2::ggplot(data = df.fit.NLS) +
  #   # NLS
  #   ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y)) + # add Re-D point pairs
  #   ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx), col = "yellow", method = "loess", size = 0.8) + # add fit of regression line (median fit)
  #   ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx_min), col = "blue", method = "loess", size = 0.7) + # add fit of regression line (median fit)
  #   ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx_max), col = "blue", method = "loess", size = 0.7) + # add fit of regression line (median fit)
  #   # LS
  #   ggplot2::geom_line(data = df.fit.LS, mapping = ggplot2::aes(x = x, y = y_fit_Lx), col = "orange") + # add fit of regression line (median fit)
  #   ggplot2::geom_line(data = df.fit.LS, mapping = ggplot2::aes(x = x, y = y_fit_Lx_min), col = "red") + # add fit of regression line (median fit)
  #   ggplot2::geom_line(data = df.fit.LS, mapping = ggplot2::aes(x = x, y = y_fit_Lx_max), col = "red") + # add fit of regression line (median fit)
  #   ggplot2::xlim(0, 6) + ggplot2::ylim(-1, 6)


  if(log10.transform)
  {
    alpha.Tx <- 10^Tx
    names(alpha.Tx) <- "alpha_Tx"

    alpha.res <- alpha.Tx


    if(bootstrapping)
    {
      alpha.Tx.min <- 10^Tx.min
      names(alpha.Tx.min) <- "alpha_Tx_min"

      alpha.Tx.max <- 10^Tx.max
      names(alpha.Tx.max) <- "alpha_Tx_max"

      alpha.res <- c(alpha.res, alpha.Tx.min, alpha.Tx.max)
    }
  } else {# end of log10 alpha check

    alpha.Tx <- Tx
    names(alpha.Tx) <- "alpha_Tx"

    alpha.res <- alpha.Tx

    if(bootstrapping)
    {
      alpha.Tx.min <- Tx.min
      names(alpha.Tx.min) <- "alpha_Tx_min"

      alpha.Tx.max <- Tx.max
      names(alpha.Tx.max) <- "alpha_Tx_max"

      alpha.res <- c(alpha.res, alpha.Tx.min, alpha.Tx.max)
    }
  }


  if(method == "LS")
  {
    funct <- linFunct
  }

  if(method == "NLS")
  {
    funct <- powerLawFunct
  }


  return(list(alpha.res, gamma, funct))

} # end of function getRainThreshFreqM
