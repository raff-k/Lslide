#' Get rain event of a landslide
#'
#' This function calcuates different precipitation characteristics for a specific time-series:
#' total precipitation, number of rainfall events, weighted mean intensitiy of rainfall events (normalized by MAP, RD or RDN),
#' cumulative critical event rainfall (normalized by MAP, RD or RDN), maximum rainfall during critical rainfall event,
#' duration of critical rainfall event, critical rainfall intensitiy (normalized by MAP, RD or RDN), rainfall at day of failure (start date),
#' rainfall intensity at day of failure (start date), maximum rainfall at day of failure (start date).
#'
#' @param E vector containing the cumulated event rainfall (in mm),
#' @param D vector containing the duration of the rainfall events
#' @param prob.threshold exceedance probability level. Default: 0.05 (5 %)
#' @param log.transform log-transformation of input vectors E and/or D. Default: c(TRUE, TRUe)
#' @param bootstrapping If TRUE bootstrapping is performed. Default: TRUE
#' @param R the number of bootstrap replicates, see boot::boot() for more information. Default: 1000
#' @param seed replicable bootstrapping. Default: 123
#' @param use.integralError for estimating x of prob.threshold, the entire function is integrated first to estimate the bias: (1 - INTEGRAL)/2. Default. TRUE
#' @param ... more options passed to boot:boot() function, i.e. parallel for paralell processing
#' @return
#' vector containing rainfall metrics (see description). If return.DataFrame is TRUE a data.frame is returned containing similar
#' rain metrics for all rain events.
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
#'
#' @keywords rainfall tresholds, rainfall event, landslide, automatic appraoch
#'
#'
#' @export
library(MASS)
getRainThreshFreqM <- function(E, D, prob.threshold = 0.05, log.transform = c(TRUE, TRUE),
                               bootstrapping = TRUE, R = 10000, seed = 123, use.integralError = TRUE, ...)
{
  E <- log.cLE.H09.sum
  E <- log.LE.H09.sum

  D <- log.cLE.H09.dur
  D <- log.LE.H09.dur

  ## log transformation of input vectors
  # cumultive precipitation of rain event
  if(log.transform[1])
  {
    E <- log(E)
  }

  # duration of rain event
  if(log.transform[1])
  {
    D <- log(D)
  }



  if(method == "LS")
  {
    ## Least Squares Method (LS)
    # ED.LS <- stats::lsfit(y = E, x = D) # delivers same result as lm()
    if(bootstrapping)
    {
      # ... create bootstrap function
      boot.LS <- function(formula, data, indices)
      {
        # subset data
        data <- data[indices,]

        # boot model
        LM.boot <- stats::lm(formula, data)

        return(coef(LM.boot))
      }

      # set seed for replication
      set.seed(seed)

      # bootstrapping with R replications
      ED.LS.boot <- boot::boot(data = data.frame(E = E, D = D), statistic = boot.LS,
                      R = R, formula = E ~ D)

      # boot::boot.ci(ED.LS.boot, conf= c(.05, 0.5, .95))

      # ... according to Rossi et al. (2017: 20)
      # compute quantiles: 5th for min, 50 as median for best fit, 95th for max
      ED.LS.boot.conf <- apply(X = ED.LS.boot$t, MARGIN = 2, FUN = quantile, probs = c(.05, .50, .95), na.rm=TRUE)
      ED.LS.boot.stat <- apply(X = ED.LS.boot$t, MARGIN = 2, FUN = median, na.rm=TRUE)

      # ... intercept alpha from model
      T50 <- ED.LS.boot.stat[1]

      # ... slope gamma from model
      gamma <- ED.LS.boot.stat[2]

      # Getting residuals of fit
      # ... defining linear function
      linFunct <- function(x, intercept, gamma){return((intercept + gamma * x))}

      ## ... fitting median model and get residuals and error estimates
      ED.Model.fit <- linFunct(x = D, intercept = ED.LS.boot.stat[1], gamma = ED.LS.boot.stat[2])
      ED.Model.Res <- E - ED.Model.fit

      # ... computing error estimates using http://pages.mtu.edu/~fmorriso/cm3215/UncertaintySlopeInterceptOfLeastSquaresFit.pdf
      n <- length(D)
      SSE <- sum(ED.Model.Res^2, na.rm = TRUE) # Error sum of squares
      SST <- sum((E - mean(E))^2, na.rm = TRUE) # Total sum of squares
      SSR <- SST - SSE # Regression sum of squares
      Sxx <- sum((D - mean(D))^2, na.rm = TRUE)
      Sxy <- sum((D - mean(D)) * (E - mean(E)), na.rm = TRUE)
      SDyx <- sqrt(SSE/(n-2)) # Standard Deviation of y(x)
      SDm <- sqrt(SDyx^2/Sxx) # Standard Deviation of Slope
      SDb <- sqrt(SDyx^2 * (1/n + ((mean(D, na.rm = TRUE)^2)/Sxx))) # Standard Deviation of Intercept

      # ... get 5th and 95th percentile of gamma (slope) and intercept
      # ... gamma (slope)
      ED.Model.gamma.quant <- ED.LS.boot.stat[2] + qt(c((.05/2), .95+(0.05/2)), df = (n-2)) * SDm

      # ... intercept
      ED.Model.interc.quant <- ED.LS.boot.stat[1] + qt(c((.05/2), .95+(0.05/2)), df = (n-2)) * SDb
      T5 <- ED.Model.interc.quant[1]
      T95 <- ED.Model.interc.quant[2]

      ## ... fitting min and max model and get residuals
      # min
      ED.Model.fit.min <- linFunct(x = D, intercept = ED.Model.interc.quant[1], gamma = ED.Model.gamma.quant[1])
      ED.Model.Res.min <- E - ED.Model.fit.min

      # max
      ED.Model.fit.max <- linFunct(x = D, intercept = ED.Model.interc.quant[2], gamma = ED.Model.gamma.quant[2])
      ED.Model.Res.max <- E - ED.Model.fit.max

      # plot(E ~ D)
      # abline(lm(E ~ D))
      # abline(a = ED.LS.boot.stat[1], b = ED.LS.boot.stat[2], col = "red")
      # abline(a = ED.LS.boot.stat[1] - ED.Model.interc.quant[1], b = ED.LS.boot.stat[2], col = "orange")

    } else {
      ED.Model <- stats::lm(E ~ D)
      ED.Model.fit <- ED.Model$fitted.values
      ED.Model.Res <- ED.Model$residuals
      # ED.LS <- stats::lsfit(y = E, x = D)

      # intercept alpha from model
      T50 <- coef(ED.Model)[[1]]

      # slope gamma from model
      gamma <- coef(ED.Model)[[2]]
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
      par.opt.NLS <- optim(x = D, y = E, par = c(mean(D), sd(D), 1),
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
        ED.NLS.boot <- boot::boot(data = data.frame(y = E, x = D), statistic = boot.NLS,
                                 R = R, formula = y ~ t + a * (x^gamma), control = list(maxiter = 500),
                                 start = list(t = par.opt.NLS[1], a = par.opt.NLS[2], gamma = par.opt.NLS[3]))


        # ... according to Rossi et al. (2017: 20)
        # compute quantiles: 5th for min, 50 as median for best fit, 95th for max
        ED.NLS.boot.conf <- apply(X = ED.NLS.boot$t, MARGIN = 2, FUN = quantile, probs = c(.05, .50, .95), na.rm=TRUE)
        T5 <- ED.NLS.boot.conf[1, 1]
        T50 <- ED.NLS.boot.conf[2, 1]
        T95 <- ED.NLS.boot.conf[3, 1]

        # Getting residuals of fit
        # ... defining linear function
        powerLawFunct <- function(x, t, a, gamma){return((t + a * (x^gamma)))}

        ## ... fitting median model and get residuals
        ED.Model.fit <- powerLawFunct(x = D, t = ED.NLS.boot.conf[2, 1], a = ED.NLS.boot.conf[2, 2], gamma = ED.NLS.boot.conf[2, 3])
        ED.Model.Res <- E - ED.Model.fit

        ## ... fitting min and max model and get residuals
        # min
        ED.Model.fit.min <- powerLawFunct(x = D, t = ED.NLS.boot.conf[1, 1], a = ED.NLS.boot.conf[1, 2], gamma = ED.NLS.boot.conf[1, 3])
        ED.Model.Res.min <- E - ED.Model.fit.min

        # max
        ED.Model.fit.max <- powerLawFunct(x = D, t = ED.NLS.boot.conf[3, 1], a = ED.NLS.boot.conf[3, 2], gamma = ED.NLS.boot.conf[3, 3])
        ED.Model.Res.max <- E - ED.Model.fit.max

      } else{

        # get model
        ED.NLS <- nls(formula = y ~ t + a * (x^gamma), data = data.frame(x = D, y = E),
                                    start = list(t = par.opt.NLS[1], a = par.opt.NLS[2], gamma = par.opt.NLS[3]),
                                    control = list(maxiter = 500))


        ED.NLS <- minpack.lm::nlsLM(formula = y ~ t + a * (x^gamma), data = data.frame(x = D, y = E),
                                    start = list(t = par.opt.NLS[1], a = par.opt.NLS[2], gamma = par.opt.NLS[3]),
                                    control = nls.lm.control(maxiter = 500))

        ED.NLS.fit <- predict(ED.NLS)
        # ED.Model.fit <- powerLawFunct(x = D, t = coef(ED.NLS)[1], a = coef(ED.NLS)[2], gamma = coef(ED.NLS)[3])

      }
  } else { # end of if method == "NLS"
    stop('Selected method is not supported. Please use "NLS" or "LS"!')
  } # end of if method

  # PLOT NLS
  # plot(y = E, x = D)
  # a<-coef(ED.NLS)[1]
  # b<-coef(ED.NLS)[2]
  # k<-coef(ED.NLS)[3]
  # x <- seq(from = 0, to = 5, by = 0.001)
  # lines(x = x, a+b*x^k,col='red')
  # ggplot2::ggplot(data.frame(D = D, ED.NLS.fit = ED.NLS.fit),
  #                 ggplot2::aes(D,ED.NLS.fit)) + ggplot2::geom_point() + ggplot2::geom_smooth()




  # IS THAT NECESAIRY?
  ## standardize residuals
  # ED.Model.SDRes <- ED.Model$residuals
  # ED.Model.SDRes <- scale(ED.Model$residuals)[, 1]


  ## check normality of residuals (limit: large sample sizes easily produce significant results from small deviations from normality)
  if(shapiro.test(ED.Model.Res)$p.value <= 0.005){warning("Residuals are highly significantly not normal distributed")}
  if(bootstrapping && shapiro.test(ED.Model.Res.min)$p.value <= 0.005){warning("Residuals are highly significantly not normal distributed")}
  if(bootstrapping &&shapiro.test(ED.Model.Res.max)$p.value <= 0.005){warning("Residuals are highly significantly not normal distributed")}



  ## Probability Densitiy Function (PDF) with Kernel Density Estimation (KDE) with a Gaussian function
  ED.PDF <- stats::density(x = ED.Model.Res, kernel = "gaussian")
  ED.PDF.df <- data.frame(x = ED.PDF$x, y = ED.PDF$y, type = "PDF")


  if(bootstrapping)
  {
    # min
    ED.PDF.min <- stats::density(x = ED.Model.Res.min, kernel = "gaussian")
    ED.PDF.df.min <- data.frame(x = ED.PDF.min$x, y = ED.PDF.min$y, type = "PDF")
    # plot(ED.PDF.min)

    # max
    ED.PDF.max <- stats::density(x = ED.Model.Res.max, kernel = "gaussian")
    ED.PDF.df.max <- data.frame(x = ED.PDF.max$x, y = ED.PDF.max$y, type = "PDF")
    # plot(ED.PDF.max)
  }



  ## Modelling PDF using a Gaussian Function
  GaussianFunction <- function(ED.PDF)
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
    par.opt.Gauss <- stats::optim(x = ED.PDF$x, y = ED.PDF$y, par = c(mean(ED.PDF$x), sd(ED.PDF$x), 1), fn = opt.Gauss,
                           method = "BFGS")


    # ... comute NLS using a Gaussian Function and fit values
    ED.PDF.NLS <- stats::nls(y ~ k*exp(-1/2*(x-mu)^2/sigma^2),
                      start = c(mu = par.opt.Gauss$par[1], sigma = par.opt.Gauss$par[2], k = par.opt.Gauss$par[3]),
                      data = data.frame(x = ED.PDF$x, y = ED.PDF$y), control = list(maxiter = 10000, reltol=1e-9))


    # ... return model and fitted values
    return(list(stats::predict(ED.PDF.NLS), ED.PDF.NLS))

  } # end of Gaussian Function


  ED.PDF.result <- GaussianFunction(ED.PDF = ED.PDF)
  ED.PDF.fit <- ED.PDF.result[[1]]
  ED.PDF.NLS <- ED.PDF.result[[2]]

  if(bootstrapping)
  {
    ED.PDF.result.min <- GaussianFunction(ED.PDF = ED.PDF.min)
    ED.PDF.fit.min <- ED.PDF.result.min[[1]]
    ED.PDF.NLS.min <- ED.PDF.result.min[[2]]

    ED.PDF.result.max <- GaussianFunction(ED.PDF = ED.PDF.max)
    ED.PDF.fit.max <- ED.PDF.result.max[[1]]
    ED.PDF.NLS.max <- ED.PDF.result.max[[2]]
  }

  # plotting density curves
  ED.PDF.fit.df <- data.frame(x = ED.PDF$x, y = ED.PDF.fit, type = "Gaussian fit")
  ED.PDF.df <- rbind(ED.PDF.df, ED.PDF.fit.df)

  ggExample <- ggplot2::ggplot(data = ED.PDF.df, ggplot2::aes(x = x, y = y, color = type, linetype = type)) +
    ggplot2::geom_line(size = 0.9) +
    ggplot2::scale_color_manual(values = c("royalblue4", "black")) +
    ggplot2::scale_linetype_manual(values = c(2, 1))

  ggExample


  ## Extracting the x procent proceeding probability
  # https://stackoverflow.com/questions/44313022/how-to-compute-confidence-interval-of-the-fitted-value-via-nls
  # https://stackoverflow.com/questions/37455512/predict-x-values-from-simple-fitting-and-annoting-it-in-the-plot?rq=1
  # https://stackoverflow.com/questions/43322568/predict-x-value-from-y-value-with-a-fitted-model

  # integrate.xy is not 1 for the whole function. Therehore, the error is estimated first and divided by 2 for both sides
  if(use.integralError)
  {
    integral.error <-  (1 - sfsmisc::integrate.xy(x = ED.PDF$x, fx = ED.PDF.fit, a = min(ED.PDF$x), b = max(ED.PDF$x)))/2
  } else{
    integral.error <- 0
  }


  integrFunctToThresh <- function(ED.PDF, ED.PDF.fit, thresh, lower, upper, par, integral.error)
  {
    # ... optimizer function for integration
    opt.Integral <- function(par, x, y, thresh)
    {
      b <- par[1]
      abs(sfsmisc::integrate.xy(x = x, fx = y, b = b) - thresh)
    }


    # ... find x parameter to area threshold of integral
    par.opt <- optim(thresh = (thresh - integral.error), x = ED.PDF$x, y = ED.PDF.fit,
                          lower = lower, upper = upper,
                          par = par, # start value
                          fn = opt.Integral, method = "L-BFGS-B",
                          control = list(pgtol = 1e-9, maxit = 10000, ndeps = 1e-9))$par

    # sfsmisc::integrate.xy(x = ED.PDF$x, fx = ED.PDF.fit, b = par.opt.Mean)
    # predict(object = ED.PDF.NLS, newdata = data.frame(x = par.opt.Mean))

    # return correspond to x which fulfills threshold
    return(par.opt)
  }

  fit.Gauss.50 <- integrFunctToThresh(ED.PDF = ED.PDF, ED.PDF.fit = ED.PDF.fit, thresh = 0.5, integral.error = integral.error,
                                      lower = (mean(ED.PDF$x) - sd(ED.PDF$x)), upper = (mean(ED.PDF$x) + sd(ED.PDF$x)), par = mean(ED.PDF$x))

  fit.Gauss.x <- integrFunctToThresh(ED.PDF = ED.PDF, ED.PDF.fit = ED.PDF.fit, thresh = prob.threshold, integral.error = integral.error,
                                     lower = (min(ED.PDF$x) + 0.0001), upper = mean(ED.PDF$x), par = mean(c(min(ED.PDF$x), mean(ED.PDF$x)))) # if b becomes the min value, then in sfsmisc::integrate.xy a = b fails!


  if(bootstrapping)
  {
    # min
    fit.Gauss.50.min <- integrFunctToThresh(ED.PDF = ED.PDF.min, ED.PDF.fit = ED.PDF.fit.min, thresh = 0.5, integral.error = integral.error,
                                        lower = (mean(ED.PDF.min$x) - sd(ED.PDF.min$x)), upper = (mean(ED.PDF.min$x) + sd(ED.PDF.min$x)), par = mean(ED.PDF.min$x))

    fit.Gauss.x.min <- integrFunctToThresh(ED.PDF = ED.PDF.min, ED.PDF.fit = ED.PDF.fit.min, thresh = prob.threshold, integral.error = integral.error,
                                       lower = (min(ED.PDF.min$x) + 0.0001), upper = mean(ED.PDF.min$x), par = mean(c(min(ED.PDF.min$x), mean(ED.PDF.min$x))))

    # max
    fit.Gauss.50.max <- integrFunctToThresh(ED.PDF = ED.PDF.max, ED.PDF.fit = ED.PDF.fit.max, thresh = 0.5, integral.error = integral.error,
                                            lower = (mean(ED.PDF.max$x) - sd(ED.PDF.max$x)), upper = (mean(ED.PDF.max$x) + sd(ED.PDF.max$x)), par = mean(ED.PDF.max$x))

    fit.Gauss.x.max <- integrFunctToThresh(ED.PDF = ED.PDF.max, ED.PDF.fit = ED.PDF.fit.max, thresh = prob.threshold, integral.error = integral.error,
                                         lower = (min(ED.PDF.max$x) + 0.0001), upper = mean(ED.PDF.max$x), par = mean(c(min(ED.PDF.max$x), mean(ED.PDF.max$x))))
  }



  par.opt.Mean <- optim(thresh = (0.5 - integral.error), x = ED.PDF$x, y = ED.PDF.fit,
                        lower = (mean(ED.PDF$x) - sd(ED.PDF$x)),
                        upper = (mean(ED.PDF$x) + sd(ED.PDF$x)),
                        par = mean(ED.PDF$x), # start value
                        fn = opt.Integral, method = "L-BFGS-B",
                        control = list(pgtol = 1e-9, maxit = 10000, ndeps = 1e-9))$par



  # ... ... (2) integration for probablilty level threshold
  fit.Gauss.x  <- optim(thresh = (prob.threshold - integral.error), x = ED.PDF$x, y = ED.PDF.fit,
                            lower = (min(ED.PDF$x) + 0.0000001), #  if b becomes the min value, then in sfsmisc::integrate.xy a = b fails!
                            upper = mean(ED.PDF$x),
                            par = mean(min(ED.PDF$x), (mean(ED.PDF$x) - sd(ED.PDF$x))), #  start value
                            fn = opt.Integral, method = "L-BFGS-B",
                            control = list(pgtol = 1e-9, maxit = 10000, ndeps = 1e-9))$par


  # sfsmisc::integrate.xy(x = ED.PDF$x, fx = ED.PDF.fit, a = min(ED.PDF$x), b = par.opt.probLevel)
  ggExample + ggplot2::geom_segment(inherit.aes = FALSE, color = "red", size = 0.8,
                                    ggplot2::aes(x = fit.Gauss.x , xend = fit.Gauss.x,
                                                 y = 0, yend = predict(object = ED.PDF.NLS, newdata = data.frame(x = fit.Gauss.x))))

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
  df.fit.LS <- data.frame(x = D, y = E, y_fit = ED.Model.fit,
                   y_fit_min = ED.Model.fit.min, y_fit_max = ED.Model.fit.max,
                   y_fit_Lx = linFunct(x = D, intercept = Tx, gamma = gamma),
                   y_fit_Lx_min = linFunct(x = D, intercept = Tx.min, gamma = ED.Model.gamma.quant[1]),
                   y_fit_Lx_max = linFunct(x = D, intercept = Tx.max, gamma = ED.Model.gamma.quant[2]))


  ggplot2::ggplot(data = df.fit.LS) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y)) + # add E-D point pairs
    # ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y_fit)) + # add fit of regression line (median fit)
    ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y_fit_Lx), col = "green") + # add fit of regression line (median fit)
    ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y_fit_Lx_min), col = "blue") + # add fit of regression line (median fit)
    ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y_fit_Lx_max), col = "red") + # add fit of regression line (median fit)
    ggplot2::xlim(0, 6) + ggplot2::ylim(-1, 6)


  # ... plot best fit with 5-95th quantile
  ggplot2::ggplot(data = df.fit.LS) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y)) + # add E-D point pairs
    ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y_fit_min), col = "blue") + # add fit of regression line (median fit)
    ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y_fit), col = "black") +
    ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y_fit_max), col = "red") +
    ggplot2::xlim(0, 6) + ggplot2::ylim(-1, 6)


  ### NLS
  ## create data frame for ggplot2
  df.fit.NLS <- data.frame(x = D, y = E, y_fit = ED.Model.fit,
                          y_fit_min = ED.Model.fit.min, y_fit_max = ED.Model.fit.max,
                          y_fit_Lx = powerLawFunct(x = D, t = Tx, a = ED.NLS.boot.conf[2, 2], gamma = ED.NLS.boot.conf[2, 3]),
                          y_fit_Lx_min = powerLawFunct(x = D, t = Tx.min, a = ED.NLS.boot.conf[1, 2], gamma = ED.NLS.boot.conf[1, 3]),
                          y_fit_Lx_max = powerLawFunct(x = D, t = Tx.max, a = ED.NLS.boot.conf[3, 2], gamma = ED.NLS.boot.conf[3, 3]))


  ggplot2::ggplot(data = df.fit.NLS) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y)) + # add E-D point pairs
    ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_min), col = "orange", method = "loess", size = 0.7) +
    ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_max), col = "yellow", method = "loess", size = 0.7) +
    ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit), col = "black", method = "loess", size = 0.7)+ # add fit of regression line (median fit)
    ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx), col = "green", method = "loess", size = 0.7) + # add fit of regression line (median fit)
    ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx_min), col = "blue", method = "loess", size = 0.7) + # add fit of regression line (median fit)
    ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx_max), col = "red", method = "loess", size = 0.7) + # add fit of regression line (median fit)
    ggplot2::xlim(0, 6) + ggplot2::ylim(-1, 6)


  ggplot2::ggplot(data = df.fit.NLS) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y)) + # add E-D point pairs
    ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit), col = "black", method = "loess", size = 0.7) +
    ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_min), col = "orange", method = "loess", size = 0.7) +
    ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_max), col = "yellow", method = "loess", size = 0.7) +
    ggplot2::xlim(0, 6) + ggplot2::ylim(-1, 6)


  ggplot2::ggplot(data = df.fit.NLS) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y)) + # add E-D point pairs
    ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx), col = "green", method = "loess", size = 0.7) + # add fit of regression line (median fit)
    ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx_min), col = "blue", method = "loess", size = 0.7) + # add fit of regression line (median fit)
    ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx_max), col = "red", method = "loess", size = 0.7) + # add fit of regression line (median fit)
    ggplot2::xlim(0, 6) + ggplot2::ylim(-1, 6)


  ## combine both LS and NSL
  ggplot2::ggplot(data = df.fit.NLS) +
    # NLS
    ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y)) + # add E-D point pairs
    ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx), col = "yellow", method = "loess", size = 0.8) + # add fit of regression line (median fit)
    ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx_min), col = "blue", method = "loess", size = 0.7) + # add fit of regression line (median fit)
    ggplot2::geom_smooth(mapping = ggplot2::aes(x = x, y = y_fit_Lx_max), col = "blue", method = "loess", size = 0.7) + # add fit of regression line (median fit)
    # LS
    ggplot2::geom_line(data = df.fit.LS, mapping = ggplot2::aes(x = x, y = y_fit_Lx), col = "orange") + # add fit of regression line (median fit)
    ggplot2::geom_line(data = df.fit.LS, mapping = ggplot2::aes(x = x, y = y_fit_Lx_min), col = "red") + # add fit of regression line (median fit)
    ggplot2::geom_line(data = df.fit.LS, mapping = ggplot2::aes(x = x, y = y_fit_Lx_max), col = "red") + # add fit of regression line (median fit)
    ggplot2::xlim(0, 6) + ggplot2::ylim(-1, 6)



  return(c(alpha, gamma))

} # end of function getRainThreshFreqM
