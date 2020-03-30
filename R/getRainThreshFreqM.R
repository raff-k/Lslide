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
#' @param use.boot_median Threshold is delineated from median of bootstrapping result. Default: FALSE
#' @param seed replicable bootstrapping. Default: 123
#' @param use.integralError for estimating x of prob.threshold, the entire function is integrated first to estimate the bias: (1 - INTEGRAL)/2. Default. TRUE
#' @param cores If cores > 1 than parallisation for bootstrapping is initialized via future backend. Default: 1
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
getRainThreshFreqM <- function(Re, D, method = c("LM", "NLS"), prob.threshold = 0.05, log10.transform = FALSE,
                               bootstrapping = TRUE, R = 1000, use.boot_median = FALSE,
                               nls.bounds = list(lower = c(0, 0, 0), upper = c(500, 100, 10)), seed = 123, cores = 1){

  # get start time of process
  process.time.start <- proc.time()


  # approximate zScores
  sigma.coef <- qnorm((1 + (1-prob.threshold*2))/2)

  # init some empty variables
  boot.result <- list()

  ## log transformation of input vectors
  if(log10.transform)
  {
    Re <- log(x = Re, base = 10) # cumultive precipitation of rain event
    D <- log(x = D, base = 10) # duration of rain event
  } else {
    if(method != "NLS"){ warning("Precipitation and duration must be in log-scale!\n")}
  }


  ## Least Squares Method (LS) ---------------------------
  if(method == "LS"){


    # ... defining linear function in logit scale!
    funct <- function(x, a, gamma){return((a + gamma * x))}

    m.ReD <- stats::lm(Re ~ D)
    m.ReD.fit <- m.ReD$fitted.values
    m.ReD.Res <- m.ReD$residuals

    # intercept alpha from model
    alpha <- coef(m.ReD)[[1]]

    # slope gamma from model
    gamma <- coef(m.ReD)[[2]]

    # boot-variables
    boot.formula <- Re ~ D
    boot.control <- NULL
    boot.start <- NULL
    boot.bounds <- NULL
    t <- NULL
  }


  ## Nonlinear Least Squares Method (NLS) ----------------
  if(method == "NLS") {

   # ... defining non-linear function
    funct <- function(x, t, a, gamma){return((t + a * (x^gamma)))}

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
    par.opt.NLS <- optim(x = D, y = Re, par = c(0, 1, 1),
                         fn = opt.NLS,
                         method = "L-BFGS-B",
                         control = list(pgtol = 1e-8,
                                        maxit = 100),
                         lower = nls.bounds$lower,
                         upper = nls.bounds$upper)$par

    # get model
    m.ReD <- stats::nls(formula = Re ~ t + a * (D^gamma),
                        # data = data.frame(x = D, y = Re),
                        start = list(t = par.opt.NLS[1],
                                     a = par.opt.NLS[2],
                                     gamma = par.opt.NLS[3]),
                        control = list(maxiter = 500),
                        algorithm = "port",
                        lower = nls.bounds$lower,
                        upper = nls.bounds$upper)


    m.ReD.fit <- predict(m.ReD)
    m.ReD.Res <- residuals(m.ReD)

    # intercept i
    t <- coef(m.ReD)[[1]]

    # alpha from model
    alpha <- coef(m.ReD)[[2]]

    # slope gamma from model
    gamma <- coef(m.ReD)[[3]]


    # boot-variables
    boot.formula <- Re ~ t + a * (D^gamma)
    boot.control <- list(maxiter = 500)
    boot.start <- list(t = par.opt.NLS[1],
                       a = par.opt.NLS[2],
                       gamma = par.opt.NLS[3])
    boot.bounds <- nls.bounds
  }


  # bootstrapping --------------
  if(bootstrapping){

    if(cores > 1)
    {
      cores <- ifelse(cores > parallel::detectCores()[[1]], parallel::detectCores(), cores)
      cl <- parallel::makeCluster(cores)
    } else {
      cl <- NULL
    }

      # set seed for replication
      set.seed(seed)


    # ... create bootstrap function
      boot.fun <- function(formula, data, method, start, control, bounds, indices){

        # subset data
        data <- data[indices,]

            # boot model
            if(method == "LS")
            {
              tryCatch({
                        m.boot <- stats::lm(formula, data)
                        return(coef(m.boot)) # alpha, gamma
                       }, error = function(e) {
                          return(c(NA, NA))
                       })
            }


            if(method == "NLS")
            {
              tryCatch({
                          m.boot <-  stats::nls(formula = formula, data = data, start = start,
                                                 control = control, algorithm = "port",
                                                 lower = bounds$lower, upper = bounds$upper)
                          return(coef(m.boot)[c(2,3,1)]) # alpha, gamma, t

              }, error = function(e) {
                return(c(NA, NA, NA))
              })
            }

      } # end of boot.fun


      # bootstrapping with R replications
      boot.ReD <- boot::boot(data = data.frame(Re = Re, D = D),
                             statistic = boot.fun,
                             R = R,
                             method = method,
                             formula = boot.formula,
                             control = boot.control,
                             start = boot.start,
                             bounds = boot.bounds,
                             cl = cl, ncpus = cores)



      # ... according to Rossi et al. (2017: 20)
      # compute quantiles: 5th for min, 50 as median for best fit, 95th for max
      boot.ReD.conf <- apply(X = boot.ReD$t, MARGIN = 2, FUN = function(x) c(quantile(x = x, probs = c(.05, .50, .95), na.rm=TRUE), mean(x, na.rm = TRUE)))
      boot.ReD.conf <- unname(boot.ReD.conf)

      # ... intercept alpha from model
      boot.result$alpha <- list("a5" = boot.ReD.conf[1, 1],
                                "a50" = boot.ReD.conf[2, 1],
                                "a95" = boot.ReD.conf[3, 1],
                                "mean" = boot.ReD.conf[4, 1])


      # ... slope gamma from model
      boot.result$gamma <- list("g5" = boot.ReD.conf[1, 2],
                             "g50" = boot.ReD.conf[2, 2],
                             "g95"= boot.ReD.conf[3, 2],
                             "mean" = boot.ReD.conf[4, 2])
      if(method == "NLS")
      {
        boot.result$t <- list("t5" = boot.ReD.conf[1, 3],
                              "t50" = boot.ReD.conf[2, 3],
                              "t95"= boot.ReD.conf[3, 3],
                              "mean" = boot.ReD.conf[4, 3])
      }


      if(use.boot_median){

        # set and overwrite alpha and gamma
        alpha <- boot.T$a50
        gamma <- boot.result$gamma50

        # ... fitting median model and get residuals
        m.ReD.fit <- powerLawFunct(x = D, t = ReD.NLS.boot.conf[2, 1], a = ReD.NLS.boot.conf[2, 2], gamma = ReD.NLS.boot.conf[2, 3])
        m.ReD.Res <- Re - m.ReD.fit
      }

      if(cores > 1)
      {
        parallel::stopCluster(cl = cl)
      }

    } # end of bootstrapping



  ## FROM: "CTRL-T_code.R"
  # Kernel density calculation -------------
  # ReD.PDF <- stats::density(x = m.ReD.Res, bw = "nrd0", kernel = "gaussian", adjust = 1, from = -1, to = 1)


  # Maximum likelihood fit -------------
  m.PDF <- MASS::fitdistr(x = m.ReD.Res, densfun = "normal")
  m.PDF.sigma <- coef(m.PDF)[2] # standard deviation from nomal approximation
  m.PDF.sigma <- unname(m.PDF.sigma)


  ## get Intercept of probablilty level ----------------

  if(method == "LS")
  {
    alpha <- alpha - sigma.coef*m.PDF.sigma # ... sigma.coef ist factor for exceedance threshold
    threshold.fit <- funct(x = D, a = alpha, gamma = gamma)

    if(bootstrapping)
    {
      # ... confidence interval: shift of intercept
      boot.result$alpha$a5_thresh <- boot.result$alpha$a5 - sigma.coef*m.PDF.sigma # ... lower percentile
      boot.result$alpha$a95_thresh <- boot.result$alpha$a95 - sigma.coef*m.PDF.sigma # ... upper percentile
    }

    if(log10.transform)
    {
      alpha <- 10^alpha

      if(bootstrapping)
      {
        boot.result$alpha <- lapply(X = boot.result$alpha, function(x) 10^x)
      }
    }

  } # end of LS

  if(method == "NLS")
  {
    t <- t - sigma.coef*m.PDF.sigma # ... sigma.coef ist factor for exceedance threshold
    threshold.fit <- funct(x = D, t = t, a = alpha, gamma = gamma)

    if(bootstrapping)
    {
      # ... confidence interval: shift of intercept
      boot.result$t$t5_thresh <- boot.result$t$a5 - sigma.coef*m.PDF.sigma # ... lower percentile
      boot.result$t$t95_thresh <- boot.result$t$t95 - sigma.coef*m.PDF.sigma # ... upper percentile
    }

    if(log10.transform)
    {
      t <- 10^t

      if(bootstrapping)
      {
        boot.result$t <- lapply(X = boot.result$t, function(x) 10^x)
      }
    }


  } # end of NLS

  # get time of process
  process.time.run <- proc.time() - process.time.start
  cat("------ Run of getRainThreshFreqM: " , round(x = process.time.run["elapsed"][[1]]/60, digits = 4), " Minutes \n")


  return(list(threshold = list(fit = threshold.fit,
                               alpha = alpha,
                               gamma = gamma,
                               t = t,
                               model = m.ReD,
                               fun = funct,
                               dens = m.PDF,
                               sigma = m.PDF.sigma),
              boot = boot.result))

} # end of function getRainThreshFreqM




# CODE SNIPPETS -------------------

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
# m.ReD.SDRes <- m.ReD$residuals
# m.ReD.SDRes <- scale(m.ReD$residuals)[, 1]


## check normality of residuals (limit: large sample sizes easily produce significant results from small deviations from normality)
# if(shapiro.test(m.ReD.Res)$p.value <= 0.005){warning("Residuals are highly significantly not normal distributed")}
# if(bootstrapping && shapiro.test(m.ReD.Res.min)$p.value <= 0.005){warning("Residuals are highly significantly not normal distributed")}
# if(bootstrapping &&shapiro.test(m.ReD.Res.max)$p.value <= 0.005){warning("Residuals are highly significantly not normal distributed")}



## Probability Densitiy Function (PDF) with Kernel Density Estimation (KDE) with a Gaussian function
# ReD.PDF <- stats::density(x = m.ReD.Res, kernel = "gaussian")
# ReD.PDF.df <- data.frame(x = ReD.PDF$x, y = ReD.PDF$y, type = "PDF")



## Modelling PDF using a Gaussian Function
# GaussianFunction <- function(ReD.PDF)
# {
#   # ... using Nonlinear regression model (NLS)
#   # # https://stats.stackexchange.com/questions/220109/fit-a-gaussian-to-data-with-r-with-optim-and-nls
#
#   # ... optimizer function
#   opt.Gauss <- function(par, x, y)
#   {
#     m <- par[1]
#     sd <- par[2]
#     k <- par[3]
#     rhat <- k * exp(-0.5 * ((x - m)/sd)^2)
#     sum((y - rhat)^2)
#   }
#
#   # ... optimize parameter
#   par.opt.Gauss <- stats::optim(x = ReD.PDF$x, y = ReD.PDF$y, par = c(mean(ReD.PDF$x), sd(ReD.PDF$x), 1), fn = opt.Gauss,
#                          method = "BFGS")
#
#
#   # ... comute NLS using a Gaussian Function and fit values
#   ReD.PDF.NLS <- stats::nls(y ~ k*exp(-1/2*(x-mu)^2/sigma^2),
#                     start = c(mu = par.opt.Gauss$par[1], sigma = par.opt.Gauss$par[2], k = par.opt.Gauss$par[3]),
#                     data = data.frame(x = ReD.PDF$x, y = ReD.PDF$y), control = list(maxiter = 10000, reltol=1e-9))
#
#
#   # ... return model and fitted values
#   return(list(stats::predict(ReD.PDF.NLS), ReD.PDF.NLS))
#
# } # end of Gaussian Function


# ReD.PDF.result <- GaussianFunction(ReD.PDF = ReD.PDF)
# ReD.PDF.fit <- ReD.PDF.result[[1]]
# ReD.PDF.NLS <- ReD.PDF.result[[2]]

# if(bootstrapping)
# {
#   ReD.PDF.result.min <- GaussianFunction(ReD.PDF = ReD.PDF.min)
#   ReD.PDF.fit.min <- ReD.PDF.result.min[[1]]
#   ReD.PDF.NLS.min <- ReD.PDF.result.min[[2]]
#
#   ReD.PDF.result.max <- GaussianFunction(ReD.PDF = ReD.PDF.max)
#   ReD.PDF.fit.max <- ReD.PDF.result.max[[1]]
#   ReD.PDF.NLS.max <- ReD.PDF.result.max[[2]]
# }

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
# if(use.integralError)
# {
#   integral.error <-  (1 - sfsmisc::integrate.xy(x = ReD.PDF$x, fx = ReD.PDF.fit, a = min(ReD.PDF$x), b = max(ReD.PDF$x)))/2
# } else{
#   integral.error <- 0
# }
#
#
# integrFunctToThresh <- function(ReD.PDF, ReD.PDF.fit, thresh, lower, upper, par, integral.error)
# {
#   # ... optimizer function for integration
#   opt.Integral <- function(par, x, y, thresh)
#   {
#     b <- par[1]
#     abs(sfsmisc::integrate.xy(x = x, fx = y, b = b) - thresh)
#   }
#
#
#   # ... find x parameter to area threshold of integral
#   par.opt <- optim(thresh = (thresh - integral.error), x = ReD.PDF$x, y = ReD.PDF.fit,
#                         lower = lower, upper = upper,
#                         par = par, # start value
#                         fn = opt.Integral, method = "L-BFGS-B",
#                         control = list(pgtol = 1e-9, maxit = 10000, ndeps = 1e-9))$par
#
#   # sfsmisc::integrate.xy(x = ReD.PDF$x, fx = ReD.PDF.fit, b = par.opt.Mean)
#   # predict(object = ReD.PDF.NLS, newdata = data.frame(x = par.opt.Mean))
#
#   # return correspond to x which fulfills threshold
#   return(par.opt)
# }

# fit.Gauss.50 <- integrFunctToThresh(ReD.PDF = ReD.PDF, ReD.PDF.fit = ReD.PDF.fit, thresh = 0.5, integral.error = integral.error,
#                                     lower = (mean(ReD.PDF$x) - sd(ReD.PDF$x)), upper = (mean(ReD.PDF$x) + sd(ReD.PDF$x)), par = mean(ReD.PDF$x))
#
# fit.Gauss.x <- integrFunctToThresh(ReD.PDF = ReD.PDF, ReD.PDF.fit = ReD.PDF.fit, thresh = prob.threshold, integral.error = integral.error,
#                                    lower = (min(ReD.PDF$x) + 0.0001), upper = mean(ReD.PDF$x), par = mean(c(min(ReD.PDF$x), mean(ReD.PDF$x)))) # if b becomes the min value, then in sfsmisc::integrate.xy a = b fails!



#if(bootstrapping)
# {
# min
# fit.Gauss.50.min <- integrFunctToThresh(ReD.PDF = ReD.PDF.min, ReD.PDF.fit = ReD.PDF.fit.min, thresh = 0.5, integral.error = integral.error,
#                                     lower = (mean(ReD.PDF.min$x) - sd(ReD.PDF.min$x)), upper = (mean(ReD.PDF.min$x) + sd(ReD.PDF.min$x)), par = mean(ReD.PDF.min$x))
#
# fit.Gauss.x.min <- integrFunctToThresh(ReD.PDF = ReD.PDF.min, ReD.PDF.fit = ReD.PDF.fit.min, thresh = prob.threshold, integral.error = integral.error,
#                                    lower = (min(ReD.PDF.min$x) + 0.0001), upper = mean(ReD.PDF.min$x), par = mean(c(min(ReD.PDF.min$x), mean(ReD.PDF.min$x))))
#
# # max
# fit.Gauss.50.max <- integrFunctToThresh(ReD.PDF = ReD.PDF.max, ReD.PDF.fit = ReD.PDF.fit.max, thresh = 0.5, integral.error = integral.error,
#                                         lower = (mean(ReD.PDF.max$x) - sd(ReD.PDF.max$x)), upper = (mean(ReD.PDF.max$x) + sd(ReD.PDF.max$x)), par = mean(ReD.PDF.max$x))
#
# fit.Gauss.x.max <- integrFunctToThresh(ReD.PDF = ReD.PDF.max, ReD.PDF.fit = ReD.PDF.fit.max, thresh = prob.threshold, integral.error = integral.error,
#                                      lower = (min(ReD.PDF.max$x) + 0.0001), upper = mean(ReD.PDF.max$x), par = mean(c(min(ReD.PDF.max$x), mean(ReD.PDF.max$x))))
# }



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


### LS
## create data frame for ggplot2
# df.fit.LS <- data.frame(x = D, y = Re, y_fit = m.ReD.fit,
#                  y_fit_min = m.ReD.fit.min, y_fit_max = m.ReD.fit.max,
#                  y_fit_Lx = linFunct(x = D, intercept = Tx, gamma = gamma),
#                  y_fit_Lx_min = linFunct(x = D, intercept = Tx.min, gamma = m.ReD.gamma.quant[1]),
#                  y_fit_Lx_max = linFunct(x = D, intercept = Tx.max, gamma = m.ReD.gamma.quant[2]))
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
# df.fit.NLS <- data.frame(x = D, y = Re, y_fit = m.ReD.fit,
#                         y_fit_min = m.ReD.fit.min, y_fit_max = m.ReD.fit.max,
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



# ReD.PDF.gaussian.mean <-coef(ReD.PDF.gaussian)[1]
# ReD.PDF.gaussian.curve <- dnorm(ReD.PDF$x,
#                                 mean = ReD.PDF.gaussian.mean,
#                                 sd = ReD.PDF.gaussian.sigma)

# ggplot() +
#   geom_line(aes(x = ReD.PDF$x, y = ReD.PDF.gaussian.curve)) +
#   geom_line(aes(x = ReD.PDF$x, y = ReD.PDF$y), color = "blue")
