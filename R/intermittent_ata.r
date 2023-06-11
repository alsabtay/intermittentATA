#' @docType package
#' @keywords package
"_PACKAGE"

#' @importFrom ATAforecasting ATA.SeasAttr ATA.BoxCoxAttr ATA
#' @importFrom tsibble measured_vars as_tibble tibble index
#' @importFrom rlang expr_text
#' @importFrom stats frequency ts start
#' @importFrom fabletools get_frequencies
#' @importFrom tsbox ts_ts ts_tsibble
#' @importFrom dplyr ungroup
#' @importFrom lubridate decimal_date
#'
train_intermittentATA <- function(.data, specials, ...){
  if(length(tsibble::measured_vars(.data)) > 1){
    stop("Only univariate responses are supported by Intermittent ATA method")
  }

  # Allow only one special of each type.
  # TODO: Add message/warning if more than one of these specials is in formula

  # Get specials
  d_level <- specials$d_level[[1]]
  d_trend <- specials$d_trend[[1]]
  d_accuracy <- specials$d_accuracy[[1]]
  d_transform <- specials$d_transform[[1]]
  d_holdout <- specials$d_holdout[[1]]

  i_level <- specials$i_level[[1]]
  i_trend <- specials$i_trend[[1]]
  i_accuracy <- specials$i_accuracy[[1]]
  i_transform <- specials$i_transform[[1]]
  i_holdout <- specials$i_holdout[[1]]

  intermittent <- specials$intermittent[[1]]

  # Prepare data for modelling
  model_data <- tsibble::as_tibble(.data)[c(rlang::expr_text(tsibble::index(.data)), tsibble::measured_vars(.data))]
  colnames(model_data) <- c("idy", "y")
  # Check data
  if (any(is.na(model_data$y))) {
    stop("ATA method does not support missing values.")
  }
  if (any(model_data$y < 0)) {
    stop("All observations must be non-negative for Intermittent ATA method.")
  }
  non_zero <- which(model_data$y != 0)

  if (length(non_zero) < 2) {
    stop("At least two non-zero values are required to use Intermittent ATA method.")
  }

  # Get response
  # Croston demand/interval decomposition
  y_demand <- model_data$y[non_zero]
  idy_demand <- model_data$idy[non_zero]
  y_interval <- c(non_zero[1], diff(non_zero))
  train_demand <- stats::ts(y_demand, start = 1, frequency = 1)
  train_interval <- stats::ts(y_interval, start = 1, frequency = 1)
  if (d_holdout$holdout == TRUE & d_accuracy$criteria == "AMSE") {
    stop("ATA Method does not support 'AMSE' for 'holdout' forecasting in the demand component.")
  }
  if (i_holdout$holdout == TRUE & i_accuracy$criteria == "AMSE") {
    stop("ATA Method does not support 'AMSE' for 'holdout' forecasting in the interval component.")
  }

  d_seas_opt_crit <- i_seas_opt_crit <- ATAforecasting::ATA.SeasAttr()

  d_bc_opt_crit <- ATAforecasting::ATA.BoxCoxAttr(bcMethod = d_transform$bcMethod, bcLower = d_transform$bcLower, bcUpper = d_transform$bcUpper, bcBiasAdj = FALSE)
  i_bc_opt_crit <- ATAforecasting::ATA.BoxCoxAttr(bcMethod = i_transform$bcMethod, bcLower = i_transform$bcLower, bcUpper = i_transform$bcUpper, bcBiasAdj = FALSE)


  # Build and train model for Demand
  d_pmdl_ATA <- safely(quietly(ATAforecasting::ATA))(
                  X = train_demand,
                  Y = NULL,
                  parP = d_level$parP,
                  parQ = d_trend$parQ,
                  parPHI = d_trend$parPHI,
                  model.type = ifelse(d_trend$type=="N","A",ifelse(d_trend$type=="Ad", "A", ifelse(d_trend$type=="Md","M",d_trend$type))),
                  seasonal.test = FALSE,
                  seasonal.model = "none",
                  seasonal.period = 1,
                  seasonal.type = "M",
                  seasonal.test.attr = d_seas_opt_crit,
                  find.period = NULL,
                  accuracy.type = d_accuracy$criteria,
                  nmse = d_accuracy$nmse,
                  level.fixed = d_level$level_fixed,
                  trend.opt = d_trend$trend_opt,
                  h = 1,
                  train_test_split = NULL,
                  holdout = d_holdout$holdout,
                  holdout.adjustedP = d_holdout$adjustment,
                  holdout.set_size = d_holdout$set_size,
                  holdout.onestep = d_holdout$onestep,
                  holdin = FALSE,
                  transform.order = ifelse(d_transform$order == "none", "before", d_transform$order),
                  transform.method = switch((d_transform$method != "none") + 1, NULL, d_transform$method),
                  transform.attr = d_bc_opt_crit,
                  lambda = d_transform$lambda,
                  shift = d_transform$shift,
                  initial.level = d_level$initial_level,
                  initial.trend = d_trend$initial_trend,
                  ci.level = 95,
                  start.phi = d_trend$parPHI_range[1],
                  end.phi = d_trend$parPHI_range[2],
                  size.phi = d_trend$parPHI_increment,
                  negative.forecast = TRUE,
                  onestep = FALSE,
                  print.out = FALSE,
                  plot.out = FALSE)
  d_mdl_ATA <- d_pmdl_ATA[["result"]]
  if(d_mdl_ATA$q==0){
    d_trend_mthd <- "N"
  }else if (d_mdl_ATA$q!=0 & d_mdl_ATA$phi!=1){
    d_trend_mthd <- paste(d_mdl_ATA$model.type, "d", sep="")
  }else{
    d_trend_mthd <- d_mdl_ATA$model.type
  }
  if(d_mdl_ATA$seasonal.model == "none"){
    d_seas_mthd <- "N"
  }else{
    d_seas_mthd <- d_mdl_ATA$seasonal.type
  }

  # Build and train model for Demand
  i_pmdl_ATA <- safely(quietly(ATAforecasting::ATA))(
                  X = train_interval,
                  Y = NULL,
                  parP = i_level$parP,
                  parQ = i_trend$parQ,
                  parPHI = i_trend$parPHI,
                  model.type = ifelse(i_trend$type=="N","A",ifelse(i_trend$type=="Ad", "A", ifelse(i_trend$type=="Md","M",i_trend$type))),
                  seasonal.test = FALSE,
                  seasonal.model = "none",
                  seasonal.period = 1,
                  seasonal.type = "M",
                  seasonal.test.attr = i_seas_opt_crit,
                  find.period = NULL,
                  accuracy.type = i_accuracy$criteria,
                  nmse = i_accuracy$nmse,
                  level.fixed = i_level$level_fixed,
                  trend.opt = i_trend$trend_opt,
                  h = 1,
                  train_test_split = NULL,
                  holdout = i_holdout$holdout,
                  holdout.adjustedP = i_holdout$adjustment,
                  holdout.set_size = i_holdout$set_size,
                  holdout.onestep = i_holdout$onestep,
                  holdin = FALSE,
                  transform.order = ifelse(i_transform$order == "none", "before", i_transform$order),
                  transform.method = switch((i_transform$method != "none") + 1, NULL, i_transform$method),
                  transform.attr = i_bc_opt_crit,
                  lambda = i_transform$lambda,
                  shift = i_transform$shift,
                  initial.level = i_level$initial_level,
                  initial.trend = i_trend$initial_trend,
                  ci.level = 95,
                  start.phi = i_trend$parPHI_range[1],
                  end.phi = i_trend$parPHI_range[2],
                  size.phi = i_trend$parPHI_increment,
                  negative.forecast = TRUE,
                  onestep = FALSE,
                  print.out = FALSE,
                  plot.out = FALSE)
  i_mdl_ATA <- i_pmdl_ATA[["result"]]
  if(i_mdl_ATA$q==0){
    i_trend_mthd <- "N"
  }else if (i_mdl_ATA$q!=0 & i_mdl_ATA$phi!=1){
    i_trend_mthd <- paste(i_mdl_ATA$model.type, "d", sep="")
  }else{
    i_trend_mthd <- i_mdl_ATA$model.type
  }
  if(i_mdl_ATA$seasonal.model == "none"){
    i_seas_mthd <- "N"
  }else{
    i_seas_mthd <- i_mdl_ATA$seasonal.type
  }

interval_param <- i_mdl_ATA$p / length(i_mdl_ATA$actual)

if (intermittent == "sba") {
    coeff <- 1 - interval_param / 2
  } else if (intermittent == "sbj") {
    coeff <- 1 - interval_param / (2 - interval_param)
  } else {
    coeff <- 1
  }
  ratio <- coeff * d_mdl_ATA$fitted / i_mdl_ATA$fitted

  # In-sample demand rate
  demand_rate <- rep(c(NA_real_, ratio), diff(c(0, non_zero, length(model_data$y))))

  demand_data <- tibble(mth = yearmonth(idy_demand), value = d_mdl_ATA$actual)
  demand_data <- as_tsibble(demand_data, index = mth)

  interval_data <- tibble(mth = yearmonth(idy_demand), value = i_mdl_ATA$actual)
  interval_data <- as_tsibble(interval_data, index = mth)

  # Return model
   structure(
     list(
         par = list("d_par" = tsibble::tibble(term = names(unlist(d_mdl_ATA$par.specs)), estimate = unlist(d_mdl_ATA$par.specs)),
                    "i_par" = tsibble::tibble(term = names(unlist(i_mdl_ATA$par.specs)), estimate = unlist(i_mdl_ATA$par.specs))
                    ),
          est = dplyr::mutate(dplyr::ungroup(.data),
                              .fitted = demand_rate,
                              .resid = model_data$y - demand_rate
                            ),
          est_comp = list(d_est = dplyr::mutate(dplyr::ungroup(demand_data),
                                                .fitted = d_mdl_ATA$fitted,
                                                .resid = d_mdl_ATA$residuals
                                                ),
                          i_est = dplyr::mutate(dplyr::ungroup(interval_data),
                                                .fitted = i_mdl_ATA$fitted,
                                                .resid = i_mdl_ATA$residuals
                                                )
                          ),
          fit = list("d_fit" = d_mdl_ATA$accuracy$fits,
                    "i_fit" = i_mdl_ATA$accuracy$fits
                    ),
          components = list(d_components = list("coefp" = d_mdl_ATA$coefp,
                                                "coefq" = d_mdl_ATA$coefp,
                                                "P" = d_mdl_ATA$p,
                                                "Q" = d_mdl_ATA$q,
                                                "PHI" = d_mdl_ATA$phi,
                                                "response" = d_mdl_ATA$actual,
                                                "level" = d_mdl_ATA$level,
                                                "trend" = d_mdl_ATA$trend,
                                                "season" = d_mdl_ATA$seasonal,
                                                "seasindex" = d_mdl_ATA$seasonal.index,
                                                "seasadj" = d_mdl_ATA$seasonal.adjusted,
                                                "remainder" = d_mdl_ATA$residuals
                                              ),
                           i_components = list("coefp" = i_mdl_ATA$coefp,
                                                "coefq" = i_mdl_ATA$coefp,
                                                "P" = i_mdl_ATA$p,
                                                "Q" = i_mdl_ATA$q,
                                                "PHI" = i_mdl_ATA$phi,
                                                "response" = i_mdl_ATA$actual,
                                                "level" = i_mdl_ATA$level,
                                                "trend" = i_mdl_ATA$trend,
                                                "season" = i_mdl_ATA$seasonal,
                                                "seasindex" = i_mdl_ATA$seasonal.index,
                                                "seasadj" = i_mdl_ATA$seasonal.adjusted,
                                                "remainder" = i_mdl_ATA$residuals
                                               )
                              ),
          transform = list(d_transform = list("method" = d_transform$method,
                                              "order" = d_transform$order,
                                               "lambda" = d_mdl_ATA$lambda,
                                                "shift" = d_mdl_ATA$shift
                                              ),
                            i_transform = list("method" = i_transform$method,
                                               "order" = i_transform$order,
                                               "lambda" = i_mdl_ATA$lambda,
                                               "shift" = i_mdl_ATA$shift
                                              )
                          ),
          holdout = list(d_holdout = list("holdout" = d_mdl_ATA$holdout,
                                          "adjustment" = d_holdout$adjustment,
                                           "onestep" = d_holdout$onestep
                                          ),
                          i_holdout = list("holdout" = i_mdl_ATA$holdout,
                                            "adjustment" = i_holdout$adjustment,
                                            "onestep" = i_holdout$onestep
                                          )
                              ),
          spec = list(d_spec = list("errortype" = "A",
                                    "trendtype" = d_trend_mthd,
                                    "seasontype" = d_seas_mthd,
                                    "damped" = ifelse(d_mdl_ATA$phi==1, FALSE, TRUE),
                                    "period" = d_mdl_ATA$seasonal.period,
                                    "method" = d_mdl_ATA$method,
                                    "onestep" = d_mdl_ATA$onestep
                      ),
                      i_spec = list("errortype" = "A",
                                    "trendtype" = i_trend_mthd,
                                    "seasontype" = i_seas_mthd,
                                    "damped" = ifelse(i_mdl_ATA$phi==1, FALSE, TRUE),
                                    "period" = i_mdl_ATA$seasonal.period,
                                    "method" = i_mdl_ATA$method,
                                    "onestep" = i_mdl_ATA$onestep
                      ),
                      intermittent = list("type" = intermittent, 
                                          "method" = paste("IntermittentATA[", intermittent, ", D", substr(d_mdl_ATA$method, 4, 22), ", I", substr(i_mdl_ATA$method, 4, 22), sep="")
                                      )
                      ),
          model_output = list("d_model" = d_mdl_ATA,
                              "i_model" = i_mdl_ATA
                            )
          ),
          class = "intermittentATA")
}

specials_intermittentATA <- fabletools::new_specials(
    d_level = function(parP= NULL, level_fixed = FALSE, initial_level = "none")
                  {
                   list("parP" = parP, "level_fixed" = level_fixed, "initial_level" = initial_level)
                  },
    d_trend = function(type = "A", parQ = NULL, initial_trend = "none", trend_opt = "none",
                    parPHI = NULL, parPHI_range = c(0.0, 1.0), parPHI_increment = 0.05,
                    uroot_test = "adf", uroot_alpha = 0.05, uroot_type = "level", uroot_maxd = 2)
                   {
                     if (type == "N"){
                       parQ = 0
                       parPHI = 1
                       warning("Q has been set 0 because of no trend option.")
                     }
                     if (type == "Ad" | type == "Md"){
                       parPHI = NULL
                       warning("PHI has been set NULL as damped trend is choosen.")
                     }
                     list("type" = type, "parQ" = parQ, "initial_trend" = initial_trend, "trend_opt" = trend_opt,
                          "parPHI" = parPHI, "parPHI_range" = parPHI_range, "parPHI_increment" = parPHI_increment,
                          "uroot_test" = uroot_test, "uroot_alpha" = uroot_alpha, "uroot_type" = uroot_type, "uroot_maxd" = uroot_maxd)
                  },
    d_accuracy = function(criteria = "sMAPE", nmse = 3, ic = "AIC")
                    {
                       if (nmse > 30 & (criteria == "AMSE" | criteria == "GAMSE")) {
                         nmse <- 30
                         warning("'nmse' must be less than 30. 'nmse' is set to 30.")
                       }else if ((is.null(nmse) | nmse <= 1) & (criteria == "AMSE" | criteria == "GAMSE")) {
                         nmse <- 3
                         warning("'nmse' must be greater than 1. 'nmse' is set to 3.")
                       }else{
                       }
                        list("criteria" = criteria, "nmse" = nmse, "ic" = ic)
                    },
    d_transform = function(method = "none", order = "none", lambda = NULL, shift = 0,
                      bcMethod = "guerrero", bcLower = 0, bcUpper = 5)
                      {
                        list("method" = method, "order" = order, "lambda" = lambda, "shift" = shift,
                        "bcMethod" = bcMethod, "bcLower" = bcLower, "bcUpper" = bcUpper)
                      },
    d_holdout = function(holdout = FALSE, adjustment = TRUE, set_size = NULL, onestep = FALSE )
                    {
                        list("holdout" = holdout, "adjustment" = adjustment, "set_size" = set_size, "onestep" = onestep)
                    },
    i_level = function(parP= NULL, level_fixed = FALSE, initial_level = "none")
                  {
                   list("parP" = parP, "level_fixed" = level_fixed, "initial_level" = initial_level)
                  },
    i_trend = function(type = "A", parQ = NULL, initial_trend = "none", trend_opt = "none",
                    parPHI = NULL, parPHI_range = c(0.0, 1.0), parPHI_increment = 0.05,
                    uroot_test = "adf", uroot_alpha = 0.05, uroot_type = "level", uroot_maxd = 2)
                   {
                     if (type == "N"){
                       parQ = 0
                       parPHI = 1
                       warning("Q has been set 0 because of no trend option.")
                     }
                     if (type == "Ad" | type == "Md"){
                       parPHI = NULL
                       warning("PHI has been set NULL as damped trend is choosen.")
                     }
                     list("type" = type, "parQ" = parQ, "initial_trend" = initial_trend, "trend_opt" = trend_opt,
                          "parPHI" = parPHI, "parPHI_range" = parPHI_range, "parPHI_increment" = parPHI_increment,
                          "uroot_test" = uroot_test, "uroot_alpha" = uroot_alpha, "uroot_type" = uroot_type, "uroot_maxd" = uroot_maxd)
                  },
    i_accuracy = function(criteria = "sMAPE", nmse = 3, ic = "AIC")
                    {
                       if (nmse > 30 & (criteria == "AMSE" | criteria == "GAMSE")) {
                         nmse <- 30
                         warning("'nmse' must be less than 30. 'nmse' is set to 30.")
                       }else if ((is.null(nmse) | nmse <= 1) & (criteria == "AMSE" | criteria == "GAMSE")) {
                         nmse <- 3
                         warning("'nmse' must be greater than 1. 'nmse' is set to 3.")
                       }else{
                       }
                        list("criteria" = criteria, "nmse" = nmse, "ic" = ic)
                    },
    i_transform = function(method = "none", order = "none", lambda = NULL, shift = 0,
                      bcMethod = "guerrero", bcLower = 0, bcUpper = 5)
                      {
                        list("method" = method, "order" = order, "lambda" = lambda, "shift" = shift,
                        "bcMethod" = bcMethod, "bcLower" = bcLower, "bcUpper" = bcUpper)
                      },
    i_holdout = function(holdout = FALSE, adjustment = TRUE, set_size = NULL, onestep = FALSE )
                    {
                        list("holdout" = holdout, "adjustment" = adjustment, "set_size" = set_size, "onestep" = onestep)
                    },
    intermittent = function(type = "croston")
                    {
                        type
                    },
   .required_specials = c("d_level", "d_trend", "d_accuracy", "d_transform", "d_holdout", "i_level", "i_trend", "i_accuracy", "i_transform", "i_holdout", "intermittent")
)

#' Intermittent demand time series analysis using the Ata Method based on the ATAforecasting package.
#'
#' Intermittent Ata Method is based on Croston's (1972) <doi:10.2307/3007885> method for intermittent demand forecasting, also described in Shenstone and Hyndman (2005) <doi:10.1002/for.963>. 
#' Croston's method involves using simple exponential smoothing (SES) on the non-zero elements of the time series and a separate application of (SES) to the times between non-zero elements of the time series.
#'
#' There are two variant methods available which apply multiplicative correction factors
#' to the forecasts that result from the original Croston's method. For the
#' Syntetos-Boylan approximation (`type = "sba"`), this factor is \eqn{1 - \alpha / 2},
#' and for the Shale-Boylan-Johnston method (`type = "sbj"`), this factor is
#' \eqn{1 - \alpha / (2 - \alpha)}, where \eqn{\alpha} is the smoothing parameter for
#' the interval SES application.
#' 
#' Returns IntermittentATA[intermittent.type, D(p,q,phi)(E,T,S), I(p,q,phi)(E,T,S)] applied to time series data.
#' 
#' The Ata method based on the modified simple exponential smoothing as described in Yapar, G. (2016) <doi:10.15672/HJMS.201614320580> ,
#' Yapar G., Capar, S., Selamlar, H. T., Yavuz, I. (2017) <doi:10.15672/HJMS.2017.493> and Yapar G., Selamlar, H. T., Capar, S., Yavuz, I. (2019)
#' <doi:10.15672/hujms.461032> is a new univariate time series forecasting method which provides innovative solutions to issues faced during
#' the initialization and optimization stages of existing methods.
#' Forecasting performance of the Ata method is superior to existing methods both in terms of easy implementation and accurate forecasting.
#' It can be applied to non-seasonal or seasonal time series which can be decomposed into four components (remainder, level, trend and seasonal).
#' This methodology performed well on the M3 and M4-competition data.
#'
#' @param formula Model specification (see "Specials" section).
#' @param ... Other arguments
#'
#' @section Specials:
#'
#' The _specials_ define the methods and parameters for the components (level, trend, accuracy, transform, holdout) of an ATA method for both of demand and interval time series.
#'
#' There are a couple of limitations to note about ATA method:
#'
#' - It supports only additive error term.
#' - It does not support exogenous regressors.
#' - It does not support missing values. You can complete missing values in the data with imputed values (e.g. with [tsibble::fill_gaps()], [tidyr::fill()], or by fitting a different model type and then calling [fabletools::interpolate()]) before fitting the model.
#'
#' \subsection{d_level}{
#' The `level` special is used to specify the form of the level term for the demand time series.
#' \preformatted{
#' d_level(parP = NULL, level_fixed = TRUE, initial_level = "none")
#' }
#'
#' \tabular{ll}{
#'   `parP`     \tab The value of the smoothing parameter for the level. If `p = 0`, the level will not change over time. Conversely, if `p = 1` the level will update similarly to a random walk process. If NULL or "opt", it is estimated. \code{p} has all integer values from 1 to \code{length(data)}. \cr
#'   `level_fixed`      \tab If TRUE, "pStarQ"  --> First, fits ATA(p,0) where p = p* is optimized for q=0. Then, fits ATA(p*,q) where q is optimized for p = p*. \cr
#'   `initial_level`     \tab If NULL, "none" is default. If "none", ATA Method calculates the pth observation in \code{data} for level. If "mean", ATA Method calculates average of first p value in \code{data}for level. If "median", ATA Method calculates median of first p value in \code{data}for level. \cr
#'
#' }
#' }
#'
#' \subsection{d_trend}{
#' The `trend` special is used to specify the form of the trend term and associated parameters for the demand time series.
#' \preformatted{
#' d_trend(type = "A", parQ = NULL, initial_trend = "none", opt_trend = "none",
#'        parPHI = NULL, parPHI_range = c(0.8, 1.0), parPHI_increment = 0.01,
#'        uroot_test = "adf", uroot_alpha = 0.05, uroot_type = "level")
#' }
#'
#' \tabular{ll}{
#'   `type`     \tab The form of the trend term: either none ("N"), additive ("A"), multiplicative ("M") or damped variants ("Ad", "Md"). \cr
#'   `parQ`      \tab The value of the smoothing parameter for the slope. If `q = 0`, the slope will not change over time. Conversely, if `q = 1` the slope will have mean of past slopes. \cr
#'   `parPHI` \tab The value of the dampening parameter for the slope. If `phi = 0`, the slope will be dampened immediately (no slope). Conversely, if `phi = 1` the slope will not be dampened. \cr
#'   `parPHI_range`       \tab If `phi=NULL`, `phi_range` provides bounds for the optimised value of `phi`.\cr
#'   `parPHI_increment`  \tab If `phi=NULL`, `parPHI_increment` provides increment step for searching `phi`. If NULL, `parPHI_increment` will be determined as the value that allows the `parPHI_range` to be divided into 20 equal parts. \cr
#'   `initial_trend`     \tab If NULL, "none" is default. If "none", ATA Method calculates the qth observation in \code{data} for trend. If "mean", ATA Method calculates average of first q value in \code{X(T)-X(T-1)} for trend. If "median", ATA Method calculates median of first q value in \code{X(T)-X(T-1)} for trend. \cr
#'   `trend_opt`        \tab Default is `none`. If `fixed` is set, "pBullet" --> Fits ATA(p,1) where p = p* is optimized for q = 1. If `search` is set "qBullet" --> Fits ATA(p,q) where p = p* is optimized for q = q* (q > 0). Then, fits ATA(p*,q) where q is optimized for p = p*. \cr
#'   `uroot_test`        \tab Type of unit root test before all type seasonality test. Possible values are "adf", "pp" and "kpss". \cr
#'   `uroot_alpha`   \tab Significant level of the unit root test, possible values range from 0.01 to 0.1. \cr
#'   `uroot_type`        \tab Specification of the deterministic component in the regression for unit root test. Possible values are "level" and "trend". \cr
#'   `uroot_maxd`       \tab Maximum number of non-seasonal differences allowed. \cr
#' }
#' }
#'
#' \subsection{d_accuracy}{
#' The `accuracy` special is used to the optimization criterion for selecting the best ATA Method forecasting for the demand time series.
#' \preformatted{
#' d_accuracy(criteria = "sMAPE", nmse = 3, ic = "AIC")
#' }
#'
#' \tabular{ll}{
#'   `criteria`     \tab Accuracy measure for optimization of the best ATA Method forecasting. IF NULL, `sMAPE` is default. \cr
#'   `nmse`     \tab If `accuracy.type == "AMSE"`, `nmse` provides the number of steps for average multistep MSE `(2<=nmse<=30)`. \cr
#'   `ic`     \tab The information criterion used in selecting the model.  \cr
#' }
#' }
#'
#' \subsection{d_transform}{
#' The `transform` special is used to provide the applicability of different types of transformation techniques for the demand data to which the ATA method will be applied.
#' \preformatted{
#' d_transform(method="none", order = "none", lambda = NULL, shift = 0,
#'           bcMethod = "guerrero", bcLower = 0, bcUpper = 5)
#' }
#'
#' \tabular{ll}{
#'   `method`     \tab Transformation method  --> "Box_Cox", "Sqrt", "Reciprocal", "Log", "NegLog", "Modulus", "BickelDoksum", "Manly", "Dual", "YeoJohnson", "GPower", "GLog". If the transformation process needs shift parameter, it will be calculated required shift parameter automatically. \cr
#'   `order`     \tab Default is "none. If "before", Box-Cox transformation family will be applied and then seasonal decomposition techniques will be applied. If "after", seasonal decomposition techniques will be applied and then Box-Cox transformation family will be applied. \cr
#'   `lambda`     \tab Box-Cox power transformation family parameter. If NULL, data transformed before model is estimated.  \cr
#'   `shift`     \tab Box-Cox power transformation family shifting parameter. If NULL, data transformed before model is estimated. \cr
#'   `bcMethod`     \tab Choose method to be used in calculating lambda. "guerrero" is default. Other method is "loglik". \cr
#'   `bcLower`     \tab Lower limit for possible lambda values. The lower value is limited by -5. Default value is 0. \cr
#'   `bcUpper`     \tab Upper limit for possible lambda values. The upper value is limited by 5. Default value is 5. \cr
#' }
#' }
#'
#' \subsection{d_holdout}{
#' The `holdout` special is used to improve the optimized parameter value obtained for the ATA Method forecasting for the demand time series.
#' \preformatted{
#' d_holdout(holdout = FALSE, adjustment = TRUE, set_size = NULL, onestep = FALSE)
#' }
#'
#' \tabular{ll}{
#'   `holdout`     \tab Default is FALSE. If TRUE, ATA Method uses the holdout forecasting for accuracy measure to select the best parameter set. In holdout forecasting, this parameter divides `data` into two parts: training set (in-sample) and validation set (holdout set). \cr
#'   `adjustment`     \tab Default is TRUE. If TRUE, `parP` will be adjusted by length of training, validation sets and main data set when the holdout forecasting is active. \cr
#'   `set_size`     \tab If `holdout` is TRUE, this parameter divides `data` into two parts: training set (in-sample) and validation set (holdout set). Also, this parameter will be same as `h` for defining holdout set.  \cr
#'   `onestep`     \tab Default is FALSE. if TRUE, the dynamic forecast strategy uses a one-step model multiple times `h` (forecast horizon) where the holdout prediction for the prior time step is used as an input for making a prediction on the following time step.
#' }
#' }
#'
#' \subsection{i_level}{
#' The `level` special is used to specify the form of the level term for the interval time series.
#' \preformatted{
#' i_level(parP = NULL, level_fixed = TRUE, initial_level = "none")
#' }
#'
#' \tabular{ll}{
#'   `parP`     \tab The value of the smoothing parameter for the level. If `p = 0`, the level will not change over time. Conversely, if `p = 1` the level will update similarly to a random walk process. If NULL or "opt", it is estimated. \code{p} has all integer values from 1 to \code{length(data)}. \cr
#'   `level_fixed`      \tab If TRUE, "pStarQ"  --> First, fits ATA(p,0) where p = p* is optimized for q=0. Then, fits ATA(p*,q) where q is optimized for p = p*. \cr
#'   `initial_level`     \tab If NULL, "none" is default. If "none", ATA Method calculates the pth observation in \code{data} for level. If "mean", ATA Method calculates average of first p value in \code{data}for level. If "median", ATA Method calculates median of first p value in \code{data}for level. \cr
#'
#' }
#' }
#'
#' \subsection{i_trend}{
#' The `trend` special is used to specify the form of the trend term and associated parameters for the interval time series.
#' \preformatted{
#' i_trend(type = "A", parQ = NULL, initial_trend = "none", opt_trend = "none",
#'        parPHI = NULL, parPHI_range = c(0.8, 1.0), parPHI_increment = 0.01,
#'        uroot_test = "adf", uroot_alpha = 0.05, uroot_type = "level")
#' }
#'
#' \tabular{ll}{
#'   `type`     \tab The form of the trend term: either none ("N"), additive ("A"), multiplicative ("M") or damped variants ("Ad", "Md"). \cr
#'   `parQ`      \tab The value of the smoothing parameter for the slope. If `q = 0`, the slope will not change over time. Conversely, if `q = 1` the slope will have mean of past slopes. \cr
#'   `parPHI` \tab The value of the dampening parameter for the slope. If `phi = 0`, the slope will be dampened immediately (no slope). Conversely, if `phi = 1` the slope will not be dampened. \cr
#'   `parPHI_range`       \tab If `phi=NULL`, `phi_range` provides bounds for the optimised value of `phi`.\cr
#'   `parPHI_increment`  \tab If `phi=NULL`, `parPHI_increment` provides increment step for searching `phi`. If NULL, `parPHI_increment` will be determined as the value that allows the `parPHI_range` to be divided into 20 equal parts. \cr
#'   `initial_trend`     \tab If NULL, "none" is default. If "none", ATA Method calculates the qth observation in \code{data} for trend. If "mean", ATA Method calculates average of first q value in \code{X(T)-X(T-1)} for trend. If "median", ATA Method calculates median of first q value in \code{X(T)-X(T-1)} for trend. \cr
#'   `trend_opt`        \tab Default is `none`. If `fixed` is set, "pBullet" --> Fits ATA(p,1) where p = p* is optimized for q = 1. If `search` is set "qBullet" --> Fits ATA(p,q) where p = p* is optimized for q = q* (q > 0). Then, fits ATA(p*,q) where q is optimized for p = p*. \cr
#'   `uroot_test`        \tab Type of unit root test before all type seasonality test. Possible values are "adf", "pp" and "kpss". \cr
#'   `uroot_alpha`   \tab Significant level of the unit root test, possible values range from 0.01 to 0.1. \cr
#'   `uroot_type`        \tab Specification of the deterministic component in the regression for unit root test. Possible values are "level" and "trend". \cr
#'   `uroot_maxd`       \tab Maximum number of non-seasonal differences allowed. \cr
#' }
#' }
#'
#' \subsection{i_accuracy}{
#' The `accuracy` special is used to the optimization criterion for selecting the best ATA Method forecasting for the interval time series.
#' \preformatted{
#' i_accuracy(criteria = "sMAPE", nmse = 3, ic = "AIC")
#' }
#'
#' \tabular{ll}{
#'   `criteria`     \tab Accuracy measure for optimization of the best ATA Method forecasting. IF NULL, `sMAPE` is default. \cr
#'   `nmse`     \tab If `accuracy.type == "AMSE"`, `nmse` provides the number of steps for average multistep MSE `(2<=nmse<=30)`. \cr
#'   `ic`     \tab The information criterion used in selecting the model.  \cr
#' }
#' }
#'
#' \subsection{i_transform}{
#' The `transform` special is used to provide the applicability of different types of transformation techniques for the interval data to which the ATA method will be applied.
#' \preformatted{
#' i_transform(method="none", order = "none", lambda = NULL, shift = 0,
#'           bcMethod = "guerrero", bcLower = 0, bcUpper = 5)
#' }
#'
#' \tabular{ll}{
#'   `method`     \tab Transformation method  --> "Box_Cox", "Sqrt", "Reciprocal", "Log", "NegLog", "Modulus", "BickelDoksum", "Manly", "Dual", "YeoJohnson", "GPower", "GLog". If the transformation process needs shift parameter, it will be calculated required shift parameter automatically. \cr
#'   `order`     \tab Default is "none. If "before", Box-Cox transformation family will be applied and then seasonal decomposition techniques will be applied. If "after", seasonal decomposition techniques will be applied and then Box-Cox transformation family will be applied. \cr
#'   `lambda`     \tab Box-Cox power transformation family parameter. If NULL, data transformed before model is estimated.  \cr
#'   `shift`     \tab Box-Cox power transformation family shifting parameter. If NULL, data transformed before model is estimated. \cr
#'   `bcMethod`     \tab Choose method to be used in calculating lambda. "guerrero" is default. Other method is "loglik". \cr
#'   `bcLower`     \tab Lower limit for possible lambda values. The lower value is limited by -5. Default value is 0. \cr
#'   `bcUpper`     \tab Upper limit for possible lambda values. The upper value is limited by 5. Default value is 5. \cr
#' }
#' }
#'
#' \subsection{i_holdout}{
#' The `holdout` special is used to improve the optimized parameter value obtained for the ATA Method forecasting for the interval time series.
#' \preformatted{
#' i_holdout(holdout = FALSE, adjustment = TRUE, set_size = NULL, onestep = FALSE)
#' }
#'
#' \tabular{ll}{
#'   `holdout`     \tab Default is FALSE. If TRUE, ATA Method uses the holdout forecasting for accuracy measure to select the best parameter set. In holdout forecasting, this parameter divides `data` into two parts: training set (in-sample) and validation set (holdout set). \cr
#'   `adjustment`     \tab Default is TRUE. If TRUE, `parP` will be adjusted by length of training, validation sets and main data set when the holdout forecasting is active. \cr
#'   `set_size`     \tab If `holdout` is TRUE, this parameter divides `data` into two parts: training set (in-sample) and validation set (holdout set). Also, this parameter will be same as `h` for defining holdout set.  \cr
#'   `onestep`     \tab Default is FALSE. if TRUE, the dynamic forecast strategy uses a one-step model multiple times `h` (forecast horizon) where the holdout prediction for the prior time step is used as an input for making a prediction on the following time step.
#' }
#' }
#'
#' \subsection{intermittent}{
#' The `intermittent` special is used to improve the optimized parameter value obtained for the ATA Method forecasting for the interval time series.
#' \preformatted{
#' intermittent(type = "croston")
#' }
#'
#' \tabular{ll}{
#'   `type`     \tab Default is "croston". For the Syntetos-Boylan approximation (`type = "sba"`), this factor is \eqn{1 - \alpha / 2}, and for the Shale-Boylan-Johnston method (`type = "sbj"`), this factor is \eqn{1 - \alpha / (2 - \alpha)}, where \eqn{\alpha} is the smoothing parameter for the interval ATA Method application. \cr
#' }
#' }
#' 
#' @return A model specification.
#'
#' @importFrom fabletools new_model_class new_model_definition
#' @importFrom rlang enquo
#'
#'@examples
#' library(intermittentATA)
#' as_tsibble(fmcgData) %>% model(crostonata = intermittentATA(value ~ d_trend(type = "M", parQ = 1) + i_trend("A") + intermittent("croston")))
#'
#' @export
intermittentATA <- function(formula, ...){
        intermittent_atam_model <- fabletools::new_model_class("intermittentATA", train = train_intermittentATA, specials = specials_intermittentATA)
        fabletools::new_model_definition(intermittent_atam_model, !!rlang::enquo(formula))
}

#' Forecast a model from the fable intermittentATA model
#'
#'
#' @param object The time series model used to produce the forecasts
#' @param new_data A `tsibble` containing future information used to forecast.
#' @param h The forecast horison (can be used instead of `new_data` for regular time series with no exogenous regressors).
#' @param ci_level Confidence Interval levels for forecasting. Default value is 95.
#' @param negative_forecast Negative values are allowed for forecasting. Default value is TRUE. If FALSE, all negative values for forecasting are set to 0.
#' @param onestep Default is FALSE. if TRUE, the dynamic forecast strategy uses a one-step model multiple times `h` forecast horizon) where the prediction for the prior time step is used as an input for making a prediction on the following time step.
#' @param ... Other arguments
#'
#' @return A vector of fitted residuals.
#'
#' @examples
#' library(intermittentATA)
#' as_tsibble(fmcgData) %>% 
#'    model(crostonata = intermittentATA(value ~ d_trend(type = "M", parQ = 1) + i_trend("A") + intermittent("croston"))) %>% forecast(h=6)
#'
#' @importFrom tsibble is_tsibble as_tibble measured_vars index
#' @importFrom rlang enquo expr_text
#' @importFrom tsbox ts_ts
#' @importFrom stats frequency ts start
#' @importFrom ATAforecasting ATA.Forecast
#' @importFrom distributional dist_degenerate
#'
#' @export
forecast.intermittentATA <- function(object, new_data, h=NULL, ci_level=95, negative_forecast=TRUE, onestep=FALSE, ...){
  mdl <- object$model_output
  spec_mdl <- object$spec
  mh <- nrow(new_data)
  if (is.null(mh)){
    if (spec_mdl$period==4){
      mh <- 8
    }else if (spec_mdl$period==5){
      mh <- 10
    }else if (spec_mdl$period==7){
      mh <- 14
    }else if (spec_mdl$period==12){
      mh <- 24
    }else if (spec_mdl$period==24){
      mh <- 48
    }else {
      mh <- 6
    }
  }

# Prepare data and forecast
 if (length(tsibble::measured_vars(new_data)) == 0){
      fc_demand_rate <- rep(object$est$.fitted[length(object$est$.fitted)], mh)
  }else {
    test_set <- tsibble::as_tibble(new_data)[c(rlang::expr_text(tsibble::index(new_data)), tsibble::measured_vars(new_data))]
    colnames(test_set) <- c("ds", "yh")

  # Check data
    if (any(is.na(test_set$yh))) {
      stop("ATA method does not support missing values.")
    }
    if (any(test_set$yh < 0)) {
      stop("All observations must be non-negative for Intermittent ATA method.")
    }
    non_zero <- which(test_set$yh != 0)

    # Get response
    # Croston demand/interval decomposition
    yh_demand <- test_set$yh[non_zero]
    yh_interval <- c(non_zero[1], diff(non_zero))

    test_demand <- stats::ts(yh_demand, start = end(mdl$actual) + 1, frequency = 1)
    test_interval <- stats::ts(yh_interval, start = end(mdl$actual) + 1, frequency = 1)

    d_pfc <- safely(quietly(ATAforecasting::ATA.Forecast))(mdl$d_model, mh, test_demand, ci_level, negative_forecast, onestep, print.out = FALSE)
    i_pfc <- safely(quietly(ATAforecasting::ATA.Forecast))(mdl$i_model, mh, test_interval, ci_level, negative_forecast, onestep, print.out = FALSE)

    d_fc <- d_pfc[["result"]]
    i_fc <- i_pfc[["result"]]
  
     interval_param <- mdl$i_model$p / length(mdl$i_model$actual)
  
     if (spec_mdl$intermittent$type == "sba") {
        coeff <- 1 - interval_param / 2
      } else if (spec_mdl$intermittent$type == "sbj") {
        coeff <- 1 - interval_param / (2 - interval_param)
      } else {
        coeff <- 1
      }
      ratio <- coeff * d_fc$forecast / i_fc$forecast
    
      # Out-of-sample demand rate
      fc_demand_rate <- rep(c(NA_real_, ratio), diff(c(0, non_zero, length(test_set$yh))))
    }
 
  # Return forecasts
  distributional::dist_degenerate(fc_demand_rate)
}

#' Extract fitted values
#'
#' Extracts the fitted values from an estimated intermittentATA model.
#'
#' @param object An estimated model.
#' @param ... Unused.
#'
#' @return A vector of fitted values.
#'
#' @export
fitted.intermittentATA <- function(object, ...){
  object$est[[".fitted"]]
}

#' Extract model residuals
#'
#' Extracts the residuals from an estimated intermittentATA model.
#'
#' @param object An estimated model.
#' @param ... Unused.
#'
#' @return A vector of residuals.
#'
#' @export
residuals.intermittentATA <- function(object, ...){
  object$est[[".resid"]]
}


#' Extract estimated states from an intermittentATA model.
#'
#' @param object An estimated model.
#' @param ... Unused.
#'
#' @return A [fabletools::dable()] containing estimated states.
#'
#' @importFrom fabletools as_dable
#' @importFrom rlang sym expr list2 ":="
#' @importFrom tsibble measured_vars index
#' @importFrom dplyr transmute mutate select
#'
#' @examples
#' library(intermittentATA)
#' as_tsibble(fmcgData) %>% 
#'    model(crostonata = intermittentATA(value ~ d_trend(type = "M", parQ = 1) + i_trend("A") + intermittent("croston"))) %>% components()
#' 
#' @export
components.intermittentATA <- function(object, ...){
  d_cmp <- object$components$d_components
  i_cmp <- object$components$i_components

  d_spec <- object$spec$d_spec
  i_spec <- object$spec$i_spec

  response <- tsibble::measured_vars(object$est)[[1]]
  est_vars <- dplyr::transmute(object$est,
                        !!rlang::sym(response),
                        fitted = !!rlang::sym(".fitted"),
                        remainder = !!rlang::sym(".resid")
                       )

    # Demand model
  d_response <- tsibble::measured_vars(object$est_comp[[1]])[[1]]
  d_est_vars <- dplyr::transmute(object$est_comp[[1]],
                        !!rlang::sym(d_response),
                        remainder = !!rlang::sym(".resid")
                       )
  d_idx <- tsibble::index(d_est_vars)
  d_eqn <- rlang::expr(!!rlang::sym("level"))

  if (d_spec$trendtype == "A") {
      d_eqn <- rlang::expr(!!d_eqn + !!rlang::sym("trend"))
  } else if (d_spec$trendtype == "M") {
      d_eqn <- rlang::expr(!!d_eqn * !!rlang::sym("trend"))
  }
  if (d_spec$seasontype == "A") {
    d_eqn <- rlang::expr(!!d_eqn + !!rlang::sym("season"))
  } else if (d_spec$seasontype == "M") {
    d_eqn <- rlang::expr((!!d_eqn) * !!rlang::sym("season"))
  }
  d_eqn <- rlang::expr(!!d_eqn + !!rlang::sym("remainder"))
  if (d_spec$trendtype == "N" & d_spec$seasontype != "N"){
    d_f_cmp <- dplyr::mutate(dplyr::ungroup(d_est_vars),
                            "level" = d_cmp$level,
                            "season" = d_cmp$season)
    d_f_cmp <- dplyr::select(d_f_cmp, intersect(c(expr_text(d_idx), d_response, "level", "season", "remainder"), colnames(d_f_cmp)))
    colnames(d_f_cmp) <- c(expr_text(d_idx), "demand", "d_level", "d_season", "d_remainder")
    seasonality <- list("season" = d_cmp$season)
  } else if (d_spec$trendtype != "N" & d_spec$seasontype == "N"){
    d_f_cmp <- dplyr::mutate(dplyr::ungroup(d_est_vars),
                            "level" = d_cmp$level,
                            "trend" = d_cmp$trend)
    d_f_cmp <- dplyr::select(d_f_cmp, intersect(c(expr_text(d_idx), d_response, "level", "trend", "remainder"), colnames(d_f_cmp)))
    colnames(d_f_cmp) <- c(expr_text(d_idx), "demand", "d_level", "d_trend", "d_remainder")
    seasonality <- list()
  }else if (d_spec$trendtype == "N" & d_spec$seasontype == "N") {
    d_f_cmp <- dplyr::mutate(dplyr::ungroup(d_est_vars),
                            "level" = d_cmp$level)
    d_f_cmp <- dplyr::select(d_f_cmp, intersect(c(expr_text(d_idx), d_response, "level", "remainder"), colnames(d_f_cmp)))
    colnames(d_f_cmp) <- c(expr_text(d_idx), "demand", "d_level", "d_remainder")
    seasonality <- list()
  }else {
    d_f_cmp <- dplyr::mutate(dplyr::ungroup(d_est_vars),
                            "level" = d_cmp$level,
                            "trend" = d_cmp$trend,
                            "season" = d_cmp$season)
    d_f_cmp <- dplyr::select(d_f_cmp, intersect(c(expr_text(d_idx), d_response, "level", "trend", "season", "remainder"), colnames(d_f_cmp)))
    colnames(d_f_cmp) <- c(expr_text(d_idx), "demand", "d_level", "d_trend", "d_season", "d_remainder")
    seasonality <- list("season" = d_cmp$season)
  }

  # Interval model
  i_response <- tsibble::measured_vars(object$est_comp[[2]])[[1]]
  i_est_vars <- dplyr::transmute(object$est_comp[[2]],
                        !!rlang::sym(i_response),
                        remainder = !!rlang::sym(".resid")
                       )
  i_idx <- tsibble::index(i_est_vars)
  i_eqn <- rlang::expr(!!rlang::sym("level"))

 if (i_spec$trendtype == "A") {
      i_eqn <- rlang::expr(!!i_eqn + !!rlang::sym("trend"))
  } else if (i_spec$trendtype == "M") {
      i_eqn <- rlang::expr(!!i_eqn * !!rlang::sym("trend"))
  }
  if (i_spec$seasontype == "A") {
    i_eqn <- rlang::expr(!!i_eqn + !!rlang::sym("season"))
  } else if (i_spec$seasontype == "M") {
    i_eqn <- rlang::expr((!!i_eqn) * !!rlang::sym("season"))
  }
  i_eqn <- rlang::expr(!!i_eqn + !!rlang::sym("remainder"))
  if (i_spec$trendtype == "N" & i_spec$seasontype != "N"){
    i_f_cmp <- dplyr::mutate(dplyr::ungroup(i_est_vars),
                            "level" = i_cmp$level,
                            "season" = i_cmp$season)
    i_f_cmp <- dplyr::select(i_f_cmp, intersect(c(expr_text(i_idx), i_response, "level", "season", "remainder"), colnames(i_f_cmp)))
    colnames(i_f_cmp) <- c(expr_text(i_idx), "interval", "i_level", "i_season", "i_remainder")
    seasonality <- list("season" = i_cmp$season)
  } else if (i_spec$trendtype != "N" & i_spec$seasontype == "N"){
    i_f_cmp <- dplyr::mutate(dplyr::ungroup(i_est_vars),
                            "level" = i_cmp$level,
                            "trend" = i_cmp$trend)
    i_f_cmp <- dplyr::select(i_f_cmp, intersect(c(expr_text(i_idx), i_response, "level", "trend", "remainder"), colnames(i_f_cmp)))
    colnames(i_f_cmp) <- c(expr_text(i_idx), "interval", "i_level", "i_trend", "i_remainder")
    seasonality <- list()
  }else if (i_spec$trendtype == "N" & i_spec$seasontype == "N") {
    i_f_cmp <- dplyr::mutate(dplyr::ungroup(i_est_vars),
                            "level" = i_cmp$level)
    i_f_cmp <- dplyr::select(i_f_cmp, intersect(c(expr_text(i_idx), i_response, "level", "remainder"), colnames(i_f_cmp)))
    colnames(i_f_cmp) <- c(expr_text(i_idx), "interval", "i_level", "i_remainder")
    seasonality <- list()
  }else {
    i_f_cmp <- dplyr::mutate(dplyr::ungroup(i_est_vars),
                            "level" = i_cmp$level,
                            "trend" = i_cmp$trend,
                            "season" = i_cmp$season)
    i_f_cmp <- dplyr::select(i_f_cmp, intersect(c(expr_text(i_idx), i_response, "level", "trend", "season", "remainder"), colnames(i_f_cmp)))
    colnames(i_f_cmp) <- c(expr_text(i_idx), "interval", "i_level", "i_trend", "i_season", "i_remainder")
    seasonality <- list("season" = i_cmp$season)
  }
  c_f_cmp <- d_f_cmp %>% left_join(i_f_cmp, by = "mth")
  colnames(est_vars)[1] <- "mth"
  all_est_vars <- est_vars %>% left_join(c_f_cmp, by = "mth")

  if (object$spec$intermittent$type == "sba") {
    eqn_type <- "coeff(sba)"
  } else if (object$spec$intermittent$type == "sbj") {
    eqn_type <- "coeff(sbj)"
  } else {
    eqn_type <- "coeff(croston)"
  }

  fabletools::as_dable(all_est_vars,
                       resp = !!rlang::sym(response),
                       method = model_sum(object),
                       seasons = list(),
                       aliases = aliases <- tibble::lst(!!rlang::sym(response) := quote(!!eqn_type * demand / interval))
                      )
}

#' Glance an intermittentATA model
#'
#' @param x An estimated model.
#' @param ... Unused.
#'
#' @return A one row tibble summarising the model's fit.
#'
#' @examples
#' library(intermittentATA)
#' as_tsibble(fmcgData) %>% 
#'    model(crostonata = intermittentATA(value ~ d_trend(type = "M", parQ = 1) + i_trend("A") + intermittent("croston"))) %>% glance()
#'
#' @importFrom tibble as_tibble
#'
#' @export
glance.intermittentATA <- function(x, ...){
  tibble::as_tibble(x$fit)
}

#' Tidy a intermittentATA model
#'
#' @param x An estimated model.
#' @param ... Unused.
#'
#' @return The model's coefficients in a `tibble`.
#'
#' @examples
#' library(intermittentATA)
#' as_tsibble(fmcgData) %>% 
#'     model(crostonata = intermittentATA(value ~ d_trend(type = "M", parQ = 1) + i_trend("A") + intermittent("croston"))) %>% tidy()
#'
#' @export
tidy.intermittentATA <- function(x, ...){
  x$par
}

#' Summary of intermittentATA model
#'
#' @param x An estimated model.
#' @param ... Unused.
#'
#' @return The model's summary specs.
#'
#' @export
model_sum.intermittentATA <- function(x, ...){
  paste("IntermittentATA[", x$spec$intermittent$type, ", D", substr(x$spec$d_spec$method, 4, 22), ", I", substr(x$spec$i_spec$method, 4, 22), "]", sep="")
  }

#' Format of intermittentATA model
#'
#' @param x An estimated model.
#' @param ... Unused.
#'
#' @return The forecasting model's name.
#'
#' @export
format.intermittentATA <- function(x, ...){
  "intermittentATA"
}

#' Specialized Screen Print Function of intermittentATA model
#'
#' @param object An estimated model.
#' @param ... Unused.
#'
#' @return a summary for the results of the ATAforecasting 
#'
#' @examples
#'  library(intermittentATA)
#'  as_tsibble(fmcgData) %>% 
#'      model(crostonata = intermittentATA(value ~ d_trend(type = "M", parQ = 1) + i_trend("A") + intermittent("croston"))) %>% report()
#'
#' @export
report.intermittentATA <- function(object, ...) {
    opscipen <- options("scipen" = 100, "digits"=7)
    on.exit(options(opscipen))

    cat("----------Intermittent ATA Model------------------------", "\n\n")
    cat(paste("Intermittent Type:", object$spec$intermittent$type, "\n\n", sep=""))

    # Demand model
	  x <- object$model_output$d_model
    cat("----------ATA Model for Demand--------------------------", "\n\n")
    cat(x$method,"\n\n")
    if (x$level.fixed==TRUE){
      cat("   level.fixed: TRUE", "\n\n")
    }
    if (x$trend.opt!="none"){
      cat(paste("   trend optimization method: trend.", x$trend.opt, "\n\n", sep=""))
    }
    if(!is.null(x$transform.method)){
      cat(paste("   '",x$transform.method, "' transformation method was selected.","\n\n", sep=""))
    }
    if(!is.null(x$lambda)){
      cat("   Box-Cox transformation: lambda=",round(x$lambda,4), "\n\n")
    }
    cat(paste("   model.type:",x$model.type, "\n\n"))
    if (x$is.season==FALSE){
      cat("   seasonal.model: no seasonality","\n\n")
    }else {
      cat(paste("   seasonal.model:",x$seasonal.model, "\n\n"))
    }
    if (x$is.season==TRUE){
      cat(paste("   seasonal.type:",x$seasonal.type, "\n\n"))
    }
    cat(paste("   forecast horizon:",x$h, "\n\n"))
    cat(paste("   accuracy.type:",x$accuracy.type, "\n\n"))

    cat("In-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$fits$sigma2, x$accuracy$fits$loglik, x$accuracy$MAE$inSample$MAE, x$accuracy$MSE$inSample$MSE, x$accuracy$MSE$inSample$RMSE, x$accuracy$MPE$inSample$MPE, x$accuracy$MAPE$inSample$MAPE, x$accuracy$sMAPE$inSample$sMAPE, x$accuracy$MASE$inSample$MASE, x$accuracy$OWA$inSample$OWA)
    names(stats) <- c("sigma2", "loglik", "MAE", "MSE", "RMSE", "MPE", "MAPE", "sMAPE", "MASE", "OWA")
    cat("\n")
    print(stats)
    cat("\n")

    cat("In-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$MAE$inSample$MdAE, x$accuracy$MSE$inSample$MdSE, x$accuracy$MSE$inSample$RMdSE, x$accuracy$MPE$inSample$MdPE, x$accuracy$MAPE$inSample$MdAPE, x$accuracy$sMAPE$inSample$sMdAPE)
    names(stats) <- c("MdAE", "MdSE", "RMdSE", "MdPE", "MdAPE", "sMdAPE")
    cat("\n")
    print(stats)
    cat("\n")

    cat("Out-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$MAE$outSample$MAE, x$accuracy$MSE$outSample$MSE, x$accuracy$MSE$outSample$RMSE, x$accuracy$MPE$outSample$MPE, x$accuracy$MAPE$outSample$MAPE, x$accuracy$sMAPE$outSample$sMAPE, x$accuracy$MASE$outSample$MASE, x$accuracy$OWA$outSample$OWA)
    names(stats) <- c("MAE", "MSE", "RMSE", "MPE", "MAPE", "sMAPE", "MASE",  "OWA")
    cat("\n")
    print(stats)
    cat("\n")

    cat("Out-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$MAE$outSample$MdAE, x$accuracy$MSE$outSample$MdSE, x$accuracy$MSE$outSample$RMdSE, x$accuracy$MPE$outSample$MdPE, x$accuracy$MAPE$outSample$MdAPE, x$accuracy$sMAPE$outSample$sMdAPE)
    names(stats) <- c("MdAE", "MdSE", "RMdSE", "MdPE", "MdAPE", "sMdAPE")
    cat("\n")
    print(stats)
    cat("\n")

    cat("Information Criteria:","\n")
    stats <- c(x$accuracy$fits$AIC, x$accuracy$fits$AICc, x$accuracy$fits$BIC)
    names(stats) <- c("AIC", "AICc", "BIC")
    cat("\n")
    print(stats)
    cat("\n")

    stats <- c(x$execution.time[1], x$execution.time[2], x$execution.time[3])
    names(stats) <- c("user","system","elapsed")
    cat("\n")
    print(stats)
    cat("\n")
    cat(paste("calculation.time:",x$calculation.time, "\n\n"))
    cat("\n")


    # Interval model
    x <- object$model_output$i_model
    cat("---------ATA Model for Interval------------------------n\n")
    cat(x$method,"\n\n")
    if (x$level.fixed==TRUE){
      cat("   level.fixed: TRUE","\n\n")
    }
    if (x$trend.opt!="none"){
      cat(paste("   trend optimization method: trend.", x$trend.opt, "\n\n", sep=""))
    }
    if(!is.null(x$transform.method)){
      cat(paste("   '",x$transform.method, "' transformation method was selected.","\n\n", sep=""))
    }
    if(!is.null(x$lambda)){
      cat("   Box-Cox transformation: lambda=",round(x$lambda,4), "\n\n")
    }
    cat(paste("   model.type:",x$model.type, "\n\n"))
    if (x$is.season==FALSE){
      cat("   seasonal.model: no seasonality","\n\n")
    }else {
      cat(paste("   seasonal.model:",x$seasonal.model, "\n\n"))
    }
    if (x$is.season==TRUE){
      cat(paste("   seasonal.type:",x$seasonal.type, "\n\n"))
    }
    cat(paste("   forecast horizon:",x$h, "\n\n"))
    cat(paste("   accuracy.type:",x$accuracy.type, "\n\n"))

    cat("In-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$fits$sigma2, x$accuracy$fits$loglik, x$accuracy$MAE$inSample$MAE, x$accuracy$MSE$inSample$MSE, x$accuracy$MSE$inSample$RMSE, x$accuracy$MPE$inSample$MPE, x$accuracy$MAPE$inSample$MAPE, x$accuracy$sMAPE$inSample$sMAPE, x$accuracy$MASE$inSample$MASE, x$accuracy$OWA$inSample$OWA)
    names(stats) <- c("sigma2", "loglik", "MAE", "MSE", "RMSE", "MPE", "MAPE", "sMAPE", "MASE", "OWA")
    cat("\n")
    print(stats)
    cat("\n")

    cat("In-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$MAE$inSample$MdAE, x$accuracy$MSE$inSample$MdSE, x$accuracy$MSE$inSample$RMdSE, x$accuracy$MPE$inSample$MdPE, x$accuracy$MAPE$inSample$MdAPE, x$accuracy$sMAPE$inSample$sMdAPE)
    names(stats) <- c("MdAE", "MdSE", "RMdSE", "MdPE", "MdAPE", "sMdAPE")
    cat("\n")
    print(stats)
    cat("\n")

    cat("Out-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$MAE$outSample$MAE, x$accuracy$MSE$outSample$MSE, x$accuracy$MSE$outSample$RMSE, x$accuracy$MPE$outSample$MPE, x$accuracy$MAPE$outSample$MAPE, x$accuracy$sMAPE$outSample$sMAPE, x$accuracy$MASE$outSample$MASE, x$accuracy$OWA$outSample$OWA)
    names(stats) <- c("MAE", "MSE", "RMSE", "MPE", "MAPE", "sMAPE", "MASE",  "OWA")
    cat("\n")
    print(stats)
    cat("\n")

    cat("Out-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$MAE$outSample$MdAE, x$accuracy$MSE$outSample$MdSE, x$accuracy$MSE$outSample$RMdSE, x$accuracy$MPE$outSample$MdPE, x$accuracy$MAPE$outSample$MdAPE, x$accuracy$sMAPE$outSample$sMdAPE)
    names(stats) <- c("MdAE", "MdSE", "RMdSE", "MdPE", "MdAPE", "sMdAPE")
    cat("\n")
    print(stats)
    cat("\n")

    cat("Information Criteria:","\n")
    stats <- c(x$accuracy$fits$AIC, x$accuracy$fits$AICc, x$accuracy$fits$BIC)
    names(stats) <- c("AIC", "AICc", "BIC")
    cat("\n")
    print(stats)
    cat("\n")

    stats <- c(x$execution.time[1], x$execution.time[2], x$execution.time[3])
    names(stats) <- c("user","system","elapsed")
    cat("\n")
    print(stats)
    cat("\n")
    cat(paste("calculation.time:",x$calculation.time, "\n\n"))
    cat("\n")
}