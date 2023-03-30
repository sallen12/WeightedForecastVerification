#' Conditional PIT values for parametric distributions
#'
#' Calculate conditional probability integral transform (PIT) values for probabilistic
#' forecasts in the form of familiar parametric distributions.
#'
#' @param y vector of observations.
#' @param mean,sd,location,scale,rate,shape,shape1,shape2,df1,df2,df,ncp,meanlog,sdlog,nmeans,nranges,min,max vector of parameters corresponding
#'  to the parametric predictive distributions.
#' @param a,b numeric; lower and upper threshold defining the interval of outcomes that are of interest.
#'
#' @details
#'
#' Probability integral transform (PIT) histograms are a well-established diagnostic
#' tool with which to assess the calibration of probabilistic forecasts. Conditional
#' PIT histograms (cPIT) allow forecast calibration to be assessed whilst focusing on a particular range
#' of outcomes [\code{a}, \code{b}]. This is achieved by restricting attention to the observed values \code{y} that are between
#' \code{a} and \code{b}, and calculating the PIT values corresponding to the conditional
#' distribution given that the observation is in this range (see references for details).
#'
#' The arguments \code{a} and \code{b} are numerical values corresponding to the
#' lower and upper bound of the interval of interest. If we wish to assess forecast calibration
#' when predicting outcomes that exceed a threshold \code{t}, then we can set \code{a = t}
#' and \code{b = Inf}.
#'
#' A NA is returned for the entries of \code{y} that are not in the specified range,
#' otherwise these functions output the cPIT value associated with \code{y} and the
#' corresponding forecast distribution.
#'
#' The cPIT values can be obtained for all continuous parametric distributions
#' that are available in base \proglang{R}: beta, Cauchy, chi-square, exponential,
#' F, gamma, logistic, log-normal, normal, Student's t, Tukey, uniform, and Weibull.
#' The parameters used in these functions are the same as in the standard R
#' distribution handles, e.g. \code{pnorm}.
#'
#'
#' @return vector of conditional PIT values.
#'
#' @author Sam Allen
#'
#' @references
#' \emph{PIT histograms}
#'
#' Dawid, A. P. (1984): `Present position and potential developments: Some personal views: Statistical theory: The prequential approach'. \emph{Journal of the Royal Statistical Society: Series A (General)} 147, 278-290. \doi{10.2307/2981683}
#'
#' \emph{Conditional PIT histograms}
#'
#' Allen, S., Bhend, J., Martius, O. and Ziegel, J (2023): `Weighted verification tools to evaluate univariate and multivariate probabilistic forecasts for high-impact weather events'. \emph{Weather and Forecasting} 38, 499–516. \doi{10.1175/WAF-D-22-0161.1}
#'
#'
#' @examples
#' mu <- rnorm(20, mean = 0, sd = 5)
#' y <- rnorm(20, mean = mu, sd = 1)
#'
#' cpit_norm(y = y, mean = mu, sd = 0.5)
#' pnorm(q = y, mean = mu, sd = 0.5) # pnorm() returns the same values if a = -Inf and b = Inf
#' cpit_norm(y = y, mean = mu, sd = 0.5, a = 0) # restrict attention to values that exceed 0
#' cpit_norm(y = y, mean = mu, sd = 0.5, a = -1, b = 1) # restrict attention to values between -1 and 1
#'
#' mu <- rnorm(10000, mean = 0, sd = 5)
#' y <- rnorm(10000, mean = mu, sd = 1)
#'
#' pit_hist(cpit_norm(y = y, mean = mu, sd = 1), pitvals = TRUE) # PIT hist
#'
#' pit_hist(cpit_norm(y = y, mean = mu, sd = 1, a = 0), pitvals = TRUE) # cPIT hist
#' pit_hist(cpit_norm(y = y, mean = mu + 1, sd = 1, a = 0), pitvals = TRUE) # positive bias
#' pit_hist(cpit_norm(y = y, mean = mu - 1, sd = 1, a = 0), pitvals = TRUE) # negative bias
#' pit_hist(cpit_t(y = y, df = 1, a = 0), pitvals = TRUE)
#'
#' @name cpit_param

#' @rdname cpit_param
#' @export
cpit_beta <- function(y, shape1 = 1, shape2 = 1, a = -Inf, b = Inf){
  pit <- (pbeta(y, shape1, shape2) - pbeta(a, shape1, shape2))/(pbeta(b, shape1, shape2) - pbeta(a, shape1, shape2))
  pit[y <= a | y >= b] <- NA
  return(pit)
}

#' @rdname cpit_param
#' @export
cpit_cauchy <- function(y, location = 0, scale = 1, a = -Inf, b = Inf){
  pit <- (pcauchy(y, location, scale) - pcauchy(a, location, scale))/(pcauchy(b, location, scale) - pcauchy(a, location, scale))
  pit[y <= a | y >= b] <- NA
  return(pit)
}

#' @rdname cpit_param
#' @export
cpit_chisq <- function(y, df, ncp = 0, a = -Inf, b = Inf){
  pit <- (pchisq(y, df, ncp) - pchisq(a, df, ncp))/(pchisq(b, df, ncp) - pchisq(a, df, ncp))
  pit[y <= a | y >= b] <- NA
  return(pit)
}

#' @rdname cpit_param
#' @export
cpit_exp <- function(y, rate = 1, a = -Inf, b = Inf){
  pit <- (pexp(y, rate) - pexp(a, rate))/(pexp(b, rate) - pexp(a, rate))
  pit[y <= a | y >= b] <- NA
  return(pit)
}

#' @rdname cpit_param
#' @export
cpit_f <- function(y, df1, df2, ncp = 0, a = -Inf, b = Inf){
  pit <- (pf(y, df1, df2, ncp) - pf(a, df1, df2, ncp))/(pf(b, df1, df2, ncp) - pf(a, df1, df2, ncp))
  pit[y <= a | y >= b] <- NA
  return(pit)
}

#' @rdname cpit_param
#' @export
cpit_gamma <- function(y, shape, rate = 1, scale = 1/rate, a = -Inf, b = Inf){
  pit <- (pgamma(y, shape, rate, scale) - pgamma(a, shape, rate, scale))/(pgamma(b, shape, rate, scale) - pgamma(a, shape, rate, scale))
  pit[y <= a | y >= b] <- NA
  return(pit)
}

#' @rdname cpit_param
#' @export
cpit_logis <- function(y, location = 0, scale = 1, a = -Inf, b = Inf){
  pit <- (plogis(y, location, scale) - plogis(a, location, scale))/(plogis(b, location, scale) - plogis(a, location, scale))
  pit[y <= a | y >= b] <- NA
  return(pit)
}

#' @rdname cpit_param
#' @export
cpit_lnorm <- function(y, meanlog = 0, sdlog = 1, ncp = 0, a = -Inf, b = Inf){
  pit <- (plnorm(y, meanlog, sdlog) - plnorm(a, meanlog, sdlog))/(plnorm(b, meanlog, sdlog) - plnorm(a, meanlog, sdlog))
  pit[y <= a | y >= b] <- NA
  return(pit)
}

#' @rdname cpit_param
#' @export
cpit_norm <- function(y, mean = 0, sd = 1, a = -Inf, b = Inf){
  pit <- (pnorm(y, mean, sd) - pnorm(a, mean, sd))/(pnorm(b, mean, sd) - pnorm(a, mean, sd))
  pit[y <= a | y >= b] <- NA
  return(pit)
}

#' @rdname cpit_param
#' @export
cpit_t <- function(y, df, ncp = 0, a = -Inf, b = Inf){
  pit <- (pt(y, df, ncp) - pt(a, df, ncp))/(pt(b, df, ncp) - pt(a, df, ncp))
  pit[y <= a | y >= b] <- NA
  return(pit)
}

#' @rdname cpit_param
#' @export
cpit_tukey <- function(y, nmeans, df, nranges = 1, a = -Inf, b = Inf){
  pit <- (ptukey(y, nmeans, df, nranges) - ptukey(a, nmeans, df, nranges))/(ptukey(b, nmeans, df, nranges) - ptukey(a, nmeans, df, nranges))
  pit[y <= a | y >= b] <- NA
  return(pit)
}

#' @rdname cpit_param
#' @export
cpit_unif <- function(y, min, max, a = -Inf, b = Inf){
  pit <- (punif(y, min, max) - punif(a, min, max))/(punif(b, min, max) - punif(a, min, max))
  pit[y <= a | y >= b] <- NA
  return(pit)
}

#' @rdname cpit_param
#' @export
cpit_weibull <- function(y, shape, scale = 1, a = -Inf, b = Inf){
  pit <- (pweibull(y, shape, scale) - pweibull(a, shape, scale))/(pweibull(b, shape, scale) - pweibull(a, shape, scale))
  pit[y <= a | y >= b] <- NA
  return(pit)
}


