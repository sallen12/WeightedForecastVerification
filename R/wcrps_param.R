#' Weighted CRPS
#'
#' Weighted versions of the continuous ranked probability score (CRPS) for forecast
#' distributions in the form of a normal, logistic, or Student's t distribution.
#' The thresold-weighted, outcome-weighted, and vertically re-scaled CRPS are
#' available.
#'
#' @param y vector of observations.
#' @param mean,sd,location,scale,df,ncp vector of parameters corresponding
#'  to the parametric predictive distributions.
#' @param a,b numeric; lower and upper threshold defining the interval of outcomes that are of interest.
#' @param x0 	numeric value at which to centre the vertically re-scaled CRPS.
#' @param BS logical specifying whether the outcome-weighted scores should be complemented with the Brier score.
#'
#' @details
#'
#' Weighted scoring rules allow forecasts to be assessed and compared when
#' predicting certain outcomes of interest. Weighted scoring rules are available
#' in the \pkg{scoringRules} package, but only for forecasts in the form of predictive samples.
#' Here, a collection of weighted scoring rules is provided for familiar parametric distributions.
#' Currently, this is available for the normal, logistic, and Student's t distributions,
#' though more could be added in the future.
#'
#' The weight function used in the weighted scores is \code{w(z) = 1{a < z < b}}, where
#' \code{a} and \code{b} are numerical values corresponding to the lower and upper
#' bound of the interval of interest. The default values are \code{a = -Inf} and
#' \code{b = Inf}, which returns the unweighted CRPS. Future work could also integrate
#' the option to employ custom weight functions, as in \pkg{scoringRules}, though
#' this is generally more difficult for parametric forecast distributions.
#'
#' The parameter arguments (\code{mean}, \code{sd}, \code{location}, \code{scale}, \code{df}, \code{ncp})
#' should be vectors of the same length as the observation vector \code{y}. The naming of these parameters
#' is the same as in the standard \proglang{R} distribution handles, e.g. \code{pnorm}.
#'
#' Threshold-weighted and outcome-weighted scoring rules are provided. The vertically re-scaled CRPS
#' is also provided for forecasts in the form of a normal distribution. This could be
#' added for other distributions in the future.
#'
#' @return vector of weighted score values.
#'
#' @author Sam Allen
#'
#' @references
#' \emph{Threshold-weighted CRPS:}
#'
#' Gneiting, T. and R. Ranjan (2011):
#' `Comparing density forecasts using threshold-and quantile-weighted scoring rules',
#' \emph{Journal of Business & Economic Statistics} 29, 411-422.
#' \doi{10.1198/jbes.2010.08110}
#'
#' \emph{Outcome-weighted CRPS:}
#'
#' Holzmann, H. and B. Klar (2017):
#' `Focusing on regions of interest in forecast evaluation',
#' \emph{Annals of Applied Statistics} 11, 2404-2431.
#' \doi{10.1214/17-AOAS1088}
#'
#' \emph{Vertically re-scaled CRPS:}
#'
#' Allen, S., Ginsbourger, D. and J. Ziegel (2023):
#' `Evaluating forecasts for high-impact events using transformed kernel scores',
#' \emph{SIAM/ASA Journal on Uncertainty Quantification} 11, 906-940.
#' \doi{10.1137/22M1532184}
#'
#' @examples
#' sig <- 1
#' mu <- rnorm(20, mean = 0, sd = 5)
#' y <- rnorm(20, mean = mu, sd = sig)
#'
#' scoringRules::crps_norm(y = y, mean = mu, sd = sig)
#' twcrps_norm(y = y, mean = mu, sd = sig)
#' owcrps_norm(y = y, mean = mu, sd = sig)
#' vrcrps_norm(y = y, mean = mu, sd = sig)
#'
#' # emphasise outcomes above 0
#' twcrps_logis(y = y, location = mu, scale = sig, a = 0)
#' owcrps_logis(y = y, location = mu, scale = sig, a = 0)
#'
#' # emphasise outcomes between -1 and 1
#' twcrps_t(y = y, df = 5, ncp = mu, a = -1, b = 1)
#' owcrps_t(y = y, df = 5, ncp = mu, a = -1, b = 1)
#'
#' @name wcrps_param

#' @rdname wcrps_param
#' @export
twcrps_logis <- function(y, location = 0, scale = 1, a = -Inf, b = Inf){
  s <- scoringRules::crps_clogis(y = pmin(pmax(y, a), b), location = location, scale = scale, lower = a, upper = b)
  return(s)
}

#' @rdname wcrps_param
#' @export
owcrps_logis <- function(y, location = 0, scale = 1, a = -Inf, b = Inf, BS = T){
  w_y <- as.numeric(y > a & y < b)
  obj <- scoringRules::crps_tlogis(y, location = location, scale = scale, lower = a, upper = b)*w_y

  if(BS){
    p <- plogis(b, location, scale) - plogis(a, location, scale)
    obj <- obj + (p - w_y)^2
  }

  return(obj)
}

#' @rdname wcrps_param
#' @export
twcrps_norm <- function(y, mean = 0, sd = 1, a = -Inf, b = Inf){
  s <- scoringRules::crps_cnorm(y = pmin(pmax(y, a), b), location = mean, scale = sd, lower = a, upper = b)
  return(s)
}

#' @rdname wcrps_param
#' @export
owcrps_norm <- function(y, mean = 0, sd = 1, a = -Inf, b = Inf, BS = T){
  w_y <- as.numeric(y > a & y < b)
  obj <- scoringRules::crps_tnorm(y, location = mean, scale = sd, lower = a, upper = b)*w_y

  if(BS){
    p <- pnorm(b, mean, sd) - pnorm(a, mean, sd)
    obj <- obj + (p - w_y)^2
  }

  return(obj)
}

#' @rdname wcrps_param
#' @export
twcrps_t <- function(y, df, ncp = 0, a = -Inf, b = Inf){
  s <- scoringRules::crps_ct(y = pmin(pmax(y, a), b), df = df, location = ncp, lower = a, upper = b)
  return(s)
}

#' @rdname wcrps_param
#' @export
owcrps_t <- function(y, df, ncp = 0, a = -Inf, b = Inf, BS = T){
  w_y <- as.numeric(y > a & y < b)
  obj <- scoringRules::crps_tt(y, df = df, location = ncp, lower = a, upper = b)*w_y

  if(BS){
    p <- pt(b, df, ncp) - pt(a, df, ncp)
    obj <- obj + (p - w_y)^2
  }

  return(obj)
}

#' @rdname wcrps_param
#' @export
vrcrps_norm <- function(y, mean = 0, sd = 1, a = -Inf, b = Inf, x0 = 0){

  b_min <- pmin(y, b)
  b_max <- pmax(y, b)
  a_min <- pmin(y, a)
  a_max <- pmax(y, a)

  b_min_x0 <- pmin(x0, b)
  b_max_x0 <- pmax(x0, b)
  a_min_x0 <- pmin(x0, a)
  a_max_x0 <- pmax(x0, a)

  w_y <- as.numeric(y > a & y < b)


  F_b <- pnorm(b, mean = mean, sd = sd)
  F_a <- pnorm(a, mean = mean, sd = sd)

  d_b <- dnorm(b, mean = mean, sd = sd)
  d_a <- dnorm(a, mean = mean, sd = sd)

  F_b_sq2 <- pnorm(b, mean = mean, sd = sd/sqrt(2))
  F_a_sq2 <- pnorm(a, mean = mean, sd = sd/sqrt(2))

  F_b_min <- pnorm(b_min, mean = mean, sd = sd)
  F_b_max <- pnorm(b_max, mean = mean, sd = sd)
  F_a_min <- pnorm(a_min, mean = mean, sd = sd)
  F_a_max <- pnorm(a_max, mean = mean, sd = sd)

  d_b_min <- dnorm(b_min, mean = mean, sd = sd)
  d_b_max <- dnorm(b_max, mean = mean, sd = sd)
  d_a_min <- dnorm(a_min, mean = mean, sd = sd)
  d_a_max <- dnorm(a_max, mean = mean, sd = sd)

  F_b_min_x0 <- pnorm(b_min_x0, mean = mean, sd = sd)
  F_b_max_x0 <- pnorm(b_max_x0, mean = mean, sd = sd)
  F_a_min_x0 <- pnorm(a_min_x0, mean = mean, sd = sd)
  F_a_max_x0 <- pnorm(a_max_x0, mean = mean, sd = sd)

  d_b_min_x0 <- dnorm(b_min_x0, mean = mean, sd = sd)
  d_b_max_x0 <- dnorm(b_max_x0, mean = mean, sd = sd)
  d_a_min_x0 <- dnorm(a_min_x0, mean = mean, sd = sd)
  d_a_max_x0 <- dnorm(a_max_x0, mean = mean, sd = sd)

  o1 <- (y - mean)*(F_b_min + F_a_max - F_b_max - F_a_min)
  o2 <- (sd^2)*(d_b_min + d_a_max - d_b_max - d_a_min)
  o3 <- sd*(F_b_sq2 - F_a_sq2)/sqrt(pi)
  o4 <- (sd^2)*(d_b + d_a)*(F_b - F_a)
  obj1 <- w_y*(o1 + o2) - o3 + o4

  o1 <- (x0 - mean)*(F_b_min_x0 + F_a_max_x0 - F_b_max_x0 - F_a_min_x0)
  o2 <- (sd^2)*(d_b_min_x0 + d_a_max_x0 - d_b_max_x0 - d_a_min_x0)
  obj2_1 <- o1 + o2

  obj2 <- (obj2_1 - abs(y - x0)*w_y)*((F_b - F_a) - w_y)

  obj <- obj1 + obj2

  return(obj)
}

#' @rdname wcrps_param
#' @export
vrcrps_logis <- function(y, location = 0, scale = 1, a = -Inf, b = Inf, x0 = 0){
  # needs checking

  b_min <- pmin(y, b)
  b_max <- pmax(y, b)
  a_min <- pmin(y, a)
  a_max <- pmax(y, a)

  b_min_x0 <- pmin(x0, b)
  b_max_x0 <- pmax(x0, b)
  a_min_x0 <- pmin(x0, a)
  a_max_x0 <- pmax(x0, a)

  w_y <- as.numeric(y > a & y < b)


  F_b <- plogis(b, location = location, scale = scale)
  F_a <- plogis(a, location = location, scale = scale)

  F_b_min <- plogis(b_min, location = location, scale = scale)
  F_b_max <- plogis(b_max, location = location, scale = scale)
  F_a_min <- plogis(a_min, location = location, scale = scale)
  F_a_max <- plogis(a_max, location = location, scale = scale)

  F_b_min_x0 <- plogis(b_min_x0, location = location, scale = scale)
  F_b_max_x0 <- plogis(b_max_x0, location = location, scale = scale)
  F_a_min_x0 <- plogis(a_min_x0, location = location, scale = scale)
  F_a_max_x0 <- plogis(a_max_x0, location = location, scale = scale)

  o1 <- (b_max - y)*F_b_max + (a_min - y)*F_a_min - (b_min - y)*F_b_min + (a_max - y)*F_a_max
  o2 <- (b_max + a_min - b_min - a_max)
  o3 <- scale*(log(F_b_max) + log(F_a_min) - log(F_b_min) - log(F_a_max))
  obj1 <- w_y*(o1 - o2 + o3)

  o1 <- b*F_b - b*F_b*F_a - b + b*F_a + scale*log(F_b) - scale*F_a*log(F_b) + scale*F_a*log(F_a)
  o2 <- a*F_a - a*F_a*F_b - a + a*F_b + scale*log(F_a) - scale*F_b*log(F_a) + scale*F_b*log(F_b)
  o3 <- scale*(F_b - F_a)
  obj2 <- o1 - o2 + o3

  o1 <- (b_max_x0 - x0)*F_b_max_x0 + (a_min_x0 - x0)*F_a_min_x0 - (b_min_x0 - x0)*F_b_min_x0 + (a_max_x0 - x0)*F_a_max_x0
  o2 <- (b_max_x0 + a_min_x0 - b_min_x0 - a_max_x0)
  o3 <- scale*(log(F_b_max_x0) + log(F_a_min_x0) - log(F_b_min_x0) - log(F_a_max_x0))
  obj3 <- o1 - o2 + o3

  obj3 <- (obj3 - abs(y - x0)*w_y)*((F_b - F_a) - w_y)

  obj <- obj1 - obj2 + obj3

  return(obj)
}
