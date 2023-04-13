#' weighted CRPS for a normal distribution
#'
#' @param y vector of observations.
#' @param mean,sd,location,scale,df,ncp vector of parameters corresponding
#'  to the parametric predictive distributions.
#' @param a,b numeric; lower and upper threshold defining the interval of outcomes that are of interest.
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
#' the option to employ custom weight functions, as in \pkg{scoringRules}.
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
#' Allen, S., Ginsbourger, D. and J. Ziegel (2022):
#' `Evaluating forecasts for high-impact events using transformed kernel scores',
#' \emph{arXiv preprint} arXiv:2202.12732.
#' \doi{10.48550/arXiv.2202.12732}
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
vrcrps_norm <- function(y, mean = 0, sd = 1, a = -Inf){
  # does not use b
  t <- a
  y_star <- pmax(y, t)
  w_y <- as.numeric(y > t)

  Ft <- pnorm(t, mean = mean, sd = sd)
  Fys <- pnorm(y_star, mean = mean, sd = sd)
  F0 <- pnorm(0, mean = mean, sd = sd)
  Ft_sq2 <- pnorm(t, mean = mean, sd = sd/sqrt(2))
  dt <- dnorm(t, mean = mean, sd = sd)
  dys <- dnorm(y_star, mean = mean, sd = sd)
  d0 <- dnorm(0, mean = mean, sd = sd)

  o1 <- (y - mean)*(2*Fys - 1 - Ft)
  o2 <- (sd^2)*(2*dys - dt)
  o3 <- sd*(1 - Ft_sq2)/sqrt(pi)
  o4 <- (sd^2)*dt*(1 - Ft)

  obj1 <- w_y*(o1 + o2) - o3 + o4
  if(t <= 0){
    obj2_1 <- mean*(1 - Ft - 2*F0) + (sd^2)*(2*d0 - dt)
  }else{
    obj2_1 <- mean*(1 - Ft) + (sd^2)*dt
  }

  obj2 <- (obj2_1 - abs(y)*w_y)*((1 - Ft) - w_y)

  obj <- obj1 + obj2

  return(obj)
}
