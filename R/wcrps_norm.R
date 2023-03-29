#' weighted CRPS for a normal distribution
#'
#' @param y vector of observations.
#' @param mean vector of location parameters.
#' @param sd vector of scale parameters.
#' @param a,b lower/upper threshold used in the weight function.
#' @param BS logical specifying whether the outcome-weighted score should be complemented with the Brier score.
#'
#' @return vector of weighted score values.
#' @export
#'
#' @examples
#' twcrps_norm(y = rnorm(10))
#' owcrps_norm(y = rnorm(10))
#' vrcrps_norm(y = rnorm(10))
twcrps_norm <- function(y, mean = 0, sd = 1, a = -Inf, b = Inf){
  s <- scoringRules::crps_cnorm(y = pmin(pmax(y, a), b), location = mean, scale = sd, lower = a, upper = b)
  return(s)
}

#' @export
#' @rdname twcrps_norm
owcrps_norm <- function(y, mean = 0, sd = 1, a = -Inf, b = Inf, BS = T){
  w_y <- as.numeric(y > a & y < b)
  obj <- scoringRules::crps_tnorm(y, location = mean, scale = sd, lower = a, upper = b)*w_y

  if(BS){
    p <- pnorm(b, mean, sd) - pnorm(a, mean, sd)
    obj <- obj + (p - w_y)^2
  }

  return(obj)
}

#' @export
#' @rdname twcrps_norm
vrcrps_norm <- function(y, mean = 0, sd = 1, a = -Inf, b = -Inf){
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
