#' conditional PIT values for a normal distribution
#'
#' @param y vector of observations.
#' @param mean vector of location parameters.
#' @param sd vector of scale parameters.
#' @param a,b lower/upper threshold used in the weight function.
#'
#' @return vector of conditional PIT values.
#' @export
#'
#' @examples
#' cpit_norm(y = rnorm(10))
#' cpit_norm(y = rnorm(10), a = 0.5)
#' cpit_norm(y = runif(10, 0, 5), mean = 2, sd = 0.5, a = 1.5, b = 2.5)
cpit_norm <- function(y, mean = 0, sd = 1, a = -Inf, b = Inf){
  pit <- (pnorm(y, mean, sd) - pnorm(a, mean, sd))/(pnorm(b, mean, sd) - pnorm(a, mean, sd))
  pit[y <= a | y >= b] <- NA
  return(pit)
}
