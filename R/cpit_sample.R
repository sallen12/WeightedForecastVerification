#' conditional PIT values for simulated forecast distributions
#'
#' Using kernel density estimation
#'
#' @param y vector of observations.
#' @param dat matrix of simulation draws from forecast distribution.
#' @param a,b lower/upper threshold used in the weight function.
#' @param bw vector of bandwidth parameters to use in kernel density estimation.
#'
#' @return vector of conditional PIT values.
#' @export
#'
#' @examples
#' cpit_sample(y = rnorm(10), dat = matrix(rnorm(50), nrow = 10))
#' cpit_sample(y = rnorm(10), dat = matrix(rnorm(50), nrow = 10), a = 0.5)
#' cpit_sample(y = runif(10, 0, 5), dat = matrix(rnorm(50, 3), nrow = 10), a = 0, b = 2.5)
cpit_sample <- function(y, dat, a = -Inf, b = Inf, bw = NULL){
  if (is.null(bw)) bw <- apply(dat, 1, bw.nrd)
  pit <- sapply(seq_along(y), function(i) p_tr_kde(y[i], m = dat[i, ], bw[i], a, b))
  pit[y <= a | y >= b] <- NA
  return(pit)
}

# distribution function for normal kernel density estimation
p_tr_kde <- function(q, m_vec, s_vec, a = -Inf, b = Inf){
  p <- (mean(pnorm(q, m_vec, s_vec)) - mean(pnorm(a, m_vec, s_vec)))/
    (mean(pnorm(b, m_vec, s_vec)) - mean(pnorm(a, m_vec, s_vec)))
  return(p)
}
