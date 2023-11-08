#' Reliability indices
#'
#' Calculate reliability indices that measure the flatness of a rank histogram.
#'
#' @param z vector of ranks or PIT values.
#' @param bins number of bins in the histogram.
#' @param ranks logical specifying whether \code{z} is ranks (\code{TRUE}) or PIT
#'  values (\code{FALSE}).
#' @param method method to measure the flatness of the histogram.
#'
#' @details
#'
#' If a forecast is probabilistically calibrated, it will generate a flat
#' rank or probability integral transform (PIT) histogram. To quantify how calibrated
#' forecasts are, several metrics (called reliability indices) have been proposed
#' to measure the distance between the observed histogram and a flat histogram.
#'
#' If \code{bins} is the number of bins in the histogram, then the relative frequency
#' associated with each bin of a flat histogram is 1/\code{bins}.
#'
#' Three methods are available here to calculate reliability indices:
#' \code{method = 'absolute'} (default) measures the sum of absolute distances between the observed
#' relative frequencies and 1/\code{bins}; \code{method = 'squared'} measures the sum of
#' squared distances between the observed relative frequencies and 1/\code{bins}, which
#' follows a chi-squared distribution when appropriately scaled; and
#' \code{method = 'entropy'} measures the entropy of the relative frequencies, which
#' will be equal to 1 for a uniform rank histogram, and between 0 and 1 otherwise.
#'
#' Note that \code{method = 'absolute'} and \code{method = 'squared'} lead to reliability
#' indices for which a lower value indicates a more calibrated forecast, whereas
#' \code{method = 'entropy'} returns an index that is higher for calibrated forecasts.
#'
#' The argument \code{z} contains a vector of ranks or PIT values. If \code{z} contains
#' PIT values, they are grouped into \code{bins} bins, and the reliability index is then
#' calculated using these discretised values, rather than the original PIT values.
#' Of course, future extensions could include alternative methods to measure the distance
#' between the sample of PIT values and the standard uniform distribution, e.g. using
#' maximum mean discrepancies or scoring rule divergences.
#'
#' @return a reliability index.
#'
#' @references
#'
#' Delle Monache, L., Hacker, J. P., Zhou, Y., Deng, X., and Stull, R. B. (2006):
#' `Probabilistic aspects of meteorological and ozone regional ensemble forecasts'.
#' \emph{Journal of Geophysical Research: Atmospheres}, 111 D24307.
#' \doi{10.1029/2005JD006917}
#'
#'
#' Wilks, D. S. (2019):
#' `Indices of rank histogram flatness and their sampling properties'.
#' \emph{Monthly Weather Review}, 147, 763-769.
#' \doi{10.1175/MWR-D-18-0369.1}
#'
#' @examples
#'
#' # Ranks
#' fc_cal <- sample(1:10, 1000, replace = TRUE)
#' fc_ud <- sample(1:10, 1000, replace = TRUE, prob = abs(seq(0, 1, length.out=10) - 0.5))
#' fc_od <- sample(1:10, 1000, replace = TRUE, prob = 1 - abs(seq(0, 1, length.out=10) - 0.5))
#'
#' pit_hist(z = fc_cal, title = "Calibrated forecasts")
#' pit_hist(z = fc_ud, title = "Under-dispersed forecasts")
#' pit_hist(z = fc_od, title = "Over-dispersed forecasts")
#'
#' rel_index(z = fc_cal, method = "absolute")
#' rel_index(z = fc_ud, method = "absolute")
#' rel_index(z = fc_od, method = "absolute")
#'
#' rel_index(z = fc_cal, method = "entropy")
#' rel_index(z = fc_ud, method = "entropy")
#' rel_index(z = fc_od, method = "entropy")
#'
#'
#' # PIT values
#' fc_cal <- runif(1000)
#' fc_ud <- rbeta(1000, shape1 = 0.5, shape2 = 0.5)
#' fc_od <- rbeta(1000, shape1 = 2, shape2 = 2)
#'
#' pit_hist(z = fc_cal, ranks = FALSE, title = "Calibrated forecasts")
#' pit_hist(z = fc_ud, ranks = FALSE, title = "Under-dispersed forecasts")
#' pit_hist(z = fc_od, ranks = FALSE, title = "Over-dispersed forecasts")
#'
#' rel_index(z = fc_cal, ranks = FALSE)
#' rel_index(z = fc_ud, ranks = FALSE)
#' rel_index(z = fc_od, ranks = FALSE)
#'
#' @export
rel_index <- function(z, bins = NULL, ranks = TRUE, method = "absolute"){

  if(is.null(bins)) {
    if(!ranks) {
      bins <- 10
    }else {
      bins <- max(z)
    }
  }

  if(!ranks) { z <- floor(z*bins) + 1; z[z > bins] <- bins}

  rank_freq <- sapply(1:bins, function(i) mean(z == i, na.rm = T))

  if(method == "absolute") {
    ri <- sum(abs(rank_freq - (1/bins)))
  } else if(method == "squared") {
    ri <- sum((rank_freq - (1/bins))^2)
  } else if(method == "entropy") {
    ri <- -sum(log(rank_freq ** rank_freq))/log(bins)
  } else {
    stop("method not recognised. method must be one of 'absolute', 'squared',
         and 'entropy'")
  }

  return(ri)
}
