#' PIT reliability diagrams
#'
#' Plot probability integral transform (PIT) reliability diagrams from a vector
#' of PIT values.
#'
#' @param z vector or list of vector of pit values.
#' @param ranks logical specifying whether the input are ranks (TRUE) or PIT values (FALSE).
#' @param resampling logical specifying whether resampling should be used to calculate consistency regions.
#' @param n_resamples number of resamples to use when calculating consistency regions.
#' @param region_level significance level of the consistency regions.
#' @param title optional title.
#'
#' @details
#'
#' This code has been adapted from https://github.com/resinj/replication_GR21.
#'
#' PIT histograms display the distribution of probability integral transform (PIT)
#' values (see \code{\link{pit_hist}} for details). While PIT histograms have become
#' well-established when assessing forecast calibration, there is generally no
#' canonical choice for the number of bins to use in the histogram, and this choice
#' can have a large impact on the histogram's interpretation, particularly when
#' the sample size is small.
#'
#' Instead, it has been argued that it is more appropriate to display the (empirical)
#' distribution function of the PIT values, rather than their (empirical) density function.
#' If the PIT values are independent samples from a standard uniform distribution,
#' as is the case for a probabilistically calibrated forecast, then this distribution
#' function should lie on the line y = x, subject to sampling variation.
#' This permits a straightforward assessment of whether or not forecasts are
#' calibrated, and deviations from the diagonal can be used to identify systematic
#' errors that occur in the forecasts. These plots of the PIT values' distribution function
#' are called PIT reliability diagrams.
#'
#' To assess whether the distribution function is significantly different from the
#' diagonal, we can calculate consistency regions around the diagonal, which would
#' contain a calibrated forecast with a certain significance level. To calculate and
#' display these regions, we can set the \code{resampling} argument of the
#' \code{pit_reldiag()} function to \code{TRUE} (the default).
#'
#' These consistency regions are obtained using parametric bootstrap resampling,
#' with the number of resamples given by the argument \code{n_resamples}, and the
#' significance level by the argument \code{region_level}.
#'
#' \code{z} is either a vector of PIT values, or a list of vectors of PIT values to
#' be plotted together. If \code{z} is a list and \code{resampling = TRUE}, resampling
#' is performed based on the number of elements in the first element of \code{z}.
#' This implicitly assumes the number of elements is the same in all elements of \code{z}.
#'
#' @return ggplot object containing the PIT reliability diagram
#'
#' @author Sam Allen
#'
#' @references
#'
#' Gneiting, T., and Resin, J. (2021):
#' `Regression diagnostics meets forecast evaluation: Conditional calibration, reliability diagrams, and coefficient of determination'.
#' \emph{arXiv preprint.}
#' \doi{arXiv:2108.03210}
#'
#' @examples
#' fc_cal <- runif(1000)
#' pit_reldiag(z = fc_cal, title = "Example 1")
#' pit_hist(z = fc_cal, ranks = FALSE, title = "Example 1")
#'
#' fc_ud <- rbeta(1000, shape1 = 0.5, shape2 = 0.5)
#' pit_reldiag(z = fc_ud, title = "Example 2")
#' pit_hist(z = fc_ud, ranks = FALSE, title = "Example 2")
#'
#' pit_reldiag(z = fc_ud, region_level = 0.5, title = "Example 2")
#' pit_reldiag(z = fc_ud, region_level = 0.99, title = "Example 2")
#'
#' fc_od <- rbeta(1000, shape1 = 2, shape2 = 2)
#' pit_reldiag(z = fc_od, resampling = FALSE, title = "Example 3")
#' pit_hist(z = fc_od, ranks = FALSE, title = "Example 3")
#'
#' z <- list(Calibrated = fc_cal, Underdispersed = fc_ud, Overdispersed = fc_od)
#' pit_reldiag(z = z, title = "Example 4")
#'
#' @export
pit_reldiag <- function(z, ranks = FALSE, resampling = TRUE, n_resamples = 1000, region_level = 0.9, title = NULL){
  if (is.list(z)) {

    if (ranks) z <- lapply(z, function(zz) (zz + runif(length(zz)) - 1)/max(zz)) # convert ranks to PIT values

    if (is.null(names(z))) names(z) <- 1:length(z)

    dist_pit <- lapply(z, ecdf)

    x <- seq(0, 1, 0.01)
    df <- data.frame(x = x,
                     Fx = as.vector(sapply(dist_pit, function(edf) edf(x))),
                     mth = rep(names(z), each = length(x)))
    plot_reldiag <- ggplot2::ggplot(df) +
      ggplot2::geom_step(ggplot2::aes(x = x, y = Fx, col = as.factor(mth))) +
      ggplot2::geom_abline(ggplot2::aes(slope = 1, intercept = 0), lty = "dotted") +
      ggplot2::scale_x_continuous(name = expression(z), limits = c(0, 1), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(name = expression("fraction of PIT values" <= z), limits = c(0, 1), expand = c(0, 0)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.title = ggplot2::element_blank(),
                     legend.justification = c(0, 1),
                     legend.position = c(0.01, 0.99),
                     panel.grid = ggplot2::element_blank()) +
      ggplot2::ggtitle(title)

  } else {

    if (ranks) z <- (z + runif(length(z)) - 1)/max(z) # convert ranks to PIT values

    dist_pit = ecdf(z)

    x <- seq(0, 1, 0.01)
    plot_reldiag <- ggplot2::ggplot(data.frame(x = x, Fx = dist_pit(x))) +
      ggplot2::geom_step(ggplot2::aes(x = x, y = Fx)) +
      ggplot2::geom_abline(ggplot2::aes(slope = 1, intercept = 0), lty = "dotted") +
      ggplot2::scale_x_continuous(name = expression(z), limits = c(0, 1), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(name = "fraction of PIT values \u2264 z", limits = c(0, 1), expand = c(0, 0)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank()) +
      ggplot2::ggtitle(title)

  }

  if (resampling) {
    low <- floor(n_resamples * (1 - region_level)/2)
    up <- n_resamples - low

    if (is.list(z)) {
      len <- lengths(z[[1]])
    } else {
      len <- length(z)
    }
    resamples <- sapply(1:n_resamples, function(i) runif(len))

    dist_resamples_x <- apply(resamples, 2, function(s) ecdf(s)(x))
    dist_resamples_x_sorted <- apply(dist_resamples_x, 1, sort)

    bounds_low <- dist_resamples_x_sorted[low, ]
    bounds_up <- dist_resamples_x_sorted[up, ]

    if (is.list(z)) {
      bounds_low <- rep(bounds_low, length(z))
      bounds_up <- rep(bounds_up, length(z))
    }

    plot_reldiag <- plot_reldiag + ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = bounds_low, ymax = bounds_up),
                                                        fill = "lightblue", alpha = 0.3)
  }

  plot_reldiag

}
