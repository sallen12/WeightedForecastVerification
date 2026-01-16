#' Threshold calibration reliability diagrams
#'
#' Diagnostic plots to assess threshold calibration. These are reliability diagrams
#' when predicting whether multiple thresholds will be exceeded.
#'
#' @param x matrix of probability forecasts, with columns corresponding to different thresholds.
#' @param y vector of observations. Length should be the same as the number of rows in \code{x}.
#' @param t_vec vector of thresholds.
#' @param xlab,y_lab optional labels of x- and y-axes.
#' @param pointSize,textSize,spaceLegend optional arguments for the plot legend.
#' @param title optional title.
#'
#' @details
#'
#' Probabilistic forecasts are threshold calibrated if they issue calibrated forecasts
#' for the outcome exceeding all possible thresholds. That is, a (random) predictive
#' distribution function $F$ is threshold calibrated if (almost surely)
#' \deqn{P(Y \leq y \mid F(y)) = F(y)}
#' for all real values \eqn{y}. See references below.
#'
#' This can be assessed by choosing multiple thresholds \eqn{t} and plotting standard
#' reliability diagrams to check the conditional calibration of \eqn{F(t)}. These
#' reliability diagrams are then overlaid on one plot, called a threshold calibration
#' plot. If the forecasts are threshold calibrated, then all curves should be close
#' to the diagonal line.
#'
#' \code{t_vec} contains the thresholds \eqn{t} at which calibration will be assessed.
#'
#' \code{x} contains the probability forecasts \eqn{F(t)} at each of these thresholds.
#' \code{x} should therefore be a matrix with the number of columns equal to the length
#' of the vector \code{t_vec}.
#'
#' \code{y} contains the corresponding observations, from which the conditional probabilities
#' \eqn{P(Y \leq y \mid F(y)) = F(y)} can be calculated. The length of \code{y} should be
#' equal to the number of rows of \code{x}.
#'
#' A plot title can be added using \code{title}, and x- and y-axis labels can be
#' specified using \code{xlab} and \code{ylab} respectively. \code{pointSize},
#' \code{textSize}, and \code{spaceLegend} can be used to change the size and spacing
#' of the legend.
#'
#'
#' @return ggplot object containing the threshold calibration plot
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
#' n <- 1000
#' mu <- rnorm(n)
#' y <- rnorm(n, mu)
#'
#' t_vec <- -2:2
#'
#' # threshold calibrated
#' x <- sapply(t_vec, pnorm, mean = mu)
#' tc_reldiag(x, y, t_vec, title = "A calibrated example")
#'
#' # not threshold calibrated
#' x <- sapply(t_vec, pnorm, mean = mu, sd = 2)
#' tc_reldiag(x, y, t_vec, title = "A miscalibrated example")
#'
#'
#' @export
tc_reldiag <- function(x, y, t_vec, title = NULL, xlab = "x", ylab = "x_rc", pointSize = NULL, textSize = NULL, spaceLegend = NULL){

  y <- sapply(t_vec, function(t) as.numeric(y <= t))

  na_ind <- is.na(x)
  x <- matrix(x[!na_ind], ncol = length(t_vec))
  y <- matrix(y[!na_ind], ncol = length(t_vec))

  x_rc <- sapply(seq_along(t_vec), function(i) isoreg(x[, i], y[, i])$yf) # values correspond to ORDERED forecast values!
  x <- apply(x, 2, sort)

  df <- data.frame(x = as.vector(x), x_rc = as.vector(x_rc), t = rep(round(t_vec, 2), each = nrow(x)))
  plt <- ggplot2::ggplot(df) + ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = 1), lty = "dotted") +
    ggplot2::geom_line(aes(x = x, y = x_rc, col = as.factor(t))) +
    ggplot2::scale_x_continuous(name = xlab, limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(name = ylab, limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                  legend.title = ggplot2::element_blank(),
                  legend.justification = c(0, 1),
                  legend.position = c(0.01, 0.99),
                  plot.margin = ggplot2::margin(5.5, 10.5, 5.5, 5.5)) +
    ggplot2::ggtitle(title)

  if (!is.null(pointSize) && !is.null(textSize) && !is.null(spaceLegend)) {
    plt <- plt + ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = pointSize)),
                                color = ggplot2::guide_legend(override.aes = list(size = pointSize))) +
      ggplot2::theme(legend.title = ggplot2::element_text(size = textSize),
                    legend.text  = ggplot2::element_text(size = textSize),
                    legend.key.size = ggplot2::unit(spaceLegend, "lines"))
  }

  return(plt)
}
