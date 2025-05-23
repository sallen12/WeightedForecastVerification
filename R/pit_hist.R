#' Rank and PIT histograms
#'
#' Plot rank and probability integral transform (PIT) histograms from a vector
#' of ranks or PIT values.
#'
#' @param z vector of ranks or PIT values.
#' @param bins number of bins to use in the histogram.
#' @param ranks logical specifying whether the input are ranks (TRUE) or PIT values (FALSE).
#' @param title optional title.
#' @param ymax optional upper limit of the y axis.
#' @param xlab,ylab optional x and y axes labels.
#' @param xticks,yticks option to remove x and y ticks.
#' @param linecol optional colour of horizontal line (default is red).
#' @param linetype optional type of horizontal line (default is dashed).
#'
#' @details
#'
#' Rank and probability integral transform (PIT) histograms are well-established
#' diagnostic tools with which to assess the calibration of probabilistic forecasts.
#'
#' If the forecasts are samples from a predictive distribution (i.e. ensemble forecasts),
#' then rank histograms display the relative frequency that the observation assumes
#' each rank among the sample members. If the observations and ensemble members are exchangeable
#' then the observation should be equally likely to assume any rank, meaning the ranks
#' should be uniformly distributed. The resulting rank histogram should therefore be
#' flat.
#'
#' If the forecasts are continuous predictive distributions, then probabilistic
#' calibration can be assessed using PIT values, defined as the forecast
#' distribution function evaluated at the observed outcome (with suitable
#' randomisation for forecast distributions that are not continuous). If the forecasts
#' are probabilistically calibrated, then these PIT values should
#' follow a standard uniform distribution. Plotting the distribution of
#' these PIT values should yield a flat histogram.
#'
#' In both cases, deviations from a flat histogram can be used to identify systematic
#' deficiencies in the forecasts.
#'
#' By default, \code{pit_hist()} assumes the input \code{z} is a vector of ranks (integers),
#' and thus that a rank histogram is to be plotted. Unless specified otherwise,
#' the number of bins to display in the rank histogram is the maximum
#' of the vector of ranks, though the number of bins can be specified manually
#' using the \code{bins} argument. For PIT histograms, the default is \code{bins = 10}.
#'
#' If a PIT histogram is desired, then the user should set \code{ranks = FALSE}, and
#' the vector \code{z} should contain values between zero and one.
#'
#' The rank and PIT histograms display a dashed red line at the height of
#' a uniform histogram. The colour of this line can be altered using the argument
#' \code{linecol}, while the type of this line can be altered using \code{linetype}.
#' \code{linecol} can be either a colour string (e.g. \code{"red"}) or a hex code
#' (e.g. \code{"#FF0000"}). \code{linetype} can be any valid ggplot linetype,
#' (e.g. \code{"solid"}, \code{"dashed"}, \code{"dotted"}, \code{"dotdash"}). The
#' uniformity line can be omitted by setting either of these arguments to \code{NULL}.
#'
#'
#' @return ggplot object containing the rank/PIT histogram
#'
#' @author Sam Allen
#'
#' @references
#'
#' Dawid, A. P. (1984):
#' `Present position and potential developments: Some personal views: Statistical theory: The prequential approach'.
#' \emph{Journal of the Royal Statistical Society: Series A (General)} 147, 278-290.
#' \doi{10.2307/2981683}
#'
#' Gneiting, T., Balabdaoui, F., and Raftery, A. E. (2007):
#' `Probabilistic forecasts, calibration and sharpness'.
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 69, 243-268.
#' \doi{10.1111/j.1467-9868.2007.00587.x}
#'
#' @examples
#'
#' # Rank histograms
#' fc_cal <- sample(1:10, 1000, replace = TRUE)
#' fc_ud <- sample(1:10, 1000, replace = TRUE, prob = abs(seq(0, 1, length.out=10) - 0.5))
#' fc_od <- sample(1:10, 1000, replace = TRUE, prob = 1 - abs(seq(0, 1, length.out=10) - 0.5))
#'
#' pit_hist(z = fc_cal, title = "Calibrated example", ymax = 0.3, linecol = "blue")
#' pit_hist(z = fc_ud, title = "Under-dispersed example", ymax = 0.3, linetype = "solid")
#' pit_hist(z = fc_od, title = "Over-dispersed example", ymax = 0.3, linecol = NULL)
#'
#' # PIT histograms (ranks = FALSE)
#' pit_hist(z = runif(1000), ranks = FALSE, title = "Calibrated example", ymax = 0.5, xlab = "")
#' pit_hist(z = rbeta(1000, 0.5, 0.5), ranks = FALSE, title = "Under-dispersed")
#' pit_hist(z = rbeta(1000, 2, 2), bins = 20, ranks = FALSE, title = "Over-dispersed")
#'
#' @export
pit_hist <- function(z, bins = NULL, ranks = TRUE, title = NULL, ymax = NULL,
                     ylab = "Rel. Freq.", xlab = "Rank",
                     yticks = TRUE, xticks = TRUE,
                     linecol = "red", linetype = "dashed") {

  if (is.null(bins)) {
    if(!ranks) {
      bins <- 10
    }else {
      bins <- max(z, na.rm = T)
    }
  }

  if (!ranks) z <- floor(z*bins) + 1; z[z > bins] <- bins

  rank_freq <- sapply(1:bins, function(i) mean(z == i, na.rm = T))
  if (is.null(ymax)) ymax <- 1.5*max(rank_freq)
  if (any(rank_freq > ymax)) {
    message(paste(sum(rank_freq > ymax), "bin(s) have been truncated at ymax"))
    rank_freq[rank_freq > ymax] <- ymax
  }

  alpha <- 1
  if (is.null(linecol) || is.null(linetype)) {linecol <- "red"; linetype <- "dashed"; alpha <- 0}

  df <- data.frame(freq = rank_freq, rank = as.factor(1:bins))
  out_plot <- ggplot2::ggplot(df, ggplot2::aes_string(x = "rank", y = "freq")) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 1/bins), col = linecol, alpha = alpha, lty = linetype) +
    ggplot2::scale_x_discrete(name = xlab) +
    ggplot2::scale_y_continuous(name = ylab, limits = c(0, ymax), expand = c(0, 0)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(c(5.5, 10.5, 5.5, 5.5))) +
    ggplot2::ggtitle(title)
  if (!yticks) out_plot <- out_plot + ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                                                     axis.text.y = ggplot2::element_blank())
  if (!xticks) out_plot <- out_plot + ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                                                     axis.text.x = ggplot2::element_blank())
  return(out_plot)
}
