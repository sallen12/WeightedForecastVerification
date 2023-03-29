#' PIT reliability diagrams (adapted from https://github.com/resinj/replication_GR21)
#'
#' @param pit vector of pit values.
#' @param resampling logical specifying whether resampling should be used to calculate consistency regions
#' @param n_resamples number of resamples to use when calculating consistency regions.
#' @param region_level significance level of the consistency regions.
#' @param title optional title of plot.
#'
#' @export
#'
#' @examples
#' pit_reldiag(pit = runif(1000), title = "Example 1")
#' pit_reldiag(pit = rbeta(1000, shape1 = 0.5, shape2 = 0.5), title = "Example 2")
#' pit_reldiag(pit = rbeta(1000, shape1 = 0.5, shape2 = 0.5), region_level = 0.5, title = "Example 3")
#' pit_reldiag(pit = rbeta(1000, shape1 = 2, shape2 = 2), resampling = FALSE, title = "Example 4")
pit_reldiag <- function(pit, resampling = TRUE, n_resamples = 1000, region_level = 0.9, title = NULL){
  dist_pit = ecdf(pit)

  x <- seq(0, 1, 0.01)
  plot_reldiag <- ggplot2::ggplot(data.frame(x = x, Fx = dist_pit(x))) +
    ggplot2::geom_step(ggplot2::aes(x = x, y = Fx)) +
    ggplot2::geom_abline(ggplot2::aes(slope = 1, intercept = 0), lty = "dotted") +
    ggplot2::scale_x_continuous(name = expression(z), limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(name = expression("fraction of PIT values" <= z), limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank()) +
    ggplot2::ggtitle(title)

  if(resampling){
    low = floor(n_resamples * (1 - region_level)/2)
    up = n_resamples - low

    resamples = sapply(1:n_resamples, function(i) runif(length(pit)))

    dist_resamples_x = apply(resamples, 2, function(s) ecdf(s)(x))
    dist_resamples_x_sorted = apply(dist_resamples_x, 1, sort)

    bounds_low <- dist_resamples_x_sorted[low,]
    bounds_up <- dist_resamples_x_sorted[up,]

    plot_reldiag <- plot_reldiag + ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = bounds_low, ymax = bounds_up),
                                                        fill = "lightblue", alpha = 0.3)
  }

  plot_reldiag

}
