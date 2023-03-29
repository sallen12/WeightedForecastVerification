#' plot rank and probability integral transform (PIT) histograms
#'
#' @param ranks vector of ranks.
#' @param bins number of bins to use in the histogram.
#' @param pitvals logical specifying whether the input are PIT values.
#' @param title optional title of plot.
#' @param ymax optional upper limit of the y axis.
#'
#' @export
#'
#' @examples
#' pit_hist(ranks = sample(1:10, 1000, replace = TRUE), title = "Example 1", ymax = 1)
#' pit_hist(ranks = runif(1000), bins = 10, pitvals = TRUE, title = "Example 2", ymax = 0.3)
#' pit_hist(ranks = rbeta(1000, shape1 = 0.5, shape2 = 0.5), pitvals = TRUE, title = "Example 3")
#' pit_hist(ranks = rbeta(1000, shape1 = 2, shape2 = 2), bins = 20, pitvals = TRUE, title = "Example 4")
pit_hist <- function(ranks, bins = NULL, pitvals = FALSE, title = NULL, ymax = NULL){

  if(is.null(bins)){
    if(pitvals){
      bins <- 10
    }else{
      bins <- max(ranks)
    }
  }

  if(pitvals){ ranks <- floor(ranks*bins) + 1; ranks[ranks > bins] <- bins}

  rank_freq <- sapply(1:bins, function(i) mean(ranks == i, na.rm = T))
  if(is.null(ymax)) ymax <- 1.5*max(rank_freq)

  df <- data.frame(freq = rank_freq, rank = as.factor(1:bins))
  ggplot2::ggplot(df, ggplot2::aes(x = rank, y = freq)) + ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 1/bins), col = "red", lty = "dashed") +
    ggplot2::scale_x_discrete(name = "Rank") +
    ggplot2::scale_y_continuous(name = "Rel. Freq.", limits = c(0, ymax), expand = c(0, 0)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank()) +
    ggplot2::ggtitle(title)
}
