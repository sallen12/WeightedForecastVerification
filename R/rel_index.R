#' function to calculate reliability index from ranks
#'
#' @param ranks vector of ranks.
#'
#' @return vector of the reliability index values.
#' @export
#'
#' @examples
#' rel_index(ranks = sample(1:10, 100, replace = TRUE))
rel_index <- function(ranks, n_bins = NULL){

  if(is.null(n_bins)) n_bins <- max(ranks)

  rank_freq <- sapply(1:n_bins, function(i) mean(ranks == i, na.rm = T))
  ri <- sum(abs(rank_freq - (1/n_bins)))
  return(ri)
}
