#' noisy_data
#'
#' @description
#' Simulated data designed to resemble 2m temperature forecasts issued by
#' MeteoSwiss's COSMO-E ensemble prediction system.
#'
#' @format
#' An object of type list containing 8 elements:
#' \describe{
#'   \item{clim_mean}{Matrix of dimension (3, 2000) containing the mean of a
#'   climatological forecast distribution at 3 consecutive lead times and for
#'   2000 forecast cases.}
#'   \item{clim_sd}{Matrix of dimension (3, 2000) containing the standard deviation of a
#'   climatological forecast distribution at 3 consecutive lead times and for
#'   2000 forecast cases.}
#'   \item{ens_clim}{Array of dimension (3, 2000, 10) containing the 10 ensemble
#'   members comprising the climatological forecast distribution for each lead time
#'   and forecast case.}
#'   \item{ens_pp}{Array of dimension (3, 2000, 10) containing the 10 ensemble
#'   members comprising the statistically post-processed forecast distribution for each
#'   lead time and forecast case.}
#'   \item{ens_raw}{Array of dimension (3, 2000, 10) containing the 10 ensemble
#'   members comprising the COSMO-E forecast distribution for each lead time
#'   and forecast case.}
#'   \item{obs_dat}{Matrix of dimension (3, 2000) containing the temperature observations
#'   for each lead time and forecast case.}
#'   \item{pp_mean}{Matrix of dimension (3, 2000) containing the mean of a
#'   statistically post-processed forecast distribution at 3 consecutive lead times and for
#'   2000 forecast cases.}
#'   \item{pp_sd}{Matrix of dimension (3, 2000) containing the standard deviation of a
#'   statistically post-processed forecast distribution at 3 consecutive lead times and for
#'   2000 forecast cases.}
#' }
#'
#' @details
#' Forecasts are available for three probabilistic prediction systems: the COSMO-E
#' ensemble prediction system, a climatological forecast distribution, and a statistically
#' post-processed forecast distribution. The climatological and statistical forecast
#' distributions are normal distributions with means \code{clim_mean} and \code{pp_mean},
#' respectively, and standard deviations \code{clim_sd} and \code{pp_sd}. The corresponding
#' ensemble forecasts \code{ens_clim} and \code{ens_pp} are 10 evenly-spaced quantiles
#' from these normal distributions.
#'
#' The three forecast methods can be compared to the observation matrix \code{obs_dat},
#' which provides temperature observations for each lead time and forecast case.
#'
#' @usage data("noisy_data")
#'
"noisy_data"
