#' set_weights
#'
#' @importFrom EpiForecastsUtils pweibull_with_mean_sd
#' @param timeseries a matrix of timeseries of cases for each geographical region [region x timestep]
#' @param shapes sf or sp object with polygon geometry of regions and population size
#' @param timestep timestep to use in fit (relative to the "timeseries"" timestep)
#' @param period_and_lag vector of 2 values in timesteps: first the period of influencial cases, second the lag between influential cases and fitted data
#' @param identifier name of the column to sort regions by - must cooporate with timeseries order
#' @param popid population column in spatial data
#' @param interaction vector with numbers to indicate the types of interaction to include in the model: 1. Gravity model. 2. Gravity model with population density. 3. Power law (no population info). 4. Adjacency model.
#' @param distrib  poisson = 0, negative binomial = 1
#'
#' @export
set_weights = function( day, days, mu, sig, FUN=pweibull_with_mean_sd){

  times = -(days-day)
  cumweib = sapply(X=times, FUN = function(X) {FUN(X, mu, sig)})

  cumweib[1:(length(times)-1)] = cumweib[1:(length(times)-1)] - cumweib[2:length(times)]

  weights = cumweib

  weights
}
