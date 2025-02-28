% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_weights.R
\name{set_weights}
\alias{set_weights}
\title{set_weights}
\usage{
set_weights(day, days, mu, sig, FUN = pweibull_with_mean_sd)
}
\arguments{
\item{timeseries}{a matrix of timeseries of cases for each geographical region [region x timestep]}

\item{shapes}{sf or sp object with polygon geometry of regions and population size}

\item{timestep}{timestep to use in fit (relative to the "timeseries"" timestep)}

\item{period_and_lag}{vector of 2 values in timesteps: first the period of influencial cases, second the lag between influential cases and fitted data}

\item{identifier}{name of the column to sort regions by - must cooporate with timeseries order}

\item{popid}{population column in spatial data}

\item{interaction}{vector with numbers to indicate the types of interaction to include in the model: 1. Gravity model. 2. Gravity model with population density. 3. Power law (no population info). 4. Adjacency model.}

\item{distrib}{poisson = 0, negative binomial = 1}
}
\description{
What is happening here
}
