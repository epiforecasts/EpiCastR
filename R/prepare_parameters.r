
#' Prepare stan inputs
#'
#' @importFrom sf st_area st_centroid st_distance
#' @importFrom spdep poly2nb nb2mat
#' @importFrom zoo rollapply
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
prepare_stan_inputs <- function(timeseries, shapes, timestep=1, period_and_lag=c(5,7), identifier="ADM2_NAME",
                                popid='totpop2019', interaction=c(1), distrib=0){




  R = dim(timeseries)[1]                                                   # Set number of areas R
  Ti = dim(timeseries)[2]/timestep                                         # Set nuber of days

  shapes_ordered = shapes[order(shapes[[identifier]]),]                    # set order of areas to alphabetical to ensure consistancy

  shapes_ordered['AREA'] = st_area(shapes_ordered)

  sf_cent = st_centroid(shapes_ordered)  # find centroids of the polygons

  distmat = st_distance(sf_cent)                                           # matrix of distance between centroids as 1D matrix for denom of grav model
  distmat = matrix(distmat, ncol=R, nrow=R)                                # format matrix to correct shape
  distmat[distmat < 1.] = 1.                                               # set minimum disance to 1. (avoiding devide by 0 (over written in stan for final interaction mat))







  # SET GRAVITY INTERACTION


        if (1 %in% interaction){
          popmat = outer(as.vector(shapes_ordered[[popid]]), as.vector(shapes_ordered[[popid]]))      # set population matrix for numerator of grav model
          popmat[is.na(popmat)] = 1.
        } else if (2 %in% interaction){
          popmat = outer(as.vector(shapes_ordered[[popid]] / shapes_ordered[["AREA"]]), as.vector(shapes_ordered[[popid]]/shapes_ordered[["AREA"]]))      # set population matrix for numerator of grav model
          popmat[is.na(popmat)] = 1.
        } else if (3 %in% interaction){
          popmat = matrix(1, ncol=R, nrow=R)      # set population matrix for numerator of grav model
          popmat[is.na(popmat)] = 1.
        } else{
          popmat = matrix(0., ncol=R, nrow=R)      # set population matrix for numerator of grav model
        }

  # SET ADJACENCY INTERACTION

        if (4 %in% interaction){
          nbsdrc = poly2nb(shapes_ordered)
          adjmat = nb2mat(nbsdrc, zero.policy=TRUE)
        } else {
          adjmat = matrix(0., ncol=R, nrow=R)
        }


  if (timestep != 1){
    cases_to = t(rollapply(t(timeseries), timestep, sum, by=timestep))
  } else {
    cases_to = timeseries
  }


#  set_weights = function( day, days, mu, sig, FUN=EpiForecastsUtils::pweibull_with_mean_sd){
#
#    times = -(days-day)
#    cumweib = sapply(X=times, FUN = function(X) {FUN(X, mu, sig)})
#
#    cumweib[1:(length(times)-1)] = cumweib[1:(length(times)-1)] - cumweib[2:length(times)]
#
#    weights = cumweib
#
#    weights
#  }

  if (period_and_lag[1] != 1000){

  rolledsums = t(rollapply(t(cases_to), period_and_lag[1], sum))

  summed_cases_offsett = cbind(matrix(0, nrow=R, ncol=sum(period_and_lag)), rolledsums)
  summed_cases_offsett = summed_cases_offsett[,1:Ti]

  cases_from = summed_cases_offsett
  }

  else if (period_and_lag[1] == 1000){
   cases_from = c()

   for(d in 1:Ti){
     weights = set_weights(d, 1:Ti, 5,3)
     cases_from = cbind(cases_from, colSums(weights * t(cases_to)))
   }

  }


  datalist = list(

    R = R,
    T = Ti,
    MIJ = distmat/1000.0,
    popmat = popmat,
    adjmat = adjmat,
    N = cases_to,
    Nsum = cases_from,
    distrib = distrib
           )

  list(datalist, shapes_ordered)


}
