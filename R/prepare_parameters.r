
#' Prepare stan inputs
#'
#' @importFrom sf st_area st_centroid st_distance
#' @importFrom spdep poly2nb nb2mat
#' @importFrom zoo rollapply
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

  rolledsums = t(rollapply(t(cases_to), period_and_lag[1], sum))

  summed_cases_offsett = cbind(matrix(0, nrow=R, ncol=sum(period_and_lag)), rolledsums)
  summed_cases_offsett = summed_cases_offsett[,1:Ti]

  cases_from = summed_cases_offsett


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
