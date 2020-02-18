
#' Prepare stan inputs
#' @description What does this do?
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
prepare_stan_inputs <- function(timeseries, shapes, timestep=1, period_and_lag = c(5,7), identifier = "ADM2_NAME",
                                popid = 'totpop2019', interaction = c(1), distrib = 0, con_mat = NULL, con_mat_times=NULL){



  ## Set number of areas R
  R = dim(timeseries)[1]
  ## Set number of days
  Ti = dim(timeseries)[2]/timestep

  ## Set order of areas to alphabetical to ensure consistancy
  shapes_ordered = shapes[order(shapes[[identifier]]),]

  shapes_ordered['AREA'] = st_area(shapes_ordered)

  ## Find centroids of the polygons
  sf_cent = st_centroid(shapes_ordered)

  ## Matrix of distance between centroids as 1D matrix for denom of grav model
  distmat = st_distance(sf_cent)
  ## Format matrix to correct shape
  distmat = matrix(distmat, ncol = R, nrow = R)
  ## Set minimum disance to 1. (avoiding devide by 0 (over written in stan for final interaction mat))
  ## Why is this 1.
  distmat[distmat < 1.] = 1.







  ## Set gravity interaction
        if (1 %in% interaction){
          ## Set population matrix for numerator of grav model
          popmat = outer(as.vector(shapes_ordered[[popid]]), as.vector(shapes_ordered[[popid]]))
          popmat[is.na(popmat)] = 1.
        } else if (2 %in% interaction){
          ## Set population matrix for numerator of grav model
          popmat = outer(as.vector(shapes_ordered[[popid]] / shapes_ordered[["AREA"]]),
                         as.vector(shapes_ordered[[popid]] / shapes_ordered[["AREA"]]))
          popmat[is.na(popmat)] = 1.
        } else if (3 %in% interaction){
          ## Set population matrix for numerator of grav model
          popmat = matrix(1, ncol=R, nrow=R)
          popmat[is.na(popmat)] = 1.
        } else{
          ## Set population matrix for numerator of grav model
          popmat = matrix(1., ncol=R, nrow=R)
        }

  ## SET ADJACENCY INTERACTION
        if (4 %in% interaction){
          nbsdrc = poly2nb(shapes_ordered)
          adjmat = nb2mat(nbsdrc, zero.policy=TRUE)
        } else {
          adjmat = matrix(0., ncol=R, nrow=R)
        }

  ## SET PRECALCULATED CONNECTIVITY MATRIX
        if (6 %in% interaction){
          if (is.list(con_mat)){
            con_mat_multi = con_mat
            con_mat = matrix(0., ncol=R, nrow=R)
          }
          else if (is.matrix(con_mat)){
            con_mat = con_mat
            con_mat_multi = NULL
          }
          else {
            print('Warning: interaction matrix must be matrix or list of matices - not using predefined interaction matrix')
            con_mat = matrix(0., ncol=R, nrow=R)
          }
        }
        else{
          con_mat = matrix(0., ncol=R, nrow=R)
        }

  ## What is happening here?
  if (timestep != 1){
    cases_to = t(rollapply(t(timeseries), timestep, sum, by=timestep))
  } else {
    cases_to = timeseries
  }

  ## What is happening here?
  ## Comments on this bit please!
  if (period_and_lag[1] != 1000){

  rolledsums = t(rollapply(t(cases_to), period_and_lag[1], sum))

  summed_cases_offsett = cbind(matrix(0, nrow=R, ncol=sum(period_and_lag)), rolledsums)
  summed_cases_offsett = summed_cases_offsett[,1:Ti]

  cases_from = summed_cases_offsett
  } else if (period_and_lag[1] == 1000){
   cases_from = c()

   for(d in 1:Ti){
     weights = set_weights(d, 1:Ti, 8.5,2.6)
     cases_from = cbind(cases_from, colSums(weights * t(cases_to)))
   }

  }

  datalist = list(
    R = R,
    T = Ti,
    MIJ = distmat/1000.0,
    popmat = popmat,
    adjmat = adjmat,
    con_mat = con_mat,
    N = cases_to,
    Nsum = cases_from,
    distrib = distrib
           )
  if (length(con_mat_times) > 1){
    timevecs = list()
    for(t in con_mat_times){
      timevecs[paste0("timevec_", s.character(i))] = rep(0, Ti)
      timevecs[paste0("timevec_", s.character(i))][con_mat_times[i]:con_mat_times[i+1]] = 1
    }
  }

  if(!is.null(con_mat_multi)){
    if(length(con_mat_multi) == lenght(con_mat_times)){
      con_mat_times = append(con_mat_times, Ti)
      for(i in 1:length(con_mat_multi)){
        datalist[paste0('con_mat_', as.character(i))] = con_mat_multi[i]
        datalist[paste0('Nsum_', as.character(i))] = t(timevecs[i] * t(cases_from))

      }
    }
    else{
      print('Incorrect no. times - using first matrix only')
      datalist['con_mat'] = con_mat_multi[1]
    }
  }

  return(list(datalist, shapes_ordered))
}
