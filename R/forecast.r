
fix_nan <- function(x){
  x[is.nan(x)] <- 0
  x
}


fix_inf <- function(x){
  x[is.infinite(x)] <- 0
  x
}

#' Simulate cases
#' @param start a matrix of timeseries of cases for each geographical region [region x timestep]
#' @param days sf or sp object with polygon geometry of regions and population size
#' @param case_mat timestep to use in fit (relative to the "timeseries"" timestep)
#' @param cases_from vector of 2 values in timesteps: first the period of influencial cases, second the lag between influential cases and fitted data
#' @param gamma name of the column to sort regions by - must cooporate with timeseries order
#' @param alpha_adj population column in spatial data
#' @param alpha_spat vector with numbers to indicate the types of interaction to include in the model: 1. Gravity model. 2. Gravity model with population density. 3. Power law (no population info). 4. Adjacency model.
#' @param k  poisson = 0, negative binomial = 1
#' @param beta path to stan template to build model
#' @param D variational bayes 'vb', Hamiltonian MC 'nuts'
#' @param Dprime number of MC chains
#' @param distrib number of cores to use
#' @param R number of iterations per chain
#' @param dist_mat number of iterations in warmup phase
#' @param popmat tbc
#' @param adjmat tbc
#' @importFrom stats rpois rnbinom rexp
#' @export
#'
simulate_cases <-  function(start=0, days=30, case_mat = NULL, cases_from = NULL, gamma = 0.14, alpha_adj = 0.03, alpha_spat = 0.03, alpha_con=0.0,
                            k=2.5, beta=1., D=5, Dprime=7, distrib=0, R = 0, dist_mat=0, popmat=0, adjmat=0, con_mat=0)

{
  W_ij = popmat / (dist_mat^k)                                                                # calculate and normalise gravity kernal matrix
  diag(W_ij) = 0.

  W_ij = t(W_ij/rowSums(W_ij))



  W_ij = fix_nan(W_ij )

  #
  # initiate simulation model from 0: random selection of HZs, 1: Man(first hz infected), 2: Reported case data (case_mat fed to function)
  if (start == 1){
    initial_vec = rep(0, R)
    initial_vec[89] = 1
    case_mat = t(matrix(initial_vec, nrow=R, ncol=D+Dprime))
  }

  else if (start == 0) {
    initial_vec = floor(rexp(R,rate=50) * 10 )
    case_mat = t(matrix(initial_vec, nrow=R, ncol=D+Dprime))
  }

  else if (start == 2) {

    case_mat = case_mat

  }

#  if (is.null(cases_from) ){
#      rolledsums = t(rollapply(case_mat,Dprime, sum))
#
#      summed_cases_offsett = cbind(matrix(0, nrow=R, ncol=D+Dprime), rolledsums)
#      summed_cases_offsett = summed_cases_offsett[,1:Ti]
#
#      cases_from = summed_cases_offsett
#  }
#

  #case_mat = t(case_mat)
  if (distrib == 0) {
    sampler = function (x) {rpois(n=1, lambda=x )[1]}                                          # generates a number of cases based on poisson distribution with rate x
  }

  if (distrib == 1) {
    sampler = function (x) {rnbinom(n=1, size=1/beta, mu=x)[1]}                                          # generates a number of cases based on negative binomial distribution with rate x
  }



  list(case_mat, W_ij)

  if (D != 1000) {

  for (i in 1:days) {                                                                         # run for "days" timesteps
    dvec = colSums(case_mat[(nrow(case_mat)-(D+Dprime)):(nrow(case_mat)-Dprime),])                              # sum over appropriate days

    case_mat = rbind(case_mat, sapply(gamma * dvec + alpha_spat * W_ij %*% dvec + alpha_adj * adjmat %*% dvec + alpha_con * con_mat %*% dvec , sampler)) # add timestep of cases sampled at rate to case_mat
  }

  }

  if (D == 1000) {

    for (i in 1:days) {         # run for "days" timesteps
      d = dim(case_mat)[2]
      weights = set_weights(d, 1:d, 8.5,2.6)
      dvec = colSums(weights * case_mat)                              # sum over appropriate days
      case_mat = rbind(case_mat, sapply(gamma * dvec + alpha_spat * W_ij %*% dvec + alpha_adj * adjmat %*% dvec + alpha_con * con_mat %*% dvec , sampler)) # add timestep of cases sampled at rate to case_mat
    }


  }

  case_mat                                                                                    # return case_mat
}


#' Pump posteriors
#' @param fit rstan fit object
#' @param data datalist - used in stan model
#' @param iters number of iterations per sample for simulated cases
#' @param time_horizons vector of timehorizons to return forecasts for
#' @importFrom rstan extract
#' @importFrom rlist list.append
#' @importFrom utils tail
pump_posteriors_multi <- function(fit, data, iters=1, time_horizons=c(7,14,28), D=5, Dprime=7, close_down=FALSE, buffer=0) {

  fitmat = rstan::extract(fit)

  if (!is.null(fitmat$gamma)){
  gammas = fitmat$gamma  }   else {
  gammas = rep(0,length(fitmat[[1]]))
    }# extract the samples for gamma

  if (!is.null(fitmat$alpha_spat)){
    alpha_spats = fitmat$alpha_spat  }    else {
      alpha_spats = rep(0,length(fitmat[[1]]))  }                      # extract the samples for alpha_spat

  if (!is.null(fitmat$k)){
    ks = fitmat$k  } else {
      ks = rep(0,length(fitmat[[1]]))
                  }                                      # extract the samples for gamma

  if (!is.null(fitmat$alpha_adj)){
    alpha_adjs = fitmat$alpha_adj  } else {
      alpha_adjs = rep(0,length(fitmat[[1]])) }

  if (!is.null(fitmat$alpha_con)){
    alpha_cons = fitmat$alpha_con  } else {
      alpha_cons = rep(0,length(fitmat[[1]])) }

  if (!is.null(fitmat$beta)){
    betas = fitmat$beta  } else {
      betas = rep(0,length(fitmat[[1]])) }

  R = data$R
  dist_mat = data$MIJ
  popmat=data$popmat
  adjmat = data$adjmat
  con_mat = data$con_mat

  case_mat_int = data$N

  cases_from = data$Nsum

  all_outs = list()
  for (time_horizon in time_horizons){
    all_outs = list.append(all_outs, as.matrix(rowSums(data$N)))
  }




  for (i in seq(iters)){

    for (n in 1:length(fitmat$gamma)){
      #for (n in 1:10){
      gamma = gammas[n]                             # set parameter values for realisation
      alpha_spat = alpha_spats[n]
      alpha_adj = alpha_adjs[n]
      alpha_con = alpha_cons[n]

      if (close_down == TRUE){
        alpha_spat = 0.0
        alpha_adj = 0.0
      }

      k = ks[n]
      beta = betas[n]
      # Simulate forecast data VVV
      case_mat_forcast = simulate_cases(start=2, days=max(time_horizons), case_mat = t(data$N), cases_from = cases_from,
                                        gamma = gamma, alpha_spat = alpha_spat, alpha_adj= alpha_adj, alpha_con= alpha_con, k = k, beta =beta, R = R,
                                        dist_mat = dist_mat, popmat=popmat, adjmat = adjmat, con_mat=con_mat, D=D, Dprime=Dprime)


      # Return binary descriptor of the presence of cases VVV

      for (i in 1:length(time_horizons)) {
        predicted_cases = colSums(case_mat_forcast[(dim(t(data$N))[1] + buffer) :(dim(t(data$N))[1]+time_horizons[i]),])

        # Add to the output object
        all_outs[[i]] = t(rbind(t(all_outs[[i]]), as.vector(predicted_cases)))

      }
    }

  }

  for (i in 1:length(time_horizons)){


    all_outs[[i]] = t(tail(t(all_outs[[i]]), -2))


  }


  all_outs


}


