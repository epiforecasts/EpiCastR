#' Fit stan model
#' @param timeseries a matrix of timeseries of cases for each geographical region [region x timestep]
#' @param shapes sf or sp object with polygon geometry of regions and population size
#' @param timestep timestep to use in fit (relative to the "timeseries"" timestep)
#' @param period_and_lag vector of 2 values in timesteps: first the period of influencial cases, second the lag between influential cases and fitted data
#' @param identifier name of the column to sort regions by - must cooporate with timeseries order
#' @param popid population column in spatial data
#' @param interaction vector with numbers to indicate the types of interaction to include in the model: 1. Gravity model. 2. Gravity model with population density. 3. Power law (no population info). 4. Adjacency model.
#' @param distrib  poisson = 0, negative binomial = 1
#' @param model_path path to stan template to build model
#' @param fit_meth variational bayes 'vb', Hamiltonian MC 'nuts'
#' @param chains number of MC chains
#' @param cores number of cores to use
#' @param iter number of iterations per chain
#' @param warmup number of iterations in warmup phase
#'
#'
#' @importFrom rstan stan_model vb sampling
#' @export



fit_model <- function(timeseries, shapes, timestep=1, period_and_lag=c(5,7), identifier="ADM2_NAME",
                      popid='totpop2019', interaction=c(1), distrib=0, model_path=system.file("extdata/template.stan",package = "EpiCastR"),
                      fit_meth = 'vb', chains=1, iter=100, warmup=50, cores=1, con_mat=NULL) {

  params = prepare_stan_inputs(timeseries, shapes, timestep, period_and_lag, identifier,
                                 popid, interaction, distrib, con_mat=NULL)

  datalist = params[[1]]
  ordered_shapes = params[[2]]

  if (distrib == 1){

    interactions = append(interaction, c(5))

  } else {

    interactions = interaction

    }

  model = construct_stan_model(model_path, interactions)

  write(model, "../R_HA/model_running.stan")

  pars = c("epsilon")

  if (1 %in% interactions){
    pars = append(pars, c("k", "gamma", "alpha_spat"))
  }

  if (2 %in% interactions){
    pars = append(pars, c("k", "gamma", "alpha_spat"))
  }

  if (3 %in% interactions){
    pars = append(pars, c("k", "gamma", "alpha_spat"))
  }

  if (4 %in% interactions){
    pars = append(pars, c("gamma", "alpha_adj"))
  }

  if (5 %in% interactions){
    pars = append(pars, c("beta"))
  }

  if (6 %in% interactions){
    pars = append(pars, c("alpha_con"))
  }

  pars = unique(pars)

  print("writing model")
  sm = rstan::stan_model(model_code = model)

  print("fitting model")

  if (fit_meth == 'nuts'){
    fit1 <- rstan::sampling(sm, data = datalist, chains = chains, iter = iter, warmup = warmup, cores = cores,pars =pars)    # fit stan model
  }
  if (fit_meth == 'vb'){
    fit1 <- rstan::vb(sm, data = datalist, pars = pars)    # fit stan model
  }

  list(fit = fit1, data = datalist, ordered_shapes)

}


