#' Fit stan model
#'
#' @importFrom rstan stan_model vb sampling
fit_model <- function(timeseries, shapes, timestep=1, period_and_lag=c(5,7), identifier="ADM2_NAME",
                      popid='totpop2019', interaction=c(1), distrib=0, model_path=system.file("extdata/template.stan",package = "EpiCastR"), fit_meth = 'vb') {

  params = prepare_stan_inputs(timeseries, shapes, timestep, period_and_lag, identifier,
                                 popid, interaction, distrib)

  datalist = params[[1]]
  ordered_shapes = params[[2]]

  if (distrib == 1){

    interactions = append(interaction, c(5))

  } else {

    interactions = interaction

    }

  model = construct_stan_model(model_path, interactions)

  write(model, "../R_HA/model_running.stan")

  pars = c()

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

  pars = unique(pars)


  sm = stan_model(model_code = model)

  if (fit_meth == 'nuts'){
    fit1 <- sampling(sm, data = datalist, chains = chains, iter = iter, warmup = warmup, cores = cores,pars =pars)    # fit stan model
  }
  if (fit_meth == 'vb'){
    fit1 <- vb(sm, data = datalist, pars = pars)    # fit stan model
  }

  list(fit = fit1, data = datalist, ordered_shapes)

}


