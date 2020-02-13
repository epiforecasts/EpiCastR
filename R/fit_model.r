#' Fit stan model
#'
#' @description What does this function do?
#' @param timeseries a matrix of timeseries of cases for each geographical region [region x timestep]
#' @param shapes sf or sp object with polygon geometry of regions and population size
#' @param timestep timestep to use in fit (relative to the "timeseries"" timestep)
#' @param period_and_lag vector of 2 values in timesteps: first the period of influencial cases, second the lag between influential cases and fitted data
#' @param identifier name of the column to sort regions by - must cooporate with timeseries order
#' @param popid population column in spatial data
#' @param interaction vector with numbers to indicate the types of interaction to include in the model: 1. Gravity model. 2. Gravity model with population density. 3. Power law (no population info). 4. Adjacency model.
#' @param distrib  poisson = 0, negative binomial = 1
#' @param base_model_path path to stan template to build model. If not supplied the default package template is used.
#' @param final_model_path path to final stan model. If not supplied a temporary directory is used.
#' @param fit_meth variational bayes 'vb', Hamiltonian MC 'nuts'
#' @param chains number of MC chains
#' @param cores number of cores to use
#' @param iter number of iterations per chain
#' @param warmup number of iterations in warmup phase
#'
#'
#' @importFrom rstan stan_model vb sampling
#' @export




fit_model <- function(timeseries, shapes, timestep = 1, period_and_lag = c(5,7),
                      identifier = "ADM2_NAME", popid = 'totpop2019', interaction = c(1), distrib = 0,
                      base_model_path = NULL, final_model_path = NULL,
                      fit_meth = 'vb', chains = 1, iter=100, warmup=50, cores = 1, con_mat = 0) {

  ## Set default model path to be within package
  if (is.null(base_model_path)) {
    base_model_path <- system.file("extdata/template.stan",package = "EpiCastR")
  }

  ## Set default final model path to be a tmp directory
  if (is.null(final_model_path)) {
    final_model_path <- tempdir()
  }


  ## Prepare stan inputs based on options
  params = prepare_stan_inputs(timeseries, shapes, timestep, period_and_lag, identifier,
                                 popid, interaction, distrib, con_mat = con_mat)

  ## Extract data list and shapefile
  datalist = params[[1]]
  ordered_shapes = params[[2]]


  # Add interactions? - this would be much better as a character string of the options rather than a very cryptic number
  if (distrib == 1){

    interactions = append(interaction, c(5))

  } else {

    interactions = interaction

    }

  ## Construct the stan model based on the template and specified interactions
  model = construct_stan_model(base_model_path, interactions)

  ## Write the model to specific directiory - tmp by default
  final_model_path <- file.path(final_model_path, "model_running.stan")
  message("Saving the model to ", final_model_path)
  write(model, final_model_path)


  ## Add additional interaction parameters
  ## Numeric structure makes this hard to understand
  ## Is this possible to vectorise?
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

  if (6 %in% interactions){
    pars = append(pars, c( "gamma", "alpha_con"))
  }

  pars = unique(pars)

  message("Writing model")
  sm = rstan::stan_model(model_code = model)

  message("Fitting model")
  ## Fit the stan model using either MCMC or variational bayes
  if (fit_meth == 'nuts'){
    fit1 <- rstan::sampling(sm, data = datalist, chains = chains,
                            iter = iter, warmup = warmup, cores = cores,
                            pars = pars)
  }
  if (fit_meth == 'vb'){
    fit1 <- rstan::vb(sm, data = datalist, pars = pars)
  }

  # Why is ordered_shapes not named?
  return(list(fit = fit1, data = datalist, ordered_shapes))
}


