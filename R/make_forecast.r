
#' make_forecast
#' @description What is happening here
#' @importFrom stats rpois rnbinom rexp
#' @param timeseries_mat a matrix of timeseries of cases for each geographical region [region x timestep]
#' @param shapes sf or sp object with polygon geometry of regions and population size
#' @param day_of_forecast day of outbreak to make forecast - default last day of data
#' @param fit_over number of timesteps to fit over.
#' @param timestep timestep to use in fit (relative to the "timeseries"" timestep)
#' @param period_and_lag vector of 2 values in timesteps: first the period of influencial cases, second the lag between influential cases and fitted data
#' @param interaction vector with numbers to indicate the types of interaction to include in the model: 1. Gravity model. 2. Gravity model with population density. 3. Power law (no population info). 4. Adjacency model.
#' @param distrib  poisson = 0, negative binomial = 1
#' @param fit_meth variational bayes 'vb', Hamiltonian MC 'nuts'
#' @param chains number of MC chains
#' @param cores number of cores to use
#' @param iter number of iterations per chain
#' @param warmup number of iterations in warmup phase
#' @inheritParams fit_model
#' @export
#'


make_forecast <- function(timeseries_mat, shapes, identifier = 'ADM2_NAME', day_of_forecast = NULL, do_score_forecast = TRUE,
                          fit_over = 1000, timestep = 1 , period_and_lag = c(5,7) ,interaction = c(1), distrib = 1, fit_meth='vb',
                          chains = 1, iter = 100, warmup = 50, cores = 1, timehorizons = c(7, 14, 28), thresholds = c(1, 2, 5, 10, 20),
                          close_down = FALSE, con_mat = NULL, compiled_model = NULL) {


  message(do_score_forecast)

  if (is.null(day_of_forecast)){
    day_of_forecast = dim(timeseries_mat)[2]
    messsage(day_of_forecast)
  }


  fit_over = min(fit_over, day_of_forecast)




  diff_cases_slice = timeseries_mat[,(day_of_forecast-fit_over):day_of_forecast]

  ## Change FitModel to match naming scheme
  FitModel = fit_model(diff_cases_slice, shapes, distrib=distrib, interaction = interaction, fit_meth = fit_meth,
                       period_and_lag = period_and_lag,  chains = chains, iter = iter, warmup = warmup, cores = cores,
                       identifier = identifier, con_mat = con_mat, compiled_model = compiled_model)

  message("Model fitted")
  ds_ordered = FitModel$ordered_shapes

  message("Got shapes")


  ## Also match naming scheme
  ForecastCases = pump_posteriors_multi(FitModel$fit, FitModel$data, D = period_and_lag[1], Dprime = period_and_lag[2],
                                        time_horizons = timehorizons, close_down = close_down)

  message("Forecasts made")

  for (i in 1:length(timehorizons)){
    for (h in thresholds){
      risks = round(rowSums(ForecastCases[[i]] >= h) / dim(ForecastCases[[i]])[2], 3)
      ds_ordered[[paste0("risk_", as.character(timehorizons[i]), '_', as.character(h))]] = risks
      }
  }

  scores = list()
  log_scores = list()


  if (do_score_forecast == TRUE)  {
    message('Scoring forecasts')
    for (i in 1:length(timehorizons)){
      t = timehorizons[i]
      for (h in thresholds){
        message(dim(timeseries_mat)[2] - day_of_forecast)
        if (dim(timeseries_mat)[2] - day_of_forecast  > t){

          message('Scoring forecasts')
          score = score_forecasts(ds_ordered[[paste0("risk_", as.character(t), '_', as.character(h))]],
                                   timeseries_mat, t, day_of_forecast, h)
          log_score = log_prob_score(ds_ordered[[paste0("risk_", as.character(t), '_', as.character(h))]],
                                      timeseries_mat, t, day_of_forecast, h)
          message(score)
          message(log_score)

          scores[[paste0("score_", as.character(t), '_', as.character(h))]] = score
          log_scores[[paste0("score_", as.character(t), '_', as.character(h))]] = log_score
        }
      }
    }
  }

  ## What is happening?
  casestodate = rowSums(timeseries_mat[,1:(day_of_forecast)])
  ds_ordered$casestodate = casestodate



  return(
    list(risks = ds_ordered, scores = scores, log_scores = log_scores,
         fit = FitModel$fit, case_counts = ForecastCases)
  )

}
