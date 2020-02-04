
#' make_forecast
#'
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
#' @export
#'

make_forecast <- function(timeseries_mat, shapes, identifier='ADM2_NAME', day_of_forecast=NULL, do_score_forecast=TRUE, fit_over=1000, timestep =1 ,  period_and_lag = c(5,7) ,interaction = c(1), distrib=1, fit_meth='vb', chains=1, iter=100, warmup=50, cores=1) {


  if (is.null(day_of_forecast)){
    day_of_forecast = dim(timeseries_mat)[2]
    print(day_of_forecast)
  }


  fit_over = min(fit_over, day_of_forecast)




  diff_cases_slice = timeseries_mat[,(day_of_forecast-fit_over):day_of_forecast]

  FitModel = fit_model(diff_cases_slice, shapes, distrib=distrib, interaction = interaction, fit_meth = fit_meth, period_and_lag = period_and_lag, chains=chains, iter=iter, warmup=warmup, cores=cores, identifier=identifier)

  print("model fitted")
  ds_ordered = FitModel$ordered_shapes


  ForecastCases = pump_posteriors_multi(FitModel$fit, FitModel$data, D=period_and_lag[1], Dprime=period_and_lag[2])

  print("forecasts made")

  outs_7 = ForecastCases[[1]]
  outs_14 = ForecastCases[[2]]
  outs_28 = ForecastCases[[3]]

  ds_ordered$risk_7 = round(rowSums(outs_7 > 0)/dim(outs_7)[2],3)
  ds_ordered$risk_14 = round(rowSums(outs_14 > 0)/dim(outs_14)[2],3)
  ds_ordered$risk_28 = round(rowSums(outs_28 > 0)/dim(outs_28)[2],3)

  ds_ordered$risk_7_2 = round(rowSums(outs_7 >= 2)/dim(outs_7)[2],3)
  ds_ordered$risk_14_2 = round(rowSums(outs_14 >= 2)/dim(outs_14)[2],3)
  ds_ordered$risk_28_2 = round(rowSums(outs_28 >= 2)/dim(outs_28)[2],3)

  ds_ordered$risk_7_5 = round(rowSums(outs_7 >= 5)/dim(outs_7)[2],3)
  ds_ordered$risk_14_5 = round(rowSums(outs_14 >= 5)/dim(outs_14)[2],3)
  ds_ordered$risk_28_5 = round(rowSums(outs_28 >= 5)/dim(outs_28)[2],3)

  ds_ordered$risk_7_6 = round(rowSums(outs_7 >= 6)/dim(outs_7)[2],3)
  ds_ordered$risk_14_6 = round(rowSums(outs_14 >= 6)/dim(outs_14)[2],3)
  ds_ordered$risk_28_6 = round(rowSums(outs_28 >= 6)/dim(outs_28)[2],3)

  ds_ordered$risk_7_10 = round(rowSums(outs_7 >= 10)/dim(outs_7)[2],3)
  ds_ordered$risk_14_10 = round(rowSums(outs_14 >= 10)/dim(outs_14)[2],3)
  ds_ordered$risk_28_10 = round(rowSums(outs_28 >= 10)/dim(outs_28)[2],3)

  ds_ordered$risk_7_20 = round(rowSums(outs_7 >= 20)/dim(outs_7)[2],3)
  ds_ordered$risk_14_20 = round(rowSums(outs_14 >= 20)/dim(outs_14)[2],3)
  ds_ordered$risk_28_20 = round(rowSums(outs_28 >= 20)/dim(outs_28)[2],3)


  scores = list()


  if (do_score_forecast == TRUE)  {

  prev_28 = 1. * (rowSums(timeseries_mat[,(day_of_forecast-28):day_of_forecast]) > 0)



  score_7 = 1000.
  score_14 = 1000.
  score_28 = 1000.
  score_7_log = 1000.
  score_14_log = 1000.
  score_28_log = 1000.

  if (dim(timeseries_mat)[2] - day_of_forecast  < 7)

  {print('no evaluation possible - not enough data')}

  if (dim(timeseries_mat)[2] - day_of_forecast  > 7)
  {print('evaluating 7 day forecast')

    scores$score_7 = score_forecasts(ds_ordered$risk_7, timeseries_mat, 7, day_of_forecast)
    scores$score_7_log = log_prob_score(ds_ordered$risk_7, timeseries_mat, 7, day_of_forecast)
    scores$null_7 = score_forecasts(prev_28, timeseries_mat, 7, day_of_forecast)
    scores$null_7_log = log_prob_score(prev_28, timeseries_mat, 7, day_of_forecast)

  }



  if (dim(timeseries_mat)[2] - day_of_forecast  > 14)
  {print('evaluating 14 day forecast')


    scores$score_14 = score_forecasts(ds_ordered$risk_14, timeseries_mat, 14, day_of_forecast)
    scores$score_14_log = log_prob_score(ds_ordered$risk_14, timeseries_mat, 14, day_of_forecast)
    scores$null_14 = score_forecasts(prev_28, timeseries_mat, 14, day_of_forecast)
    scores$null_14_log = log_prob_score(prev_28, timeseries_mat, 14, day_of_forecast)



  }

  if (dim(timeseries_mat)[2] - day_of_forecast  > 28)
  {print('evaluating 28 day forecast')

    scores$score_28 = score_forecasts(ds_ordered$risk_28, timeseries_mat, 28, day_of_forecast)
    scores$score_28_log = log_prob_score(ds_ordered$risk_28, timeseries_mat, 28, day_of_forecast)
    scores$null_28 = score_forecasts(prev_28, timeseries_mat, 28, day_of_forecast)
    scores$null_28_log = log_prob_score(prev_28, timeseries_mat, 28, day_of_forecast)
  }

  }


  casestodate = rowSums(timeseries_mat[,1:(day_of_forecast)])
  ds_ordered$casestodate = casestodate



  list(risks = ds_ordered, scores = scores, fit=FitModel$fit)


}
