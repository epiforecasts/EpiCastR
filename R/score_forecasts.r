#' Score forecasts
#'
#' @importFrom DescTools BrierScore
score_forecasts= function(risk_values, diff_cases, timehorizon, day_of_forecast) {

  wherecaseswere = 1*(rowSums(diff_cases[,day_of_forecast:(day_of_forecast+timehorizon)]) > 0)
  score = BrierScore(risk_values, wherecaseswere)

  score
}



log_prob_score = function(risk_values, diff_cases, timehorizon, day_of_forecast) {


  wherecaseswere = 1*(rowSums(diff_cases[,day_of_forecast:(day_of_forecast+timehorizon)]) > 0)

  log_probs = wherecaseswere * log(abs(risk_values-1e-10)) + ((1 - wherecaseswere ) * log(1 - (abs(risk_values-1e-10))))
  sum_log_probs = sum(log_probs)

  sum_log_probs

}


find_nulls = function(diff_cases, timehorizon, day_of_forecast) {

  prev_28 = 1*(rowSums(diff_cases[,(day_of_forecast - 28):(day_of_forecast)]) > 0)

  score_briar = score_forecasts(prev_28, diff_cases, timehorizon, day_of_forecast)
  score_log = log_prob_score(prev_28, diff_cases, timehorizon, day_of_forecast)

  c(day_of_forecast, score_briar, score_log)

  }
