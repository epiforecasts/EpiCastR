
#' Run historical ha
#'
#' @importFrom data.table data.table
#' @importFrom utils write.table
run_historical_ha = function(days, fit_over, iters=10, warmup=5, cores=1, chains=1) {
  seven = t(rep(0,8))
  three = t(rep(0,4))
  dir_name = paste0(paste0("..//R_HA/output/",Sys.Date()), '_reports')

  risk_dir_name = paste0(dir_name, '/risks')
  score_dir_name = paste0(dir_name, '/scores')
  param_dir_name = paste0(dir_name, '/params')
  dir.create(dir_name)
  dir.create(risk_dir_name)
  dir.create(score_dir_name)
  dir.create(param_dir_name)

  #write.csv(seven, file=paste0(param_dir_name,"/kvals.csv"))
  #write.csv(seven, file=paste0(param_dir_name,"/avals.csv"))
  #write.csv(seven, file=paste0(param_dir_name,"/gvals.csv"))
  #write.csv(three, file=paste0(score_dir_name,"/scores_briar.csv"))
  #write.csv(three, file=paste0(score_dir_name,"/scores_log.csv"))

  for(day in days) {
          outs = fit_and_forecast_ha(day, fit_over=fit_over, popdens=TRUE)

          fit_summary = data.table(rstan::summary(outs[[1]]$fit1, pars = c("k", "alpha_spat", "gamma"), probs = c(0.1, 0.9))$summary)
          fit_summary$date = rep(day, 3)

          ds_ordered = data.frame(outs[[2]])

          for (ha in outs[[2]]$adm3_pc){

            risks = ds_ordered[ds_ordered$adm3_pc == ha,colnames(ds_ordered)[grepl("^risk", colnames(ds_ordered))]]
            risks$day = day

            write.table(risks, file=paste0(risk_dir_name, '/' , ha ,'_risks.csv'), sep=',', append = T, quote=F, col.names = F, row.names = F)





          }

          write.table(fit_summary[1, ], file=paste0(param_dir_name, '/kvals.csv'), sep=',', append = T, quote=F, col.names = F, row.names = T)
          write.table(fit_summary[2, ], file=paste0(param_dir_name, '/avals.csv'), sep=',', append = T, quote=F, col.names = F, row.names = T)
          write.table(fit_summary[3, ], file=paste0(param_dir_name, '/gvals.csv'), sep=',', append = T, quote=F, col.names = F, row.names = T)

          write.table(data.table(t(c(day, outs[[3]], outs[[4]], outs[[5]]))), file=paste0(score_dir_name, '/scores_briar.csv'), sep=',', append = T, quote=F, col.names = F, row.names = T)
          write.table(data.table(t(c(day, outs[[6]], outs[[7]], outs[[8]]))), file=paste0(score_dir_name, '/scores_log.csv'), sep=',', append = T, quote=F, col.names = F, row.names = T)
       }
}


historical_nulls = function(diff_cases, days, timehorizons) {

  dir_name = paste0(paste0('../R_HA/output/', Sys.Date()), '_nulls')
  dir.create(dir_name)

  for (timehorizon in timehorizons){
    print(timehorizon)
    nulls = c()
    for (day in days) {

      nulls = cbind(nulls, find_nulls(diff_cases, timehorizon, day))
    }

    filename = paste0(paste0(paste0(dir_name, '/nulls_'), as.character(timehorizon) ),'.csv')
    print(filename)
    write.table(t(nulls), file=filename, sep=',', quote=F, col.names = F)
  }
  nulls
}

