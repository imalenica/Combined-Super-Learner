make_online_sl <- function(adapt_fit, outcome, past_online_risk, 
                           past_onlineSL){
  forecast_data <- adapt_fit$forecast_task$data
  truth <- forecast_data[[outcome]]
  Time <- head(forecast_data[["time_and_date"]], 1)
  Min <- head(forecast_data[["min"]], 1)
  
  all_forecasts <- data.table(cbind(adapt_fit$sl_forecasts, 
                                    adapt_fit$lrnr_forecasts))
  
  # isolate historical/individual forecasts
  lrnr_forecasts <- data.table(adapt_fit$lrnr_forecasts)
  lrnr_cols <- colnames(lrnr_forecasts)
  historical_lrnrs <- lrnr_cols[grep("historical_", lrnr_cols)]
  historical_forecasts <- lrnr_forecasts[, historical_lrnrs, with=FALSE]
  individual_lrnrs <- lrnr_cols[grep("individual_", lrnr_cols)]
  individual_forecasts <- lrnr_forecasts[, individual_lrnrs, with=FALSE]
  
  if(!is.na(adapt_fit$sl_ind_forecasts)){
    all_forecasts <- data.table(cbind(adapt_fit$sl_ind_forecasts,
                                      all_forecasts))
    individual_forecasts <- data.table(cbind(adapt_fit$sl_ind_forecasts,
                                             individual_forecasts))
  }
  forecast_data_list <- list(Individual = individual_forecasts,
                             Historical = historical_forecasts,
                             All = all_forecasts)
  
  new_onlineSL <- lapply(seq_along(forecast_data_list), function(i){
    
    # select new onlineSL
    forecasts <- forecast_data_list[[i]]
    Type <- names(forecast_data_list)[i]
    past_risk <- past_online_risk[[Type]]
    current_risk <- apply(forecasts, 2, function(pred) mean((pred-truth)^2))
    online_risk <- as.numeric(current_risk + past_risk)
    
    Name <- colnames(forecasts)[which.min(online_risk)]
    onlineSL <- c(Time, Min, Type, Name)
    
    return(list(onlineSL = onlineSL, online_risk = online_risk))
  })
  
  onlineSL <- data.table(do.call(rbind, lapply(new_onlineSL, '[[', 'onlineSL')))
  
  online_risk <- lapply(new_onlineSL, '[[', 'online_risk')
  names(online_risk) <- names(forecast_data_list)
  
  # evaluate performance of prior onlineSL
  if(evaluate_onlineSL_performance){
    onlineSL_MSE <- lapply(seq_along(forecast_data_list), function(i){
      forecasts <- data.table(forecast_data_list[[i]])
      onlineSL <- past_onlineSL[i, ]
      onlineSL_forecasts <- forecasts[[onlineSL$Name]]
      onlineSL$MSE <- mean((onlineSL_forecasts-truth)^2)
      return(onlineSL)
    })
    onlineSL_performance <- data.table(do.call(rbind, onlineSL_MSE))
  } else {
    onlineSL_performance <- NA
  }
  return(list(onlineSL = new_onlineSL, 
              online_risk = online_risk, 
              onlineSL_performance = onlineSL_performance))
}