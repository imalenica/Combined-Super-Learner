### Functions to support simulations

run_posl <- function(data, covs, outcome, learners, splits){
  
  data_full <- data$full
  data_historical <- data$historical
  data_individual <- data$individual
  
  #Save all results
  osl_fit                 <- list()
  result_list             <- list()
  
  result_weights          <- list()
  result_ind_weights      <- list()
  
  result_forecasts        <- list()
  result_ind_forecasts    <- list()
  result_Vfold_forecasts   <- list()
  result_hist_forecasts   <- list()
  result_online_forecasts <- list()
  
  loss                    <- list()
  
  ### Learn the full fit with V-fold CV
  full_fit <- make_historical_sl(
    historical_data = data_full,
    outcome = outcome, 
    covariates = covs, 
    id = "id", 
    historical_stack = learners
  )
  
  ### Learn the historical fit
  historical_fit <- make_historical_sl(
    historical_data = data_historical,
    outcome = outcome, 
    covariates = covs, 
    id = "id", 
    historical_stack = learners
  )
  
  split_data <- lapply(splits, function(x) data.table(data_individual[1:x,]))
  split_datas <- lapply(splits, function(x){
    data_full %>%
      dplyr::group_by(id) %>%
      dplyr::top_n(x)
  })
    
  for(i in 1:length(splits)){
    
    forecast_data <- data_individual[c((nrow(split_data[[i]])+1):(nrow(split_data[[i]])+5)),]
    
    #Train Online SL:
    osl_fit[[i]] <- make_online_sl(
      online_data=split_datas[[i]],
      outcome=outcome,
      covariates=covs,
      id = "id",
      online_stack=learners
      )
    
    #Train Personalized Online SL:
    result_list[[i]] <- make_adapt_sl(
      individual_training_data = split_data[[i]], 
      indiviual_forecast_data = forecast_data,
      outcome = outcome, 
      covariates = covs, 
      id = "id",
      historical_fit = historical_fit,
      individual_stack = learners,
      return_individual_fit=TRUE
    )
    
    ##################################################
    ###  Evaluate all performance on forecast task
    ##################################################
    forecast_task <- result_list[[i]]$forecast_task
    
    #Combined model fit on forecast data:
    result_weights[[i]] <- t(result_list[[i]]$sl_weights)
    result_forecasts[[i]] <- t(result_list[[i]]$sl_forecasts)
    
    #Individual model fit on forecast data:
    result_ind_weights[[i]] <- t(result_list[[i]]$sl_ind_weights)
    result_ind_forecasts[[i]] <- t(result_list[[i]]$sl_ind_forecasts)
    
    #Historical model fit on forecast data:
    result_hist_forecasts[[i]] <- t(cbind.data.frame(discreteSL=historical_fit$sl_discrete$predict(forecast_task),
                                                     nnls_convexSL=historical_fit$sl_nnls_convex$predict(forecast_task),
                                                     nnls_SL=historical_fit$sl_nnls$predict(forecast_task)))
    #Full model fit on forecast data (V-fold):
    result_Vfold_forecasts[[i]] <- t(cbind.data.frame(discreteSL=full_fit$sl_discrete$predict(forecast_task),
                                                     nnls_convexSL=full_fit$sl_nnls_convex$predict(forecast_task),
                                                     nnls_SL=full_fit$sl_nnls$predict(forecast_task)))
    #Online SL fit on forecast data:
    result_online_forecasts[[i]] <- t(cbind.data.frame(discreteSL=osl_fit[[i]]$sl_discrete$predict(forecast_task),
                                                       nnls_convexSL=osl_fit[[i]]$sl_nnls_convex$predict(forecast_task),
                                                       nnls_SL=osl_fit[[i]]$sl_nnls$predict(forecast_task)))
    
    #Evaluate overall loss on the forecasts:
    truth <- forecast_data$Y
    mse   <- function(pred) {mean((pred-truth)^2)}
    loss[[i]] <- cbind.data.frame(VfoldSL = apply(result_Vfold_forecasts[[i]], 1, mse),
                                  histSL = apply(result_hist_forecasts[[i]], 1, mse),
                                  indSL  = apply(result_ind_forecasts[[i]], 1, mse),
                                  onlSL  = apply(result_online_forecasts[[i]], 1, mse),
                                  persSL = apply(result_forecasts[[i]], 1, mse))
    }
  
  #Save only results for POSL
  res_discrete <- lapply(result_weights, function(x) x[,1])
  res_discrete <- do.call(cbind, res_discrete)
  colnames(res_discrete) <- paste0("t=", splits)
  
  res_nnls_convex <- lapply(result_weights, function(x) x[,2])
  res_nnls_convex <- data.frame(do.call(cbind, res_nnls_convex))
  colnames(res_nnls_convex) <- paste0("t=", splits)
  
  res_nnls <- lapply(result_weights, function(x) x[,3])
  res_nnls <- do.call(cbind, res_nnls)
  colnames(res_nnls) <- paste0("t=", splits)
  
  return(list(res_discrete=res_discrete,
              res_nnls_convex=res_nnls_convex,
              res_nnls=res_nnls,
              loss=loss))
 
}