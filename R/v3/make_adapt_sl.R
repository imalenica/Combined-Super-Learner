################################################################################
# adaptive (individual + historical) super learner
################################################################################

# individual_training_data: data.table with outcome observed
# indiviual_forecast_data: data.table with outcome unobserved and covariates observed
# outcome: name of outcome variable which is found in both data.tables 
# covariates: vector of covariate names which are columns in both data.tables  
# id: column name indicating unique subject identifier 
# historical_fit: trained model on historical data
# individual_stack: optional stack of sl3 learners to model individual data  
# past_individual_fit: optional past sl3 fit to update with new data
# acf: autocorrelation in the residuals
# fit_individualSL: returns SL for just individual sl

make_adapt_sl <- function(individual_training_data, 
                          indiviual_forecast_data,
                          outcome, 
                          covariates, 
                          time = "min", 
                          first_window = 30, 
                          batch = 5, 
                          gap = 0,
                          validation_size = batch, 
                          fold_fun = folds_rolling_origin,
                          individual_stack, 
                          historical_fit,
                          past_individual_fit = NULL,
                          fit_individualSL = TRUE, 
                          acf = FALSE,
                          fit_onlineSL = NULL,
                          names_onlineSL = NULL,
                          return_historical_fit_object = FALSE){
  
  ################################ check arguments #############################
  
  # need one of these two args
  if(is.null(past_individual_fit) & is.null(individual_stack)){
    stop("Individual fit and individual stack both missing. Provide ",
         "individual_stack to initialize a new fit, and past_individual_fit ",
         "to update with new data.")
  }

  ################################## train SLs #################################
  
  individual_training_data <- data.table::data.table(individual_training_data)
  
  # make folds with training data
  folds <- origami::make_folds(individual_training_data, 
                               fold_fun = fold_fun,
                               first_window = first_window, 
                               validation_size = validation_size, 
                               gap = gap,
                               batch = batch
                               )

  individual_training_task <- sl3::make_sl3_Task(
    data = individual_training_data, 
    covariates = covariates,
    outcome = outcome, 
    folds = folds,
    time = time
    )
  
  # issue with mismatch between historical and individual delta column:
  if(is.list(historical_fit)){
    historical_task <- historical_fit$cv_fit$fit_object$full_fit$training_task
  } else {
    historical_task <- historical_fit$training_task
  }
  ind_training_task <- process_task(individual_training_task, historical_task)
  
  # fit super learner if individual_stack is provided
  if(!is.null(individual_stack)){
    # is individual_stack already a cv stack?
    if(grepl("Lrnr_cv", class(individual_stack)[1])){ 
      cv_stack <- individual_stack 
    } else {
      cv_stack <- Lrnr_cv$new(individual_stack)
    }
    ind_fit <- cv_stack$train(ind_training_task)
  } else {
    # update individualized learners
    ind_fit <- past_individual_fit$update(ind_training_task)
  }
  
  # predict with individualized learners and historical learners
  ind_preds <- ind_fit$predict()
  
  if(is.list(historical_fit)){
    # omit the task returned by make_historical_sl when fit_sl = TRUE
    historical_fits <- historical_fit[-(which(names(historical_fit) == "task"))]
    hist_preds <- dplyr::bind_rows(lapply(folds, function(fold){
      test_set_in_training_task <- validation_task(ind_training_task, fold)
      cbind.data.frame(lapply(historical_fits, function(fit){
        fit$predict_fold(test_set_in_training_task, "full")
      }))
    }))
  } else {
    hist_preds <- dplyr::bind_rows(lapply(folds, function(fold){
      test_set_in_training_task <- validation_task(ind_training_task, fold)
      historical_fit$predict_fold(test_set_in_training_task, "full")
      }))
  }

  # combine predictions
  learners <- c(paste0("historical_", colnames(hist_preds)),
                paste0("individual_", colnames(ind_preds)))
  training_preds <- cbind.data.frame(hist_preds, ind_preds)
  colnames(training_preds) <- learners

  # get true Y
  truth <- lapply(folds, function(fold){
    data.frame(individual_training_data[fold$validation, get(outcome)])
    })
  truth <- dplyr::bind_rows(truth)[,1]
  
  # evaluate empirical loss for training
  loss <- apply(training_preds, 2, function(pred) mean((pred-truth)^2))
  # establish various super learners
  sl_weights <- sl_weights_fit(training_preds, truth, loss)
  # order as original training prediction columns
  sl_weights <- sl_weights[, colnames(training_preds)]
  # use sl weights for prediction with sl
  sl_pred <- sl_weights_predict(sl_weights, training_preds)
  
  if(fit_individualSL){
    # evaluate empirical loss and weights just for individual learners
    ind_loss <- apply(ind_preds, 2, function(pred) mean((pred-truth)^2))
    sl_ind_weights <- sl_weights_fit(ind_preds, truth, ind_loss)
    sl_ind_weights <- sl_ind_weights[, colnames(ind_preds)] 
    sl_ind_pred <- sl_weights_predict(sl_ind_weights, ind_preds)
  } else {
    ind_loss <- NA
    sl_ind_weights <- NA
    sl_ind_pred <- NA
  }
  
  if(acf){
    # get the residuals -- check dependence (option to use get_acf function)
    residuals <- apply(sl_pred, 2, function(pred) (pred-truth))
    lags <- apply(residuals, 2, get_acf)
  } else {
    lags <- NA
  }
  
  ############################ forecast with SLs ###############################
  
  indiviual_forecast_data <- data.table(indiviual_forecast_data)
  
  # avoid issues with NA as outcome when outcome hasn't yet been observed
  unknown_outcome <- rep(0, nrow(indiviual_forecast_data))
  forecast_data <- data.table(unknown_outcome, indiviual_forecast_data)
  
  # obtain forecast with individual fit and historical fit
  forecast_fold <- origami::make_folds(forecast_data, fold_fun=folds_vfold, V=1)
  forecast_task <- make_sl3_Task(
    data = data.table(forecast_data), 
    covariates = covariates,
    outcome = "unknown_outcome", 
    folds = forecast_fold,
    outcome_type = individual_training_task$outcome_type$type,
    time = time
  )
  forecast_task <- process_task(forecast_task, historical_task)
  forecast_task <- process_task(forecast_task, ind_training_task)
  
  ind_forecast <- ind_fit$predict(forecast_task)
  
  if(is.list(historical_fit)){
    # omit the task returned by make_historical_sl when fit_sl = TRUE
    historical_fits <- historical_fit[-(which(names(historical_fit) == "task"))]
    hist_forecast <- cbind.data.frame(lapply(historical_fits, function(fit){
      fit$predict_fold(forecast_task, "full")
    }))
  } else {
    hist_forecast <- historical_fit$predict_fold(forecast_task, "full")
  }
  
  # combine forecasts
  learners <- c(paste0("historical_", colnames(hist_forecast)),
                paste0("individual_", colnames(ind_forecast)))
  forecasts <- cbind.data.frame(hist_forecast, ind_forecast)
  names(forecasts) <- learners
  ordered_forecasts <- forecasts[, colnames(sl_weights)]
  
  # use sl weights for forecast with sl
  sl_forecasts <- sl_weights_predict(sl_weights, ordered_forecasts)

  # use ind_sl weights for forecast with ind_sl
  if(fit_individualSL){
    ordered_ind_forecast <- ind_forecast[, colnames(sl_ind_weights), with=F]
    sl_ind_forecasts <- sl_weights_predict(sl_ind_weights, ordered_ind_forecast)
  } else {
    sl_ind_forecasts <- NA
  }
  
  # prep forecast tables for variations of online SL 
  all_forecasts <- data.table(cbind(sl_forecasts, ordered_forecasts))
  if(fit_individualSL){
    colnames(sl_ind_forecasts) <- paste0("indSL_", colnames(sl_ind_forecasts))
    all_forecasts <- data.table(cbind(sl_ind_forecasts, all_forecasts))
  }
  
  # get forecasts from the current online SL fit
  if(!is.null(fit_onlineSL) && !is.null(names_onlineSL)){
    onlineSL_forecasts <- get_onlineSL_forecasts(
      fit_onlineSL, names_onlineSL, forecast_task, hist_forecast
    )
  } else {
    onlineSL_forecasts <- NA
  }
  
  if(!return_historical_fit_object){
    historical_fit <- NA
  }
  
  return(list(
    # sl3 fit objects
    individual_fit = ind_fit,
    historical_fit = historical_fit,
    # objects used for creating adaptive/combined SL:
    training_preds = training_preds,
    sl_truth = truth, 
    loss = loss, 
    sl_weights = sl_weights, 
    sl_pred = sl_pred,
    # objects used for training individualized SL:
    fit_individualSL = fit_individualSL,
    ind_loss = ind_loss,
    sl_ind_weights = sl_ind_weights,
    sl_ind_pred = sl_ind_pred,
    # optional acf info:
    lags = lags,
    # forecast-related items (not used for training):
    indiviual_forecast_data = indiviual_forecast_data,
    forecast_task = forecast_task,
    sl_ind_forecasts = sl_ind_forecasts,
    sl_forecasts = sl_forecasts,
    lrnr_forecasts = ordered_forecasts,
    all_forecasts = all_forecasts,
    onlineSL_forecasts = onlineSL_forecasts
  ))
}

get_onlineSL_forecasts <- function(fit_onlineSL, names_onlineSL, 
                                   forecast_task, historical_forecasts){
  past_fit_weights <- fit_onlineSL$sl_weights
  past_fit_ind_fit <- fit_onlineSL$individual_fit
  past_fit_ind_forecast <- past_fit_ind_fit$predict(forecast_task)
  past_fit_learners <- c(paste0("historical_", colnames(historical_forecasts)),
                         paste0("individual_", colnames(past_fit_ind_forecast)))
  past_fit_forecasts <- cbind.data.frame(historical_forecasts, 
                                         past_fit_ind_forecast)
  names(past_fit_forecasts) <- past_fit_learners
  past_fit_ordered_forecasts <- past_fit_forecasts[, colnames(past_fit_weights)]
  past_fit_sl_forecasts <- sl_weights_predict(past_fit_weights,
                                              past_fit_ordered_forecasts)
  past_fit_all_forecasts <- data.table(cbind(past_fit_sl_forecasts, 
                                             past_fit_ordered_forecasts))
  if(fit_onlineSL$fit_individualSL){
    past_fit_ind_weights <- fit_onlineSL$sl_ind_weights
    past_fit_ord <- past_fit_ind_forecast[, colnames(past_fit_ind_weights), with=F]
    past_fit_sl_ind_forecasts <- sl_weights_predict(past_fit_ind_weights, 
                                                    past_fit_ord)
    colnames(past_fit_sl_ind_forecasts) <- paste0(
      "indSL_", colnames(past_fit_sl_ind_forecasts)
    )
    past_fit_all_forecasts <- data.table(cbind(past_fit_sl_ind_forecasts, 
                                               past_fit_all_forecasts))
  }
  onlineSL_forecasts <- lapply(names_onlineSL, function(x){
    past_fit_all_forecasts[[as.character(x)]]
  })
  names(onlineSL_forecasts) <- names(names_onlineSL)
  onlineSL_forecasts <- c(list(forecast_time = forecast_task$time), 
                          onlineSL_forecasts)
  return(onlineSL_forecasts)
}
