# TODO: ensure individual task nodes = historical task nodes, 
#       update loss and weights w/o recalculating,
#       inference, nested super learners, new covariates, 
#       correspondence, update coefficients / add offset with new batch,
#       online learners -- offset as previous fit, 
#       linreg updates second derivative, rf random node update
#       loss of historical fit as a prior -- model historical fit learner 
#       performance as a function of the covariates, sl fit of historical loss, 
#       stratified metalearner to pair current patient with historical matches 
#       based on historical pattern of patient, weights heavily people with similar 
#       patterns.

################################################################################
# adaptive (individual + historical) individualized online super learner
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
# return_individual_fit: returns SL for just individual sl

make_adapt_sl <- function(individual_training_data, indiviual_forecast_data, 
                          outcome, covariates, historical_fit, id, time = "min",
                          burn_in = 30, batch = NULL, past_individual_fit = NULL, 
                          individual_stack = NULL, return_individual_fit = TRUE,
                          acf = FALSE, parallelize = FALSE, ncores = NULL) {
  
  # check arguments
  if( is.null(past_individual_fit) & is.null(individual_stack) ) {
    stop("Individual fit and individual stack both missing. Provide ",
         "individual_stack to initialize a new fit, and past_individual_fit ",
         "to update with new data.")
  }
  
  if(is.null(batch)){
    batch <- nrow(indiviual_forecast_data)
  }
  
  if(parallelize & is.null(ncores)){
    print("Cannot parallelize when no `ncores` provided")
    parallelize <- FALSE
  }
  
  ################################## train SLs #################################
  
  # make folds with training data
  folds <- origami::make_folds(data.table(individual_training_data), 
                               fold_fun = folds_rolling_origin,
                               first_window = burn_in, 
                               validation_size = batch, 
                               gap = 0,
                               batch = batch
                               )

  individual_training_task <- make_sl3_Task(
    data = data.table(individual_training_data), 
    covariates = covariates,
    outcome = outcome, 
    folds = folds,
    time = time
    )
  
  # issue with mismatch between historical and individual delta column:
  if(class(historical_fit) == "list"){
    historical_task <- historical_fit$cv_fit$fit_object$full_fit$training_task
  } else {
    historical_task <- historical_fit$fit_object$full_fit$training_task
  }
  training_task <- process_task(individual_training_task, historical_task)
  
  # fit super learner if individual_stack is provided
  if(!is.null(individual_stack)){
    
    if(grepl("Lrnr_cv", class(individual_stack)[1])){ # already a cv stack
      cv_stack <- individual_stack 
    } else {
      cv_stack <- Lrnr_cv$new(individual_stack)
    }
    
    if(parallelize){
      plan(multicore, workers = ncores)
      ind_delayed <- delayed_learner_train(cv_stack, training_task)
      sched <- Scheduler$new(ind_delayed, FutureJob, nworkers = ncores)
      ind_fit <- sched$compute()
      } else {
      ind_fit <- cv_stack$train(training_task)
      }
    
    } else {
      # update individualized learners
      ind_fit <- past_individual_fit$update(training_task)
    }
  
  # predict with individualized learners and historical learners
  ind_preds <- ind_fit$predict(training_task)
  
  if(class(historical_fit) == "list"){
    # -1 to omit the task returned by make_historical_sl when fit_sl = TRUE
    historical_fits <- historical_fit[1:(length(historical_fit)-1)]
    hist_preds <- bind_rows(lapply(folds, function(fold){
      test_set_in_training_task <- validation_task(training_task, fold)
      cbind.data.frame(lapply(historical_fits, function(fit){
        fit$predict_fold(test_set_in_training_task, "full")
      }))
    }))
  } else {
    hist_preds <- bind_rows(lapply(folds, function(fold){
      test_set_in_training_task <- validation_task(training_task, fold)
      historical_fit$predict_fold(test_set_in_training_task, "full")
      }))
  }

  # combine predictions
  learners <- c(paste0("historical_", colnames(hist_preds)),
                paste0("individual_", colnames(ind_preds)))
  training_preds <- cbind.data.frame(hist_preds, ind_preds)
  names(training_preds) <- learners

  # get true Y
  truth <- lapply(folds, function(v){
    data.frame(individual_training_data[v$validation, get(outcome)])
    })
  truth <- bind_rows(truth)[,1]
  
  # evaluate empirical loss for training
  loss <- apply(training_preds, 2, function(pred) mean((pred-truth)^2))
  # establish various super learners
  sl_weights <- sl_weights_fit(training_preds, truth, loss)
  # use sl weights for prediction with sl
  sl_pred <- sl_weights_predict(sl_weights, training_preds)
  
  if(return_individual_fit){
    # evaluate empirical loss and weights just for individual learners
    ind_loss <- apply(ind_preds, 2, function(pred) mean((pred-truth)^2))
    sl_ind_weights <- sl_weights_fit(ind_preds, truth, ind_loss)
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
  
  # avoid issues with NA as outcome when outcome hasn't yet been observed
  unknown_outcome <- rep(0, nrow(indiviual_forecast_data))
  forecast_data <- data.table(unknown_outcome, indiviual_forecast_data)
  
  # obtain forecast with individual fit and historical fit
  forecast_fold <- make_folds(forecast_data, fold_fun=folds_vfold, V=1)
  forecast_task <- make_sl3_Task(
    data = data.table(forecast_data), 
    covariates = covariates,
    outcome = "unknown_outcome", 
    folds = forecast_fold
  )

  ind_forecast <- ind_fit$predict(forecast_task)
  
  if(class(historical_fit) == "list"){
    # -1 to omit the task returned by make_historical_sl when fit_sl = TRUE
    historical_fits <- historical_fit[1:(length(historical_fit)-1)]
    hist_forecast <- cbind.data.frame(lapply(historical_fits, function(fit){
      fit$predict(forecast_task)
    }))
  } else {
    hist_forecast <- historical_fit$predict(forecast_task, "full")
  }
  
  # combine forecasts
  learners <- c(paste0("historical_", colnames(hist_forecast)),
                paste0("individual_", colnames(ind_forecast)))
  forecast <- cbind.data.frame(hist_forecast, ind_forecast)
  names(forecast) <- learners
  
  # use sl weights for forecast with sl
  sl_forecasts <- sl_weights_predict(sl_weights, forecast)

  # use ind_sl weights for forecast with ind_sl
  if(return_individual_fit){
    sl_ind_forecasts <- sl_weights_predict(sl_ind_weights, ind_forecast)
  } else {
    sl_ind_forecasts <- NA
  }
  
  return(list(
    # return objects used for training adaptive/combined SL:
    individual_fit = ind_fit,
    historical_fit = historical_fit,
    training_preds = training_preds, 
    sl_truth = truth, 
    loss = loss, 
    sl_weights = sl_weights, 
    sl_pred = sl_pred,
    # return objects used for training individualized SL:
    ind_loss = ind_loss,
    sl_ind_weights = sl_ind_weights,
    sl_ind_pred = sl_ind_pred,
    # return optional acf info:
    lags = lags,
    # return forecast-related items (not used for training):
    forecast_task = forecast_task,
    sl_ind_forecasts = sl_ind_forecasts,
    sl_forecasts = sl_forecasts
    ))
}