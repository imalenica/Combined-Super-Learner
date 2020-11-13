slstream <- function(training_data, outcome, covariates, id, time, fold_fun, 
                     batch, window_size = NULL, first_window = NULL, 
                     individual_stack = NULL, previous_individual_fit = NULL, 
                     historical_fit, first_fit = FALSE, loss_table = NULL, 
                     weights_control, forecast_data){
  
  ################################ check arguments #############################
  
  # need one of these two args
  if(is.null(previous_individual_fit) & is.null(individual_stack)){
    stop("Previous individual fit and individual stack both missing. Provide ",
         "individual_stack to initialize a new fit, and previous_individual_fit ",
         "to update with new data.")
  }
  
  if(!first_fit){
    if(is.null(loss_table)){
      stop("Table of all learner losses must be provided if not first fit.")
    }
  }
  
  ################################## train SLs #################################
  
  training_data <- data.table::data.table(training_data)
  
  # make folds with training data
  if(fold_fun == "folds_rolling_origin"){
    folds <- origami::make_folds(training_data, fold_fun = folds_rolling_origin,
                                 first_window = first_window, gap = 0, 
                                 validation_size = batch, batch = batch)
  } else if(fold_fun == "folds_rolling_window"){
    folds <- origami::make_folds(training_data, fold_fun = folds_rolling_window, 
                                 window_size = window_size, gap = 0,
                                 validation_size = batch, batch = batch)
  }
  
  if(length(folds) < 2){
    stop("Insufficient data: ", batch, " more time points needed.")
  }
  
  train_task <- sl3::make_sl3_Task(training_data, covariates, outcome, 
                                   time = time, folds = folds)
  
  # issue with mismatch between historical and individual delta column:
  if(is.list(historical_fit)){
    # historical_fit is a list when fit with SL, see fun make_historical_sl
    historical_task <- historical_fit$cv_fit$fit_object$full_fit$training_task
  } else {
    historical_task <- historical_fit$training_task
  }
  train_task <- process_task(train_task, historical_task)
  
  # fit stack if individual_stack is provided
  if(!is.null(individual_stack)){
    # is individual_stack already a cv stack?
    if(grepl("Lrnr_cv", class(individual_stack)[1])){ 
      cv_stack <- individual_stack 
    } else {
      cv_stack <- Lrnr_cv$new(individual_stack, full_fit = T)
    }
    ind_fit <- cv_stack$train(train_task)
  } else {
    ind_fit <- previous_individual_fit$update(train_task) # update learners
  }
  
  ########################### get fold-specific predictions ####################
  fold_fits <- lapply(folds, function(fold){
    # retain validation task fold v
    v_task <- validation_task(train_task, fold)
    
    # get historical fit predictions for fold v
    if(is.list(historical_fit)){
      hist_pred <- data.table::data.table(cbind.data.frame(
        lapply(historical_fit, function(fit) fit$predict_fold(v_task, "full"))
      ))
    } else {
      hist_pred <- historical_fit$predict(v_task)
    }
    colnames(hist_pred) <- paste0("historical_", colnames(hist_pred))
    
    # get individualized fit predictions for fold v
    ind_pred <- ind_fit$predict_fold(v_task, fold$v)
    colnames(ind_pred) <- paste0("individual_", colnames(ind_pred))
    
    # put all fold v predictions together
    pred <- cbind(hist_pred, ind_pred)
    
    # get truth for fold v
    valid_data <- train_task$data[fold$validation, ]
    truth <- valid_data[[outcome]]
    time <- valid_data[[time]]
    
    return(list(pred = pred, truth = truth, time = time))
  })
  
  ############## use v-1 fold-specific predictions to train SLs ################
  
  # train SLs on every fold fit except the most recent
  # SL lrnr weights are wrt to 3 lrnr sets: all, ind, & hist lrnrs
  trainSL_fold_fits <- fold_fits[-length(fold_fits)]
  t_truth <- Reduce(c, lapply(trainSL_fold_fits, '[[', 'truth'))
  t_pred <- rbindlist(lapply(trainSL_fold_fits, '[[', 'pred'), fill=T)
  t_pred_ind_lrnrs <- t_pred[, grep("individual_", colnames(t_pred)), with=F]
  t_pred_hist_lrnrs <- t_pred[, grep("historical_", colnames(t_pred)), with=F]
  
  wts_SLall <- sl_weights_fit(t_pred, t_truth)
  wts_SLind <- sl_weights_fit(t_pred_ind_lrnrs, t_truth)
  wts_SLhist <- sl_weights_fit(t_pred_hist_lrnrs, t_truth)
  
  ########### get ALL learner (SL and base learner) predictions ################
  # predictions correspond to validation data that was not used in training. 
  # Note this only works when the validation data never overlaps across 
  # timeseries folds, which is why we force equivalence of validation_size & 
  # batch in the timeseries cv specification.
  valid_fold_fit <- fold_fits[[length(fold_fits)]]
  v_pred <- valid_fold_fit$pred
  
  # get SL learner predictions
  v_pred_SLall <- sl_weights_predict(wts_SLall, v_pred)
  
  v_pred_ind_lrnrs <- v_pred[, grep("individual_", colnames(v_pred)), with=F]
  v_pred_SLind <- sl_weights_predict(wts_SLind, v_pred_ind_lrnrs)
  colnames(v_pred_SLind) <- paste0("individualSL_", colnames(v_pred_SLind))
  
  v_pred_hist_lrnrs <- v_pred[, grep("historical_", colnames(v_pred)), with=F]
  v_pred_SLhist <- sl_weights_predict(wts_SLhist, v_pred_hist_lrnrs)
  colnames(v_pred_SLhist) <- paste0("historicalSL_", colnames(v_pred_SLhist))
  
  # put SLs predictions together with individual learner predictions
  v_pred <- cbind(v_pred_SLall, v_pred_SLind, v_pred_SLhist, v_pred)
  
  # retain other info that will be relevant
  v_time <- valid_fold_fit$time
  v_truth <- valid_fold_fit$truth
  
  # calculate the loss under these honest predictions
  v_loss <- data.table::data.table(apply(v_pred, 2, function(p) (p-v_truth)^2))
  loss_tbl <- data.table::data.table(time=v_time, v_loss)
  
  # combine with previous losses 
  if(!first_fit){
    # add these current_loss rows to existing table of losses that were also 
    # obtained from predictions that were never used in training 
    loss_tbl <- data.table::data.table(rbind(loss_table, loss_tbl, fill=T))
  }
  
  ########################### select the discrete SL ###########################
  
  # create weights based on lags, ie. distance from current the time
  times <- loss_tbl[["time"]]
  
  # intialize weights of 1 for all losses
  weights <- rep(1, length(times))
  
  # update weights based on weights_control list
  if (!is.null(weights_control$window)) {
    window <- max(times) - weights_control$window
    weights <- weights * ifelse(times <= window, 0, 1)
  }
  
  if (!is.null(weights_control$rate_decay)) {
    lags <- max(times) - times
    if (!is.null(weights_control$delay_decay)) {
      lags_delayed <- lags - weights_control$delay_decay
      lags <- ifelse(lags_delayed < 0, 0, lags_delayed)
    }
    weights <- weights * (1 - weights_control$rate_decay)^lags
  }
  
  # calculate the cv_risk, ie. weighted mean loss for each learner
  losses <- loss_tbl[, -"time", with=F]
  cv_risk <- t(apply(losses, 2, function(loss) weighted.mean(loss, weights)))
  # subset to set of individually trained lrnrs, and historically trained lrnrs
  cv_risk_ind <- cv_risk[, grep("individual", colnames(cv_risk))]
  cv_risk_hist <- cv_risk[, grep("historical", colnames(cv_risk))]
  
  # select the online SL under various learner subsets
  dSL_all <- colnames(cv_risk)[which.min(cv_risk)]
  dSL_ind <- names(cv_risk_ind)[which.min(cv_risk_ind)]
  dSL_hist <- names(cv_risk_hist)[which.min(cv_risk_hist)]
  
  SL_tbl <- data.table::data.table(id, min(v_time), max(v_time), dSL_all, 
                                   dSL_ind, dSL_hist)
  colnames(SL_tbl) <- c("ID", "min_start", "min_end", "SL_adapt", 
                        "SL_individual", "SL_historical")
  
  ################# forecast with online SL and other learners #################
  forecast_data <- data.table::data.table(forecast_data)
  no_outcome <- rep(0, nrow(forecast_data))
  forecast_data <- data.table::data.table(no_outcome, forecast_data)
  
  forecast_fold <- origami::make_folds(forecast_data, folds_vfold, V = 1)
  forecast_task <- sl3::make_sl3_Task(forecast_data, covariates, "no_outcome", 
                                      folds = forecast_fold, time = time, 
                                      outcome_type = train_task$outcome_type$type)
  forecast_task <- process_task(forecast_task, historical_task)
  forecast_task <- process_task(forecast_task, train_task)
  
  # get forecasts from all learners
  if(is.list(historical_fit)){
    hist_forecast <- data.table::data.table(cbind.data.frame(
      lapply(historical_fit, function(fit) fit$predict_fold(forecast_task, "full"))
    ))
  } else {
    hist_forecast <- historical_fit$predict(forecast_task)
  }
  colnames(hist_forecast) <- paste0("historical_", colnames(hist_forecast))
  
  ind_forecast <- ind_fit$predict_fold(forecast_task, "full")
  colnames(ind_forecast) <- paste0("individual_", colnames(ind_forecast))
  
  forecast <- cbind(hist_forecast, ind_forecast)
  
  # get SL learner forecasts
  forecast_SLall <- sl_weights_predict(wts_SLall, forecast)
  forecast_SLind <- sl_weights_predict(wts_SLind, ind_forecast)
  colnames(forecast_SLind) <- paste0("individualSL_", colnames(forecast_SLind))
  forecast_SLhist <- sl_weights_predict(wts_SLhist, hist_forecast)
  colnames(forecast_SLhist) <- paste0("historicalSL_", colnames(forecast_SLhist))
  
  # retain the forecasts that correspond to the online SL
  forecast <- cbind(forecast, forecast_SLall, forecast_SLind, forecast_SLhist)
  SL_adapt <- forecast[[as.character(dSL_all)]]
  SL_individual <- forecast[[as.character(dSL_ind)]]
  SL_historical <- forecast[[as.character(dSL_hist)]]
  forecast <- cbind(SL_adapt, SL_individual, SL_historical, forecast)
  forecast <- data.table::data.table(time = forecast_task$time, forecast)
  
  return(list(loss_table = loss_tbl, SL_table = SL_tbl, 
              forecast_table = forecast, individual_fit = ind_fit))
}