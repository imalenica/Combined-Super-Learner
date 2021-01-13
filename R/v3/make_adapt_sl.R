slstream <- function(training_data, outcome, covariates, id, time, fold_fun, 
                     batch, horizon, window_size = NULL, first_window = NULL, 
                     individual_stack = NULL, previous_individual_fit = NULL, 
                     previous_fold_fits = NULL, historical_fit, first_fit = FALSE, 
                     loss_table = NULL, weights_control, forecast_data,
                     print_timers = TRUE, forecast_with_full_fit = FALSE,
                     full_fit_threshold = NULL, outcome_bounds = NULL){
  
  # start timer
  start_time <- proc.time()
  
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
  # gap = horizon-batch ensures outcomes are horizon-away from being a covariate
  if(fold_fun == "folds_rolling_origin"){
    folds <- origami::make_folds(training_data, fold_fun = folds_rolling_origin,
                                 first_window = first_window, gap = horizon-batch, 
                                 validation_size = batch, batch = batch)
  } else if(fold_fun == "folds_rolling_window"){
    folds <- origami::make_folds(training_data, fold_fun = folds_rolling_window, 
                                 window_size = window_size, gap = horizon-batch,
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
  
  if(!is.null(previous_individual_fit)){
    # update learners with new folds
    ind_fit <- previous_individual_fit$update(train_task)
  } else {
    # is individual_stack already a cv stack?
    if(grepl("Lrnr_cv", class(individual_stack)[1])){ 
      cv_stack <- individual_stack 
    } else {
      cv_stack <- Lrnr_cv$new(individual_stack)
    }
    ind_fit <- cv_stack$train(train_task)
  }
  # timer for training learners
  ind_train_time <- proc.time()
  
  ########################### get fold-specific predictions ####################
  if(!is.null(previous_fold_fits) & !is.null(previous_individual_fit)){
    # check correspondence between previous fold-specific predictions & truths
    valid_fold_fit <- unlist(lapply(seq_len(length(folds)), function(x) {
      if (x > length(previous_individual_fit$training_task$folds)) {
        return(FALSE)
      } else {
        prev_fold <- previous_individual_fit$training_task$folds[[x]]
        eq_train <- all.equal(prev_fold$training_set, folds[[x]]$training_set)
        eq_valid <- all.equal(prev_fold$validation_set, folds[[x]]$validation_set)
        train_check <- eq_train & eq_valid
        
        d <- train_task$data[folds[[x]]$validation_set, ]
        eq_time <- all.equal(previous_fold_fits[[x]]$time, d[[time]])
        eq_obs <- all.equal(previous_fold_fits[[x]]$obs, d[[outcome]])
        pred_check <- eq_time & eq_obs
        
        return(train_check & pred_check)
      }
    }))
    cat("\nUpdating", length(valid_fold_fit) - sum(valid_fold_fit), "fold(s)\n")
    previous_fold_fits <- previous_fold_fits[which(valid_fold_fit)]
    new_folds <- folds[which(!valid_fold_fit)]
  } else {
    previous_fold_fits <- NULL
    new_folds <- folds
  }
  
  fold_fits <- lapply(new_folds, function(fold){
    
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
    if(!is.null(outcome_bounds)){
      pred <- data.table(apply(pred, 2, truncate, outcome_bounds))
    }
    
    
    # get truth for fold v
    valid_data <- train_task$data[fold$validation, ]
    obs <- valid_data[[outcome]]
    time <- valid_data[[time]]
    
    return(list(pred = pred, obs = obs, time = time))
  })
  
  # recombine with previous fold-specific list of pred & obs
  if(length(previous_fold_fits) > 0){
    fold_fits <- c(previous_fold_fits, fold_fits)
  } 
  fit_folds_time <- proc.time()
  
  ############## use v-1 fold-specific predictions to train SLs ################
  
  # train SLs on every fold fit except the most recent
  # SL lrnr weights are wrt to 3 lrnr sets: all, ind, & hist lrnrs
  trainSL_fold_fits <- fold_fits[-length(fold_fits)]
  t_obs <- Reduce(c, lapply(trainSL_fold_fits, '[[', 'obs'))
  t_pred <- rbindlist(lapply(trainSL_fold_fits, '[[', 'pred'), fill=T)
  t_pred_ind_lrnrs <- t_pred[, grep("individual_", colnames(t_pred)), with=F]
  t_pred_hist_lrnrs <- t_pred[, grep("historical_", colnames(t_pred)), with=F]
  
  wts_SLall <- sl_weights_fit(t_pred, t_obs)
  wts_SLind <- sl_weights_fit(t_pred_ind_lrnrs, t_obs)
  wts_SLhist <- sl_weights_fit(t_pred_hist_lrnrs, t_obs)
  
  ########### get ALL learner (SL and base learner) predictions ################
  # predictions correspond to validation data that was not used in training. 
  # Note this only works when the validation data never overlaps across 
  # timeseries folds, which is why we force equivalence of validation_size & 
  # batch in the timeseries cv specification.
  valid_fold_fit <- fold_fits[[length(fold_fits)]]
  v_pred <- valid_fold_fit$pred
  
  # get SL learner predictions
  v_pred_SLall <- sl_weights_predict(wts_SLall, v_pred)
  colnames(v_pred_SLall) <- paste0("adaptSL_", colnames(v_pred_SLall))
  
  v_pred_ind_lrnrs <- v_pred[, grep("individual_", colnames(v_pred)), with=F]
  v_pred_SLind <- sl_weights_predict(wts_SLind, v_pred_ind_lrnrs)
  colnames(v_pred_SLind) <- paste0("individualSL_", colnames(v_pred_SLind))
  
  v_pred_hist_lrnrs <- v_pred[, grep("historical_", colnames(v_pred)), with=F]
  v_pred_SLhist <- sl_weights_predict(wts_SLhist, v_pred_hist_lrnrs)
  colnames(v_pred_SLhist) <- paste0("historicalSL_", colnames(v_pred_SLhist))
  
  # put SLs predictions together with individual learner predictions
  v_preds <- cbind(v_pred_SLall, v_pred_SLind, v_pred_SLhist, v_pred)
  
  # retain other info that will be relevant
  v_time <- valid_fold_fit$time
  v_obs <- valid_fold_fit$obs
  
  # calculate the loss under these honest predictions
  v_loss <- data.table::data.table(apply(v_preds, 2, function(p) (p-v_obs)^2))
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
  
  # select the online SL (a double SL) 
  dSL <- colnames(cv_risk)[which.min(cv_risk)]
  
  # get weights for the online SL (incase it selected a SL) 
  solo_lrnrs <- colnames(losses)[-grep("SL", colnames(losses))]
  if(grepl("SL", dSL)){ # the dSL selected a SL base learner
    if(grepl("adaptSL_", dSL)){
      dSL_wts <- retain_sl_weights(dSL, wts_SLall, v_pred)
    } else if(grepl("individualSL_", dSL)){
      dSL_wts <- retain_sl_weights(dSL, wts_SLind, v_pred_ind_lrnrs)
    } else if(grepl("historicalSL_", dSL)){
      dSL_wts <- retain_sl_weights(dSL, wts_SLhist, v_pred_hist_lrnrs)
    }
    # incorporate all learners in losses table wrt to losses table order
    dSL_wts <- dSL_wts[match(solo_lrnrs, names(dSL_wts))]
    names(dSL_wts) <- solo_lrnrs
    dSL_wts[is.na(dSL_wts)] <- 0 # assign weight 0 to unused learners
  } else { # the dSL selected a non-SL base learner
    # assign weight 0 to all learners that are not the non-SL dSL
    dSL_wts <- rep(0, length(solo_lrnrs))
    names(dSL_wts) <- solo_lrnrs
    dSL_wts[which(solo_lrnrs == dSL)] <- 1
  }
  dSL_tbl <- data.table::data.table(id, min(v_time), max(v_time), dSL, t(dSL_wts))
  colnames(dSL_tbl)[1:4] <- c("id", "min_start", "min_end", "onlineSL")
  
  sl_time <- proc.time()
  
  ################# fit relevant learners for full_fit onlineSL ################
  if(forecast_with_full_fit){
    nonzero_dSL_wts <- dSL_wts[which(dSL_wts != 0)]
    dSL_lrnrs <- names(nonzero_dSL_wts)
    full_fits <- NULL
    # do the full fit if individually trained learners were selected
    if(any(grepl("individual_", dSL_lrnrs))){ 
      ind_dSL_lrnrs <- dSL_lrnrs[grep("individual_", dSL_lrnrs)]
      # make fold for a full-ish fit trained on at most full_fit_threshold obs
      if(is.null(full_fit_threshold) | nrow(training_data) <= full_fit_threshold){
        full_fit_data <- data.table(training_data)
      } else {
        start_idx <- nrow(training_data) - full_fit_threshold + 1
        full_fit_data <- data.table(training_data[start_idx:(nrow(training_data)), ])
      }
      full_fit_task <- sl3::make_sl3_Task(full_fit_data, covariates, outcome, 
                                          time = time)
      full_fit_task <- process_task(full_fit_task, historical_task)
      
      # check if there's a pipeline embedded in the list of learners
      cv_stack_learners <- ind_fit$params$learner$params$learners
      lrnr_classes <- sapply(cv_stack_learners, function(x) class(x)[1])
      pipe <- names(lrnr_classes)[grep("Pipeline", lrnr_classes)]
      if(length(pipe) > 0){
        if(any(grepl(pipe, ind_dSL_lrnrs))){
          # train all learners in the pipeline (couldn't figure out a workaround)
          nonpipe_ind_dSL_lrnrs <- ind_dSL_lrnrs[-grep(pipe, ind_dSL_lrnrs)]
          pipe_full_fit <- cv_stack_learners[[pipe]]$train(full_fit_task)
        } else {
          nonpipe_ind_dSL_lrnrs <- ind_dSL_lrnrs
          pipe_full_fit <- NULL
        }
      } else {
        nonpipe_ind_dSL_lrnrs <- ind_dSL_lrnrs
        pipe_full_fit <- NULL
      }
      
      # train non-pipeline learners selected by the online SL
      if(length(nonpipe_ind_dSL_lrnrs) > 0){
        cat(nonpipe_ind_dSL_lrnrs)
        lrnr_full_fits <- lapply(nonpipe_ind_dSL_lrnrs, function(lrnr){
          # subset to learner as it was named in Stack, and train it
          lrnr <- gsub("individual_", "", lrnr)
          full_fit <- cv_stack_learners[[lrnr]]$train(full_fit_task)
          return(full_fit)
        })
        names(lrnr_full_fits) <- nonpipe_ind_dSL_lrnrs
      } else {
        lrnr_full_fits <- NULL
      }
      
      # augment full_fits list w the relevant learner full fits 
      stopifnot(!is.null(lrnr_full_fits) | !is.null(pipe_full_fit))
      if(!is.null(pipe_full_fit) & !is.null(lrnr_full_fits)){
        full_fits <- c(pipeline = list(pipe_full_fit), lrnr_full_fits)
        names(full_fits)[1] <- pipe
      } else if (is.null(pipe_full_fit) & !is.null(lrnr_full_fits)){
        full_fits <- lrnr_full_fits
      } else if (!is.null(pipe_full_fit) & is.null(lrnr_full_fits)){
        full_fits <- list(pipe_full_fit)
        names(full_fits)[1] <- pipe
      }
    } 
  }
  full_fit_time <- proc.time()
  
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
  
  ind_forecast <- ind_fit$predict_fold(forecast_task, length(folds))
  colnames(ind_forecast) <- paste0("individual_", colnames(ind_forecast))
  
  forecast <- cbind(hist_forecast, ind_forecast)
  if(!is.null(outcome_bounds)){
    forecast <- data.table(apply(forecast, 2, truncate, outcome_bounds))
  }
  
  # get SL learner forecasts
  forecast_SLall <- sl_weights_predict(wts_SLall, forecast)
  colnames(forecast_SLall) <- paste0("adaptSL_", colnames(forecast_SLall))
  forecast_SLind <- sl_weights_predict(wts_SLind, ind_forecast)
  colnames(forecast_SLind) <- paste0("individualSL_", colnames(forecast_SLind))
  forecast_SLhist <- sl_weights_predict(wts_SLhist, hist_forecast)
  colnames(forecast_SLhist) <- paste0("historicalSL_", colnames(forecast_SLhist))
  
  forecast <- cbind(forecast, forecast_SLall, forecast_SLind, forecast_SLhist)
  
  # retain the forecasts that correspond to the online SL
  onlineSL <- forecast[[as.character(dSL)]]
  forecast <- cbind(onlineSL, forecast)
  
  # retain the forecasts that correspond to the full_fit online SL
  if(forecast_with_full_fit){
    if(!is.null(full_fits)){
      ff_forecast_list <- lapply(full_fits, function(fit) fit$predict(forecast_task))
      ff_forecast <- do.call(cbind, lapply(seq_along(ff_forecast_list), function(i){
        forecast <- ff_forecast_list[[i]]
        if(names(ff_forecast_list)[i] == pipe){
          colnames(forecast) <- paste0("individual_", pipe, "_", colnames(forecast))
          forecast
        } else {
          forecast <- data.table(forecast)
          colnames(forecast) <- names(ff_forecast_list)[i]
          forecast
        }
      }))
      to_rm <- which(colnames(ind_forecast) %in% colnames(ff_forecast))
      ind_forecast_reduced <- ind_forecast[, -to_rm, with=F]
      dSL_forecast_tbl <- data.table(hist_forecast, ind_forecast_reduced, ff_forecast)
    } else {
      dSL_forecast_tbl <- data.table(hist_forecast, ind_forecast)
    }
    if(!is.null(outcome_bounds)){
      dSL_forecast_tbl <- data.table(apply(dSL_forecast_tbl, 2, truncate, outcome_bounds))
    }
    onlineSL_full_fit <- osl_weights_predict(dSL_wts, dSL_forecast_tbl)
    forecast <- cbind(onlineSL_full_fit, forecast)
  }
  forecast_timer <- proc.time()
  #################################### timers ##################################
  # add the time corresponding to the forecasts and the current time
  forecast_time <- rep(max(forecast_task$time), length(forecast_task$time))
  # we can be more precise with the real time, and add the fit time
  end_time <- proc.time() - start_time
  if(time == "min"){
    forecast_time_precise <- forecast_time + end_time["elapsed"]/60
    tbl <- data.table::data.table(time = forecast_task$time, forecast_time, 
                                  forecast_time_precise, forecast)
  } else {
    tbl <- data.table::data.table(time = forecast_task$time, forecast_time, forecast)
  }
  
  timers <- list(
    run_timer = end_time, 
    ind_train_timer = ind_train_time-start_time, 
    fit_folds_timer = fit_folds_time-ind_train_time, 
    SL_timer = sl_time-fit_folds_time, 
    full_fit_timer = full_fit_time-sl_time, 
    forecast_timer = forecast_timer-full_fit_time 
  )
  
  if(print_timers){
    pretty <- lapply(timers, function(x) round(as.numeric(x["elapsed"]/60), 4))
    names(pretty) <- paste0(names(pretty), "_min")
    print(pretty)
  }
  
  return(list(loss_table = loss_tbl, sl_table = dSL_tbl, forecast_table = tbl, 
              individual_fit = ind_fit, fold_fits = fold_fits, timers = timers))
}

truncate <- function(x, bounds){
  x[x < min(bounds)] <- min(bounds)
  x[x > max(bounds)] <- max(bounds)
  return(x)
}

