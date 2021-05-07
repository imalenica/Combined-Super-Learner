slstream <- function(
  training_data, outcome, outcome_bounds = NULL, covariates, id, time, 
  fold_fun, batch, horizon, window_size = NULL, first_window = NULL, 
  historical_fit, 
  individual_stack, train_individual_stack,
  first_fit = FALSE, loss_table = NULL, weights_control, 
  previous_individual_fit = NULL, previous_fold_fits = NULL,  
  forecast_with_full_fit = FALSE, full_fit_threshold = NULL, forecast_data, 
  print_timers = TRUE
  ){
  
  # start timer
  start_time <- proc.time()
  
  ################################ check arguments #############################
  if(!first_fit){
    if(is.null(loss_table)){
      stop("Table of all learner losses must be provided if not first fit.")
    }
  }
  
  ################################## format task ###############################
  training_data <- data.table::data.table(training_data)
  
  # make folds with training data
  # gap = horizon-batch ensures outcomes are horizon-away from being a covariate
  if(fold_fun == "folds_rolling_origin"){
    folds <- origami::make_folds(
      training_data, fold_fun = folds_rolling_origin, first_window = first_window, 
      gap = horizon-batch, validation_size = batch, batch = batch
    )
  } else if(fold_fun == "folds_rolling_window"){
    folds <- origami::make_folds(
      training_data, fold_fun = folds_rolling_window, window_size = window_size, 
      gap = horizon-batch, validation_size = batch, batch = batch
    )
  }
  
  train_task <- sl3::make_sl3_Task(
    training_data, covariates, outcome, time = time, folds = folds
  )
  
  # issue with mismatch between historical and individual delta column:
  if(is.list(historical_fit)){
    # historical_fit is a list when fit with SL, see fun make_historical_sl
    historical_task <- historical_fit$cv_fit$fit_object$full_fit$training_task
  } else {
    historical_task <- historical_fit$training_task
  }
  train_task <- process_task(train_task, historical_task)
  task_time <- proc.time()
  
  ################################## train individual ##########################
  if(train_individual_stack){
    # get fold-specific predictions 
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
      if(any(valid_fold_fit)){
        cat("\nUpdating individualized training of", 
        length(valid_fold_fit) - sum(valid_fold_fit), "fold(s)\n")
        ind_fit <- previous_individual_fit$update(train_task)
        previous_fold_fits <- previous_fold_fits[which(valid_fold_fit)]
        new_folds <- folds[which(!valid_fold_fit)]
      } else {
        cat("\nIndividualized training across all folds, supplied previous fits do not correspond\n")
        if(!grepl("Lrnr_cv", class(individual_stack)[1])){ 
          individual_stack <- Lrnr_cv$new(individual_stack)
        }
        ind_fit <- individual_stack$train(train_task)
        previous_fold_fits <- NULL
        new_folds <- folds
      }
    } else {
      cat("\nBeginning individualized training\n")
      if(!(grepl("Lrnr_cv", class(individual_stack)[1]))){ 
        individual_stack <- Lrnr_cv$new(individual_stack)
      }
      ind_fit <- individual_stack$train(train_task)
      previous_fold_fits <- NULL
      new_folds <- folds
    }
  } else {
    previous_fold_fits <- NULL
    new_folds <- folds
    ind_fit <- NULL
  }
  # timer for training learners
  ind_train_time <- proc.time()
  
  fold_fits <- lapply(new_folds, function(fold){
    
    # retain validation task fold v
    v_task <- sl3::validation_task(train_task, fold)
    
    # get historical fit predictions for fold v
    if(is.list(historical_fit)){
      hist_pred <- data.table::data.table(cbind.data.frame(
        lapply(historical_fit, function(fit) fit$predict_fold(v_task, "full"))
      ))
    } else {
      hist_pred <- try(historical_fit$predict(v_task), silent=TRUE)
    }
    colnames(hist_pred) <- paste0("historical_", colnames(hist_pred))
    
    if(train_individual_stack){
      # get individualized fit predictions for fold v
      ind_pred <- ind_fit$predict_fold(v_task, fold$v)
      colnames(ind_pred) <- paste0("individual_", colnames(ind_pred))
      # put all fold v predictions together
      pred <- cbind(hist_pred, ind_pred)
    } else {
      pred <- hist_pred
    }
    
    if(!is.null(outcome_bounds)){
      pred <- data.table::data.table(apply(pred, 2, bound, outcome_bounds))
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
  preds_time <- proc.time()
  
  ########### get ALL learner predictions ################
  # predictions correspond to validation data that was not used in training. 
  # Note this only works when the validation data never overlaps across 
  # timeseries folds, which is why we force equivalence of validation_size & 
  # batch in the timeseries cv specification.
  valid_fold_fit <- fold_fits[[length(fold_fits)]]
  v_pred <- valid_fold_fit$pred
  
  if(length(folds) >= 2){
    ############## use v-1 fold-specific predictions to train SLs ##############
    # train SLs on every fold fit except the most recent
    # SL lrnr weights are wrt to 3 lrnr sets: all, ind, & hist lrnrs
    trainSL_fold_fits <- fold_fits[-length(fold_fits)]
    t_obs <- Reduce(c, lapply(trainSL_fold_fits, '[[', 'obs'))
    t_pred <- rbindlist(lapply(trainSL_fold_fits, '[[', 'pred'), fill=T)
    wts_SL <- sl_weights_fit(t_pred, t_obs, train_task$outcome_type)
    
    # get SL learner predictions
    v_pred_SL <- sl_weights_predict(wts_SL, v_pred, outcome_bounds)
    colnames(v_pred_SL) <- paste0("SL_", colnames(v_pred_SL))
    
    # put SLs predictions together with individual learner predictions
    v_pred <- cbind(v_pred, v_pred_SL)
  }
  
  # retain other info that will be relevant
  v_time <- valid_fold_fit$time
  v_obs <- valid_fold_fit$obs
  
  # calculate the loss under these honest predictions
  if(train_task$outcome_type$type == "continuous"){
    v_loss <- data.table::data.table(apply(v_pred, 2, function(p) (p-v_obs)^2))
  } else if(train_task$outcome_type$type %in% c("binomial", "constant")){
    v_loss <- data.table::data.table(apply(v_pred, 2, function(p){
      -1 * ifelse(v_obs == 1, log(bound(p)), log(1 - bound(p)))
    }))
  }
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
  solo_lrnrs <- colnames(losses)[!grepl("SL", colnames(losses))]
  if(grepl("SL", dSL)){ # the dSL selected a SL base learner
    dSL_wts <- retain_sl_weights(dSL, wts_SL, v_pred)
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
  dSL_tbl <- data.table::data.table(
    id = id, time = max(v_time), onlineSL = dSL, t(dSL_wts)
  )
  
  cv_risk_tbl <- data.table::data.table(
    id = id, time = max(v_time), cv_risk
  )
  
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
        full_fit_data <- data.table::data.table(training_data)
      } else {
        start_idx <- nrow(training_data) - full_fit_threshold + 1
        full_fit_data <- data.table::data.table(
          training_data[start_idx:(nrow(training_data)), ]
        )
      }
      full_fit_task <- sl3::make_sl3_Task(
        full_fit_data, covariates, outcome, time = time
      )
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
  forecast_task <- sl3::make_sl3_Task(
    forecast_data, covariates, "no_outcome", folds = forecast_fold, time = time, 
    outcome_type = train_task$outcome_type$type
  )
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
  
  if(train_individual_stack){
    ind_forecast <- ind_fit$predict_fold(forecast_task, length(folds))
    colnames(ind_forecast) <- paste0("individual_", colnames(ind_forecast))
    forecast <- cbind(hist_forecast, ind_forecast)
  } else {
    forecast <- hist_forecast
  }
  
  if(!is.null(outcome_bounds)){
    forecast <- data.table::data.table(apply(forecast, 2, bound, outcome_bounds))
  }
  
  if(length(folds) >= 2){
    # get SL learner forecasts
    forecast_SL <- sl_weights_predict(wts_SL, forecast, outcome_bounds)
    colnames(forecast_SL) <- paste0("SL_", colnames(forecast_SL))
    forecast <- cbind(forecast, forecast_SL)
  }
  
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
          forecast <- data.table::data.table(forecast)
          colnames(forecast) <- names(ff_forecast_list)[i]
          forecast
        }
      }))
      if(train_individual_stack){
        to_rm <- which(colnames(ind_forecast) %in% colnames(ff_forecast))
        ind_forecast_reduced <- ind_forecast[, -to_rm, with=F]
        dSL_forecast_tbl <- data.table::data.table(
          hist_forecast, ind_forecast_reduced, ff_forecast
        )
      } else {
        dSL_forecast_tbl <- data.table::data.table(hist_forecast, ff_forecast)
      }
    } else {
      if(train_individual_stack){
        dSL_forecast_tbl <- data.table::data.table(hist_forecast, ind_forecast)
      } else {
        dSL_forecast_tbl <- data.table::data.table(hist_forecast)
      }
    }
    if(!is.null(outcome_bounds)){
      dSL_forecast_tbl <- data.table::data.table(apply(
        dSL_forecast_tbl, 2, bound, outcome_bounds
      ))
    }
    onlineSL_full_fit <- osl_weights_predict(
      dSL_wts, dSL_forecast_tbl, outcome_bounds
    )
    forecast <- cbind(onlineSL_full_fit, forecast)
  }
  forecast_timer <- proc.time()
  
  #################################### timers ##################################
  # add the time corresponding to the forecasts and the current time
  forecast_time <- rep(max(forecast_task$time), length(forecast_task$time))
  # we can be more precise with the real time, and add the fit time
  end_time <- proc.time() - start_time
  
  # forecast_time_precise assumes time is in minutes
  forecast_time_precise <- forecast_time + end_time["elapsed"]/60
  tbl <- data.table::data.table(
    time = forecast_task$time, forecast_time, forecast_time_precise, forecast
  )
  
  timers <- list(run_timer = end_time, task_timer = task_time-start_time,
                 ind_train_timer = ind_train_time-task_time, 
                 fit_folds_timer = preds_time-ind_train_time, 
                 SL_timer = sl_time-preds_time, 
                 full_fit_timer = full_fit_time-sl_time, 
                 forecast_timer = forecast_timer-full_fit_time)
  timers <- as.matrix(t(
    sapply(timers, function(x) round(as.numeric(x["elapsed"]/60), 4))
  ))
  colnames(timers) <- paste0(colnames(timers), "_min")
  
  if(print_timers){
    print(timers)
  }
  
  return(list(loss_table = loss_tbl, risk_table = cv_risk_tbl, 
              sl_table = dSL_tbl, forecast_table = tbl, 
              individual_fit = ind_fit, fold_fits = fold_fits, timers = timers))
}

