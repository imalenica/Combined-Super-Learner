make_adapt_and_online_sl <- function(individual_training_data, 
                                     outcome, 
                                     covariates, 
                                     id, 
                                     time = "min", 
                                     first_window, 
                                     batch, 
                                     individual_stack, 
                                     historical_fit, 
                                     individual_fit = NULL, 
                                     first_fit = FALSE, 
                                     loss_table = NULL, 
                                     previous_onlineSL_forecasts = NULL,
                                     onlineSL_table = NULL, 
                                     individual_forecast_data, 
                                     evaluate_onlineSL_window = list("all", 10, 5, 1),
                                     trainSL_window = list("all", 10, 5, 1)){
  
  ################################ check arguments #############################
  
  # check trainsl list window is not null and no duplicates
  
  # need one of these two args
  if(is.null(individual_fit) & is.null(individual_stack)){
    stop("Individual fit and individual stack both missing. Provide ",
         "individual_stack to initialize a new fit, and past_individual_fit ",
         "to update with new data.")
  }
  
  if(!first_fit){
    if(is.null(loss_table)){
      stop("Table of all learner losses must be provided if not first fit")
    }
    if(is.null(previous_onlineSL_forecasts)){
      stop("Previous onlineSL forecasts must be provided if not first fit")
    }
    if(is.null(onlineSL_table)){
      stop("Table with previous onlineSLs must be provided if not first fit.")
    }
  }
  
  ################################## train SLs #################################
  
  individual_training_data <- data.table::data.table(individual_training_data)
  
  # make folds with training data
  folds <- origami::make_folds(individual_training_data, 
                               fold_fun = folds_rolling_origin,
                               first_window = first_window, 
                               validation_size = batch, 
                               gap = 0,
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
    ind_fit <- individual_fit$update(ind_training_task) # update learners
  }
  
  ########################### get fold-specific predictions ####################
  fold_fits <- lapply(folds, function(fold){
    # retain validation task fold v
    v_task <- validation_task(ind_training_task, fold)
    
    # get predictions for fold v
    if(is.list(historical_fit)){
      hist_preds <- cbind.data.frame(lapply(historical_fit, function(fit){
        fit$predict_fold(v_task, "full")
      }))
    } else {
      hist_preds <- historical_fit$predict_fold(v_task, "full")
    }
    colnames(hist_preds) <- paste0("historical_", colnames(hist_preds))
    
    ind_preds <- ind_fit$predict_fold(v_task, fold_number = fold$v)
    colnames(ind_preds) <- paste0("individual_", colnames(ind_preds))
    
    all_preds <- data.table(cbind.data.frame(hist_preds, ind_preds))
    
    # get truth for fold v
    valid_data <- data.table(ind_training_task$data[fold$validation, ])
    truth <- valid_data[[outcome]]
    time <- valid_data[[time]]
    
    # get loss fold v
    all_loss <- data.table(t(apply(all_preds, 2, function(p) mean((p-truth)^2))))
    ind_loss <- data.table(t(apply(ind_preds, 2, function(p) mean((p-truth)^2))))
    hist_loss <- data.table(t(apply(hist_preds, 2, function(p) mean((p-truth)^2))))
    
    return(list(all_preds=all_preds, ind_preds=ind_preds, hist_preds=hist_preds,
                all_loss=all_loss, ind_loss=ind_loss, hist_loss=hist_loss,
                truth=truth, time=time))
  })
  
  ###################### fit SL to v-1 fold-specific predictions ###############
  if(first_fit){ 
    # do not train SL
    preds <- data.table(fold_fits[[1]]$all_preds)
    
  } else { 
    # train SL on every fold fit except the last one
    trainSL_fold_fits <- fold_fits[-length(fold_fits)]
    validSL_fold_fit <- fold_fits[[length(fold_fits)]]
    
    # omit windows that are not less than length(trainSL_fold_fits)
    trainSL_window <- process_window(trainSL_window, length(trainSL_fold_fits))
    
    SLfits <- lapply(trainSL_window, fitSL, trainSL_fold_fits, validSL_fold_fit)
    SLweights <- lapply(SLfits, '[[', 'weights')
    names(SLweights) <- trainSL_window
    
    SLpreds <- do.call(cbind, lapply(SLfits, '[[', 'preds'))
    preds <- data.table(cbind(SLpreds, validSL_fold_fit$all_preds))
  }
  
  ######################### save last v-fold for online SL #####################
  if(first_fit){
    times <- fold_fits[[1]]$time
    truth <- fold_fits[[1]]$truth
    
  } else {
    times <- validSL_fold_fit$time
    truth <- validSL_fold_fit$truth
    
    # retain truth, ensure correspondence to the last forecast
    previous_times <- previous_onlineSL_forecasts$forecast_time
    if(!all.equal(times, previous_times)){
      stop("Last online SL forecasts correspond to times that are not equal ",
           "these evaluation data times. Later, support will be added.")
    }
    
    #################### evaluate forecasts from last time #####################
    row_idx <- which(as.numeric(onlineSL_table$minEnd)+1 == min(previous_times))
    rm <- which(names(previous_onlineSL_forecasts) == "forecast_time")
    previous_onlineSL_forecasts <- previous_onlineSL_forecasts[-rm]
    
    previous_onlineSL <- cbind(
      onlineSL_table[row_idx, c("ID", "minStart", "minEnd"), with=F], 
      do.call(cbind, lapply(seq_along(previous_onlineSL_forecasts), function(j){
        MSE <- mean((previous_onlineSL_forecasts[[j]]-truth)^2)
        type <- names(previous_onlineSL_forecasts)[j]
        col_idx <- which(colnames(onlineSL_table) == type)
        name <- onlineSL_table[row_idx, col_idx, with=F]
        d <- data.table(name, MSE)
        colnames(d)[2] <- paste0("MSE_", type)
        return(d)
      }))
    )
    if(row_idx > 1){ 
      # previous online SLs exist before this one
      past_onlineSL_table <- onlineSL_table[-row_idx, ]
      onlineSL_table <- rbind(past_onlineSL_table, previous_onlineSL, fill=T)
    } else { 
      # previous online SLs do not exist before this one
      onlineSL_table <- previous_onlineSL
    }
  }
  
  ############################# make new online SLs ############################
  current_loss <- data.table(t(apply(preds, 2, function(p) mean((p-truth)^2))))
  
  if(first_fit){
    loss_tbl <- current_loss
  } else { # add this row to existing risk table
    loss_tbl <- data.table(rbind(loss_table, current_loss, fill=T))
  }
  
  oSL_window <- process_window(evaluate_onlineSL_window, nrow(loss_tbl))
  
  new_oSL <- do.call(
    cbind, lapply(oSL_window, fit_onlineSL, loss_tbl, current_loss, first_fit)
  )
  names_new_onlineSL <- c(new_oSL[, -grepl("MSE", colnames(new_oSL)), with=F])
  new_onlineSL <- data.table(cbind(
    data.table(t(c(ID=id, minStart=min(times), minEnd=max(times)))), new_oSL
  ))
  
  if(first_fit){
    onlineSL_table <- new_onlineSL
  } else {
    onlineSL_table <- rbind(onlineSL_table, new_onlineSL, fill=T)
  }
  
  ######################### forecast with new online SLs #######################
  individual_forecast_data <- data.table(individual_forecast_data)
  unknown_outcome <- rep(0, nrow(individual_forecast_data))
  forecast_data <- data.table(unknown_outcome, individual_forecast_data)
  
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
  
  # get forecasts from all learners
  if(is.list(historical_fit)){
    hist_forecasts <- cbind.data.frame(lapply(historical_fit, function(fit){
      fit$predict_fold(forecast_task, "full")
    }))
  } else {
    hist_forecasts <- historical_fit$predict_fold(forecast_task, "full")
  }
  colnames(hist_forecasts) <- paste0("historical_", colnames(hist_forecasts))
  
  ind_forecasts <- ind_fit$predict_fold(forecast_task, "full")
  colnames(ind_forecasts) <- paste0("individual_", colnames(ind_forecasts))
  
  forecasts <- data.table(cbind.data.frame(hist_forecasts, ind_forecasts))
  
  if(!first_fit){
    # use learner forecasts to create SL forecasts with SL weights
    SLforecasts <- do.call(cbind, lapply(seq_along(SLweights), function(j){
      
      window_specific_weights_list <- SLweights[[j]]
      awts <- window_specific_weights_list$allSL_weights
      iwts <- window_specific_weights_list$indSL_weights
      hwts <- window_specific_weights_list$histSL_weights
      
      acol <- colnames(forecasts)[which(colnames(forecasts) %in% colnames(awts))]
      aSL <- sl_weights_predict(awts[,acol], forecasts[, acol, with=F])
      
      icol <- colnames(ind_forecasts)[which(colnames(ind_forecasts) %in% colnames(iwts))]
      iSL <- sl_weights_predict(iwts[,icol], forecasts[, icol, with=F])
      
      hcol <- colnames(hist_forecasts)[which(colnames(hist_forecasts) %in% colnames(hwts))]
      hSL <- sl_weights_predict(hwts[,hcol], forecasts[, hcol, with=F])
      
      w <- names(SLweights)[j]
      if(w != "all"){
        colnames(aSL) <- paste0("window", w, "_", colnames(aSL))
        colnames(iSL) <- paste0("individualSL", "window", w, "_", colnames(iSL))
        colnames(hSL) <- paste0("historicalSL", "window", w, "_", colnames(hSL))
      } else {
        colnames(iSL) <- paste0("individualSL_", colnames(iSL))
        colnames(hSL) <- paste0("historicalSL_", colnames(hSL))
      }
      
      return(cbind(aSL, iSL, hSL))
    }))
    
    # augment forecasts with SLforecasts
    forecasts <- data.table(cbind(forecasts, SLforecasts))
  }
  
  # retain the forecasts that correspond to onlineSLs
  onlineSL_forecasts <- lapply(names_new_onlineSL, function(x){
    forecasts[[as.character(x)]]
  })
  names(onlineSL_forecasts) <- names(names_new_onlineSL)
  onlineSL_forecasts <- c(
    list(forecast_time = forecast_task$time), onlineSL_forecasts
  )
  
  # save all forecasts for book-keeping, potential results tables/figures
  all_forecasts <- data.table(do.call(cbind, onlineSL_forecasts), forecasts)
  
  return(list(
    loss_table = loss_tbl,
    onlineSL_table = onlineSL_table,
    onlineSL_forecasts = onlineSL_forecasts,
    individual_fit = ind_fit,
    all_forecasts = all_forecasts
  ))
}
