run_slstream <- function(individual_data, individual_id, outcome, covariates, 
                         time = "time", cv_stack, slstream_weights_control, 
                         individual_burn_in = 30, batch = 5, max_stop_time = 1440, 
                         batch_multiplier_for_window = 100, historical_fit,
                         drop_missing_outcome = TRUE, print_timers = TRUE,
                         forecast_with_full_fit = FALSE, 
                         full_fit_threshold = NULL, outcome_bounds, 
                         horizon, inner_nfolds){
  
  require(data.table)
  
  ############################### preliminary checks ###########################
  individual_data <- data.table(individual_data)
  if(length(unique(individual_data[["id"]])) > 1){
    individual_data$id <- as.character(individual_data$id)
    individual_data <- individual_data[id == individual_id,]
  }
  setorder(individual_data, time_and_date)
  
  missing_outcomes <- is.na(individual_data[[outcome]])
  if(sum(missing_outcomes) > 0){
    print(paste0("Dropping ", sum(missing_outcomes), " missing outcomes."))
    individual_data <- individual_data[-which(missing_outcomes), ]
  }
  
  outcomes <- individual_data[[outcome]]
  if(length(unique(outcomes)) == 2){
    min_burn_in1 <- min(which(cumsum(outcomes) == inner_nfolds))
    min_burn_in0 <- min(which(cumsum(abs(outcomes-1)) == inner_nfolds))
    min_burn_in <- max(min_burn_in1, min_burn_in0)
    min_batches <- ceiling(min_burn_in/batch)
    min_burn_in <- min_batches*batch
    if(individual_burn_in < min_burn_in){
      print(paste0("Resetting individualized training burn-in to ", min_burn_in))
      individual_burn_in <- min_burn_in
    }
  }
  
  ########################### set up pseudo-online data ########################
  max_times <- c(max_stop_time, nrow(individual_data))
  max_time <- max_times[which.min(max_times)]
  if(max_time < 100) {
    stop(paste0("Insufficient individual training time for subject :", individual_id))
  }
  
  gap <- horizon-batch

  splits <- seq((3*batch + gap), max_time, batch)
  train_data_list <- lapply(splits, function(x) individual_data[1:x, ])
  pred_data_list <- lapply(splits, function(x) individual_data[(x+gap+1):(x+batch+gap), ])
  
  na <- unlist(lapply(pred_data_list, function(x) !any(is.na(x[[outcome]]))))
  pred_data_list <- pred_data_list[na]
  train_data_list <- train_data_list[na]
  stopifnot(length(pred_data_list) == length(train_data_list))
  N <- length(pred_data_list)
  
  # no. obs for 1st rw fit
  window <- batch * batch_multiplier_for_window
  n_first_rolling_window_cv <- window + 2*batch + gap 
  
  # no. obs for 1st individualized fit
  n_first_individual_training <- individual_burn_in + 2*batch + gap
  start_individual_training <- min(which(splits >= n_first_individual_training))
  if(start_individual_training > N){
    print(paste0("Cannot individually train subject :", individual_id))
  }
  
  # function to increase burn_in, switch folds, 
  get_burn_in <- function(n){
    if (n_first_individual_training <= n_first_rolling_window_cv){
      if(n < n_first_individual_training){
        return(batch)
      } else if (n >= n_first_individual_training & n <= n_first_rolling_window_cv){
        return(individual_burn_in)
      } else if (n > n_first_rolling_window_cv){
        return(window)
      }
    } else if (n_first_individual_training > n_first_rolling_window_cv){
      if(n < n_first_rolling_window_cv){
        return(batch)
      } else if (n >= n_first_rolling_window_cv & n < n_first_individual_training){
        return(window)
      } else if (n >= n_first_individual_training){
        return(individual_burn_in)
      }
    }
  }
  
  if(forecast_with_full_fit){
    if(is.null(full_fit_threshold)){
      full_fit_threshold <- batch + n_first_rolling_window_cv
    }
  }
  
  ################################## run #########################################
  
  loss_tables <- list()
  individual_fits <- list()
  forecast_tables <- list()
  sl_tables <- list()
  risk_tables <- list()
  fold_fits <- list()
  timers <- matrix(ncol = 7, nrow = N)
  individualized_training <- c()
  fold_fun <- c()
  
  for(i in 1:N){
    
    cat("\n_______ starting run", i, "of", N, "for id", individual_id, "_______\n")
    if(i > 2){
      loss_tables[[i-2]] <- NA
      individual_fits[[i-2]] <- NA
      fold_fits[[i-2]] <- NA
    }
    
    # initiate "online" training and forecast data
    training_data <- data.table(train_data_list[[i]])
    forecast_data <- data.table(pred_data_list[[i]])
    burn_in <- get_burn_in(splits[[i]])
    
    # control to switch from folds_rolling_origin to folds_rolling_window
    if(burn_in < window){ 
      fold_fun[i] <- "folds_rolling_origin"
      first_window <- burn_in
      window_size <- NULL
    } else {
      fold_fun[i] <- "folds_rolling_window"
      first_window <- NULL
      window_size <- burn_in
    }
    
    if(i == 1){
      first_fit <- TRUE      
      loss_table <- NULL
      previous_fold_fits <- NULL
      previous_individual_fit <- NULL
      individualized_training[i] <- FALSE
    } else {
      first_fit <- FALSE
      loss_table <- loss_tables[[i-1]]
      previous_fold_fits <- fold_fits[[i-1]]
      previous_individual_fit <- individual_fits[[i-1]]
      # control to begin individualized training
      if(i == start_individual_training){
        cat("\n....... starting individualized training .......\n")
        individualized_training[i] <- TRUE
      } else {
        individualized_training[i] <- individualized_training[i-1]
      }
      
      if(fold_fun[i] != fold_fun[i-1]){
        cat("\n....... switching fold type .......\n")
        previous_fold_fits <- NULL
      }
    }
    
    fit <- slstream(
      training_data = training_data, outcome = outcome, id = individual_id, 
      covariates = covariates, time = time, historical_fit = historical_fit, 
      batch = batch, fold_fun = fold_fun[i], first_window = first_window, 
      window_size = window_size, individual_stack = cv_stack, 
      forecast_data = forecast_data, print_timers = print_timers, 
      first_fit = first_fit, weights_control = slstream_weights_control,
      loss_table = loss_table, 
      previous_fold_fits = previous_fold_fits, 
      previous_individual_fit = previous_individual_fit,
      forecast_with_full_fit = forecast_with_full_fit,
      full_fit_threshold = full_fit_threshold, outcome_bounds = outcome_bounds, 
      horizon = horizon, train_individual_stack = individualized_training[i]
    )
    loss_tables[[i]] <- fit$loss_table
    sl_tables[[i]] <- fit$sl_table
    risk_tables[[i]] <- fit$risk_table 
    forecast_tables[[i]] <- fit$forecast_table
    if(is.null(fit$individual_fit)){
      individual_fits[i] <- list(fit$individual_fit) 
    } else {
      individual_fits[[i]] <- fit$individual_fit
    }
    fold_fits[[i]] <- fit$fold_fits
    timers[i,] <- fit$timers
  }
  
  ################################ compile result ##############################
  colnames(timers) <- c("total", "task", "individualized_train", "prediction",
                        "OSL", "full_fit", "forecast")
  colnames(timers) <- paste0(colnames(timers), "_min")
  
  ##### summarize prop of individual vs. historical learner weights
  SLtbl <- rbindlist(sl_tables, fill = T)
  wts <- SLtbl[, -c("id", "time", "onlineSL"), with = F]
  ind_lrnrs <- colnames(wts)[grep("individual_", colnames(wts))]
  hist_lrnrs <- colnames(wts)[grep("historical_", colnames(wts))]
  p_ind <- rowSums(wts[, ind_lrnrs, with = F], na.rm = T) / rowSums(wts, na.rm = T)
  p_hist <- rowSums(wts[, hist_lrnrs, with = F], na.rm = T) / rowSums(wts, na.rm = T)
  SLtbl <- data.table(prop_individual = p_ind, prop_historical = p_hist, SLtbl)
  
  ##### summarize cv_risk
  risk_tables <- rbindlist(risk_tables, fill = T)
  
  ##### make table of true horizon-ahead outcomes & forecasts
  time <- individual_data[[time]] # time in which forecast was made 
  obs <- individual_data[[outcome]] # the true horizon-ahead outcome
  outcome_time <- time + horizon # time of horizon-ahead outcome
  obs_tbl <- data.table(obs, time, outcome_time)
  # merge with forecasts
  forecast_tbl <- rbindlist(forecast_tables, fill = T)
  tbl <- merge(obs_tbl, forecast_tbl, by = "time", all.x = F, all.y = T)
  
  ##### for each learner, calculate time-specific & summarized performance
  time_cols <- colnames(forecast_tbl)[grep("time", colnames(forecast_tbl))]
  preds_tbl <- forecast_tbl[, -time_cols, with = F]
  
  # calculate the loss under these honest predictions
  if(length(unique(outcomes)) == 2){
    loss_function <- function(pred, observed) {
      out <- -1 * ifelse(observed == 1, log(bound(pred)), log(1 - bound(pred)))
      return(out)
    }
  } else {
    loss_function <- function(pred, observed){ 
      (bound(pred) - observed)^2 
    }
  }
  perf_tbl <- data.table(apply(preds_tbl, 2, loss_function, tbl[["obs"]]))
  perf_summ <- c(id = individual_id, outcome = outcome, colMeans(perf_tbl, na.rm=T))
  
  ##### for each learner, calculate the number of NA
  missing_summ <- apply(perf_tbl, 2, function(result) sum(is.na(result)))
  
  # re-add time cols to MSE_tbl
  time_tbl <- forecast_tbl[, time_cols, with = F]
  perf_tbl <- data.table(time_tbl, perf_tbl)
  
  return_list <- list(
    sl_table = SLtbl, risk_table = risk_tables, forecast_table = tbl, 
    performance_summary = perf_summ, loss_table = loss_tables[[N]], 
    performance_table = perf_tbl, missing_summary = missing_summ, 
    timer_matrix = timers
  )
              
  if(length(unique(outcomes)) == 2){
    AHE_prevalence <- sum(tbl[["obs"]] == 1)/length(tbl[["obs"]])
    
    mat <- matrix(nrow = 3, ncol = ncol(preds_tbl))
    for(i in 1:ncol(preds_tbl)){
      na.pred <- is.na(preds_tbl[[i]])
      if(all(na.pred)){
        mat[1,i] <- NA
        mat[2,i] <- NA
        mat[3,i] <- NA
      } else {
        obs <- tbl[["obs"]][!na.pred]
        pred <- preds_tbl[[i]][!na.pred]
        obj <- ROCR::prediction(pred, obs)
        # baseline aucpr = 0.5
        mat[1,i] <- as.numeric(ROCR::performance(obj, "auc")@y.values)
        # baseline aucpr = fraction of positives (# positives / total # preds)
        mat[2,i] <- as.numeric(ROCR::performance(obj, "aucpr")@y.values)
        mat[3,i] <- sum(obs == 1)/length(obs)
      }
    }
    colnames(mat) <- colnames(preds_tbl)
    mat <- data.table(cbind(Metric = c("AUC", "AUCPR", "nonNA_PREV"), mat))
    return_list <- c(
      return_list, list(AHE_prevalence = AHE_prevalence, AUC_tbl = mat)
    )
  } 
  return(return_list)
}

run_slstream_allY <- function(outcomes, file_path, individual_data, 
                              individual_id, covariates, cv_stack, 
                              slstream_weights_control, time = "time",  
                              batch = 5, max_stop_time = 1440, 
                              individual_burn_in = 30, print_timers = TRUE,
                              batch_multiplier_for_window = 100,  
                              save_output = TRUE, parallel = TRUE,
                              drop_missing_outcome = TRUE, 
                              forecast_with_full_fit = FALSE,
                              full_fit_threshold = NULL, inner_nfolds = 10, 
                              cores = NULL){
  
  len <- length(outcomes)
  
  if(parallel){
    # set up parallelization
    if(is.null(cores)){
      cores <- parallel::detectCores()-1
    }
    `%dopar%` <- foreach::`%dopar%`
    doParallel::registerDoParallel(cores = cores)
    foreach::getDoParWorkers() 
    
    results <- foreach::foreach(i = 1:len) %dopar% {
      
      outcome <- outcomes[i]
      
      cat("\n_____________ STARTING OUTCOME", outcome, "FOR ID", individual_id,"_____________\n")
      
      if(grepl("abpmean_lag5_mean", colnames(individual_data))){
        smooth_type <- "mean"
      } else if(grepl("abpmean_lag5_median", colnames(individual_data))){
        smooth_type <- "median"
      }
      
      load(paste0(file_path, outcome, "_historical_", smooth_type, ".Rdata"))
      
      if(grepl("AHE", outcome)){
        outcome_bounds <- c(0.001, 0.999)
      } else {
        outcome_bounds <- c(10, 160)
      }
      
      if(grepl("Y5", outcome)){
        horizon <- 5
      } else if(grepl("10", outcome)){
        horizon <- 10
      } else if(grepl("15", outcome)){
        horizon <- 15
      } else if(grepl("20", outcome)){
        horizon <- 20
      } else if(grepl("25", outcome)){
        horizon <- 25
      } else if(grepl("30", outcome)){
        horizon <- 30
      } 
      
      result <- run_slstream(
        individual_data = individual_data, individual_id = individual_id, 
        outcome = outcome, covariates = covariates, time = time, 
        cv_stack = cv_stack, slstream_weights_control = slstream_weights_control, 
        batch = batch, max_stop_time = max_stop_time, 
        individual_burn_in = individual_burn_in, historical_fit = historical_fit,
        batch_multiplier_for_window = batch_multiplier_for_window, 
        drop_missing_outcome = drop_missing_outcome, print_timers = print_timers,
        forecast_with_full_fit = forecast_with_full_fit,
        full_fit_threshold = full_fit_threshold, outcome_bounds = outcome_bounds,
        horizon = horizon, inner_nfolds = inner_nfolds
      )
      if(save_output){
        if(grepl("AHE", outcome)){
          outcome_clear <- paste0(outcome, "_", smooth_type)
        } else {
          outcome_clear <- outcome
        }
        results_path <- paste0(outcome_clear, "_id", individual_id, ".Rdata")
        save(result, file = paste0(file_path, results_path), compress = TRUE)
        result <- paste0(outcome, " FOR ID ", individual_id, " RESULTS SAVED!")
      }
      return(result)
    }
  } else {
    results <- list()
    for(i in 1:len){
      
      outcome <- outcomes[i]
      
      cat("\n_____________ STARTING OUTCOME", outcome, "FOR ID", individual_id,"_____________\n")
      
      if(any(grepl("abpmean_lag5_mean", colnames(individual_data)))){
        smooth_type <- "mean"
      } else if(any(grepl("abpmean_lag5_median", colnames(individual_data)))){
        smooth_type <- "median"
      }
      
      load(paste0(file_path, outcome, "_historical_", smooth_type, ".Rdata"))
      
      if(grepl("AHE", outcome)){
        outcome_bounds <- c(0.001, 0.999)
      } else {
        outcome_bounds <- c(10, 160)
      }
      
      if(grepl("Y5", outcome)){
        horizon <- 5
      } else if(grepl("10", outcome)){
        horizon <- 10
      } else if(grepl("15", outcome)){
        horizon <- 15
      } else if(grepl("20", outcome)){
        horizon <- 20
      } else if(grepl("25", outcome)){
        horizon <- 25
      } else if(grepl("30", outcome)){
        horizon <- 30
      } 
      
      result <- run_slstream(
        individual_data = individual_data, individual_id = individual_id, 
        outcome = outcome, covariates = covariates, time = time, 
        cv_stack = cv_stack, slstream_weights_control = slstream_weights_control, 
        batch = batch, max_stop_time = max_stop_time, 
        individual_burn_in = individual_burn_in, historical_fit = historical_fit,
        batch_multiplier_for_window = batch_multiplier_for_window, 
        drop_missing_outcome = drop_missing_outcome, print_timers = print_timers,
        forecast_with_full_fit = forecast_with_full_fit,
        full_fit_threshold = full_fit_threshold, outcome_bounds = outcome_bounds,
        horizon = horizon, inner_nfolds = inner_nfolds
      )
      if(save_output){
        if(grepl("AHE", outcome)){
          outcome_clear <- paste0(outcome, "_", smooth_type)
        } else {
          outcome_clear <- outcome
        }
        results_path <- paste0(outcome_clear, "_id", individual_id, ".Rdata")
        save(result, file = paste0(file_path, results_path), compress = TRUE)
        result <- paste0(outcome, " FOR ID ", individual_id, " RESULTS SAVED!")
      }
      results[[i]] <- result
    }
  }
  names(results) <- outcomes
  return(results)
}

run_slstream_allID <- function(multi_individual_data, ids, covariates, 
                               outcome, horizon, historical_fit = NULL, file_path, 
                               time = "time", cv_stack, slstream_weights_control, 
                               batch = 5, max_stop_time = 1440, 
                               individual_burn_in = 30, print_timers = TRUE,
                               batch_multiplier_for_window = 100, 
                               save_output = TRUE, parallel = TRUE,
                               drop_missing_outcome = TRUE, 
                               forecast_with_full_fit = FALSE,
                               full_fit_threshold = NULL, inner_nfolds = 10,
                               cores = NULL){
  
  if(any(grepl("abpmean_lag5_mean", colnames(multi_individual_data)))){
    smooth_type <- "mean"
  } else if(any(grepl("abpmean_lag5_median", colnames(multi_individual_data)))){
    smooth_type <- "median"
  }
  
  if(grepl("AHE", outcome)){
    outcome_clear <- paste0(outcome, smooth_type)
    outcome_bounds <- c(0.001, 0.999)
  } else {
    outcome_clear <- outcome
    outcome_bounds <- c(10, 160)
  }
  
  if(is.null(historical_fit)){
    historical_path <- paste0("historical_", smooth_type, ".Rdata")
    historical_fit_path <- paste0(file_path, outcome, "_", historical_path)
    load(historical_fit_path)
  }
  
  len <- length(ids)
  if(parallel){
    # set up parallelization
    if(is.null(cores)){
      cores <- parallel::detectCores()-1
    }
    `%dopar%` <- foreach::`%dopar%`
    doParallel::registerDoParallel(cores = cores)
    foreach::getDoParWorkers() 
    
    results <- foreach::foreach(i = 1:len) %dopar% {
      
      cat("\n_____________ STARTING ID", ids[i], "FOR OUTCOME", outcome,"_____________\n")
      
      individual_data <- multi_individual_data[id == ids[i],]
      setorder(individual_data, time_and_date)
      
      result <- run_slstream(
        individual_data = individual_data, individual_id = ids[i], 
        outcome = outcome, covariates = covariates, time = time, 
        cv_stack = cv_stack, slstream_weights_control = slstream_weights_control, 
        batch = batch, max_stop_time = max_stop_time, 
        individual_burn_in = individual_burn_in, historical_fit = historical_fit,
        batch_multiplier_for_window = batch_multiplier_for_window, 
        drop_missing_outcome = drop_missing_outcome, print_timers = print_timers,
        forecast_with_full_fit = forecast_with_full_fit,
        full_fit_threshold = full_fit_threshold, outcome_bounds = outcome_bounds,
        horizon = horizon, inner_nfolds = inner_nfolds
      )
      if(save_output){
        results_path <- paste0(outcome_clear, "_id", ids[i], ".Rdata")
        save(result, file=paste0(file_path, results_path), compress=T)
        result <- paste0(outcome, " FOR ID ", ids[i], " RESULTS SAVED!")
      } 
      return(result)
    }
  } else {
    results <- list()
    for(i in 1:len){
      
      cat("\n_____________ STARTING ID", ids[i], "FOR OUTCOME", outcome,"_____________\n")
      
      individual_data <- multi_individual_data[id == ids[i],]
      setorder(individual_data, time_and_date)
      
      result <- run_slstream(
        individual_data = individual_data, individual_id = ids[i], 
        outcome = outcome, covariates = covariates, time = time, 
        cv_stack = cv_stack, slstream_weights_control = slstream_weights_control, 
        batch = batch, max_stop_time = max_stop_time, 
        individual_burn_in = individual_burn_in, historical_fit = historical_fit,
        batch_multiplier_for_window = batch_multiplier_for_window, 
        drop_missing_outcome = drop_missing_outcome, print_timers = print_timers,
        forecast_with_full_fit = forecast_with_full_fit,
        full_fit_threshold = full_fit_threshold, outcome_bounds = outcome_bounds,
        horizon = horizon, inner_nfolds = inner_nfolds
      )
      if(save_output){
        results_path <- paste0(outcome_clear, "_id", ids[i], ".Rdata")
        save(result, file=paste0(file_path, results_path), compress=T)
        result <- paste0(outcome, " FOR ID ", ids[i], " RESULTS SAVED!")
      }
      results[[i]] <- result
    }
  }
  names(results) <- ids
  return(results)
}

run_slstream_catch_error <- function(individual_data, individual_id, outcome, 
                                     covariates, time, cv_stack, 
                                     slstream_weights_control, batch, 
                                     max_stop_time, individual_burn_in, 
                                     batch_multiplier_for_window, 
                                     historical_fit, drop_missing_outcome, 
                                     print_timers, forecast_with_full_fit, 
                                     full_fit_threshold, outcome_bounds, 
                                     horizon, inner_nfolds){
  tryCatch(
    run_slstream(
      individual_data = individual_data, individual_id = individual_id, 
      outcome = outcome, covariates = covariates, time = time, 
      cv_stack = cv_stack, batch = batch, max_stop_time = max_stop_time,
      slstream_weights_control =  slstream_weights_control, 
      individual_burn_in = individual_burn_in, historical_fit = historical_fit,
      batch_multiplier_for_window = batch_multiplier_for_window, 
      drop_missing_outcome = drop_missing_outcome, print_timers = print_timers,
      forecast_with_full_fit = forecast_with_full_fit,
      full_fit_threshold = full_fit_threshold, outcome_bounds = outcome_bounds,
      horizon = horizon, inner_nfolds = inner_nfolds
    ), 
    error = function(error_message) {
      message("Caught error.")
      message("Error message from R:")
      message(error_message)
      return(NULL)
    }
  )
}
