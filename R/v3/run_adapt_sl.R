run_slstream <- function(individual_data, individual_id, outcome, covariates, 
                         time = "min", cv_stack, slstream_weights_control, 
                         batch = 5, max_stop_time = 1440, initial_burn_in = 30, 
                         burn_in_multiplier_for_window_size = 10, historical_fit,
                         drop_missing_outcome = TRUE, print_timers = TRUE,
                         forecast_with_full_fit = FALSE, 
                         full_fit_threshold = NULL, outcome_bounds = c(10,160), 
                         horizon){
  
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
  
  ########################### set up pseudo-online data ########################
  max_times <- c(max_stop_time, nrow(individual_data))
  max_time <- max_times[which.min(max_times)]
  if(max_time < 100) {
    stop(paste0("Insufficient time for subject :", individual_id))
  }
  
  gap <- horizon-batch

  splits <- seq((initial_burn_in + 2*batch + gap), max_time, batch)
  train_data_list <- lapply(splits, function(x) individual_data[1:x, ])
  pred_data_list <- lapply(splits, function(x) individual_data[(x+gap+1):(x+batch+gap), ])
  
  not_na <- unlist(lapply(pred_data_list, function(x) !any(is.na(x[[outcome]]))))
  pred_data_list <- pred_data_list[not_na]
  train_data_list <- train_data_list[not_na]
  stopifnot(length(pred_data_list) == length(train_data_list))
  N <- length(pred_data_list)
  
  # function to increase burn_in, & switch to rolling window (rw) folds
  window_size <- initial_burn_in * burn_in_multiplier_for_window_size
  n_first_rolling_window_cv <- 2 * batch + window_size + gap # no. obs for 1st rw fit
  get_burn_in <- function(n){
    ifelse(n < n_first_rolling_window_cv, initial_burn_in, window_size)
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
  fold_fits <- list()
  screener_selected_covs <- list()
  for(i in 1:N){
    
    cat("\n_____ starting run", i, "of", N, "for id", individual_id, "_____\n")
    
    # quick clean of results we no longer need
    if(i > 2){
      loss_tables[[i-2]] <- NA
      individual_fits[[i-2]] <- NA
      fold_fits[[i-2]] <- NA
    }
    
    # initiate "online" training and forecast data
    training_data <- data.table(train_data_list[[i]])
    forecast_data <- data.table(pred_data_list[[i]])
    burn_in <- get_burn_in(splits[[i]])
    
    if(i == 1){
      fit <- slstream(training_data = training_data, outcome = outcome, 
                      id = individual_id, covariates = covariates, time = time, 
                      historical_fit = historical_fit, batch = batch, 
                      fold_fun = "folds_rolling_origin", first_window = burn_in, 
                      individual_stack = cv_stack, forecast_data = forecast_data, 
                      print_timers = print_timers, first_fit = TRUE, 
                      weights_control = slstream_weights_control,
                      forecast_with_full_fit = forecast_with_full_fit,
                      full_fit_threshold = full_fit_threshold, 
                      outcome_bounds = outcome_bounds, horizon = horizon)
    } else {
      if(burn_in == initial_burn_in) { 
        # continue with rolling origin cross-validation & update fits
        fit <- slstream(training_data = training_data, outcome = outcome, 
                        id = individual_id, covariates = covariates, time = time, 
                        historical_fit = historical_fit, batch = batch,
                        fold_fun = "folds_rolling_origin", first_window = burn_in, 
                        individual_stack = NULL, forecast_data = forecast_data, 
                        print_timers = print_timers, loss_table = loss_tables[[i-1]],
                        weights_control = slstream_weights_control, 
                        previous_fold_fits = fold_fits[[i-1]],
                        previous_individual_fit = individual_fits[[i-1]],
                        forecast_with_full_fit = forecast_with_full_fit,
                        full_fit_threshold = full_fit_threshold, 
                        outcome_bounds = outcome_bounds, horizon = horizon)
      } else {
        # switch to rolling window cross-validation
        if(nrow(training_data) <= n_first_rolling_window_cv){
          # instantiate new non-updated individual fits by providing stack
          fit <- slstream(training_data = training_data, outcome = outcome, 
                          id = individual_id, covariates = covariates, time = time, 
                          historical_fit = historical_fit, batch = batch, 
                          fold_fun = "folds_rolling_window", window_size = burn_in, 
                          individual_stack = cv_stack, forecast_data = forecast_data, 
                          print_timers = print_timers, loss_table = loss_tables[[i-1]], 
                          weights_control = slstream_weights_control,
                          forecast_with_full_fit = forecast_with_full_fit,
                          full_fit_threshold = full_fit_threshold,
                          outcome_bounds = outcome_bounds, horizon = horizon)
        } else {
          # now we can update the rolling window fits by providing NULL stack
          fit <- slstream(training_data = training_data, outcome = outcome, 
                          id = individual_id, covariates = covariates, time = time, 
                          historical_fit = historical_fit, batch = batch, 
                          fold_fun = "folds_rolling_window", window_size = burn_in, 
                          individual_stack = NULL, forecast_data = forecast_data, 
                          print_timers = print_timers, loss_table = loss_tables[[i-1]],
                          weights_control = slstream_weights_control,
                          previous_fold_fits = fold_fits[[i-1]],
                          previous_individual_fit = individual_fits[[i-1]],
                          forecast_with_full_fit = forecast_with_full_fit,
                          full_fit_threshold = full_fit_threshold,
                          outcome_bounds = outcome_bounds, horizon = horizon)
        }
      }
    }
    loss_tables[[i]] <- fit$loss_table
    sl_tables[[i]] <- fit$sl_table
    forecast_tables[[i]] <- fit$forecast_table
    individual_fits[[i]] <- fit$individual_fit
    fold_fits[[i]] <- fit$fold_fits
    oSL_fit <- fit$individual_fit$fit_object$fold_fits[[length(fit$individual_fit$fit_object$folds)]]
    screener_fit <- oSL_fit$fit_object$learner_fits$lasso_screen$fit_object$learner_fits$Lrnr_screener_coefs_0_NULL
    screener_selected_covs[[i]] <- screener_fit$fit_object$selected
  }
  
  ############################## assess performance ############################
  
  ##### summarize prop of individual vs. historical learner weights
  SLtbl <- rbindlist(sl_tables, fill = T)
  wts <- SLtbl[, -c("id", "min_start", "min_end", "onlineSL"), with = F]
  ind_lrnrs <- colnames(wts)[grep("individual_", colnames(wts))]
  hist_lrnrs <- colnames(wts)[grep("historical_", colnames(wts))]
  p_ind <- rowSums(wts[, ind_lrnrs, with = F], na.rm = T) / rowSums(wts, na.rm = T)
  p_hist <- rowSums(wts[, hist_lrnrs, with = F], na.rm = T) / rowSums(wts, na.rm = T)
  SLtbl <- data.table(prop_individual = p_ind, prop_historical = p_hist, SLtbl)
  
  ##### make table of true horizon-ahead outcomes & forecasts
  time <- individual_data[[time]] # time in which forecast was made 
  obs <- individual_data[[outcome]] # the true horizon-ahead outcome
  outcome_time <- time + horizon # time of horizon-ahead outcome
  obs_tbl <- data.table(obs, time, outcome_time)
  # merge with forecasts
  forecast_tbl <- rbindlist(forecast_tables, fill = T)
  tbl <- merge(obs_tbl, forecast_tbl, by = "time", all.x = F, all.y = T)
  
  ##### for each learner, calculate time-specific & summarized MSE
  time_cols <- colnames(forecast_tbl)[grep("time", colnames(forecast_tbl))]
  preds_tbl <- forecast_tbl[, -time_cols, with = F]
  MSE_tbl <- data.table(apply(preds_tbl, 2, function(pred) (pred-tbl[["obs"]])^2))
  MSEsummary <- c(id = individual_id, outcome = outcome, colMeans(MSE_tbl, na.rm=T))
  
  ##### for each learner, calculate the number of non-NA
  NAsummary <- apply(MSE_tbl, 2, function(res) sum(is.na(res)))
  
  # re-add time cols to MSE_tbl
  time_tbl <- forecast_tbl[, time_cols, with = F]
  MSE_tbl <- data.table(time_tbl, MSE_tbl)
  
  return(list(sl_table = SLtbl, loss_table = loss_tables[[N]], 
              forecast_table = forecast_tbl, NA_summary = NAsummary,
              MSE_table = MSE_tbl, MSE_summary = MSEsummary, 
              screener_selected_covs = screener_selected_covs))
}

run_slstream_allY <- function(outcomes, historical_fit_list, file_path,
                              individual_data = NULL, individual_id, covariates, 
                              time = "min", cv_stack, slstream_weights_control, 
                              batch = 5, max_stop_time = 1440, 
                              initial_burn_in = 30, print_timers = TRUE,
                              burn_in_multiplier_for_window_size = 10,  
                              save_output = TRUE, parallel = TRUE,
                              drop_missing_outcome = TRUE, 
                              forecast_with_full_fit = FALSE,
                              full_fit_threshold = 315,  
                              outcome_bounds = c(10,160), horizon){
  
  if(is.null(individual_data)){
    load(paste0(file_path, "individual_mean.Rdata"))
    individual_mean <- individual
    individual_mean$id <- as.character(individual_mean$id)
    individual_mean <- individual_mean[id == individual_id,]
    
    load(paste0(file_path, "individual_median.Rdata"))
    individual_median <- individual
    individual_median$id <- as.character(individual_median$id)
    individual_median <- individual_median[id == individual_id,]
    
    rm(individual)
  }
  
  len <- length(outcomes)
  if(parallel){
    # set up parallelization
    cores <- parallel::detectCores()-1
    `%dopar%` <- foreach::`%dopar%`
    doParallel::registerDoParallel(cores = cores)
    foreach::getDoParWorkers() 
    
    results <- foreach::foreach(i = 1:len) %dopar% {
      
      outcome <- outcomes[i]
      cat("\n_____________ STARTING OUTCOME", outcome, "FOR ID", individual_id,"_____________\n")
      
      if(is.null(individual_data)){
        if(grepl("mean", outcome)){
          individual_data <- individual_mean
        } else if(grepl("median", outcome)){
          individual_data <- individual_median
        }
        rm(list = c("individual_mean", "individual_median"))
        setorder(individual_data, min)
      }
      
      result <- run_slstream(
        individual_data = individual_data, individual_id = individual_id, 
        outcome = outcome, covariates = covariates, time = time, 
        cv_stack = cv_stack, slstream_weights_control = slstream_weights_control, 
        batch = batch, max_stop_time = max_stop_time, 
        initial_burn_in = initial_burn_in, historical_fit = historical_fit_list[[i]],
        burn_in_multiplier_for_window_size = burn_in_multiplier_for_window_size, 
        drop_missing_outcome = drop_missing_outcome, print_timers = print_timers,
        forecast_with_full_fit = forecast_with_full_fit,
        full_fit_threshold = full_fit_threshold, outcome_bounds = outcome_bounds,
        horizon = horizon
      )
      if(save_output){
        results_path <- paste0(outcome, "_id", individual_id, ".Rdata")
        save(result, file=paste0(file_path, results_path), compress=T)
        result <- paste0(outcome, " id ", individual_id, " results saved")
      }
      return(result)
    }
  } else {
    results <- list()
    for(i in 1:len){
      
      outcome <- outcomes[i]
      cat("\n_____________ STARTING OUTCOME", outcome, "FOR ID", individual_id,"_____________\n")
      
      if(is.null(individual_data)){
        if(grepl("mean", outcome)){
          individual_data <- individual_mean
        } else if(grepl("median", outcome)){
          individual_data <- individual_median
        }
        rm(list = c("individual_mean", "individual_median"))
        setorder(individual_data, min)
      }
      
      result <- run_slstream(
        individual_data = individual_data, individual_id = individual_id, 
        outcome = outcome, covariates = covariates, time = time, 
        cv_stack = cv_stack, slstream_weights_control = slstream_weights_control, 
        batch = batch, max_stop_time = max_stop_time, 
        initial_burn_in = initial_burn_in, historical_fit = historical_fit_list[[i]],
        burn_in_multiplier_for_window_size = burn_in_multiplier_for_window_size, 
        drop_missing_outcome = drop_missing_outcome, print_timers = print_timers,
        forecast_with_full_fit = forecast_with_full_fit,
        full_fit_threshold = full_fit_threshold, outcome_bounds = outcome_bounds,
        horizon = horizon
      )
      if(save_output){
        results_path <- paste0(outcome, "_id", individual_id, ".Rdata")
        save(result, file=paste0(file_path, results_path), compress=T)
        msg <- paste0(outcome, " id ", individual_id, " results saved")
        results[[i]] <- msg
      } else {
        results[[i]] <- result
      }
    }
  }
  
  if(is.null(names(outcome_list))){
    list_names <- unlist(outcome_list)
  } else {
    list_names <- names(outcome_list)
  }
  names(results) <- list_names
  return(results)
}

run_slstream_allID <- function(multi_individual_data, ids, covariates, 
                               outcome, historical_fit = NULL, file_path, 
                               time = "min", cv_stack, slstream_weights_control, 
                               batch = 5, max_stop_time = 1440, 
                               initial_burn_in = 30, print_timers = TRUE,
                               burn_in_multiplier_for_window_size = 10, 
                               save_output = TRUE, parallel = TRUE,
                               drop_missing_outcome = TRUE, 
                               forecast_with_full_fit = FALSE,
                               full_fit_threshold = 315,  
                               outcome_bounds = c(10,160), horizon){
  
  if(is.null(historical_fit)){
    if(grepl("mean", outcome)){
      smooth_type <- "mean"
    } else if(grepl("median", outcome)){
      smooth_type <- "median"
    }
    historical_path <- paste0("historical_", smooth_type, ".Rdata")
    historical_fit_path <- paste0(file_path, outcome, "_", historical_path)
    load(historical_fit_path)
    historical_fit <- fit
    rm(fit)
  }
  
  len <- length(ids)
  if(parallel){
    # set up parallelization
    cores <- parallel::detectCores()-1
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
        initial_burn_in = initial_burn_in, historical_fit = historical_fit,
        burn_in_multiplier_for_window_size = burn_in_multiplier_for_window_size, 
        drop_missing_outcome = drop_missing_outcome, print_timers = print_timers,
        forecast_with_full_fit = forecast_with_full_fit,
        full_fit_threshold = full_fit_threshold, outcome_bounds = outcome_bounds,
        horizon = horizon
      )
      if(save_output){
        results_path <- paste0(outcome, "_id", ids[i], ".Rdata")
        save(result, file=paste0(file_path, results_path), compress=T)
        result <- paste0(outcome, " id ", ids[i], " results saved")
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
        initial_burn_in = initial_burn_in, historical_fit = historical_fit,
        burn_in_multiplier_for_window_size = burn_in_multiplier_for_window_size, 
        drop_missing_outcome = drop_missing_outcome, print_timers = print_timers,
        forecast_with_full_fit = forecast_with_full_fit,
        full_fit_threshold = full_fit_threshold, outcome_bounds = outcome_bounds,
        horizon = horizon
      )
      if(save_output){
        results_path <- paste0(outcome, "_id", ids[i], ".Rdata")
        save(result, file=paste0(file_path, results_path), compress=T)
        result <- paste0(outcome, " id ", ids[i], " results saved")
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
                                     max_stop_time, initial_burn_in, 
                                     burn_in_multiplier_for_window_size, 
                                     print_timers, historical_fit, 
                                     drop_missing_outcome, forecast_with_full_fit,
                                     full_fit_threshold, outcome_bounds, horizon){
  out <- tryCatch(
    run_slstream(
      individual_data = individual_data, individual_id = individual_id, 
      outcome = outcome, covariates = covariates, time = time, 
      cv_stack = cv_stack, batch = batch, max_stop_time = max_stop_time,
      slstream_weights_control =  slstream_weights_control, 
      initial_burn_in = initial_burn_in, historical_fit = historical_fit,
      burn_in_multiplier_for_window_size = burn_in_multiplier_for_window_size, 
      drop_missing_outcome = drop_missing_outcome, print_timers = print_timers,
      forecast_with_full_fit = forecast_with_full_fit,
      full_fit_threshold = full_fit_threshold, outcome_bounds = outcome_bounds,
      horizon = horizon
    ), error = function(e) NULL
  )
  return(out)
}
