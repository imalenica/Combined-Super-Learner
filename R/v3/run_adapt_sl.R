run_slstream <- function(individual_data, individual_id, outcome, covariates, 
                         time = "min", cv_stack, slstream_weights_control, 
                         batch = 5, max_stop_time = 2160, initial_burn_in = 30, 
                         burn_in_multiplier_for_window_size = 10, 
                         check_baseline_covariates = TRUE, historical_fit,
                         drop_missing_outcome = TRUE){
  
  ############################### preliminary checks ###########################
  individual_data <- data.table::data.table(individual_data)
  
  if(length(unique(individual_data[["id"]])) > 1){
    individual_data$id <- as.character(individual_data$id)
    individual_data <- individual_data[id == individual_id,]
    setorder(individual_data, time_and_date)
  }
  
  if(check_baseline_covariates){
    W <- c("rank_icu", "sex", "age", "care_unit", "sapsi_first", "sofa_first", 
           "bmi", "admission_type_descr")
    if(any(is.na(individual_data[, W, with=FALSE]))){
      stop("NA baseline covariates")
    }
  }
  
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
  
  splits <- seq((initial_burn_in+(2*batch)), max_time, batch)
  train_data_list <- lapply(splits, function(x) individual_data[1:x, ])
  pred_data_list <- lapply(splits, function(x) individual_data[(x+1):(x+batch), ])
  
  # function to increase burn_in, & switch to rolling window (rw) folds
  window_size <- initial_burn_in*burn_in_multiplier_for_window_size
  n_first_rw <- (2*batch)+window_size # no. obs under first rw fit
  get_burn_in <- function(n) ifelse(n < n_first_rw, initial_burn_in, window_size)
  
  ################################## run #########################################
  
  loss_tables <- list()
  SL_tables <- list()
  individual_fits <- list()
  forecast_tables <- list()
  
  sl3_debug_mode()
  N <- length(splits)-1
  for(i in 1:N){
    
    # quick clean of results we no longer need
    if(i > 2){
      loss_tables[[i-2]] <- NA
      individual_fits[[i-2]] <- NA
    }
    
    # initiate "online" training and forecast data
    training_data <- data.table::data.table(train_data_list[[i]])
    forecast_data <- data.table::data.table(pred_data_list[[i]])
    burn_in <- get_burn_in(splits[[i]])
    
    if(i == 1){
      fit <- slstream(training_data = training_data, outcome = outcome, 
                      id = individual_id, covariates = covariates, time = time,
                      batch = batch, historical_fit = historical_fit, 
                      fold_fun = "folds_rolling_origin", first_window = burn_in, 
                      individual_stack = cv_stack, forecast_data = forecast_data, 
                      first_fit = TRUE, weights_control = slstream_weights_control)
    } else {
      if(burn_in == initial_burn_in) {
        fit <- slstream(training_data = training_data, outcome = outcome, 
                        id = individual_id, covariates = covariates, time = time,
                        batch = batch, historical_fit = historical_fit, 
                        fold_fun = "folds_rolling_origin", 
                        first_window = burn_in, forecast_data = forecast_data, 
                        weights_control = slstream_weights_control,
                        loss_table = loss_tables[[i-1]], 
                        previous_individual_fit = individual_fits[[i-1]])
      } else {
        if(nrow(training_data) > n_first_rw){
          # we can start updating the fits under rolling window cv after 1st one
          cv_stack <- NULL
        }
        fit <- slstream(training_data = training_data, outcome = outcome, 
                        id = individual_id, covariates = covariates, time = time,
                        historical_fit = historical_fit, batch = batch, 
                        fold_fun = "folds_rolling_window", window_size = burn_in,
                        individual_stack = cv_stack, forecast_data = forecast_data, 
                        first_fit = FALSE, loss_table = loss_tables[[i-1]], 
                        previous_individual_fit = individual_fits[[i-1]],
                        weights_control = slstream_weights_control)
      }
    }
    loss_tables[[i]] <- fit$loss_table
    SL_tables[[i]] <- fit$SL_table
    forecast_tables[[i]] <- fit$forecast_table
    individual_fits[[i]] <- fit$individual_fit
  }
  
  return_list <- list(SL_table = rbindlist(SL_tables, fill=T), 
                      loss_table = loss_tables[[N]],
                      forecast_table = rbindlist(forecast_tables, fill=T))
  
  return(return_list)
}

run_slstream_catch_error <- function(individual_data, individual_id, outcome, 
                                     covariates, time, cv_stack, 
                                     slstream_weights_control, batch, 
                                     max_stop_time, initial_burn_in, 
                                     burn_in_multiplier_for_window_size, 
                                     check_baseline_covariates, 
                                     historical_fit){
  out <- tryCatch(
    run_slstream(
      individual_data = individual_data, individual_id = individual_id, 
      outcome = outcome, covariates = covariates, time = time, 
      cv_stack = cv_stack, batch = batch, max_stop_time = max_stop_time,
      slstream_weights_control =  slstream_weights_control, 
      initial_burn_in = initial_burn_in, historical_fit = historical_fit,
      burn_in_multiplier_for_window_size = burn_in_multiplier_for_window_size, 
      check_baseline_covariates = check_baseline_covariates), 
    error = function(e) NULL
  )
  return(out)
}

run_slstream_allY <- function(outcomes, historical_fit_list, file_path,
                              individual_data = NULL, individual_id, covariates, 
                              time = "min", cv_stack, slstream_weights_control, 
                              batch = 5, max_stop_time = 2160, 
                              initial_burn_in = 30, 
                              burn_in_multiplier_for_window_size = 10, 
                              check_baseline_covariates = TRUE, cores = NULL,
                              save_output = FALSE){
  
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
  
  # set up parallelization
  if(is.null(cores)){
    cores <- parallel::detectCores()-1
  }
  
  len <- length(outcomes)
  `%dopar%` <- foreach::`%dopar%`
  doParallel::registerDoParallel(cores = cores)
  foreach::getDoParWorkers() 
  id_specific_results <- foreach::foreach(i = 1:len) %dopar% {
    
    outcome <- outcomes[i]
    print(paste0("Starting outcome ", outcome, " for id ", individual_id))
    
    if(is.null(individual_data)){
      if(grepl("mean", outcome)){
        individual_data <- individual_mean
      } else if(grepl("median", outcome)){
        individual_data <- individual_median
      }
      rm(list = c("individual_mean", "individual_median"))
      setorder(individual_data, min)
    }
    
    run_slstream_catch_error(
      individual_data = individual_data, individual_id = individual_id, 
      outcome = outcome, covariates = covariates, time = time, 
      cv_stack = cv_stack, slstream_weights_control = slstream_weights_control, 
      batch = batch, max_stop_time = max_stop_time, 
      initial_burn_in = initial_burn_in, historical_fit = historical_fit_list[[i]],
      burn_in_multiplier_for_window_size = burn_in_multiplier_for_window_size, 
      check_baseline_covariates = check_baseline_covariates
    )
  }
  if(is.null(names(outcome_list))){
    list_names <- unlist(outcome_list)
  } else {
    list_names <- names(outcome_list)
  }
  names(id_specific_results) <- list_names
  
  if(save_output){
    results_path <- paste0("id", individual_id, ".Rdata")
    save(id_specific_results, file=paste0(file_path, results_path), compress=T)
    msg <- paste0("id ", id, " results saved")
    return(msg)
  } else {
    return(id_specific_results)
  }
}

run_slstream_allID <- function(multi_individual_data, ids, covariates, 
                               outcome, historical_fit = NULL, file_path, 
                               time = "min", cv_stack, slstream_weights_control, 
                               batch = 5, max_stop_time = 2160, 
                               initial_burn_in = 30, 
                               burn_in_multiplier_for_window_size = 10, 
                               check_baseline_covariates = TRUE,
                               cores = NULL, save_output = FALSE){
  
  if(is.null(historical_fit)){
    if(grepl("mean", outcome)){
      smooth <- "mean"
    } else if(grepl("median", outcome)){
      smooth <- "median"
    }
    historical_path <- paste0("historical_", smooth, ".Rdata")
    historical_fit_path <- paste0(file_path, outcome, "_", historical_path)
    load(historical_fit_path)
    historical_fit <- fit
    rm(fit)
  }
  
  # set up parallelization
  if(is.null(cores)){
    cores <- parallel::detectCores()-1
  }
  
  len <- length(ids)
  `%dopar%` <- foreach::`%dopar%`
  doParallel::registerDoParallel(cores = cores)
  foreach::getDoParWorkers() 
  Y_specific_results <- foreach::foreach(i = 1:len) %dopar% {
    individual_id <- ids[i]
    
    print(paste0("Starting id ", individual_id, " for outcome ", outcome))
    
    individual_data <- multi_individual_data[id == individual_id,]
    setorder(individual_data, time_and_date)
    
    run_slstream_catch_error(
      individual_data = individual_data, individual_id = individual_id, 
      outcome = outcome, covariates = covariates, time = time, 
      cv_stack = cv_stack, slstream_weights_control = slstream_weights_control, 
      batch = batch, max_stop_time = max_stop_time, 
      initial_burn_in = initial_burn_in, historical_fit = historical_fit,
      burn_in_multiplier_for_window_size = burn_in_multiplier_for_window_size, 
      check_baseline_covariates = check_baseline_covariates
    )
  }
  names(Y_specific_results) <- ids
  
  if(save_output){
    results_path <- paste0(outcome, ".Rdata")
    save(Y_specific_results, file=paste0(file_path, results_path), compress=T)
    msg <- paste0("Y ", outcome, " results saved")
    return(msg)
  } else {
    return(Y_specific_results)
  }
}