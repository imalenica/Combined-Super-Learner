################################################################################
# utility functions
################################################################################

get_weights <- function(pred, observed, loss, convex = FALSE, discrete = FALSE){
  
  if(discrete == TRUE) {
    fit_coef <- as.numeric(loss)
    fit_coef[which.min(loss)] <- 1
    fit_coef[-(which.min(loss))] <- 0
  } else if(discrete == FALSE & convex == TRUE) {
    fit_coef <- lsei::pnnls(as.matrix(pred), as.matrix(observed), sum = 1)
    fit_coef <- fit_coef$x
  } else if(discrete == FALSE & convex == FALSE) {
    fit_coef <- nnls::nnls(as.matrix(pred), as.numeric(observed))
    fit_coef <- fit_coef$x
  }
  return(fit_coef)
}

################################################################################
# historical SLs
################################################################################

# historical_data: data.table of all observed data on all previously seen subjects
# outcome: name of outcome variable which is found in historical data.table
# covariates: vector of covariate names which are column names in the data.table  
# historical_stack: stack of sl3 learners to model historical data  
# id: column name indicating unique subject identifier for appropraite V-fold cv

make_historical_fit <- function(historical_data, outcome, covariates, id,
                                historical_stack){
  task <- make_sl3_Task(
    data = historical_data, 
    covariates = covariates,
    outcome = outcome,
    id = id,
    drop_missing_outcome = T
  )
  cv_stack <- Lrnr_cv$new(historical_stack, full_fit = TRUE)
  fit <- cv_stack$train(task)
  chained_task <- fit$chain(task)
  
  metalearner_nnls <- make_learner(Lrnr_nnls)
  nnls_fit <- metalearner_nnls$train(chained_task)
  sl_nnls <- make_learner(Pipeline, fit, nnls_fit)
  
  metalearner_nnls_convex <- make_learner(Lrnr_nnls, convex = TRUE)
  nnls_fit_convex <- metalearner_nnls_convex$train(chained_task)
  sl_nnls_convex <- make_learner(Pipeline, fit, nnls_fit_convex)
  
  metalearner_discrete <- make_learner(Lrnr_cv_selector)
  discrete_fit <- metalearner_discrete$train(chained_task)
  sl_discrete <- make_learner(Pipeline, fit, discrete_fit)
  
  fits <- list(cv_fit = fit,
               sl_nnls = sl_nnls,
               sl_nnls_convex = sl_nnls_convex,
               sl_discrete = sl_discrete)
  return(fits)
}


################################################################################
# adaptive (individual + historical) individualized online super learner
################################################################################

# individual_training_data: data.table with outcome observed
# indiviual_forecast_data: data.table with outcome unobserved and covariates observed
# outcome: name of outcome variable which is found in both data.tables 
# covariates: vector of covariate names which are columns in both data.tables  
# subject_id: unique subject identifier for patient 
# historical_fit: trained model on historical data
# individual_stack: optional stack of sl3 learners to model individual data  
# past_individual_fit: optional past sl3 fit to update with new data

make_adapt_sl <- function(individual_training_data, indiviual_forecast_data, 
                          outcome, covariates, subject_id, historical_fit, 
                          individual_stack = NULL, past_individual_fit = NULL){
  
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
  
  ################################## train SLs #################################
  
  # make folds with training data
  folds <- origami::make_folds(individual_training_data, 
                               fold_fun = folds_rolling_origin,
                               first_window = 5, 
                               validation_size = 5, 
                               gap = 0,
                               batch = 5
                               )

  
  training_task <- make_sl3_Task(
    data = individual_training_data, 
    covariates = covariates,
    outcome = outcome, 
    folds = folds
    )
  
  # fit initial superlearner if past_individual_fit is not provided
  if(is.null(past_individual_fit) & is.null(individual_stack)) {
    print("Error: Individual fit and individual stack both missing. Provide       
           individual_stack to initialize a new fit, and past_individual_fit
           to update with new data.")
    stop()
 
    } else if(is.null(past_individual_fit) & !(is.null(individual_stack))) {
    
    print(paste0("Training learners for subject ", subject_id, 
                 " with ", nrow(individual_training_data), " observations."))
    
    # fit individualized learners
    cv_stack <- Lrnr_cv$new(individual_stack)
    ind_fit <- cv_stack$train(training_task)
    
    } else if(!is.null(past_individual_fit) & is.null(individual_stack)) {
    print(paste0("Updating past fit for subject ", subject_id, " with ", 
                 nrow(individual_training_data), " observations."))
    
    # update individualized learners
    ind_fit <- past_individual_fit$update(training_task)
    
    } else if(!is.null(past_individual_fit) & !is.null(individual_stack)) {
    print(paste0("Ignoring past fit and training individual_stack for subject ",
                 subject_id, " with ", nrow(individual_training_data), 
                 " observations."))
  
    # fit individualized learners
    cv_stack <- Lrnr_cv$new(individual_stack)
    ind_fit <- cv_stack$train(training_task)
    }
  
  # predict with individualized learners and historical learners
  ind_preds <- ind_fit$predict(training_task)
  
  hist_preds <- bind_rows(lapply(folds, function(fold) {
    test_set_in_training_task <- validation_task(training_task, fold)
    hist_fits <- cbind.data.frame(lapply(historical_fit, function(fit){
      fit$predict_fold(test_set_in_training_task, "full")
    }))
  }))
  
  # combine predictions
  learners <- c(paste0("individual_", colnames(ind_preds)), 
                paste0("historical_", colnames(hist_preds)))
  training_preds <- cbind.data.frame(ind_preds, hist_preds)
  names(training_preds) <- learners

  # get true Y
  truth <- lapply(folds, function(i) data.frame(individual_training_data[i$validation,
                                                                         get(outcome)]))
  truth <- bind_rows(truth)[,1]
  
  # evaluate empirical loss for training
  loss <- apply(training_preds, 2, function(pred) mean((pred-truth)^2))
  
  # establish various super learners
  weights_discrete <- get_weights(training_preds, truth, loss, convex, discrete = T)
  weights_nnls_convex <- get_weights(training_preds, truth, loss, convex = T, discrete = F)
  weights_nnls <- get_weights(training_preds, truth, loss)
  sl_weights <- rbind(weights_discrete, weights_nnls_convex, weights_nnls)
  colnames(sl_weights) <- names(training_preds)
  
  ############################ forecast with SLs ###############################
  
  # avoid issues with NA as outcome when outcome hasn't yet been observed
  unknown_outcome <- rep(0, nrow(indiviual_forecast_data))
  forecast_data <- data.table(unknown_outcome, indiviual_forecast_data)
  
  # obtain forecast with individual fit and historical fit
  forecast_fold <- origami:: make_folds(forecast_data, fold_fun = folds_vfold, 
                                        V = 1)
  forecast_task <- make_sl3_Task(
    data = forecast_data, 
    covariates = covariates,
    outcome = "unknown_outcome", 
    folds = forecast_fold
  )
  ind_forecast <- ind_fit$predict(forecast_task)
  hist_forecast <- cbind.data.frame(lapply(historical_fit, function(fit){
      fit$predict(forecast_task)
    }))
  
  # combine forecasts
  learners <- c(paste0("individual_", colnames(ind_forecast)), 
                paste0("historical_", colnames(hist_forecast)))
  forecast <- cbind.data.frame(ind_forecast, hist_forecast)
  names(forecast) <- learners
  
  # use sl weights for forecast with sl
  forecast_discrete <- as.matrix(forecast) %*% weights_discrete
  forecast_nnls_convex <- as.matrix(forecast) %*% weights_nnls_convex
  forecast_nnls <- as.matrix(forecast) %*% weights_nnls
  sl_forecasts <- data.frame(forecast_discrete, forecast_nnls_convex, forecast_nnls)
  colnames(sl_forecasts) <- c("discreteSL", "nnls_convexSL", "nnls_SL")
  
  return_list <- list(sl_forecasts = sl_forecasts,
                      sl_weights = sl_weights, 
                      individual_fit = ind_fit, 
                      historical_fit = historical_fit,
                      training_preds = training_preds,
                      loss = loss)
  
  return(return_list)
}
