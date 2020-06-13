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
                          outcome, covariates, historical_fit, id, 
                          return_individual_fit=FALSE,
                          individual_stack = NULL, past_individual_fit = NULL, 
                          acf = FALSE, parallelize = FALSE, 
                          cpus_logical = NULL) {
  
  if(parallelize & is.null(cpus_logical)){
    print("Cannot parallelize when no `cpus_logical` provided")
    parallelize <- FALSE
  }
  
  subject_id <- unique(individual_training_data[,id])
  
  ################################## train SLs #################################
  
  # make folds with training data
  folds <- origami::make_folds(data.table(individual_training_data), 
                               fold_fun = folds_rolling_origin,
                               first_window = 5, 
                               validation_size = 5, 
                               gap = 0,
                               batch = 5
                               )

  training_task <- make_sl3_Task(
    data = data.table(individual_training_data), 
    covariates = covariates,
    outcome = outcome, 
    folds = folds,
    id = id
    )
  # issue with mismatch between historical and individual delta column:
  training_task <- process_task(
    individual_training_task = training_task, 
    historical_task = historical_fit$cv_fit$fit_object$full_fit$training_task)

  # fit initial superlearner if past_individual_fit is not provided
  if(is.null(past_individual_fit) & is.null(individual_stack)) {
    print("Error: Individual fit and individual stack both missing. Provide       
           individual_stack to initialize a new fit, and past_individual_fit
           to update with new data.")
    stop()
 
    }else if(is.null(past_individual_fit) & !(is.null(individual_stack))) {
    
    print(paste0("Training learners for subject ", subject_id, 
                 " with ", nrow(individual_training_data), " observations."))
    
    # fit individualized learners
    cv_stack <- Lrnr_cv$new(individual_stack)
    if(parallelize){
      plan(multicore, workers = cpus_logical)
      test <- delayed_learner_train(cv_stack, training_task)
      sched <- Scheduler$new(test, FutureJob, nworkers = cpus_logical,
                             verbose = FALSE)
      ind_fit <- sched$compute()
    } else {
      ind_fit <- cv_stack$train(training_task)
    }
    
    }else if(!is.null(past_individual_fit) & is.null(individual_stack)) {
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
      if(parallelize){
        plan(multicore, workers = cpus_logical)
        test <- delayed_learner_train(cv_stack, training_task)
        sched <- Scheduler$new(test, FutureJob, nworkers = cpus_logical,
                               verbose = FALSE)
        ind_fit <- sched$compute()
      } else {
        ind_fit <- cv_stack$train(training_task)
      }
    }
  
  # predict with individualized learners and historical learners
  ind_preds <- ind_fit$predict(training_task)
  
  if(class(historical_fit) == "list"){
    #Omit the task
    historical_fits <- historical_fit[1:(length(historical_fit)-1)]
    
    hist_preds <- bind_rows(lapply(folds, function(fold) {
      test_set_in_training_task <- validation_task(training_task, fold)
      hist_fits <- cbind.data.frame(lapply(historical_fits, function(fit){
        fit$predict_fold(test_set_in_training_task, "full")
      }))
    }))
  } else {
    hist_preds <- bind_rows(lapply(folds, function(fold) {
      test_set_in_training_task <- validation_task(training_task, fold)
      historical_fit$predict_fold(test_set_in_training_task, "full")
      }))
  }
  # historical_fit$fit_object$full_fit$learner_fits$Lrnr_screener_coefs_0.1$fit_object$selected
  # combine predictions
  learners <- c(paste0("historical_", colnames(hist_preds)),
                paste0("individual_", colnames(ind_preds)))
  training_preds <- cbind.data.frame(hist_preds, ind_preds)
  names(training_preds) <- learners

  # get true Y
  truth <- lapply(folds, function(i) data.frame(individual_training_data[i$validation,
                                                                         get(outcome)]))
  truth <- bind_rows(truth)[,1]
  
  # evaluate empirical loss for training
  loss <- apply(training_preds, 2, function(pred) mean((pred-truth)^2))
  
  # establish various super learners
  weights_discrete <- suppressWarnings(get_weights(training_preds, truth, 
                                                   loss, convex, discrete = T))
  weights_nnls_convex <- suppressWarnings(get_weights(training_preds, truth, 
                                                      loss, convex = T, discrete = F))
  weights_nnls <- suppressWarnings(get_weights(training_preds, truth, loss))
  sl_weights <- suppressWarnings(rbind(weights_discrete, weights_nnls_convex, weights_nnls))
  colnames(sl_weights) <- names(training_preds)
  
  # use sl weights for prediction with sl
  pred_discrete <- as.matrix(training_preds) %*% weights_discrete
  pred_nnls_convex <- as.matrix(training_preds) %*% weights_nnls_convex
  pred_nnls <- as.matrix(training_preds) %*% weights_nnls
  sl_pred <- data.frame(pred_discrete, pred_nnls_convex, pred_nnls)
  colnames(sl_pred) <- c("discreteSL", "nnls_convexSL", "nnls_SL")
  
  # check dependence (option to use get_acf function)
  if(acf){
    # get the residuals
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
  forecast_fold <- origami:: make_folds(forecast_data, fold_fun = folds_vfold, 
                                        V = 1)
  forecast_task <- make_sl3_Task(
    data = data.table(forecast_data), 
    covariates = covariates,
    outcome = "unknown_outcome", 
    folds = forecast_fold,
    id = id
  )
  
  ind_forecast <- ind_fit$predict(forecast_task)
  historical_fits <- historical_fit[1:4]
  hist_forecast <- cbind.data.frame(lapply(historical_fits, function(fit){
      fit$predict(forecast_task)
    }))
  
  # combine forecasts
  learners <- c(paste0("historical_", colnames(hist_forecast)),
                paste0("individual_", colnames(ind_forecast)))
  forecast <- cbind.data.frame(hist_forecast, ind_forecast)
  names(forecast) <- learners
  
  # use sl weights for forecast with sl
  forecast_discrete <- as.matrix(forecast) %*% weights_discrete
  forecast_nnls_convex <- as.matrix(forecast) %*% weights_nnls_convex
  forecast_nnls <- as.matrix(forecast) %*% weights_nnls
  sl_forecasts <- data.frame(forecast_discrete, forecast_nnls_convex, forecast_nnls)
  colnames(sl_forecasts) <- c("discreteSL", "nnls_convexSL", "nnls_SL")
  
  ############################ just individual SL ###############################
  if(return_individual_fit){
    # evaluate empirical loss just for individual learners
    ind_loss <- apply(ind_preds, 2, function(pred) mean((pred-truth)^2))
    
    # establish various super learners
    ind_weights_discrete <- suppressWarnings(get_weights(ind_preds, truth, 
                                                     ind_loss, convex, discrete = T))
    ind_weights_nnls_convex <- suppressWarnings(get_weights(ind_preds, truth, 
                                                        ind_loss, convex = T, discrete = F))
    ind_weights_nnls <- suppressWarnings(get_weights(ind_preds, truth, ind_loss))
    sl_ind_weights <- suppressWarnings(rbind(ind_weights_discrete, ind_weights_nnls_convex, 
                                             ind_weights_nnls))
    colnames(sl_ind_weights) <- names(ind_preds)
    
    # use sl weights for prediction with sl
    pred_ind_discrete <- as.matrix(ind_preds) %*% ind_weights_discrete
    pred_ind_nnls_convex <- as.matrix(ind_preds) %*% ind_weights_nnls_convex
    pred_ind_nnls <- as.matrix(ind_preds) %*% ind_weights_nnls
    sl_ind_pred <- data.frame(pred_ind_discrete, pred_ind_nnls_convex, pred_ind_nnls)
    colnames(sl_ind_pred) <- c("discreteSL", "nnls_convexSL", "nnls_SL")
    
    # forecast:
    ind_forecast_discrete <- as.matrix(ind_forecast) %*% ind_weights_discrete
    ind_forecast_nnls_convex <- as.matrix(ind_forecast) %*% ind_weights_nnls_convex
    ind_forecast_nnls <- as.matrix(ind_forecast) %*% ind_weights_nnls
    sl_ind_forecasts <- data.frame(ind_forecast_discrete, ind_forecast_nnls_convex, ind_forecast_nnls)
    colnames(sl_ind_forecasts) <- c("discreteSL", "nnls_convexSL", "nnls_SL")
    
  }else{
    ind_loss <- NULL
    sl_ind_weights <- NULL
    sl_ind_forecasts <- NULL
  }
  
  return_list <- list(sl_forecasts = sl_forecasts,
                      sl_weights = sl_weights,
                      individual_fit = ind_fit, 
                      historical_fit = historical_fit,
                      training_preds = training_preds,
                      loss = loss, 
                      lags=lags,
                      #Return just individual predictions:
                      sl_ind_forecasts=sl_ind_forecasts,
                      sl_ind_weights=sl_ind_weights,
                      ind_loss=ind_loss,
                      #Forecast task:
                      forecast_task = forecast_task)
  
  return(return_list)
}