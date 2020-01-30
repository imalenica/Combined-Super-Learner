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
make_historical_fit <- function(historical_data, outcome, covariates, id,
                                historical_stack){
  task <- make_sl3_Task(
    data = historical_data, 
    covariates = covariates,
    outcome = outcome,
    id = id,
    drop_missing_outcome = T,
  )
  cv_stack <- Lrnr_cv$new(historical_stack)
  fit <- cv_stack$train(task)
  return(fit)
}


################################################################################
make_combined_sl <- function(individual_data, outcome, covariates, subject_id,
                             historical_fit, individual_stack = NULL, 
                             cv_type = NULL, past_individual_fit = NULL, 
                             past_cv_type = NULL, past_sl_weights = NULL){

  # TODO: ensure individual task nodes = historical task nodes, 
  #       ind learner changes, update loss and weights w/o recalculating,
  #       inference, nested super learners, new covariates, check cv type 
  #       correspondence 
  
  if(!is.null(past_individual_fit)){
    if(is.null(past_sl_weights) | is.null(past_cv_type)){
      print("Error: Individual fit provided but previous SL weights and/or
            cross-validation type are not. Provide all three to update previous 
            fit and calculate online cross-validated risk.")
      stop()
    }
  }
  
  if(!is.null(past_cv_type)){
    cv_type <- past_cv_type
  }
  
  # set up folds with specified time-series cv sheme
  if(cv_type == "folds_rolling_origin") {
    folds <- origami::make_folds(individual_data, 
                                 fold_fun = folds_rolling_origin,
                                 first_window = 5, 
                                 validation_size = 5, 
                                 gap = 0,
                                 batch = 5
    )
    
  } else if(cv_type == "folds_rolling_window") {
    folds <- origami::make_folds(individual_data,
                                 fold_fun = folds_rolling_window,
                                 window_size = 5,
                                 validation_size = 5, 
                                 gap = 0,
                                 batch = 5
    )
  }
  
  ind_task <- make_sl3_Task(
    data = individual_data, 
    covariates = covariates,
    outcome = outcome, 
    folds = folds
  )
  
  # fit initial superlearner if ind_fit is not provided
  if(is.null(past_individual_fit) & !(is.null(individual_stack))) {
    
    print(paste0("Training initial super learner for subject ", subject_id, 
                 " with ", nrow(individual_data), " observations."))
    
    # fit individualized learners
    cv_stack <- Lrnr_cv$new(individual_stack)
    ind_fit <- cv_stack$train(ind_task)
  } else if(!is.null(past_individual_fit)) {
    print(paste0("Updating previous fit with previous learners for subject ",
                 subject_id, " with ", nrow(individual_data), "observations."))
    
    # update individualized learners
    ind_fit <- past_individual_fit$update(ind_task)
  }
  
  # predict with individualized learners and historical learners
  ind_preds <- ind_fit$predict(ind_task)
  hist_preds <- historical_fit$predict(ind_task)
  
  # combine predictions
  learners <- c(paste0("individual_", colnames(ind_preds)), 
                paste0("historical_", colnames(hist_preds)))
  preds <- cbind.data.frame(ind_preds, hist_preds)
  names(preds) <- learners
  
  # get true Y
  truth <- lapply(folds, function(i) data.frame(individual_data[i$validation,
                                                                get(outcome)]))
  truth <- bind_rows(truth)[,1]
  
  # evaluate loss 
  loss <- apply(preds, 2, function(pred) sum((pred-truth)^2))
  
  # establish various super learners
  weights_discrete <- get_weights(preds, truth, loss, convex, discrete = T)
  weights_convex <- get_weights(preds, truth, loss, convex = T, discrete = F)
  weights_regular <- get_weights(preds, truth, loss)
  sl_weights <- rbind(weights_discrete, weights_convex, weights_regular)
  colnames(sl_weights) <- names(preds)
    
  # obtain predictions for super learner
  preds_discrete <- as.matrix(preds) %*% weights_discrete
  preds_convex <- as.matrix(preds) %*% weights_convex
  preds_regular <- as.matrix(preds) %*% weights_regular
  
  # put it all together
  sl_preds <- data.frame(preds_discrete, preds_convex, preds_regular)
  colnames(sl_preds) <- c("discreteSL", "convexSL", "regularSL")
  all_preds <- data.frame(sl_preds, preds)
  
  return_list <- list(predictions = all_preds, 
                      sl_weights = sl_weights, 
                      individual_fit = ind_fit, 
                      historical_fit = historical_fit,
                      cv_type = cv_type,
                      loss = loss)
  
  # evaluate loss of super learner
  if(!is.null(past_sl_weights)) {
    print(paste0("Evaluating SL performace with online cross-validated risk for subject ", 
    subject_id, "."))
    past_sl <- apply(past_sl_weights, 1, function(w) as.matrix(preds) %*% w)
    sl_loss <- apply(past_sl, 2, function(pred) sum((pred-truth)^2))
    names(sl_loss) <- c("discreteSL", "convexSL", "regularSL")
    return_list[["sl_loss"]] <- sl_loss
  }
  
  return(return_list)
}



  
  metalearner <- make_learner(Lrnr_nnls)
  sl_task <- historical_fit$chain(task)
  ml_fit <- metalearner$train(sl_task)
  sl_pipeline <- make_learner(Pipeline, ind_fit, ml_fit)
  sl_preds <- sl_pipeline$predict(ind_task)