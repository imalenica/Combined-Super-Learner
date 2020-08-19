get_covariates <- function(data) {
  outcomes <- colnames(data)[grepl("Y", colnames(data))]
  not_covariates <- c(outcomes, "time_and_date", "init_time_and_date",
                      "subject_id", "icustay_id", "id")
  covariates <- colnames(data)[-which(colnames(data) %in% not_covariates)]
  return(covariates)
}

make_historical_task <- function(historical_data, outcome, folds, 
                                 covariates = NULL, id = "id", time = "min") {
  if(is.null(covariates)) {
    covariates <- get_covariates(data)
  }
  
  if(is.null(folds)) {
    task <- make_sl3_Task(data = data, outcome = outcome, 
                          covariates = covariates, id = id, time = time)
  } else {
    task <- make_sl3_Task(data = data, outcome = outcome, 
                          covariates = covariates, id = id, time = time, 
                          folds = folds)
  }
  
  return(task)
}

make_historical_cv_stack <- function() {
  # regular learners
  grid_params <- list(max_depth = c(3, 5), eta = c(0.05, 0.2), 
                      nrounds = c(50, 100))
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
  params_default <- list(nthread = getOption("sl.cores.learners", 1))
  xgb_learners <- apply(grid, MARGIN = 1, function(params_tune) {
    do.call(Lrnr_xgboost$new, c(params_default, as.list(params_tune)))
  })
  
  grid_params <- list(alpha = seq(0, 1, 0.1), nfolds = 5)
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
  params_default <- list(nthread = getOption("sl.cores.learners", 1))
  glmnet_learners <- apply(grid, MARGIN = 1, function(params_tune) {
    do.call(Lrnr_glmnet$new, c(params_default, as.list(params_tune)))
  })
  
  lrnr_mean <- make_learner(Lrnr_mean)
  lrnr_dbarts <- make_learner(Lrnr_dbarts, sigdf = 3, sigquant = 0.90, k = 2,
                              power = 2.0, base = 0.95, binaryOffset = 0.0,
                              ntree = 200, ndpost = 1000,  nskip = 100,
                              printevery = 100,  keepevery = 1,
                              keeptrainfits = TRUE, usequants = FALSE)
  lrnr_glm <- make_learner(Lrnr_glm)
  lrnr_ranger <- make_learner(Lrnr_ranger)
  lrnr_polspline <- make_learner(Lrnr_polspline)
  
  stack <- make_learner(Stack, unlist(list(xgb_learners, glmnet_learners,
                                           lrnr_dbarts, lrnr_glm, lrnr_ranger,
                                           lrnr_polspline), recursive = TRUE))
  
  # screener & pipeline
  screener_rf <- make_learner(Lrnr_screener_randomForest, nVar = 30,
                              ntree = 500, nodesize = 50000,
                              keep.forest = FALSE, mtry = 50)
  screen_rf_pipe <- make_learner(Pipeline, screener_rf, stack)
  
  # final stack
  screen_stack <- make_learner(Stack, screen_rf_pipe)
  cv_stack <- Lrnr_cv$new(screen_stack, full_fit = TRUE)
  return(cv_stack)
}


make_individual_cv_stack <- function(n) {
  
  lrnr_mean <- sl3::make_learner(Lrnr_mean)
  lrnr_glm <- sl3::make_learner(Lrnr_glm)
  lrnr_arima <- sl3::make_learner(Lrnr_arima)
  lrnr_gam <- sl3::make_learner(Lrnr_gam)
  lrnr_ranger <- sl3::make_learner(Lrnr_ranger)
  lrnr_spline <- sl3::make_learner(Lrnr_polspline)
  lrnr_lasso <- sl3::make_learner(Lrnr_glmnet)
  lrnr_lstm <- make_learner(Lrnr_lstm, epochs = 500, batch_size = 5,
                            early_stop = TRUE)
  
  grid_params <- list(max_depth = c(3,4), eta = c(0.05,0.2),nrounds = c(50,100))
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
  params_default <- list(nthread = getOption("sl.cores.learners", 1))
  xgb_learners <- apply(grid, MARGIN = 1, function(params_tune) {
    do.call(Lrnr_xgboost$new, c(params_default, as.list(params_tune)))
  })
  
  grid_params <- list(alpha = seq(0, 1, 0.1))
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
  params_default <- list(nthread = getOption("sl.cores.learners", 1))
  glmnet_learners <- apply(grid, MARGIN = 1, function(params_tune) {
    do.call(Lrnr_glmnet$new, c(params_default, as.list(params_tune)))
  })
  
  tmle.SL.dbarts2 <- sl3::make_learner(
    Lrnr_dbarts,
    sigdf = 3, sigquant = 0.90, k = 2, power = 2.0, base = 0.95,
    binaryOffset = 0.0, ntree = 200, ndpost = 1000, nskip = 100,
    keepevery = 1, printevery = 100, keeptrainfits = T, usequants = F,
    numcut = 100, printcutoffs = 0, nthread = 1, keepcall = T
  )
  tmle.SL.dbarts.k.5 <- sl3::make_learner(
    Lrnr_dbarts,
    sigdf = 3, sigquant = 0.90, k = 0.5, power = 2.0,
    base = 0.95, binaryOffset = 0.0, ntree = 200, ndpost = 1000, nskip = 100,
    keepevery = 1, printevery = 100, keeptrainfits = T, usequants = F,
    numcut = 100, printcutoffs = 0, nthread = 1, keepcall = T
  )
  
  
  if(n < 50){
    stack <- make_learner(Stack, lrnr_mean, lrnr_glm, lrnr_arima, lrnr_gam)
    no.covs <- 5
  } else if (n >= 50 && n < 100) {
    stack <- make_learner(Stack, lrnr_mean, lrnr_glm, lrnr_arima, lrnr_gam, 
                          lrnr_spline, lrnr_lasso, lrnr_lstm)
    no.covs <- 7
  } else if (n >= 100 && n < 200) {
    stack <- make_learner(
      Stack, 
      unlist(list(glmnet_learners, lrnr_mean, lrnr_glm, lrnr_arima, lrnr_gam, 
                  lrnr_spline, lrnr_lstm, lrnr_ranger), recursive = TRUE)
    )
    no.covs <- 10
  } else if (n >= 200 && n < 500) {
    stack <- make_learner(
      Stack, 
      unlist(list(glmnet_learners, tmle.SL.dbarts.k.5, tmle.SL.dbarts2, 
                  lrnr_mean, lrnr_glm, lrnr_arima, lrnr_gam, lrnr_spline, 
                  lrnr_lstm, lrnr_ranger), recursive = TRUE)
    )
    no.covs <- 20
  } else {
    stack <- make_learner(
      Stack, 
      unlist(list(glmnet_learners, tmle.SL.dbarts.k.5, tmle.SL.dbarts2, 
                  lrnr_mean, lrnr_glm, lrnr_arima, lrnr_gam, lrnr_spline, 
                  lrnr_lstm, lrnr_ranger), recursive = TRUE)
    )
    no.covs <- 30
  }
  
  # screener & pipeline
  screener_rf <- make_learner(Lrnr_screener_randomForest, nVar = no.covs,
                              ntree = 500, nodesize = 5, keep.forest = FALSE, 
                              mtry = 50)
  screen_rf_pipe <- make_learner(Pipeline, screener_rf, stack)
  
  # final stack
  stack <- make_learner(Stack, screen_rf_pipe)
  cv_stack <- Lrnr_cv$new(stack, full_fit = TRUE)
  return(cv_stack)
}

make_individual_args <- function(individual_data, burn_in = 30, batch = 5, 
                                 max_stop_time = 1440){
  ids <- as.character(unique(individual_data$id))
  individual_args_list <- lapply(ids, function(sub){
    d <- droplevels(individual_data[id == sub,])
    d <- setorder(d, time_and_date)
    
    max_times <- c(max_stop_time, nrow(d))
    max_time <- max_times[which.min(max_times)]
    
    if(max_time < 100) {
      return(NULL)
    } else {
      # create pseudo-online data 
      splits <- seq(burn_in, max_time, batch)
      split_data <- lapply(splits, function(x) data.table(d[1:x, ]))
      split_group <- function(n){
        ifelse(n < 50, 1, 
               ifelse(n >= 50 && n < 100, 2,
                      ifelse(n >= 100 && n < 200, 3,
                             ifelse(n >= 200 && n < 500, 4, 5))))
      }
      # make non-identical stacks
      cv_stack <- list()
      training_data <- list()
      forecast_data <- list()
      for(i in 1:length(splits)){
        if(i == 1){
          cv_stack[[i]] <- make_ind_stack(splits[[i]])
        } else {
          previous_split_grp <- split_group(splits[[i-1]])
          current_split_grp <- split_group(splits[[i]])
          cv_stack[[i]] <- ifelse(previous_split_grp == current_split_grp, 
                                  NA, make_ind_stack(splits[[i]]))
        }
        training_data[[i]] <- split_data[[i]]
        forecast_data[[i]] <- head(split_data[[i+1]], batch)
      }
      return(list(training_data = training_data, 
                  forecast_data = forecast_data, 
                  cv_stack = cv_stack))
    }
  })
  names(individual_args_list) <- ids
  populated_list <- individual_args_list[!sapply(individual_args_list, is.null)]
  return(populated_list)
}



