get_covariates <- function(data){
  outcomes <- colnames(data)[grepl("Y", colnames(data))]
  not_covs <- c(outcomes, "time_and_date", "init_time_and_date", "subject_id", 
                "icustay_id", "id")
  covs <- colnames(data)[-which(colnames(data) %in% not_covs)]
  return(covs)
}

make_historical_stack <- function(){
  
  mean <- make_learner(Lrnr_mean)
  dbarts2 <- Lrnr_dbarts$new(k = 2, ntree = 200, ndpost = 1000, nskip = 100, 
                             keeptrees = TRUE)
  glm <- make_learner(Lrnr_glm)
  ranger <- make_learner(Lrnr_ranger)
  spline <- make_learner(Lrnr_polspline)

  grid_params <- list(max_depth=c(3,5), eta=c(0.05,0.2), nrounds=c(50,100))
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
  params_default <- list(nthread = getOption("sl.cores.learners", 1))
  xgboosts <- apply(grid, MARGIN = 1, function(params_tune) {
    do.call(Lrnr_xgboost$new, c(params_default, as.list(params_tune)))
  })
  
  grid_params <- list(alpha = seq(0, 1, 0.1))
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
  params_default <- list(nthread = getOption("sl.cores.learners", 1))
  glmnets <- apply(grid, MARGIN = 1, function(params_tune) {
    do.call(Lrnr_glmnet$new, c(params_default, as.list(params_tune)))
  })
  
  stack <- make_learner_stack(
    unlist(list(xgboosts, glmnets, dbarts2, spline, ranger, mean, glm), 
           recursive = TRUE)
  )
  
  # screener & pipeline ranger
  ranger_fast <- make_learner(Lrnr_ranger, write.forest = FALSE, 
                              importance = "impurity_corrected",
                              min.node.size = 20000)
  screener_ranger <- Lrnr_screener_importance$new(ranger_fast, num_screen = 50)
  W <- c("rank_icu", "sex", "age", "care_unit", "sapsi_first", "sofa_first", 
         "bmi", "admission_type_descr", "delta_sapsi_first", "delta_bmi",
         "delta_sofa_first")
  screener_ranger_augmented <- Lrnr_screener_augment$new(screener_importance, W)
  ranger_pipe <- make_learner(Pipeline, screener_ranger_augmented, stack)
  
  return(make_learner(Stack, ranger_pipe))
}

make_individual_cv_stack <- function(n = NULL, n_adaptive = TRUE) {
  
  mean <- Lrnr_mean$new()
  glm <- Lrnr_glm$new()
  ranger <- Lrnr_ranger$new()
  spline <- Lrnr_polspline$new()
  earth <- Lrnr_earth$new()
  lasso <- Lrnr_glmnet$new()
  
  arima_aicc <- Lrnr_arima$new(stepwise=FALSE, approximation=FALSE, ic="aicc")
  arima_bic <- Lrnr_arima$new(stepwise=FALSE, approximation=FALSE, ic="bic")
  nlts_lstar <- Lrnr_tsDyn$new(learner="lstar", m=25, mTh=25)
  nlts_setar <- Lrnr_tsDyn$new(learner="setar", m=25, model="MTAR", mTh=rep(1,25))
  nlts_nnet <- Lrnr_tsDyn$new(learner="nnetTs", m=25, size=5)
  
  grid_params <- list(max_depth=c(3,5), eta=c(0.05,0.2), nrounds=c(50,100))
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
  params_default <- list(nthread = getOption("sl.cores.learners", 1))
  xgboosts <- apply(grid, MARGIN = 1, function(params_tune) {
    do.call(Lrnr_xgboost$new, c(params_default, as.list(params_tune)))
  })
  
  grid_params <- list(alpha=seq(0.1, 1, 0.1))
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
  params_default <- list(nthread = getOption("sl.cores.learners", 1))
  glmnets <- apply(grid, MARGIN = 1, function(params_tune) {
    do.call(Lrnr_glmnet$new, c(params_default, as.list(params_tune)))
  })
  
  dbarts2 <- Lrnr_dbarts$new(k = 2, ntree = 200, ndpost = 1000, nskip = 100, 
                             keeptrees = TRUE)
  dbarts.5 <- Lrnr_dbarts$new(k = 0.5, ntree = 200, ndpost = 1000, 
                              nskip = 100, keeptrees = TRUE)
  
  # create the initial base_stack
  stack_init <- make_learner(Stack, spline, earth, glm, arima_aicc, arima_bic)
  lasso_screener <- Lrnr_screener_coefs$new(lasso, threshold = 0)
  lasso_pipe <- make_learner(Pipeline, lasso_screener, stack_init)
  
  # build off base_stack as n increases
  if(n_adaptive && !is.null(n)){
    if(n < 50){
      stack <- make_learner(Stack, lasso_pipe, mean, nlts_lstar, nlts_setar,
                            nlts_nnet, lasso)
    } else if (n >= 50 && n < 200) {
      stack <- make_learner(
        Stack,
        unlist(list(glmnets, ranger, dbarts2, dbarts.5, lasso_pipe, mean, 
                    nlts_lstar, nlts_setar, nlts_nnet), 
               recursive = TRUE)
      )
    } else {
      stack <- make_learner(
        Stack,
        unlist(list(xgboosts, glmnets, ranger, dbarts2, dbarts.5, lasso_pipe, 
                    mean, nlts_lstar, nlts_setar, nlts_nnet), 
               recursive = TRUE)
      )
    }
  } else {
    stack <- make_learner(
      Stack,
      unlist(list(xgboosts, glmnets, ranger, dbarts2, dbarts.5, lasso_pipe, 
                  mean, nlts_lstar, nlts_setar, nlts_nnet), 
             recursive = TRUE)
    )
  }
  return(Lrnr_cv$new(stack, full_fit = TRUE)) 
}

# not currently using this
make_individual_args <- function(individual_data, burn_in = 30, batch = 5, 
                                 max_stop_time = 1440){
  ids <- as.character(unique(individual_data$id))
  individual_args_list <- lapply(ids, function(sub){
    d <- droplevels(individual_data[id == sub,])
    d <- setorder(d, time_and_date)
    
    W <- c("rank_icu", "sex", "age", "care_unit", "sapsi_first", "sofa_first", 
           "bmi", "admission_type_descr")
    allW <- unname(unlist(sapply(W, function(x) colnames(d)[grep(x, colnames(d))])))
    
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
                      ifelse(n >= 100 && n < 200, 3, 4)))
      }
      # make non-identical stacks
      cv_stack <- list()
      training_data <- list()
      forecast_data <- list()
      for(i in 1:length(splits)){
        if(i == 1){
          cv_stack[[i]] <- make_ind_stack(n=splits[[i]], W=allW)
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