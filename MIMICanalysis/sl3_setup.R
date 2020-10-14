get_covariates <- function(data){
  outcomes <- colnames(data)[grepl("Y", colnames(data))]
  not_covs <- c(outcomes, "time_and_date", "init_time_and_date", "subject_id", 
                "icustay_id", "id")
  covs <- colnames(data)[-which(colnames(data) %in% not_covs)]
  return(covs)
}


load_and_prep_historical_data <- function(file, 
                                          historical_covs_omit_pattern = NULL,
                                          location = "savio"){
  if(location == "savio"){
    data_path <- paste0("/global/scratch/rachelvphillips/symphony-data/", file)
    load(data_path)
  } else if (location == "bluevelvet"){
    load(here::here("data", file))
  }
  
  if(!is.null(historical_covs_omit_pattern)){
    historical_covs <- c("trend_strength", "stl_e_acf1", "stl_e_acf10", 
                         "spectral_entropy", "_lag_",  "spikiness", "linearity",
                         "curvature", "n_flat_spots", "coef_hurst", "_Min", 
                         "_1stQ", "_Median", "_Mean", "_3rdQ", "_Max")
    all_historical_covs <- unname(unlist(sapply(historical_covs, function(x){
      colnames(historical)[grep(x, colnames(historical))]
    })))
    to_rm <- as.character(sapply(historical_covs_omit_pattern, function(x){
      all_historical_covs[grep(x, all_historical_covs)]
    }))
    historical <- historical[, -to_rm, with = FALSE]
  }
  return(historical)
}

make_historical_stack <- function(){
  
  # learners with no internal screening
  lrn_glm <- Lrnr_glm$new(name = "historical_glm")
  lrn_mean <- Lrnr_mean$new(name = "historical_mean")
  lrn_spline <- Lrnr_polspline$new(name = "historical_spline")
  lrn_earth <- Lrnr_earth$new(name = "historical_earth")
  stack_no_screen <- make_learner(Stack, lrn_glm, lrn_spline, lrn_earth, lrn_mean)
  
  # learners with internal screening
  grid_params <- list(max_depth=c(3, 5), eta=c(0.05, 0.2), nrounds=c(50, 100))
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS=F)
  params_default <- list(nthread = getOption("sl.cores.learners", 1))
  xgbs <- apply(grid, 1, function(params_tune) {
    do.call(Lrnr_xgboost$new, c(params_default, as.list(params_tune)))
  })
  
  grid_params <- list(alpha = seq(0, 1, 0.1))
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS=F)
  params_default <- list(nthread = getOption("sl.cores.learners", 1))
  glmnets <- apply(grid, 1, function(params_tune) {
    do.call(Lrnr_glmnet$new, c(params_default, as.list(params_tune)))
  })
  
  db2 <- Lrnr_dbarts$new(k=2, ntree=200, ndpost=1000, nskip=100, keeptrees=T)
  db.5 <- Lrnr_dbarts$new(k=0.5, ntree=200, ndpost=1000, nskip=100, keeptrees=T)
  rf <- Lrnr_ranger$new()
  stack_screen <- make_learner(Stack, unlist(list(xgbs, glmnets, db2, db.5, rf), 
                                             recursive=T))
  lrnrs <- c(xgbs, glmnets, db2, db.5, rf)
  names(lrnrs) <- c(
    "xgb_50_3_0.05", "xgb_50_5_0.05", "xgb_50_3_0.2", "xgb_50_5_0.2", 
    "xgb_100_3_0.05", "xgb_100_5_0.05", "xgb_100_1_3_0.2", "xgb_100_5_0.2", 
    "glmnet_0", "glmnet_0.1", "glmnet_0.2", "glmnet_0.3", "glmnet_0.4",
    "glmnet_0.5", "glmnet_0.6", "glmnet_0.7", "glmnet_0.8", "glmnet_0.9", 
    "glmnet_1", "dbarts_2", "dbarts_0.5", "ranger"
  )
  names(lrnrs) <- paste0("historical_", names(lrnrs))
  stack_screen <- make_learner(Stack, lrnrs)
  
  # screener & pipeline ranger for learners with no internal screening
  rf_fast <- Lrnr_ranger$new(min.node.size=20000, write.forest=F,
                             importance="impurity_corrected")
  screener_importance <- Lrnr_screener_importance$new(rf_fast, num_screen=50)
  W <- c("rank_icu", "sex", "age", "care_unit", "sapsi_first", "sofa_first", 
         "bmi", "admission_type_descr", "delta_bmi", "delta_sapsi_first", 
         "delta_sofa_first")
  screener_augmented_rf <- Lrnr_screener_augment$new(screener_importance, W, 
                                                     name="rf_screen")
  pipe_screen <- make_learner(Pipeline, screener_augmented_rf, stack_no_screen)
  
  # final stack
  return(make_learner(Stack, pipe_screen, stack_screen))
}

make_individual_cv_stack <- function(n=NULL, n_adaptive=T) {
  
  lrn_glm <- Lrnr_glm$new(name = "individual_glm")
  lrn_spline <- Lrnr_polspline$new(name = "individual_spline")
  lrn_mars <- Lrnr_earth$new(name = "individual_earth")
  lrn_lasso <- Lrnr_glmnet$new("individual_lasso")
  lrn_mean <- Lrnr_mean$new(name="individual_mean")
  
  arima_aicc <- Lrnr_arima$new(stepwise=F, approximation=F, ic="aicc", 
                               name="individual_arima_aicc")
  arima_bic <- Lrnr_arima$new(stepwise=F, approximation=F, ic="bic", 
                              name="individual_arima_bic")
  nlts_lstar15 <- Lrnr_tsDyn$new(learner="lstar", m=15, mTh=15,
                                 name="individual_nlts_lstar15")
  nlts_lstar <- Lrnr_tsDyn$new(learner="lstar", m=25, mTh=25,
                               name="individual_nlts_lstar")
  nlts_setar <- Lrnr_tsDyn$new(learner="setar", model="MTAR",
                               name="individual_nlts_setar")
  nlts_nnet15 <- Lrnr_tsDyn$new(learner="nnetTs", m=15, size=5,
                                name="individual_nlts_nnet15")
  nlts_nnet <- Lrnr_tsDyn$new(learner="nnetTs", m=25, size=5,
                              name="individual_nlts_nnet")
  
  lrn_ranger <- Lrnr_ranger$new(name="individual_ranger")
  
  grid_params <- list(max_depth=c(3,5), eta=c(0.05,0.2), nrounds=c(50,100))
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS=F)
  params_default <- list(nthread = getOption("sl.cores.learners", 1))
  xgbs <- apply(grid, 1, function(params_tune) {
    do.call(Lrnr_xgboost$new, c(params_default, as.list(params_tune)))
  })
  
  grid_params <- list(alpha=seq(0.1, 1, 0.1))
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS=F)
  params_default <- list(nthread = getOption("sl.cores.learners", 1))
  glmnets <- apply(grid, 1, function(params_tune) {
    do.call(Lrnr_glmnet$new, c(params_default, as.list(params_tune)))
  })
  
  db2 <- Lrnr_dbarts$new(k=2, ntree=200, ndpost=1000, nskip=100, keeptrees=T)
  db.5 <- Lrnr_dbarts$new(k=0.5, ntree=200, ndpost=1000, nskip=100, keeptrees=T)
  
  # create the initial stack, no internal screening
  stack_no_screen <- make_learner(Stack, lrn_spline, lrn_mars, lrn_glm, 
                                  arima_aicc, arima_bic)
  screener_lasso <- Lrnr_screener_coefs$new(lrn_lasso, threshold=0, 
                                            name="lasso_screen")
  pipe_screen <- make_learner(Pipeline, screener_lasso, stack_no_screen)
  
  # build off base_stack as n increases
  if(n_adaptive && !is.null(n)){
    if(n < 50){
      stack <- make_learner(Stack, pipe_screen, lrn_mean, nlts_lstar, nlts_setar, 
                            nlts_nnet, nlts_lstar15, nlts_nnet15, lrn_lasso)
    } else if (n >= 50 && n < 200) {
      lrnrs <- c(glmnets, lrn_ranger, db2, db.5, lrn_mean, nlts_lstar, 
                 nlts_setar, nlts_nnet, nlts_lstar15, nlts_nnet15)
      stack_screen <- make_learner(Stack, lrnrs)
      stack <- make_learner(Stack, pipe_screen, stack_screen)
    } else {
      stack_screen <- make_learner(
        Stack, unlist(list(xgbs, glmnets, lrn_ranger, db2, db.5, lrn_mean, 
                           nlts_lstar, nlts_setar, nlts_nnet, nlts_lstar15,
                           nlts_nnet15), recursive=T)
      )
      stack <- make_learner(Stack, pipe_screen, stack_screen)
    }
  } else {
    stack <- make_learner(
      Stack, unlist(list(xgbs, glmnets,  lrn_ranger, db2, db.5, lrn_mean, 
                         nlts_lstar, nlts_setar, nlts_nnet, nlts_lstar15,
                         nlts_nnet15), recursive=T)
    )
    stack <- make_learner(Stack, pipe_screen, stack_screen)
  }
  return(Lrnr_cv$new(stack, full_fit=T)) 
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