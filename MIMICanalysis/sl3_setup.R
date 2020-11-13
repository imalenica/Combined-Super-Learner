make_historical_stack <- function(){
  
  # learners with no internal screening
  lrn_glm <- Lrnr_glm$new()
  lrn_mean <- Lrnr_mean$new()
  lrn_spline <- Lrnr_polspline$new()
  
  # learners with internal screening
  grid_params <- list(max_depth=c(3, 6), eta=c(0.05, 0.2), nrounds=c(50, 100),
                      colsample_bytree = 0.8, early_stopping_rounds = 25)
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
  rf <- Lrnr_ranger$new(min.node.size=500)
  lrnrs <- c(xgbs, glmnets, db2, db.5, rf, lrn_glm, lrn_mean, lrn_spline)
  names(lrnrs) <- c(
    "xgb_50_3_0.05", "xgb_50_6_0.05", "xgb_50_3_0.2", "xgb_50_6_0.2", 
    "xgb_100_3_0.05", "xgb_100_6_0.05", "xgb_100_1_3_0.2", "xgb_100_6_0.2", 
    "glmnet_0", "glmnet_0.1", "glmnet_0.2", "glmnet_0.3", "glmnet_0.4",
    "glmnet_0.5", "glmnet_0.6", "glmnet_0.7", "glmnet_0.8", "glmnet_0.9", 
    "glmnet_1", "dbarts_2", "dbarts_0.5", "ranger", "glm", "mean", "spline"
  )
  stack <- make_learner(Stack, lrnrs)
  
  # screener & pipeline ranger for learners with no internal screening
  rf_fast <- Lrnr_ranger$new(min.node.size=10000, write.forest=F,
                             importance="impurity_corrected")
  screener_importance <- Lrnr_screener_importance$new(rf_fast, num_screen=100)
  W <- c("rank_icu", "sex", "age", "care_unit", "sapsi_first", "sofa_first", 
         "bmi", "admission_type_descr", "delta_bmi", "delta_sapsi_first", 
         "delta_sofa_first")
  screener_augmented_rf <- Lrnr_screener_augment$new(screener_importance, W)
  
  # final stack
  return(make_learner(Pipeline, screener_augmented_rf, stack))
}

get_covariates <- function(data){
  outcomes <- colnames(data)[grepl("Y", colnames(data))]
  not_covs <- c(outcomes, "time_and_date", "init_time_and_date", "subject_id", 
                "icustay_id", "id")
  covs <- colnames(data)[-which(colnames(data) %in% not_covs)]
  return(covs)
}

make_individual_cv_stack <- function() {
  
  # no internal screening stack
  lrn_glm <- Lrnr_glm$new()
  lrn_mean <- Lrnr_mean$new()
  arima_aicc <- Lrnr_arima$new(stepwise=F, approximation=F, ic="aicc")
  arima_bic <- Lrnr_arima$new(stepwise=F, approximation=F, ic="bic")
  lrnrs_no_screen <- c(lrn_glm, lrn_mean, arima_aicc, arima_bic)
  names(lrnrs_no_screen) <- c("glm", "mean", "arima_aicc", "arima_bic")
  stack_no_screen <- make_learner(Stack, lrnrs_no_screen)
  
  lasso <- Lrnr_glmnet$new()
  screener_lasso <- Lrnr_screener_coefs$new(lasso, threshold=0)
  pipe_screen <- make_learner(Pipeline, screener_lasso, stack_no_screen)
  
  # internal screening stack
  grid_params <- list(alpha=seq(0, 1, 0.2))
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS=F)
  params_default <- list(nthread = getOption("sl.cores.learners", 1))
  glmnets <- apply(grid, 1, function(params_tune) {
    do.call(Lrnr_glmnet$new, c(params_default, as.list(params_tune)))
  })
  
  rf <- Lrnr_ranger$new()
  xgb <- Lrnr_xgboost$new(max_depth=5, eta=0.1, nrounds=100, 
                          colsample_bytree=0.8, early_stopping_rounds=25)
  nlts_setar <- Lrnr_tsDyn$new(learner="setar", model="MTAR")
  nlts_nnet15 <- Lrnr_tsDyn$new(learner="nnetTs", m=15, size=5)
  nlts_nnet25 <- Lrnr_tsDyn$new(learner="nnetTs", m=25, size=5)
  lrnrs <- c(xgb, glmnets, rf, nlts_lstar15, nlts_lstar25, nlts_setar, 
             nlts_nnet15, nlts_nnet25, pipe_screen)
  names(lrnrs) <- c(
    "xgb", "glmnet_0", "glmnet_0.2",  "glmnet_0.4", "glmnet_0.6", "glmnet_0.8", 
    "glmnet_1", "ranger", "nlts_lstar15", "nlts_lstar25", "nlts_setar", 
    "nlts_nnet15", "nlts_nnet25", "lasso_screen"
  )
  stack <- make_learner(Stack, lrnrs)
  cv_stack <- Lrnr_cv$new(stack, full_fit = T)
  # all_lrnrs <- c(cv_stack, historical_fit)
  # names(all_lrnrs) <- c("individual", "historical")
  # final_stack <- make_learner(Stack, all_lrnrs)
  return(cv_stack) 
}
