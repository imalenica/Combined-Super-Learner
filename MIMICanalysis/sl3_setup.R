get_covariates <- function(data){
  outcomes <- colnames(data)[grepl("Y", colnames(data))]
  not_covs <- c(outcomes, "time_and_date", "init_time_and_date", "subject_id", 
                "icustay_id", "id", "min_elapsed")
  covs <- colnames(data)[-which(colnames(data) %in% not_covs)]
  return(covs)
}

make_historical_stack <- function(outcome_type, screen = TRUE, 
                                  num_screen = 200){
  
  if(outcome_type == "binomial"){
    scale_pos_weight <- 7 # boosted-tree args
    k <- 0.5 # dbarts args
  } else {
    scale_pos_weight <- 1 # boosted-tree args
    k <- 2 # dbarts args
  }
  xgb_grid <- expand.grid(
    list(eta = 0.1, nrounds = 1000, min_child_weight = 5000, verbose = 2,
         subsample = c(0.4, 0.7, 1), early_stopping_rounds = 50, 
         max_depth = 10, scale_pos_weight = scale_pos_weight), 
    KEEP.OUT.ATTRS = FALSE
  )
  xgbs <- apply(xgb_grid, 1, function(args) do.call(Lrnr_xgboost$new, as.list(args)))
    
  rf_grid <- expand.grid(
    list(min.node.size = 5000, oob.error = FALSE, sample.fraction = c(0.4, 0.7, 1),
         num.trees = 1000, save.memory = TRUE), KEEP.OUT.ATTRS = FALSE
  )
  rfs <- apply(rf_grid, 1, function(args) do.call(Lrnr_ranger$new, as.list(args)))
  
  enet_grid <- expand.grid(
    list(alpha = seq(0, 1, 0.2)), KEEP.OUT.ATTRS = FALSE
  )
  enets <- apply(enet_grid, 1, function(args) do.call(Lrnr_glmnet$new, as.list(args)))
  
  dbarts <- Lrnr_dbarts$new(
    k = k, ntree = 1000, ndpost = 1000, nskip = 100, keeptrees = TRUE
  )
  
  lrn_glm <- Lrnr_glm$new()
  bayesglm <- Lrnr_bayesglm$new()
  lrnrs_linear <- c(lrn_glm, bayesglm)
  names(lrnrs_linear) <- c("glm", "bayesglm")
  stack_linear <- make_learner(Stack, lrnrs_linear)
  lasso <- Lrnr_glmnet$new()
  screener_lasso <- Lrnr_screener_coefs$new(lasso, threshold=0)
  pipe_lasso <- make_learner(Pipeline, screener_lasso, stack_linear)
  
  lrnrs <- c(xgbs, rfs, dbarts, enets, pipe_lasso)
  names(lrnrs) <- c(
    "xgboost_subsample.4", "xgboost_subsample.7", "xgboost_subsample1", 
    "ranger_subsample.4", "ranger_subsample.7", "ranger_subsample1",
    "dbarts", "enet_alpha0", "enet_alpha.2", "enet_alpha.4", "enet_alpha.6", 
    "enet_alpha.8", "enet_alpha1", "screen_lasso"
  )
  
  stack <- make_learner(Stack, lrnrs)
  
  if(screen){
    # screener & pipeline ranger 
    rf_fast <- Lrnr_ranger$new(
      min.node.size = 10000, write.forest = FALSE, importance = "impurity_corrected", 
      oob.error = FALSE, save.memory = TRUE, sample.fraction = 0.7
    )
    screener_importance <- Lrnr_screener_importance$new(rf_fast, num_screen = num_screen)
    W <- c("sex", "age", "care_unit", "sapsi_first", "sofa_first", "bmi", 
           "admission_type_descr", "delta_bmi", "delta_sapsi_first", 
           "delta_sofa_first", "min")
    screener_augmented_rf <- Lrnr_screener_augment$new(screener_importance, W)
    
    # final stack
    final_stack <- make_learner(Pipeline, screener_augmented_rf, stack)
    return(final_stack)
  } else {
    return(stack)
  }  
}

make_individual_cv_stack <- function(){
  
  ##### learners that accommodate external regressors
  # no internal covariate screening
  lrn_glm <- Lrnr_glm$new()
  lrn_bayesglm <- Lrnr_bayesglm$new()
  lrn_mean <- Lrnr_mean$new() #fix-ish naming issue when all except 1 in pipe fail
  arima_aicc_ord15 <- Lrnr_arima$new(stepwise=F, approximation=F, ic="aicc", 
                                     max.order=15)
  arima_bic_ord15 <- Lrnr_arima$new(stepwise=F, approximation=F, ic="bic", 
                                    max.order=15)
  arima_aicc_ord5 <- Lrnr_arima$new(stepwise=F, approximation=F, ic="aicc", 
                                    max.order=5)
  arima_bic_ord5 <- Lrnr_arima$new(stepwise=F, approximation=F, ic="bic", 
                                   max.order=5)
  
  lrnrs_no_screen <- c(lrn_glm, lrn_bayesglm, lrn_mean, arima_aicc_ord5, 
                       arima_bic_ord5, arima_aicc_ord15, arima_bic_ord15)
  names(lrnrs_no_screen) <- c(
    "glm", "bayesglm", "mean", "arima_aicc_order5", "arima_bic_order5",  
    "arima_aicc_order15", "arima_bic_order15"
  )
  stack_no_screen <- make_learner(Stack, lrnrs_no_screen)
  lasso <- Lrnr_glmnet$new(stratify_cv = TRUE)
  screener_lasso <- Lrnr_screener_coefs$new(lasso, threshold=0)
  pipe_screen <- make_learner(Pipeline, screener_lasso, stack_no_screen)

  # internal screening 
  grid_params <- list(alpha = seq(0, 1, 0.2))
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS=F)
  glmnets <- apply(grid, 1, function(par) do.call(Lrnr_glmnet$new, as.list(par)))
  
  rf <- Lrnr_ranger$new(oob.error = FALSE)
  xgb <- Lrnr_xgboost$new(verbose=0, eta=0.1, nrounds=500, 
                          early_stopping_rounds=25)
  # tensorflow: Early stopping conditioned on metric `val_loss` which is 
  # not available. Available metrics are: loss
  # earlystop <- keras::callback_early_stopping(patience = 10, monitor = "loss")
  # callbacks <- list(earlystop)
  # grid_params <- list(batch_size=c(5,10), epochs=c(200,1000), lr=c(0.0005,0.01), 
  #                     units=c(32,64), layers=2, callbacks=callbacks)
  # grid <- expand.grid(grid_params, KEEP.OUT.ATTRS=F)
  # gru_keras <- apply(grid, 1, function(par) do.call(Lrnr_gru_keras$new, as.list(par)))
  # lstm_keras <- apply(grid, 1, function(par) do.call(Lrnr_lstm_keras$new, as.list(par)))
  # 
  ##### learners that DO NOT accommodate external regressors   
  nlts_setar <- Lrnr_tsDyn$new(learner="setar", m=1)
  nlts_nnet5 <- Lrnr_tsDyn$new(learner="nnetTs", m=5, size=2, 
                               control=list(trace=F))
  nlts_nnet15 <- Lrnr_tsDyn$new(learner="nnetTs", m=15, size=2, 
                                control=list(trace=F))
  
  ##### put it all together        
  # lrnrs <- c(glmnets, xgb, rf, gru_keras, lstm_keras, nlts_setar, nlts_nnet5, 
  #            nlts_nnet15, pipe_screen)
  # names(lrnrs) <- c(
  #   "glmnet_0", "glmnet_0.2", "glmnet_0.4", "glmnet_0.6",  "glmnet_0.8", 
  #   "glmnet_1", "xgboost", "ranger", "gru_5_32_200_5e-04", "gru_10_32_200_5e-04", 
  #   "gru_5_32_1000_5e-04", "gru_10_32_1000_5e-04", "gru_5_32_200_0.01", 
  #   "gru_10_32_200_0.01", "gru_5_32_1000_0.01", "gru_10_32_1000_0.01", 
  #   "gru_5_64_200_5e-04", "gru_10_64_200_5e-04", "gru_5_64_1000_5e-04", 
  #   "gru_10_64_1000_5e-04", "gru_5_64_200_0.01", "gru_10_64_200_0.01", 
  #   "gru_5_64_1000_0.01", "gru_10_64_1000_0.01", "lstm_5_32_200_5e-04", 
  #   "lstm_10_32_200_5e-04", "lstm_5_32_1000_5e-04", "lstm_10_32_1000_5e-04", 
  #   "lstm_5_32_200_0.01", "lstm_10_32_200_0.01", "lstm_5_32_1000_0.01", 
  #   "lstm_10_32_1000_0.01", "lstm_5_64_200_5e-04", "lstm_10_64_200_5e-04", 
  #   "lstm_5_64_1000_5e-04", "lstm_10_64_1000_5e-04", "lstm_5_64_200_0.01", 
  #   "lstm_10_64_200_0.01", "lstm_5_64_1000_0.01", "lstm_10_64_1000_0.01", 
  #   "nlts_setar", "nlts_nnet5", "nlts_nnet15", "lasso_screen"
  # )
  # lrnrs <- lrnrs[-grep("200_5e-04", names(lrnrs))] # rm keras with epoch=200 & lr=5e-4 
  
  lrnrs <- c(glmnets, xgb, rf, nlts_setar, nlts_nnet5, nlts_nnet15, pipe_screen)
  names(lrnrs) <- c(
    "enet_alpha0", "enet_alpha.2", "enet_alpha.4", "enet_alpha.6",  
    "enet_alpha.8", "enet_alpha1", "xgboost", "ranger", "nlts_setar", 
    "nlts_nnet5", "nlts_nnet15", "lasso_screen"
  )
  stack <- make_learner(Stack, lrnrs)
  cv_stack <- Lrnr_cv$new(stack)
  # all_lrnrs <- c(cv_stack, historical_fit)
  # names(all_lrnrs) <- c("individual", "historical")
  # final_stack <- make_learner(Stack, all_lrnrs)
  return(cv_stack) 
}

