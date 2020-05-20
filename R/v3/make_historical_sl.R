# TODO: add acf option?

################################################################################
# historical SLs
################################################################################

# historical_data: data.table of all observed data on all previously seen subjects
# outcome: name of outcome variable which is found in historical data.table
# covariates: vector of covariate names which are column names in the data.table  
# historical_stack: stack of sl3 learners to model historical data  
# id: column name indicating unique subject identifier for appropraite V-fold cv

make_historical_sl <- function(historical_data, outcome, covariates, id,
                               historical_stack, V = 5, parallelize = FALSE, 
                               cpus_logical = NULL, fit_sl = TRUE){
  
  task <- make_sl3_Task(
    data = historical_data, 
    covariates = covariates,
    outcome = outcome,
    id = id,
    drop_missing_outcome = T,
    folds = origami::make_folds(historical_data, fold_fun = folds_vfold, V = V)
  )
  cv_stack <- Lrnr_cv$new(historical_stack, full_fit = TRUE)
  if(parallelize){
    plan(multicore, workers = cpus_logical)
    test <- delayed_learner_train(cv_stack, task)
    sched <- Scheduler$new(test, FutureJob, nworkers = cpus_logical,
                           verbose = FALSE)
    fit <- sched$compute()
  } else {
    fit <- cv_stack$train(task)
  }
  
  if(fit_sl){
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
    
    fit_return_obj <- list(cv_fit = fit, sl_discrete = sl_discrete, 
                           sl_nnls = sl_nnls, sl_nnls_convex = sl_nnls_convex)
  } else {
    fit_return_obj <- fit
  }
  return(fit_return_obj)
}