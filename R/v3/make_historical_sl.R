################################################################################
# historical SLs
################################################################################

# historical_data: data.table of all observed data on all previously seen subjects
# outcome: name of outcome variable which is found in historical data.table
# covariates: vector of covariate names which are column names in the data.table  
# historical_stack: stack of sl3 learners to model historical data  
# id: column name indicating unique subject identifier for appropraite V-fold cv

make_historical_sl <- function(historical_data, outcome, covariates, id,
                               historical_stack, drop_missing_outcome = TRUE, 
                               fit_sl = TRUE, folds = NULL, V = 5){
 
  dt <- data.table::data.table(historical_data)
  
  task <- sl3::make_sl3_Task(dt, covariates, outcome, id = id,
                             drop_missing_outcome = drop_missing_outcome)
  
  if(fit_sl){
    if(!grepl("Lrnr_cv", class(historical_stack)[1])){ 
      historical_stack <- Lrnr_cv$new(historical_stack, full_fit = T)
    }
  }
  
  fit <- historical_stack$train(historical_stack)
  
  if(fit_sl){
    chained_task <- fit$chain(historical_task)
    
    metalearner_nnls <- make_learner(Lrnr_nnls)
    nnls_fit <- metalearner_nnls$train(chained_task)
    sl_nnls <- make_learner(Pipeline, fit, nnls_fit)
    
    metalearner_nnls_convex <- make_learner(Lrnr_nnls, convex = TRUE)
    nnls_fit_convex <- metalearner_nnls_convex$train(chained_task)
    sl_nnls_convex <- make_learner(Pipeline, fit, nnls_fit_convex)
    
    metalearner_discrete <- make_learner(Lrnr_cv_selector)
    discrete_fit <- metalearner_discrete$train(chained_task)
    sl_discrete <- make_learner(Pipeline, fit, discrete_fit)
    
    return(list(cv_fit = fit, sl_discrete = sl_discrete, sl_nnls = sl_nnls, 
                sl_nnls_convex = sl_nnls_convex, task = historical_task))
  } else {
    return(fit)
  }
}
