# solves issue with mismatch between historical and individual delta column
process_task <- function(individual_training_task, historical_fit){
  
  hist_cols <- historical_fit$task$nodes$covariates
  ind_cols <- individual_training_task$nodes$covariates
  ind_cols_miss <- hist_cols[!(hist_cols %in% ind_cols)]
  delta_cols <- grep("delta_", ind_cols_miss)
  
  if(length(delta_cols) > 0){
  individual_training_data <- individual_training_task$data
  ind_cols_miss_extra <- data.frame(matrix(nrow=nrow(individual_training_data), 
                                           ncol=length(delta_cols)))
  ind_cols_miss_extra[, delta_cols] <- rep(0, nrow(individual_training_data))
  names(ind_cols_miss_extra) <- ind_cols_miss[delta_cols]
  #training_task$add_columns(ind_cols_miss_extra)
  ind_train_data_new <- cbind.data.frame(individual_training_data, 
                                         ind_cols_miss_extra)
  processed_task <- make_sl3_Task(
    data = ind_train_data_new, 
    covariates = c(ind_cols, names(ind_cols_miss_extra)),
    outcome = outcome, 
    folds = folds
  )
  } else {
    processed_task <- individual_training_task
  }
  return(processed_task)
}