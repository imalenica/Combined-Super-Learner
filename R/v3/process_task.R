# solves issue with mismatch between historical and individual delta column
process_task <- function(individual_training_task, historical_task){
  
  hist_cols <- names(historical_task$data)
  ind_cols <- names(individual_training_task$data)
  ind_cols_miss <- hist_cols[!(hist_cols %in% ind_cols)]

  if(length(ind_cols_miss) > 0){
  individual_training_data <- individual_training_task$data
  ind_cols_miss_extra <- data.frame(matrix(nrow=nrow(individual_training_data), 
                                           ncol=length(ind_cols_miss)))
  names(ind_cols_miss_extra) <- ind_cols_miss[seq(ind_cols_miss)]
  ind_cols_miss_extra[,ind_cols_miss] <- rep(0, nrow(individual_training_data))
  ind_train_data_new <- cbind.data.frame(individual_training_data, 
                                         ind_cols_miss_extra)
  
  processed_task <- make_sl3_Task(
    data = data.table(ind_train_data_new), 
    covariates = historical_task$nodes$covariates,
    outcome = individual_training_task$nodes$outcome,
    folds = individual_training_task$folds,
    id = "subject_id"
  )
  } else {
    processed_task <- individual_training_task
  }
  return(processed_task)
}