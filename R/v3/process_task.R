# solves issue with mismatch between historical and individual delta column
process_task <- function(individual_training_task, historical_fit){
  
  hist_cols <- names(historical_fit$task$data)
  ind_cols <- names(individual_training_task$data)
  ind_cols_miss <- hist_cols[!(hist_cols %in% ind_cols)]
  delta_cols <- grep("delta_", ind_cols_miss)
  
  if(length(ind_cols_miss) > 0){
  individual_training_data <- individual_training_task$data
  ind_cols_miss_extra <- data.frame(matrix(nrow=nrow(individual_training_data), 
                                           ncol=length(ind_cols_miss)))
  names(ind_cols_miss_extra) <- ind_cols_miss[seq(ind_cols_miss)]
  ind_cols_miss_extra[,ind_cols_miss] <- rep(0, nrow(individual_training_data))
  
  #training_task$add_columns(ind_cols_miss_extra)
  ind_train_data_new <- cbind.data.frame(individual_training_data, 
                                         ind_cols_miss_extra)
  ind_train_data_new <- as.data.frame(ind_train_data_new)
  
  #Arrange to match historical:
  output <- ind_train_data_new[,hist_cols]
  
  processed_task <- make_sl3_Task(
    data =  output, 
    covariates = historical_fit$task$nodes$covariates,
    outcome = outcome, 
    folds = individual_training_task$folds
  )
  } else {
    processed_task <- individual_training_task
  }
  return(processed_task)
}