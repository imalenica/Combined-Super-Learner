# solves issue with mismatch between historical and individual delta column
process_task <- function(training_task, historical_task, verbose = FALSE){
  
  historical_cols <- colnames(historical_task$data)
  training_cols <- colnames(training_task$data)
  training_cols_missing <- historical_cols[!(historical_cols %in% training_cols)]
  
  if(length(training_cols_missing) > 0){
    
    # print missing columns in training data
    if(verbose){
      print(paste(c("Columns to be added as 0 vectors to training data, ",
                    "to avoid mismatch with historically-trained learners:", 
                    training_cols_missing), collapse=" "))
    }
    
    # make new data with additional columns that were previously missing
    training_data <- training_task$data
    missing_data <- data.frame(
      matrix(nrow = nrow(training_data), ncol = length(training_cols_missing))
    )
    names(missing_data) <- training_cols_missing[seq(training_cols_missing)]
    missing_data[,training_cols_missing] <- rep(0, nrow(training_data))
    
    # put together with training data
    training_data_processed <- cbind.data.frame(training_data, missing_data)
    
  } else {
    training_data_processed <- training_task$data
  }
  
  # make sure the training task has the nodes from the historical task
  historical_covs <- historical_task$nodes$covariates
  training_covs <- training_task$nodes$covariates
  if(!all(historical_covs %in% training_covs)){
    # print missing covariates in training task
    if(verbose){
      missing <- historical_covs[-which(historical_covs %in% training_covs)]
      print(paste(c("Covariates to be added to training task covariates, ",
                    "to avoid mismatch with covariates used for training ",
                    "historically-trained learners:", missing), collapse=" "))
    }
    processed_covs <- unique(c(training_covs, historical_covs))
  } else {
    processed_covs <- training_covs
  }
  
  # make new task
  return(sl3::make_sl3_Task(
    data = data.table(training_data_processed), 
    covariates = processed_covs,
    outcome = training_task$nodes$outcome,
    id = training_task$nodes$id,
    time = training_task$nodes$time,
    weights = training_task$nodes$weights,
    offset = training_task$nodes$offset,
    folds = training_task$folds,
    outcome_type = training_task$outcome_type$type
  ))
}
