### Functions to support simulations

run_posl <- function(data, covs, outcome, learners, splits){
  
  data_full <- data$full
  data_historical <- data$historical
  data_individual <- data$individual
  
  #Save all results
  osl_fit                 <- list()
  result_list             <- list()
  
  result_weights          <- list()
  result_ind_weights      <- list()
  
  result_forecasts        <- list()
  result_ind_forecasts    <- list()
  result_Vfold_forecasts   <- list()
  result_hist_forecasts   <- list()
  result_online_forecasts <- list()
  
  loss                    <- list()
  
  ### Learn the full fit with V-fold CV
  full_fit <- make_historical_sl(
    historical_data = data_full,
    outcome = outcome, 
    covariates = covs, 
    id = "id", 
    historical_stack = learners
  )
  
  ### Learn the historical fit
  historical_fit <- make_historical_sl(
    historical_data = data_historical,
    outcome = outcome, 
    covariates = covs, 
    id = "id", 
    historical_stack = learners
  )

  split_data <- lapply(splits, function(x) data.table(data_individual[1:x,]))
  split_datas <- lapply(splits, function(x){
    data_full %>%
      dplyr::group_by(id) %>%
      dplyr::slice(c(1:x))
  })
    
  for(i in 1:length(splits)){
    
    forecast_data <- data_individual[c((nrow(split_data[[i]])+1):(nrow(split_data[[i]])+5)),]
    
    #Train Online SL:
    osl_fit[[i]] <- make_online_sl(
      online_data=split_datas[[i]],
      outcome=outcome,
      covariates=covs,
      id = "id",
      online_stack=learners
      )
    
    #Train Personalized Online SL:
    result_list[[i]] <- make_adapt_sl(
      individual_training_data = split_data[[i]], 
      indiviual_forecast_data = forecast_data,
      outcome = outcome, 
      covariates = covs, 
      id = "id",
      historical_fit = historical_fit,
      individual_stack = learners,
      return_individual_fit=TRUE
    )
    
    ##################################################
    ###  Evaluate all performance on forecast task
    ##################################################
    forecast_task <- result_list[[i]]$forecast_task
    
    #Combined model fit on forecast data:
    result_weights[[i]] <- t(result_list[[i]]$sl_weights)
    result_forecasts[[i]] <- t(result_list[[i]]$sl_forecasts)
    
    #Individual model fit on forecast data:
    result_ind_weights[[i]] <- t(result_list[[i]]$sl_ind_weights)
    result_ind_forecasts[[i]] <- t(result_list[[i]]$sl_ind_forecasts)
    
    #Historical model fit on forecast data:
    result_hist_forecasts[[i]] <- t(cbind.data.frame(discreteSL=historical_fit$sl_discrete$predict(forecast_task),
                                                     nnls_convexSL=historical_fit$sl_nnls_convex$predict(forecast_task),
                                                     nnls_SL=historical_fit$sl_nnls$predict(forecast_task)))
    #Full model fit on forecast data (V-fold):
    result_Vfold_forecasts[[i]] <- t(cbind.data.frame(discreteSL=full_fit$sl_discrete$predict(forecast_task),
                                                     nnls_convexSL=full_fit$sl_nnls_convex$predict(forecast_task),
                                                     nnls_SL=full_fit$sl_nnls$predict(forecast_task)))
    #Online SL fit on forecast data:
    result_online_forecasts[[i]] <- t(cbind.data.frame(discreteSL=osl_fit[[i]]$sl_discrete$predict(forecast_task),
                                                       nnls_convexSL=osl_fit[[i]]$sl_nnls_convex$predict(forecast_task),
                                                       nnls_SL=osl_fit[[i]]$sl_nnls$predict(forecast_task)))
    
    #Evaluate overall loss on the forecasts:
    truth <- forecast_data$Y
    mse   <- function(pred) {mean((pred-truth)^2)}
    loss[[i]] <- cbind.data.frame(VfoldSL = apply(result_Vfold_forecasts[[i]], 1, mse),
                                  histSL = apply(result_hist_forecasts[[i]], 1, mse),
                                  indSL  = apply(result_ind_forecasts[[i]], 1, mse),
                                  onlSL  = apply(result_online_forecasts[[i]], 1, mse),
                                  persSL = apply(result_forecasts[[i]], 1, mse))
    }
  
  #Save only results for POSL
  res_discrete <- lapply(result_weights, function(x) x[,1])
  res_discrete <- do.call(cbind, res_discrete)
  colnames(res_discrete) <- paste0("t=", splits)
  
  res_nnls_convex <- lapply(result_weights, function(x) x[,2])
  res_nnls_convex <- data.frame(do.call(cbind, res_nnls_convex))
  colnames(res_nnls_convex) <- paste0("t=", splits)
  
  res_nnls <- lapply(result_weights, function(x) x[,3])
  res_nnls <- do.call(cbind, res_nnls)
  colnames(res_nnls) <- paste0("t=", splits)
  
  return(list(res_discrete=res_discrete,
              res_nnls_convex=res_nnls_convex,
              res_nnls=res_nnls,
              loss=loss))
 
}

folds_rolling_origin_pooled_edit <- function(n, t, id = NULL, time = NULL,
                                        first_window, validation_size,
                                        gap = 0, batch = 1) {
  if ((!is.null(id) & is.null(time)) | (is.null(id) & !is.null(time))) {
    stop("Cannot create flexible folds (allow for variability in the amount of 
         time observed for each id) unless both `time` and `id` argments are 
         provided. Either provide both `time` and `id` or neither.")
  }
  if ((length(id) != length(time)) & (length(id) > 1)) {
    stop("Cannot create flexible folds (allow for variability in the amount of 
         `time` observed for each `id`) unless `time` vector is of same length 
         as `id` vector. `time` is a vector of integers of time points observed 
         for each subject, and `id` is a vector of unique identifiers which
         correspond to the time vector. The `id` vector is used to subset the 
         `time` vector.")
  }
  if (length(id) == 1) {
    id <- rep(id, length(time))
  }
  
  if (is.null(id) & is.null(time)) {
    dat <- cbind.data.frame(
      time = rep(seq(t), n / t),
      id = rep(seq(n / t), each = t)
    )
  } else {
    # Index times by id (allows variability in time observed for each subject)
    dat <- cbind.data.frame(time = time, id = id)
  }
  
  ids <- unique(dat$id)
  times <- unique(dat$time)
  
  message(paste("Processing", length(ids), "samples with", t, "time points."))
  
  # establish rolling origin forecast for time-series cross-validation
  rolling_origin_skeleton <- folds_rolling_origin(
    t, first_window,
    validation_size, gap, batch
  )
  
  folds_rolling_origin <- lapply(rolling_origin_skeleton, function(fold) {
    train_times <- training(times)
    valid_times <- validation(times)
    train_idx <- which(dat$time%in%train_times)
    valid_idx <- which(dat$time%in%valid_times)
    fold <- make_fold(fold_index(), train_idx, valid_idx)
  })
  return(folds_rolling_origin)
}
