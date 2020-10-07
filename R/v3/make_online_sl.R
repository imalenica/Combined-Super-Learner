make_online_sl <- function(adapt_fit, complete_forecast_data, outcome, 
                           time = "min", risk_window = NULL, first_fit = FALSE, 
                           loss_table = NULL, onlineSL_table = NULL,
                           centered_sqloss_table = NULL, n = NULL){
  
  if(!first_fit){
    if(is.null(loss_table)){
      stop("Past risks must be provided if not first fit")
    }
    if(is.null(centered_sqloss_table)){
      stop("Sq. loss centered must be provided if not first fit")
    }
    if(is.null(onlineSL_table)){
      stop("Past online SL must be provided if not first fit.")
    }
    if(is.null(adapt_fit)){
      stop("Adapt fit with onlineSL forecasts must be provided if not first fit.")
    }
  }
  # retain truth, make sure it corresponds to the what we were forecasting
  complete_forecast_data <- data.table::data.table(complete_forecast_data)
  if(!all.equal(adapt_fit$indiviual_forecast_data[[time]], 
                complete_forecast_data[[time]])){
    stop("Adapt fit forecasts correspond to times that is are not equal to ",
         "complete forecast data times.")
  }
  
  if(!first_fit){
    past_onlineSL_forecasts <- adapt_fit$onlineSL_forecasts
    if(!all.equal(past_onlineSL_forecasts$forecast_time, 
                  complete_forecast_data[[time]])){
      stop("Past online SL forecasts correspond to times that is are not equal ",
           "to complete forecast data times.")
    }
  }
  
  truth <- complete_forecast_data[[outcome]]
  
  # data table with forecast data from all learners
  forecasts <- data.table(adapt_fit$all_forecasts)
  loss <- data.table(t(apply(forecasts, 2, function(pred) mean((pred-truth)^2))))
  
  # add to old risk table
  if(!first_fit){
    loss_tbl <- data.table(rbind(loss_table, loss, fill = TRUE))
  } else {
    loss_tbl <- loss
  }
  
  ########################### select new onlineSLs #############################
  
  # initiate output
  minStart <- min(complete_forecast_data[["min"]])
  minEnd <- max(complete_forecast_data[["min"]])
  ID <- unique(complete_forecast_data[["id"]])
  summ <- data.table(t(c(ID = ID, minStart = minStart, minEnd = minEnd)))
  
  # what if we only care about choosing the onlineSL wrt recent previous risks
  if(!all(is.null(risk_window))){
    # only include windows that are less than the current loss table
    risk_window <- risk_window[-which(sapply(risk_window, as.numeric) >=
                                        nrow(loss_tbl))]
    
    if(length(risk_window) > 1){
      if(length(risk_window) > length(centered_sqloss_table)){
        id_rep <- which.max(lapply(centered_sqloss_table, nrow))
        len <- length(centered_sqloss_table)
        name <- names(risk_window)[-which(names(risk_window) %in% 
                                            names(centered_sqloss_table))]
        centered_sqloss_table <- c(centered_sqloss_table, 
                                   centered_sqloss_table[id_rep])
        names(centered_sqloss_table)[(len+1)] <- name
      }
      new <- lapply(seq_along(risk_window), function(j){
        window <- risk_window[[j]]
        tbl_idx <- which(names(centered_sqloss_table) == names(risk_window)[j])
        tbl <- centered_sqloss_table[[tbl_idx]]
        
        new <- make_new_onlineSL(risk_window = window, loss_tbl = loss_tbl, 
                                 loss = loss, n = n, first_fit = first_fit,
                                 centered_sqloss_table = tbl)
        if(is.null(window)){
          return(new)
        }
        colnames(new$new_onlineSL) <- paste0(colnames(new$new_onlineSL), 
                                             "_window", window)
        return(new)
      })
      new_oSL <- do.call(cbind, lapply(new, '[[', 'new_onlineSL'))
      new_onlineSL <- data.table(summ, new_oSL)
      online_CIrisk_upper <- do.call(cbind, 
                                     lapply(new, '[[', 'online_CIrisk_upper'))
      centered_sqloss_table <- lapply(new, '[[', 'centered_sqloss_table')
      names(centered_sqloss_table) <- names(risk_window)
       
    } else {
      if(is.list(centered_sqloss_table)){
        centered_sqloss_table <- unlist(centered_sqloss_table)
      }
      new <- make_new_onlineSL(risk_window = unlist(risk_window), 
                               loss_tbl = loss_tbl, loss = loss, n = n, 
                               first_fit = first_fit, 
                               centered_sqloss_table = centered_sqloss_table)
      if(!is.null(risk_window)){
        colnames(new$new_onlineSL) <- paste0(colnames(new$new_onlineSL),
                                             "_window", risk_window)
      }
      new_oSL <- new$new_onlineSL
      new_onlineSL <- data.table(summ, new_oSL)
      online_CIrisk_upper <- new$online_CIrisk_upper
      centered_sqloss_table <- list(new$centered_sqloss_table)
      names(centered_sqloss_table) <- names(risk_window)
    }
  } else {
    new <- make_new_onlineSL(risk_window = NULL, loss_tbl = loss_tbl, 
                             loss = loss, n = n, first_fit)
    new_oSL <- new$new_onlineSL
    new_onlineSL <- data.table(summ, new_oSL)
    online_CIrisk_upper <- new$online_CIrisk_upper
    centered_sqloss_table <- list(new$centered_sqloss_table)
    names(centered_sqloss_table) <- "NULL"
  }
  
  names_new_onlineSL <- new_oSL[, -grepl("MSE", colnames(new_oSL)), with=F]

  
  #################### evaluate MSE of onlineSL forecasts ######################
  if(!first_fit){
    # which row is relevant here?
    row_idx <- which(
      as.numeric(onlineSL_table$minEnd)+1 == min(past_onlineSL_forecasts$forecast_time)
    )
    rm <- which(names(adapt_fit$onlineSL_forecasts) == "forecast_time")
    past_onlineSL_forecasts <- past_onlineSL_forecasts[-rm]
    
    # update recent online SL with MSE 
    recent <- cbind(
      onlineSL_table[row_idx, c("ID", "minStart", "minEnd"), with=F], 
      do.call(cbind, lapply(seq_along(past_onlineSL_forecasts), function(j){
        mse <- mean((past_onlineSL_forecasts[[j]]-truth)^2)
        forecast_name <- names(past_onlineSL_forecasts)[j]
        name_idx <- which(colnames(onlineSL_table) == forecast_name)
        name <- onlineSL_table[row_idx, name_idx, with = F]
        d <- data.table(name, mse)
        colnames(d)[2] <- paste0("MSE_", forecast_name)
        return(d)
      }))
    )
    
    if(row_idx > 1){
      # retain previous online SLs
      past <- onlineSL_table[-row_idx, ]
      onlineSL_table <- rbind(past, recent, new_onlineSL, fill = T)
    } else {
      onlineSL_table <- rbind(recent, new_onlineSL, fill = T)
    }
    
  } else {
    onlineSL_table <- new_onlineSL
  }
  
  fit_onlineSL <- list(
    sl_ind_weights = adapt_fit$sl_ind_weights,
    sl_weights = adapt_fit$sl_weights,
    individual_fit = adapt_fit$individual_fit,
    fit_individualSL = adapt_fit$fit_individualSL
  )
  
  return(list(
    onlineSL_table = onlineSL_table, 
    loss_table = loss_tbl,
    centered_sqloss_table = centered_sqloss_table, 
    fit_onlineSL = fit_onlineSL,
    online_CIrisk_upper = online_CIrisk_upper, 
    names_onlineSL = names_new_onlineSL
  ))
}

# what if we only care about choosing the onlineSL wrt recent previous risks
make_new_onlineSL <- function(risk_window = NULL, loss_tbl, loss, n, 
                              first_fit, centered_sqloss_table){

  if(!is.null(risk_window) && risk_window < nrow(loss_tbl)){
    loss_tbl_relevant <- tail(loss_tbl, risk_window)
  } else {
    loss_tbl_relevant <- loss_tbl
  }
  
  if(!is.null(risk_window) && risk_window < nrow(centered_sqloss_table)){
    centered_sqloss_table <- tail(centered_sqloss_table, risk_window)
  } else {
    centered_sqloss_table <- centered_sqloss_table
  }
  
  
  if(is.null(n)){
    n <- apply(loss_tbl_relevant, 2, function(rows) sum(!is.na(rows)))
  }
  
  online_cv_risk <- apply(loss_tbl_relevant, 2, sum, na.rm = TRUE)
  online_cv_risk_standardized <- data.table(t(online_cv_risk/n))
  
  if(!first_fit){
    if(any(n == 1)){
      idx <- which(n == 1)
      online_cv_risk_standardized[1,idx] <- 0
    }
    online_cv_risk_standardized <- online_cv_risk_standardized[,which(
      colnames(online_cv_risk_standardized) %in% colnames(loss)
    ), with=F]
    current_centered_sqloss <- (loss-online_cv_risk_standardized)^2
    centered_sqloss_table <- data.table(
      rbind(centered_sqloss_table, current_centered_sqloss, fill = TRUE)
    )
  } else {
    current_centered_sqloss <- (loss-0)^2
    centered_sqloss_table <- data.table(current_centered_sqloss)
  }
  online_centered_sqloss <- apply(centered_sqloss_table, 2, sum, na.rm = TRUE)
  online_se <- sqrt(online_centered_sqloss/n)
  online_CIrisk_upper <- data.table(loss + qnorm(.975)*online_se/sqrt(n))
  
  # isolate risks from historical and individual learners
  historical_lrnrs <- colnames(loss)[grep("historical_", colnames(loss))]
  risk_hist <- online_cv_risk_standardized[, historical_lrnrs, with=F]
  CIrisk_hist <- online_CIrisk_upper[, historical_lrnrs, with=F]
  
  individual_lrnrs <- c(colnames(loss)[grep("individual_", colnames(loss))], 
                        colnames(loss)[grep("indSL_", colnames(loss))])
  risk_ind <- online_cv_risk_standardized[, individual_lrnrs, with=F]
  CIrisk_ind <- online_CIrisk_upper[, individual_lrnrs, with=F]
  
  risk_adapt <- online_cv_risk_standardized # shorten name for convenience
  CIrisk_adapt <- online_CIrisk_upper
  onlineSL_adapt <- colnames(risk_adapt)[which.min(risk_adapt)]
  onlineSL_ind <- colnames(risk_ind)[which.min(risk_ind)]
  onlineSL_hist <- colnames(risk_hist)[which.min(risk_hist)]
  onlineSL.CI_adapt <- colnames(CIrisk_adapt)[which.min(CIrisk_adapt)]
  onlineSL.CI_ind <- colnames(CIrisk_ind)[which.min(CIrisk_ind)]
  onlineSL.CI_hist <- colnames(CIrisk_hist)[which.min(CIrisk_hist)]
  new_onlineSL <- data.table(t(c(
    onlineSL_adapt = onlineSL_adapt, MSE_onlineSL_adapt = NA,
    onlineSL_individual = onlineSL_ind, MSE_onlineSL_individual = NA,
    onlineSL_historical = onlineSL_hist, MSE_onlineSL_historical = NA,
    onlineSL.CI_adapt = onlineSL.CI_adapt, MSE_onlineSL.CI_adapt = NA,
    onlineSL.CI_individual = onlineSL.CI_ind, MSE_onlineSL.CI_individual = NA, 
    onlineSL.CI_historical = onlineSL.CI_hist, MSE_onlineSL.CI_historical = NA
  )))
  return(list(new_onlineSL = new_onlineSL, 
              centered_sqloss_table = centered_sqloss_table,
              online_CIrisk_upper = online_CIrisk_upper))
}
