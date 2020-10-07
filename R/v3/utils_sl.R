trainSL <- function(trainSL_fold_fits, truth = NULL,
                    type = c("all", "individual", "historical")){
  
  if(type == "all"){
    name <- 'all_preds'
  } else if(type == "individual"){
    name <- 'ind_preds'
  } else if(type == "historical"){
    name <- 'hist_preds'
  }
  
  if(is.null(truth)){
    truth <- Reduce(c, lapply(trainSL_fold_fits, '[[', 'truth'))
  }
  
  preds <- rbindlist(lapply(trainSL_fold_fits, '[[', name), fill = T)
  weights <- sl_weights_fit(preds, truth)
  return(weights[, colnames(preds)])
}

process_window <- function(window, total_length){
  if(length(window) == 1){
    if(window != "all" && window >= total_length){
      window <- "all"
    }
  } else {
    if(any(window == "all")){
      window <- window[-which(window == "all")]
      is_all <- TRUE
    }
    retain_idx <- which(window < total_length)
    if(length(retain_idx) > 0){
      window <- window[c(retain_idx)]
      if(is_all){
        window <- c(list("all"), window)
      }
    } else {
      window <- list("all")
    }
  }
  return(window)
}

fitSL <- function(window, trainSL_fold_fits, validSL_fold_fit){
  
  if(window != "all" && length(trainSL_fold_fits) > window){
    cutoff <- length(trainSL_fold_fits)-window+1
    trainSL_fold_fits <- trainSL_fold_fits[cutoff:length(trainSL_fold_fits)]
  }
  
  ################################## get weights ###############################
  truth <- Reduce(c, lapply(trainSL_fold_fits, '[[', 'truth'))
  awts <- trainSL(trainSL_fold_fits, truth, "all")
  iwts <- trainSL(trainSL_fold_fits, truth, "individual")
  hwts <- trainSL(trainSL_fold_fits, truth, "historical")
  
  ######################### get predictions using weights ######################
  all_preds <- validSL_fold_fit$all_preds
  acol <- colnames(all_preds)[which(colnames(all_preds) %in% colnames(awts))]
  apred <- sl_weights_predict(awts[,acol], all_preds[, acol, with=F])
  
  ind_preds <- validSL_fold_fit$ind_preds
  icol <- colnames(ind_preds)[which(colnames(ind_preds) %in% colnames(iwts))]
  ipred <- sl_weights_predict(iwts[,icol], ind_preds[, icol, with=F])
  
  hist_preds <- validSL_fold_fit$hist_preds
  hcol <- colnames(hist_preds)[which(colnames(hist_preds) %in% colnames(hwts))]
  hpred <- sl_weights_predict(hwts[,hcol], hist_preds[, hcol, with=F])
  
  if(window != "all"){
    colnames(apred) <- paste0("window", window, "_", colnames(apred))
    colnames(ipred) <- paste0("individualSL", "window", window, "_", colnames(ipred))
    colnames(hpred) <- paste0("historicalSL", "window", window, "_", colnames(hpred))
  } else {
    colnames(ipred) <- paste0("individualSL_", colnames(ipred))
    colnames(hpred) <- paste0("historicalSL_", colnames(hpred))
  }
  
  wts <- list(allSL_weights = awts, indSL_weights = iwts, histSL_weights = hwts)
  preds <- cbind(apred, ipred, hpred)
  
  return(list(weights = wts, preds = preds))
}

fit_onlineSL <- function(window, loss_tbl, current_loss, first_fit){
  
  if(window != "all" && window < nrow(loss_tbl)){
    loss_tbl <- tail(loss_tbl, window)
  }
  
  n <- apply(loss_tbl, 2, function(rows) sum(!is.na(rows)))
  
  online_cv_risk <- apply(loss_tbl, 2, sum, na.rm = TRUE)
  online_cv_risk_standardized <- data.table(t(online_cv_risk/n))
  
  # isolate risks from historical and individual learners
  hist_lrnrs <- colnames(current_loss)[grep("historical", colnames(current_loss))]
  risk_hist <- online_cv_risk_standardized[, hist_lrnrs, with=F]
  
  ind_lrnrs <- c(colnames(current_loss)[grep("individual", colnames(current_loss))])
  risk_ind <- online_cv_risk_standardized[, ind_lrnrs, with=F]
  
  risk_adapt <- online_cv_risk_standardized # shorten name for convenience
  
  onlineSL_adapt <- colnames(risk_adapt)[which.min(risk_adapt)]
  onlineSL_ind <- colnames(risk_ind)[which.min(risk_ind)]
  onlineSL_hist <- colnames(risk_hist)[which.min(risk_hist)]
  
  dt <- data.table(t(c(
    onlineSL_adapt = onlineSL_adapt, MSE_onlineSL_adapt = NA,
    onlineSL_individual = onlineSL_ind, MSE_onlineSL_individual = NA,
    onlineSL_historical = onlineSL_hist, MSE_onlineSL_historical = NA
  )))
  
  if(window != "all"){
    colnames(dt) <- paste0(colnames(dt), "_window", window)
  }
  return(dt)
}
