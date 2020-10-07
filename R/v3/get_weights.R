get_weights <- function(pred, observed, loss=NULL, convex=F, discrete=F){
  
  if(discrete && !is.null(loss)) {
    fit_coef <- as.numeric(loss)
    fit_coef[which.min(loss)] <- 1
    fit_coef[-(which.min(loss))] <- 0
  } else {
    if(convex) {
      fit_coef <- lsei::pnnls(as.matrix(pred), as.matrix(observed), sum = 1)
      fit_coef <- fit_coef$x
    } else {
      fit_coef <- nnls::nnls(as.matrix(pred), as.numeric(observed))
      fit_coef <- fit_coef$x
    }
  }
  return(fit_coef)
}
################################################################################
sl_weights_fit <- function(pred, obs){
  
  # remove NA predictions, loss
  na_idx <- which(colSums(is.na(pred)) > 0)
  if(length(na_idx) > 0){
    weight_NA <- matrix(data = 0, nrow = 2, ncol = length(na_idx))
    pred_noNA <- pred[, -na_idx, with = F]
  } else {
    pred_noNA <- pred
  }
  
  # get sl weights
  weights_nnls_convex <- get_weights(pred_noNA, obs, discrete=F, convex=T)
  weights_nnls <- get_weights(pred_noNA, obs, discrete=F, convex=F)
  sl_weights <- suppressWarnings(rbind(weights_nnls_convex, weights_nnls))
  
  if(length(na_idx) > 0){
    sl_weights <- cbind(sl_weights, weight_NA)
    colnames(sl_weights) <- c(names(pred_noNA), names(na_idx))
  } else {
    colnames(sl_weights) <- names(pred_noNA)
  }
  return(sl_weights)
}
################################################################################
sl_weights_predict <- function(sl_weights, pred){
  
  # remove NA predictions, weights
  na_idx <- which(colSums(is.na(pred)) > 0)
  if(length(na_idx) > 0){
    pred_noNA <- pred[,-na_idx,with=F]
    wts_noNA <- sl_weights[,-na_idx]
  } else {
    pred_noNA <- pred
    wts_noNA <- sl_weights
  }
  
  # use non-NA sl weights for prediction with sl
  pred_nnls_convex <- as.matrix(pred_noNA) %*% wts_noNA["weights_nnls_convex",]
  pred_nnls <- as.matrix(pred_noNA) %*% wts_noNA["weights_nnls",]
  sl_pred <- data.table(pred_nnls_convex, pred_nnls)
  colnames(sl_pred) <- c("nnls_convexSL", "nnls_SL")
  return(sl_pred)
}