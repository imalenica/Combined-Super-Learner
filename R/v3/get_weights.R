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


get_weights_solnp <- function(pred, observed, make_sparse = TRUE, 
                              convex_combination = TRUE, init_0 = FALSE, 
                              tol = 1e-5){
  
  risk <- function(alphas) {
    preds <- pred %*% alphas
    losses <- (preds - observed)^2
    return(mean(losses))
  }
  
  if (convex_combination) {
    eq_fun <- function(alphas) {
      sum(alphas)
    }
    eqB <- 1
    LB <- rep(0L, ncol(pred))
  } else {
    LB <- eqB <- eq_fun <- NULL
  }
  
  p <- ncol(pred)
  if (init_0) {
    init_alphas <- rep(0, p)
  } else {
    init_alphas <- rep(1 / p, p)
  }
  
  fit_object <- Rsolnp::solnp(init_alphas, risk, eqfun = eq_fun, eqB = eqB,
                              LB = LB, control = list(trace = 0, tol = tol))
  coefs <- fit_object$pars
  names(coefs) <- colnames(pred)
  
  if (make_sparse){
    threshold <- (max(coefs)) / 1000
    coefs[coefs < threshold] <- 0
    coefs <- coefs / sum(coefs)
  }
  return(as.numeric(coefs))
}
################################################################################
sl_weights_fit <- function(pred, obs){
  
  # remove NA predictions, loss
  na_idx <- which(colSums(is.na(pred)) > 0)
  if(length(na_idx) > 0){
    weight_NA <- matrix(data = 0, nrow = 4, ncol = length(na_idx))
    pred_noNA <- pred[, -na_idx, with = F]
  } else {
    pred_noNA <- pred
  }
  
  # get sl weights
  weights_nnls_convex <- get_weights(pred_noNA, obs, discrete=F, convex=T)
  weights_nnls <- get_weights(pred_noNA, obs, discrete=F, convex=F)
  weights_solnp_convex <- get_weights_solnp(as.matrix(pred_noNA), obs, convex_combination=T)
  weights_solnp <- get_weights_solnp(as.matrix(pred_noNA), obs, convex_combination=F)
  sl_weights <- suppressWarnings(rbind(weights_nnls_convex, weights_nnls,
                                       weights_solnp_convex, weights_solnp))
  
  if(length(na_idx) > 0){
    sl_weights <- cbind(sl_weights, weight_NA)
    colnames(sl_weights) <- c(names(pred_noNA), names(na_idx))
  } else {
    colnames(sl_weights) <- names(pred_noNA)
  }
  return(sl_weights[, colnames(pred)])
}
################################################################################
sl_weights_predict <- function(sl_weights, pred){
  
  # ensure correspondence between weights and predictions
  cols <- colnames(pred)[which(colnames(pred) %in% colnames(sl_weights))]
  sl_weights <- sl_weights[, cols]
  pred <- pred[, cols, with=F]
  
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
  pred_solnp_convex <- as.matrix(pred_noNA) %*% wts_noNA["weights_solnp_convex",]
  pred_solnp <- as.matrix(pred_noNA) %*% wts_noNA["weights_solnp",]
  sl_pred <- data.table::data.table(pred_nnls_convex, pred_nnls, 
                                    pred_solnp_convex, pred_solnp)
  colnames(sl_pred) <- c("nnls_convexSL", "nnls_SL", "solnp_convexSL", "solnp_SL")
  return(sl_pred)
}
