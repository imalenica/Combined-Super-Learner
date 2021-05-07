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


get_weights_solnp <- function(pred, observed, outcome_type, make_sparse = TRUE, 
                              convex_combination = TRUE, init_0 = FALSE, 
                              tol = 1e-5){
                              
  if(outcome_type$type %in% c("constant", "binomial")){
    risk <- function(alphas) {
      preds <- plogis(trim_logit(pred) %*% alphas)
      losses <- (-1 * ifelse(observed == 1, log(bound(preds)), log(1 - bound(preds))))
      return(mean(losses))
    }
  } else if (outcome_type$type == "continuous"){
    risk <- function(alphas) {
      preds <- pred %*% alphas
      losses <- (observed - preds)^2
      return(mean(losses))
    }
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
sl_weights_fit <- function(pred, obs, outcome_type){
  
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
  weights_solnp_convex <- get_weights_solnp(
    pred = as.matrix(pred_noNA), observed = obs, convex_combination=T, 
    outcome_type = outcome_type
  )
  weights_solnp <- get_weights_solnp(
    as.matrix(pred_noNA),  observed = obs, convex_combination=F, 
    outcome_type = outcome_type
  )
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
sl_weights_predict <- function(sl_weights, pred, outcome_bounds){
  
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
  sl_pred <- data.table::data.table(apply(sl_pred, 2, bound, outcome_bounds))
  colnames(sl_pred) <- c("nnls_convex", "nnls", "solnp_convex", "solnp")
  return(sl_pred)
}

# we include pred, since oSL selected an SL wrt to those preds
retain_sl_weights <- function(online_sl_name, learner_sl_weights, pred){
  
  # which metalearner did the online sl select
  if(grepl("nnls_convex", online_sl_name)){
    online_sl_metalearner <- "weights_nnls_convex"
  } else if (grepl("nnls", online_sl_name)){
    online_sl_metalearner <- "weights_nnls"
  } else if (grepl("solnp_convex", online_sl_name)){
    online_sl_metalearner <- "weights_solnp_convex"
  } else if (grepl("solnp", online_sl_name)){
    online_sl_metalearner <- "weights_solnp"
  }
  
  # ensure correspondence between weights and predictions
  cols <- colnames(pred)[which(colnames(pred) %in% colnames(learner_sl_weights))]
  learner_sl_weights <- learner_sl_weights[, cols]
  pred <- pred[, cols, with=F]
  
  # remove NA predictions, weights
  na_idx <- which(colSums(is.na(pred)) > 0)
  if(length(na_idx) > 0){
    wts_noNA <- learner_sl_weights[,-na_idx]
  } else {
    wts_noNA <- learner_sl_weights
  }
  
  return(wts_noNA[online_sl_metalearner,])
}

################################################################################
osl_weights_predict <- function(osl_weights, pred, outcome_bounds){
  
  # ensure correspondence between weights and predictions
  cols <- colnames(pred)[which(colnames(pred) %in% names(osl_weights))]
  osl_weights <- osl_weights[cols]
  pred <- pred[, cols, with=F]
  
  # remove NA predictions, weights
  na_idx <- which(colSums(is.na(pred)) > 0)
  if(length(na_idx) > 0){
    pred_noNA <- pred[,-na_idx,with=F]
    wts_noNA <- osl_weights[-na_idx]
  } else {
    pred_noNA <- pred
    wts_noNA <- osl_weights
  }
  
  # use non-NA osl weights for prediction with osl
  bound(as.numeric(as.matrix(pred_noNA) %*% wts_noNA), outcome_bounds)
}

# utility function
bound <- function(preds, bounds = 0.001) {
  lower <- bounds[[1]]
  if (length(bounds) > 1) {
    upper <- bounds[[2]]
  } else {
    upper <- 1 - lower
  }
  preds_bounded <- pmin(pmax(preds, lower), upper)
  return(preds_bounded)
}
trim_logit <- function(x, trim = 1e-05) {
  x[x < trim] <- trim
  x[x > (1 - trim)] <- (1 - trim)
  foo <- log(x / (1 - x))
  return(foo)
}

