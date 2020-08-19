get_weights <- function(pred, observed, loss, convex = FALSE, 
                        discrete = FALSE){
  
  if(discrete) {
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

sl_weights_fit <- function(pred, observed, loss){
  
  # get sl weights
  weights_discrete <- get_weights(pred, observed, loss, discrete=T) 
  weights_nnls_convex <- get_weights(pred, observed, loss, discrete=F, convex=T)
  weights_nnls <- get_weights(pred, observed, loss, discrete=F, convex=F)
  
  sl_weights <- suppressWarnings(rbind(weights_discrete, weights_nnls_convex, 
                                       weights_nnls))
  colnames(sl_weights) <- names(pred)
  return(sl_weights)
}

sl_weights_predict <- function(sl_weights, pred)
  # use sl weights for prediction with sl
  pred_discrete <- as.matrix(pred) %*% sl_weights["weights_discrete",]
  pred_nnls_convex <- as.matrix(pred) %*% sl_weights["weights_nnls_convex",]
  pred_nnls <- as.matrix(pred) %*% sl_weights["weights_nnls",]
  sl_pred <- data.frame(pred_discrete, pred_nnls_convex, pred_nnls)
  colnames(sl_pred) <- c("discreteSL", "nnls_convexSL", "nnls_SL")
  return(sl_pred)
}
  