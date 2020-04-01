get_weights <- function(pred, observed, loss, convex = FALSE, discrete = FALSE){
  
  if(discrete == TRUE) {
    fit_coef <- as.numeric(loss)
    fit_coef[which.min(loss)] <- 1
    fit_coef[-(which.min(loss))] <- 0
  } else if(discrete == FALSE & convex == TRUE) {
    fit_coef <- lsei::pnnls(as.matrix(pred), as.matrix(observed), sum = 1)
    fit_coef <- fit_coef$x
  } else if(discrete == FALSE & convex == FALSE) {
    fit_coef <- nnls::nnls(as.matrix(pred), as.numeric(observed))
    fit_coef <- fit_coef$x
  }
  return(fit_coef)
}