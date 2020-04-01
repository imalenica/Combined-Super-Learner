# residuals: vector of residuals
get_acf <- function(residuals){
  residuals <- data.frame(residuals, time = seq(1:length(residuals)))
  tbl <- tsibble::as_tsibble(residuals, index = time)
  acfs <- tbl %>% feasts::ACF(residuals, lag_max = nrow(residuals))
  sig_acfs <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(tbl$residuals)))
  return(min.true(acfs$acf > sig_acfs))
}

min.true <- function(x){
  min(which(x==!TRUE))-1
}
