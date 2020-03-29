process_data <- function(X) {
  
  # give all columns missing indicators to avoid mismatch in covariate dim
  all <- names(X)
  missing_indicators <- X[, lapply(.SD, function(x) as.numeric(!is.na(x))),
                          .SDcols = all
                          ]
  delta_names <- sprintf("delta_%s", all)
  setnames(missing_indicators, delta_names)
  
  # impute the covariates with missingness
  p_missing <- sapply(X, function(x) mean(is.na(x)))
  
  # nodes that are already complete still have missingness indicators
  no_missing <- names(p_missing[p_missing == 0])
  processed <- X[, no_missing, with = FALSE]
  processed <- cbind(processed, missing_indicators)
  
  # nodes to impute have missingness indicators and are imputed
  to_impute <- names(p_missing[(0 < p_missing)])
  if (length(to_impute) > 0) {
    imputed <- X[, lapply(.SD, impute_by_type), .SDcols = to_impute]
    processed <- cbind(processed, imputed)
  } 
  return(processed)
}

impute_median <- function(x) {
  value <- median(as.numeric(x[!is.na(x)]))
  x[is.na(x)] <- value
  x
}

#' @importFrom stats aggregate
impute_mode <- function(x) {
  count_df <- aggregate(count ~ x, data = data.frame(count = 1, x = x), sum)
  value <- count_df$x[which.max(count_df$count)]
  x[is.na(x)] <- value
  return(x)
}

impute_by_type <- function(x) {
  if (is.factor(x) || is.character(x)) {
    return(impute_mode(x))
  } else {
    return(impute_median(x))
  }
}
