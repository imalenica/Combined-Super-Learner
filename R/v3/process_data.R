process_data <- function(X, missingness_name = "delta", strata = NULL,
                         threshold_impute = 0.5, missing_indicator_all = FALSE,
                         impute_factor_level = TRUE, new_level_name = "Unknown") {

  # identify covariates with missingness
  p_missing <- sapply(X, function(x) mean(is.na(x)))

  # variables that have no missingness
  no_missing <- names(p_missing[p_missing == 0])
  processed <- X[, no_missing, with = FALSE]

  # impute nodes with missingness < threshold_impute
  to_impute <- names(p_missing[(0 < p_missing & p_missing <= threshold_impute)])

  if (missing_indicator_all) {
    # add missingness indicators for all columns with missingness
    missing <- names(p_missing[(0 < p_missing)])
  } else {
    # add missingness indicators only when missingness < threshold_impute
    missing <- to_impute
  }

  if (length(to_impute) > 0) {
    if (impute_factor_level) {
      X1 <- X[, to_impute, with = FALSE]
      facs <- colnames(X1[, unlist(lapply(X1, is.factor)), with = FALSE])
      imputed_facs <- X[, lapply(.SD, function(x){
        impute_new_level(x, new_level_name)
      }), .SDcols = facs]

      not_facs <- colnames(X1[, !which(colnames(X1) %in% facs), with = FALSE])
      if (is.null(strata)) {
        imputed_not_facs <- X[, lapply(.SD, impute_by_type), .SDcols = not_facs]
      } else {
        X$index <- 1:nrow(X)
        not_facs2 <- append(not_facs, "index")
        not_facs_imputed <- X[, lapply(.SD, impute_by_type),.SDcols = not_facs2,
                              by = strata]
        imputed_not_facs <- setorder(not_facs_imputed,index)[,not_facs,with = F]
      }

      imputed <- cbind(imputed_facs, imputed_not_facs)

    } else {
      if (is.null(strata)) {
        imputed <- X[, lapply(.SD, impute_by_type), .SDcols = to_impute]
      } else {
        X$index <- 1:nrow(X)
        to_impute2 <- append(to_impute, "index")
        imputed <- X[, lapply(.SD, impute_by_type), .SDcols = to_impute2,
                     by = strata]
        imputed <- setorder(imputed, index)[, to_impute, with = FALSE]
      }
    }
    processed <- cbind(processed, imputed)
  }

  if (length(missing) > 0) {
    missing_indicators <- X[, lapply(.SD, function(x) as.numeric(!is.na(x))),
                            .SDcols = missing]
    missing_names <- sprintf(paste0(missingness_name,"_%s"), missing)
    setnames(missing_indicators, missing_names)
    processed <- cbind(processed, missing_indicators)
  }
  return(processed)
}

impute_median <- function(x) {
  value <- median(as.numeric(x[!is.na(x)]))
  x[is.na(x)] <- value
  x
}

impute_mode <- function(x){
  count_df <- aggregate(count ~ x, data = data.frame(count = 1, x = x), sum)
  value <- count_df$x[which.max(count_df$count)]
  x[is.na(x)] <- value
  return(x)
}

impute_new_level <- function(x, new_level_name) {
  existing_levels <- levels(x)
  # if(new_level_name %in% existing_levels){
  #   stop("New level name detected in existing factor levels, modify `new_level_name`")
  # }
  setattr(x, "levels", c(existing_levels, new_level_name))
  x[is.na(x)] <- new_level_name
  return(x)
}

impute_by_type <- function(x) {
  if (is.factor(x) || is.character(x)) {
    return(impute_mode(x))
  } else {
    return(impute_median(x))
  }
}
