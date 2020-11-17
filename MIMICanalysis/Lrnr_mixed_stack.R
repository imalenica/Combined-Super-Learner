Lrnr_mixed_stack <- R6Class(
  classname = "Lrnr_mixed_stack", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(mixed_stack, name_ts_stack, name_non_cv_stack, ...) {
      # learner is an already initialized learner
      params <- list(
        mixed_stack, ...
      )
      super$initialize(params = params, ...)
    }
  ),
  
  private = list(

    .train = function(task) {
      mixed_stack <- self$params$mixed_stack
      fit_object <- mixed_stack$train(task)
      return(fit_object)
    },
    
    .predict = function(task = NULL) {
      mixed_fits <- self$fit_object$learner_fits
      
      # seperate the fits 
      ts_fit <- mixed_fits[[self$params$name_ts_stack]]
      non_cv_fit <- mixed_fits[[self$params$name_non_cv_stack]]

      pred_list <- lapply(task$folds, function(fold){
        v_task <- validation_task(training_task, fold)
        ts_preds <- ts_fit$predict_fold(v_task, fold_number = fold$v)
        non_cv_preds <- non_cv_fit$predict_fold(v_task, fold_number = "full")
        cbind(ts_preds, non_cv_preds)
      })
      
      
      
      learner_dict <- self$fit_object
      variable_stratify_stratas <- as.numeric(names(learner_dict))
      variable_stratify <- self$params$variable_stratify
      
      X_new <- as.matrix(task$X)
      variable_stratify_stratas_new <- unique(X_new[, variable_stratify])
      if (
        length(
          setdiff(variable_stratify_stratas_new, variable_stratify_stratas)
        ) > 0
      ) {
        stop("There is new strata in the prediction data that is not present in
             training data!")
      }
      
      prediction_df_dict <- list()
      # predictions <- aorder(results$predictions, order(results$index))
      
      for (strata in variable_stratify_stratas_new) {
        index_subtask <- which(X_new[, variable_stratify] == strata)
        # construct subtask
        sub_task <- task$subset_task(row_index = index_subtask)
        sub_task <- sub_task$next_in_chain(
          covariates = sub_task$nodes$covariates[
            sub_task$nodes$covariates != variable_stratify
            ]
        )
        # predict on the subtask
        prediction_subtask <- learner_fit_predict(
          learner_dict[[as.character(strata)]],
          sub_task
        )
        result <- list(
          prediction = prediction_subtask,
          original_index = index_subtask
        )
        prediction_df_dict[[as.character(strata)]] <- result
      }
      results <- apply(do.call(rbind, prediction_df_dict), 2, as.list)
      results <- origami::combine_results(results)
      
      predictions <- aorder(results$prediction, order(results$original_index))
      return(predictions)
      },
    .required_packages = NULL
  )
)