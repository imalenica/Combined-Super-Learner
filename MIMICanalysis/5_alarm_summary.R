library(openxlsx)
library(data.table)
library(readxl)

# ==============================================================================
process_alarm <- function(alarm_idx, alarm_reset_idx){
  na_idx <- list()
  for(i in 1:length(alarm_idx)){
    if(alarm_idx[i] %in% unlist(na_idx)){
      na_idx[[i]] <- NA
    } else {
      alarm_silence_idx <- alarm_idx[i] + 1:6
      if(any(alarm_idx %in% alarm_silence_idx)){
        if(any(alarm_silence_idx %in% alarm_reset_idx)){
          reset_idx <- min(which(alarm_silence_idx %in% alarm_reset_idx))
          if(reset_idx > 1){
            alarm_silence_idx <- alarm_silence_idx[c(1:(reset_idx-1))]
            na_idx[[i]] <- alarm_idx[which(alarm_idx %in% alarm_silence_idx)]
          } else {
            na_idx[[i]] <- NA
          }
        } else {
          na_idx[[i]] <- alarm_idx[which(alarm_idx %in% alarm_silence_idx)]
        }
      } else {
        na_idx[[i]] <- NA
      }
    }
  }
  na_idx <- unique(unlist(na_idx[!sapply(na_idx, function(x) any(is.na(x)))]))
  if(length(na_idx) > 0){
    alarm_idx <- alarm_idx[!alarm_idx %in% na_idx]
  }
  return(alarm_idx)
}

eval_AUC <- function(pred_prob, labels, cutoffs){
  stopifnot(length(pred_prob) == length(labels))

  AUC_mat <- matrix(nrow = length(cutoffs), ncol = 2)
  AUCPR_mat <- matrix(nrow = length(cutoffs), ncol = 2)
  for(i in 1:length(cutoffs)){
    alarm <- rep(0, length(pred_prob))
    alarm_idx <- which(pred_prob <= cutoffs[i])
    if(length(alarm_idx) >= 1){
      alarm_reset_idx <- which(pred_prob > cutoffs[i])
      alarm_idx <- process_alarm(alarm_idx, alarm_reset_idx)
      if(length(alarm_idx) >= 1){
        alarm[alarm_idx] <- 1
        # label alarm as correct if it occurred no more than 10 min before label
        for(h in 1:(length(alarm)-1)){
          if(alarm[h] == 1 & labels[h] == 0){
            if(labels[h+1] == 1){
              alarm[h] <- 0
              alarm[h+1] <- 1
            }
            if(labels[h+2] == 1 & h < length(alarm)-1){
              alarm[h] <- 0
              alarm[h+2] <- 1
            }
          }
        }
      }
    }

    tp <- sum(which(alarm == 1) %in% which(labels == 1))
    fp <- sum(!which(alarm == 1) %in% which(labels == 1))
    tn <- sum(which(alarm == 0) %in% which(labels == 0))
    fn <- sum(!which(alarm == 0) %in% which(labels == 0))

    recall <- sensitivity <- tpr <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    specificity <- ifelse(tn + fp == 0, 0, tn / (tn + fp))
    precision <- ppv <- ifelse(tp + fp == 0, 0, tp / (tp + fp))

    AUC_mat[i,1] <- fpr <- 1-specificity
    AUC_mat[i,2] <- sensitivity
    AUCPR_mat[i,1] <- recall
    AUCPR_mat[i,2] <- precision
  }
  ord <- order(AUC_mat[,1])
  AUC <- sum(diff(AUC_mat[ord,1]) * zoo::rollmean(AUC_mat[ord,2], 2))

  ord <- order(AUCPR_mat[,1])
  AUCPR <- sum(diff(AUCPR_mat[ord,1]) * zoo::rollmean(AUCPR_mat[ord,2], 2))
  return(c(AUC = AUC, AUCPR = AUCPR))
}

get_alarm_performance <- function(outcome, id){
  if(grepl("Y5", outcome)){
    horizon <- 5
  } else if(grepl("10", outcome)){
    horizon <- 10
  } else if(grepl("15", outcome)){
    horizon <- 15
  } else if(grepl("20", outcome)){
    horizon <- 20
  } else if(grepl("25", outcome)){
    horizon <- 25
  } else if(grepl("30", outcome)){
    horizon <- 30
  }
  if(grepl("mean", outcome)){
    smooth_type <- "mean"
  } else if(grepl("median", outcome)){
    smooth_type <- "median"
  }
  file <- paste0(results_path, "id", id, "/tables/", outcome, "_long_tables.xlsx")
  dt <- data.table(read_excel(file, sheet = "forecasts"))
  not_forecasts <- c("time", "currentY", "futureY", "futureY_time",
                     "outcome_time","forecast_time_precise")
  learners <- colnames(dt[, -c(not_forecasts, "forecast_time"), with=FALSE])
  dt <- dt[, -not_forecasts, with = FALSE]
  dt[, forecast_time := factor(forecast_time)]

  # dt with binarized predictions, predicted AHE prevalence for each 5-min batch
  dt_score <- data.table(
    apply(dt[,-"forecast_time", with = F], 2, function(x) ifelse(x <= 65, 1, 0))
  )
  dt_score <- cbind("forecast_time" = factor(dt[["forecast_time"]]), dt_score)
  dt_score <- dt_score[, lapply(.SD, mean), by=forecast_time, .SDcols=learners]

  # mean MAP for each 5-min batch
  dt <- dt[, lapply(.SD, mean), by=forecast_time, .SDcols=learners]

  learners <- learners[-which(learners == "obs")]

  # labels (true alarm pattern) based on observed MAP average over 5-min batches
  obs <- dt[["obs"]]
  alarm0_idx <- which(obs <= 65)
  if (length(alarm0_idx) > 0) {
    alarm0_reset_idx <- which(obs > 65)
    alarm0_idx <- process_alarm(alarm0_idx, alarm0_reset_idx)

    if (length(alarm0_idx) > 0) {
      alarm0 <- rep(0, length(obs))
      alarm0[alarm0_idx] <- 1

      # mat <- matrix(nrow = 6, ncol = length(learners))
      mat <- matrix(nrow = 4, ncol = length(learners))
      for(j in 1:length(learners)){
        learner <- learners[j]
        pred <- dt[[learner]]
        alarm <- rep(0, length(pred))
        alarm_idx <- which(pred <= 65)
        alarm_reset_idx <- which(pred > 65)
        if(length(alarm_idx) >= 1){
          alarm_idx <- process_alarm(alarm_idx, alarm_reset_idx)
          if(length(alarm_idx) >= 1){
            alarm[alarm_idx] <- 1
            # label alarm as correct if it occurred no more than 10 min before label
            for(h in 1:(length(alarm)-1)){
              if(alarm[h] == 1 & alarm0[h] == 0){
                if(alarm0[h+1] == 1){
                  alarm[h] <- 0
                  alarm[h+1] <- 1
                }
                if(alarm0[h+2] == 1 & h < length(alarm)-1){
                  alarm[h] <- 0
                  alarm[h+2] <- 1
                }
              }
            }
          }
        }

        tp <- sum(which(alarm == 1) %in% which(alarm0 == 1))
        fp <- sum(!which(alarm == 1) %in% which(alarm0 == 1))
        tn <- sum(which(alarm == 0) %in% which(alarm0 == 0))
        fn <- sum(!which(alarm == 0) %in% which(alarm0 == 0))
        mat[1,j] <- recall <- sensitivity <- tpr <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
        specificity <- ifelse(tn + fp == 0, 0, tn / (tn + fp))
        mat[2,j] <-  precision <- ppv <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
        misclassification <- (fn + fp) / (tn + fp + fn + tp)
        mat[3,j] <- accuracy <- 1 - misclassification
        mat[4,j] <- f1 <- ifelse(tpr == 0 | ppv == 0, 0, 2 / ((1/ppv) + (1/tpr)))

        # AUCs <- eval_AUC(dt_score[[learner]], alarm0, seq(0,1,0.1))
        # mat[5,j] <- AUCs["AUC"]
        # mat[6,j] <- AUCs["AUCPR"]
      }
      colnames(mat) <- learners

      # return(cbind(
      #   "id" = rep(id, 6), "horizon" = rep(horizon, 6),
      #   "smooth_type" = rep(smooth_type, 6),
      #   "metric" = c("TPR", "PPV", "Accuracy", "F1score", "AUC", "AUCPR"),
      #   data.table(mat)
      # ))
      return(cbind(
        "id" = rep(id, 4), "horizon" = rep(horizon, 4),
        "smooth_type" = rep(smooth_type, 4),
        "metric" = c("TPR", "PPV", "Accuracy", "F1score"),
        data.table(mat)
      ))
    } else {
      cat("\n No true alarms for id", id, "with outcome", outcome, "\n")
      return(NULL)
    }
  } else {
    cat("\n No true alarms for id", id, "with outcome", outcome, "\n")
    return(NULL)
  }
}
# ==============================================================================

results_path <- "~/Google Drive/My Drive/symphony-results/"
# get ids
ids <- list.files(results_path)[grep("^id", list.files(results_path))]
ids <- gsub("id", "", ids)
horizons <- c("5", "10", "15", "20", "30")
result <- lapply(horizons, function(horizon){
  smooth_types <- c("mean", "median")
  outcomes <- paste0("Y", horizon, "_lag5_", smooth_types)
  result <- lapply(outcomes, function(outcome){
    result <- lapply(ids, function(id){
      get_alarm_performance(outcome, id)
    })
    result <- do.call(rbind, result[!sapply(result, is.null)])
    n <- length(unique(result$id))
    group_cols <- c("smooth_type", "metric", "horizon")
    learners <- colnames(result[,-c(group_cols, "id"), with = FALSE])
    result_summary <- result[, lapply(.SD, mean, na.rm = T), by=group_cols, .SDcols=learners]
    result_summary <- cbind("n" = rep(n, nrow(result_summary)), result_summary)
    return(list(result = result, result_summary = result_summary))
  })
  result_summary <- do.call(rbind, lapply(result, '[[', 'result_summary'))
  result <- do.call(rbind, lapply(result, '[[', 'result'))
  setorder(result, horizon, smooth_type, id)
  return(list(result = result, result_summary = result_summary))
})

result_summary <- do.call(rbind, lapply(result, '[[', 'result_summary'))
# metric_levs <- c("TPR", "PPV", "F1score", "Accuracy", "AUC", "AUCPR")
metric_levs <- c("TPR", "PPV", "F1score", "Accuracy")
result_summary[, metric := factor(metric, levels = metric_levs)]
setorder(result_summary, metric, horizon, smooth_type)
write.csv(result_summary, file = paste0(results_path, "alarm_allID.csv"),
          row.names = FALSE)

result <- lapply(result, '[[', 'result')
names(result) <- horizons
wb <- createWorkbook("Workbook")
lapply(seq_along(result), function(x){
  addWorksheet(wb, paste0("horizon", names(result)[x]))
  writeDataTable(wb, paste0("horizon", names(result)[x]), result[[x]])
})
saveWorkbook(wb, file = paste0(results_path, "alarm_byID.xlsx"), overwrite = T)
