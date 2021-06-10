library(openxlsx)
library(ROCR)
library(data.table)
library(readxl)

# ==============================================================================
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
  dt <- dt[, lapply(.SD, mean), by=forecast_time, .SDcols=learners]
  learners <- learners[-which(learners == "obs")]
  obs <- dt[["obs"]]

  # when should the alarming occur?
  true_alarm_idx <- which(obs <= 65)
  true_alarm_reset_idx <- which(obs > 65)
  if(length(true_alarm_idx) > 0){
    na_idx <- lapply(seq_along(true_alarm_idx), function(i){
      true_alarm_silence_idx <- true_alarm_idx[i] + 1:6
      if(any(true_alarm_idx %in% true_alarm_silence_idx)){
        if(any(true_alarm_silence_idx %in% true_alarm_reset_idx)){
          reset_idx <- min(which(true_alarm_silence_idx %in% true_alarm_reset_idx))
          if(reset_idx > 1){
            true_alarm_silence_idx <- true_alarm_silence_idx[c(1:(reset_idx-1))]
            return(which(true_alarm_idx %in% true_alarm_silence_idx))
          } else {
            return(NULL)
          }
        } else {
          return(which(true_alarm_idx %in% true_alarm_silence_idx))
        }
      } else {
        return(NULL)
      }
    })
    na_idx <- unique(unlist(na_idx[!sapply(na_idx, is.null)]))
    true_alarm_idx <- ifelse(length(na_idx) > 0, true_alarm_idx[-na_idx], true_alarm_idx)
    true_alarm_obs_binary <- rep(0, length(obs))

    if(length(true_alarm_idx) > 0){
      true_alarm_obs_binary[true_alarm_idx] <- 1

      mat <- matrix(nrow = 4, ncol = length(learners))
      for(j in 1:length(learners)){
        learner <- learners[j]
        pred <- dt[[learner]]
        pred_binary <- rep(0, length(pred))
        alarm_idx <- which(pred <= 65)
        alarm_reset_idx <- which(pred > 65)
        if(length(alarm_idx) >= 1){
          na_idx <- lapply(seq_along(alarm_idx), function(i){
            alarm_silence_idx <- alarm_idx[i] + 1:6
            if(any(alarm_idx %in% alarm_silence_idx)){
              if(any(alarm_silence_idx %in% alarm_reset_idx)){
                reset_idx <- min(which(alarm_silence_idx %in% alarm_reset_idx))
                if(reset_idx > 1){
                  alarm_silence_idx <- alarm_silence_idx[c(1:(reset_idx-1))]
                  return(which(alarm_idx %in% alarm_silence_idx))
                } else {
                  return(NULL)
                }
              } else {
                return(which(alarm_idx %in% alarm_silence_idx))
              }
            } else {
              return(NULL)
            }
          })
          na_idx <- unique(unlist(na_idx[!sapply(na_idx, is.null)]))
          alarm_idx <- ifelse(length(na_idx) > 0, alarm_idx[-na_idx], alarm_idx)
          if(length(alarm_idx) >= 1){
            pred_binary[alarm_idx] <- 1
            alarm_pred_binary <- rep(1, length(alarm_idx))
            alarm_obs_binary <- ifelse(obs[alarm_idx] <= 65, 1, 0)
            if(length(unique(alarm_obs_binary)) > 1){
              obj <- ROCR::prediction(alarm_pred_binary, alarm_obs_binary)
              mat[1,j] <- as.numeric(ROCR::performance(obj, "auc")@y.values)
              mat[2,j] <- as.numeric(ROCR::performance(obj, "aucpr")@y.values)
            } else {
              if(sum(alarm_obs_binary) == length(alarm_obs_binary)){
                mat[1,j] <- mat[2,j] <- 1
              } else {
                mat[1,j] <- mat[2,j] <- 0
              }
            }
          } else {
            mat[1,j] <- mat[2,j] <- 0
          }
        } else {
          mat[1,j] <- mat[2,j] <- 0
        }

        obj <- ROCR::prediction(pred_binary, true_alarm_obs_binary)
        mat[3,j] <- as.numeric(ROCR::performance(obj, "auc")@y.values)
        mat[4,j] <- as.numeric(ROCR::performance(obj, "aucpr")@y.values)
      }
      colnames(mat) <- learners

      return(cbind(
        "id" = rep(id, 4), "horizon" = rep(horizon, 4),
        "smooth_type" = rep(smooth_type, 4),
        "metric" = c("AUC_event", "AUCPR_event", "AUC_alarm", "AUCPR_alarm"),
        data.table(mat)
      ))
    } else {
      cat("\n No true alarms for id ", id, "outcome ", outcome, "\n")
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
    result_summary <- result[, lapply(.SD, mean), by=group_cols, .SDcols=learners]
    result_summary <- cbind("n" = rep(n, nrow(result_summary)), result_summary)
    return(list(result = result, result_summary = result_summary))
  })
  result_summary <- do.call(rbind, lapply(result, '[[', 'result_summary'))
  result <- do.call(rbind, lapply(result, '[[', 'result'))
  setorder(result, horizon, smooth_type, id)
  return(list(result = result, result_summary = result_summary))
})

result_summary <- do.call(rbind, lapply(result, '[[', 'result_summary'))
metric_levels <- c("AUCPR_event","AUC_event", "AUCPR_alarm", "AUC_alarm")
result_summary[, metric := factor(metric, levels = metric_levels)]
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
saveWorkbook(wb, file = paste0(results_path, "alarm_byID.xlsx"))
