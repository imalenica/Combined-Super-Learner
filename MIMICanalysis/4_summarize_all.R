library(XLConnect)
library(writexl)
library(data.table)

summarize_all_ids <- function(results_path, continuousYonly = TRUE){
  # get ids
  ids <- list.files(results_path)[grep("^id", list.files(results_path))]
  ids <- gsub("id", "", ids)

  # result-type specific tables
  if(!continuousYonly){
    AUC <- list()
    AUCPR <- list()
    mean_performance_binaryY <- list()
    median_performance_binaryY <- list()
  }
  mean_performance_continuousY <- list()
  median_performance_continuousY <- list()
  for(i in 1:length(ids)){
    id <- ids[i]
    path <- paste0(results_path, "id", ids[i], "/summarized_performance_id", ids[i], ".xlsx")
    wb <- loadWorkbook(path)
    wb <- readWorksheet(wb, sheet = getSheets(wb))
    if(!continuousYonly){
      AUC[[i]] <- wb$AUC
      AUCPR[[i]] <- wb$AUCPR
      mean_performance_binaryY[[i]] <- wb$mean_performance_binaryY
      median_performance_binaryY[[i]] <- wb$median_performance_binaryY
    }
    mean_performance_continuousY[[i]] <- wb$mean_performance_continuousY
    median_performance_continuousY[[i]] <- wb$median_performance_continuousY
  }
  summarize_list <- function(list_object, ignore = "id"){
    dt <- data.table(do.call(rbind, list_object))
    num_ids <- length(unique(dt[["id"]]))
    groups <- c("horizon", "smooth_type", "outcome_type")
    cols <- colnames(dt)[-which(colnames(dt) %in% c(groups, ignore))]
    dt <- dt[, lapply(.SD, function(x) mean(x, na.rm=TRUE)), by=groups, .SDcols=cols]
    return(cbind(num_ids = rep(num_ids, nrow(dt)), dt))
  }
  mean_performance_continuousY <- summarize_list(mean_performance_continuousY)
  median_performance_continuousY <- summarize_list(median_performance_continuousY)
  sheets <- list(
    "mean_performance_continuousY" = mean_performance_continuousY,
    "median_performance_continuousY" = median_performance_continuousY
  )
  if(!continuousYonly){
    AUC <- summarize_list(AUC, ignore = c("Metric", "id"))
    AUCPR <- summarize_list(AUCPR, ignore = c("Metric", "id"))
    mean_performance_binaryY <- summarize_list(mean_performance_binaryY)
    median_performance_binaryY <- summarize_list(median_performance_binaryY)
    sheets <- c(
      sheets, list("AUC" = AUC, "AUCPR" = AUCPR,
                   "mean_performance_binaryY" = mean_performance_binaryY,
                   "median_performance_binaryY" = median_performance_binaryY)
    )
  }
  write_xlsx(sheets, paste0(results_path, "summary_all.xlsx"))
  cat("Summarized ids: ", stringr::str_c(ids))
}


summarize_all_ids("~/Google Drive/My Drive/symphony-results/")

