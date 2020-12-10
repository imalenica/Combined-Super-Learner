library(here)
library(tidyverse)
library(data.table)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(xlsx)
library(animation)

source(here::here("MIMICanalysis", "utils_plotting.R"))
source(here::here("MIMICanalysis", "utils_animation.R"))

args <- list(outcome = "Y15_lag5_mean", ids = "9275", data_path = "~/Downloads/") 

# get forecast horizon from outcome name
if(grepl("15", args$outcome)){
  horizon <- 15
} else if(grepl("20", args$outcome)){
  horizon <- 20
} else if(grepl("25", args$outcome)){
  horizon <- 25
} else if(grepl("30", args$outcome)){
  horizon <- 30
} else if(grepl("35", args$outcome)){
  horizon <- 35
} else if(grepl("40", args$outcome)){
  horizon <- 40
} else if(grepl("45", args$outcome)){
  horizon <- 45
}

# load individual data
if(grepl("mean", args$outcome)){
  smooth_type <- "mean"
  current_outcome <- "abpmean_lag5_mean"
} else if(grepl("median", args$outcome)){
  smooth_type <- "median"
  current_outcome <- "abpmean_lag5_median"
}
load(paste0(args$data_path, "individual_", smooth_type, ".Rdata"))
individual$id <- as.character(individual$id)

################################## save results ################################

for(i in 1:length(ids)){
  # subset data
  individual_id <- ids[i]
  d <- individual[id == individual_id,]
  setorder(d, time_and_date)
  
  # retain true outcomes
  truth_tbl <- d[, c(current_outcome, args$outcome, "min"), with=F]
  colnames(truth_tbl) <- c("current_outcome","future_outcome","forecast_time")
  truth_tbl$future_outcome_time <- truth_tbl$forecast_time + horizon
  
  # make main folder and subfolders
  newfolder <- paste0("id", individual_id)
  dir <- file.path(dirname(args$data_path), newfolder)
  if(!dir.exists(dir)){
    dir.create(dir)
  }
  movie_folder <- file.path(dir, "movies")
  if(!dir.exists(movie_folder)){
    dir.create(movie_folder)
  }
  plots_folder <- file.path(dir, "plots")
  if(!dir.exists(plots_folder)){
    dir.create(plots_folder)
  }
  chart_folder <- file.path(dir, "patient_charts")
  if(!dir.exists(chart_folder)){
    dir.create(chart_folder)
  }
  
  # patient chart
  chart_name <- paste0(chart_folder, "/", smooth_type, "_smooth_chart.pdf")
  if(!file.exists(chart_name)){
    chart <- make_patient_chart(individual_data = d, smooth_type = smooth_type)
    pdf(file = chart_name)
    chart
    dev.off()
  }
  
  # load results
  results_path <- paste0(args$outcome, "_", "id", individual_id,  ".Rdata")
  load(file = paste0(args$data_path, results_path))
  
  # movie
  forecast_tbl <- result$forecast_table
  movie_tbl <- suppressMessages(full_join(truth_tbl, forecast_tbl))
  osl_forecasts <- data.frame(time=movie_tbl$future_outcome_time, 
                              real_time=movie_tbl$real_time_precise,
                              forecast=movie_tbl$onlineSL)
  truth <- data.frame(time=movie_tbl$forecast_time, truth=movie_tbl$current_outcome)
  movie_tbl <- suppressMessages(data.table(full_join(truth, osl_forecasts)))
  movie_name <- paste0("/", args$outcome, "_", "id", individual_id, ".mp4")
  make_forecast_movie(movie_tbl, name=paste0(movie_folder, movie_name))
  
  # forecast plots
  forecast_tbl <- merge(truth_tbl, forecast_tbl, by="forecast_time", all.x=F, all.y=T)
  forecast_plot <- make_forecast_plot(forecast_tbl, type = "overlapping")
  pdf(file=paste0(plots_folder, "/", args$outcome, "forecasts_overlapping.pdf"), 
      height=8, width=10)
  forecast_plot
  dev.off()
  
  # individual/historical learner proportion plots
  prop_tbl <- data.frame(Individual = result$sl_table$prop_individual,
                         Historical = result$sl_table$prop_historical, 
                         time = result$sl_table$min_start)
  pdf(file=paste0(plots_folder, "/", args$outcome, "_props.pdf"), height=8, width=10)
  make_forecast_plot(forecast_tbl, type = "overlapping")
  dev.off()
  
  # tables
  tbl_path <- paste0(dir, "/", args$outcome, "_tables.xlsx")
  write.xlsx(result$sl_table, tbl_path, sheetName = "onlineSL", row.names=F, showNA=F)
  write.xlsx(forecast_tbl, path, sheetName="forecasts", append=T, row.names=F, showNA=F)
  MSE_tbl <- result$MSE_table
  MSE_tbl <- data.table(merge(truth_tbl, MSE_tbl, by="forecast_time", all.x=F, all.y=T))
  write.xlsx(MSE_tbl, path, sheetName = "MSE", row.names=F, showNA=F,append=T)
  write.xlsx(result$loss_table, path, sheetName="learner_losses", append=T, 
             row.names=F, showNA=F)
}

make_visuals(args = list(outcome="Y15_lag5_mean", id="9275"))
make_visuals(args = list(outcome="Y15_lag5_mean", id="3094"))
make_visuals(args = list(outcome="Y15_lag5_mean", id="26451"))
make_visuals(args = list(outcome="Y15_lag5_mean", id="791"))
make_visuals(args = list(outcome="Y15_lag5_mean", id="20934"))
# onlineSL_fullfitF <- forecast_tbl_fullfitF[["onlineSL"]]
# onlineSL_fullfitT <- forecast_tbl_fullfitT[["onlineSL"]]
# ### comparator plot
# #opts <- c("onlineSL_fullfitTRUE", "onlineSL_fullfitFALSE") 
# learners <- data.table(onlineSL_fullfitF, onlineSL_fullfitT)
# #colnames(learners)  <- c("individual_glm", "individual_arima", "onlineSL", "historical_glm")
# 
# 
# df <- cbind.data.frame(time=forecast_tbl_fullfitF$future_outcome_time, 
#                        truth=forecast_tbl_fullfitF$future_outcome,learners)
# df <- melt(df, id.vars = "time")
# df$mysize <- rep(0.4, nrow(df))
# df$mysize[df$variable=="truth"] <- 1
# pdf(file=paste0(plots_folder, "/oSL_forecast_", args$outcome, ".pdf"), 
#     height=8, width=10)
# make_plot(df, title="Forecasts overlapped with True BP")
# dev.off()
# 
# 
# df <- cbind.data.frame(time=forecast_tbl$min, 
#                        truth=forecast_tbl$current_outcome,
#                        onlineSL=forecast_tbl$SL_adapt)
# df <- melt(df, id.vars = "time")
# df$mysize <- rep(0.4, nrow(df))
# df$mysize[df$variable=="truth"] <- 1
# pdf(file=paste0(plots_folder, "/doSL_realtime_", args$outcome, ".pdf"), 
#     height=8, width=10)
# make_plot(df)
# dev.off()
# 
# df <- cbind.data.frame(time=forecast_tbl$future_outcome_time, 
#                        truth=forecast_tbl$future_outcome, 
#                        onlineSL=forecast_tbl$SL_adapt)
# df <- melt(df, id.vars = "time")
# df$mysize <- rep(0.4, nrow(df))
# df$mysize[df$variable=="truth"] <- 1
# pdf(file=paste0(plots_folder, "/doSL_forecast_", args$outcome, ".pdf"), 
#     height=8, width=10)
# make_plot(df, title="Forecasts overlapped with True BP")
# dev.off()
# 
# 
# SLtbl_folder <- file.path(dir, "SL_tables")
# if(!dir.exists(SLtbl_folder)){
#   dir.create(SLtbl_folder)
# }
# SL_tbl_path <- paste0(SLtbl_folder, "/", args$outcome, ".csv")
# write.csv(fit$SL_table, file=SL_tbl_path, row.names=F)
# 
# loss_tbl_folder <- file.path(dir, "loss_tables")
# if(!dir.exists(loss_tbl_folder)){
#   dir.create(loss_tbl_folder)
# }
# loss_tbl_path <- paste0(loss_tbl_folder, "/", args$outcome, ".csv")
# write.csv(fit$loss_table, file=loss_tbl_path, row.names=F)
# 
# 
# process_tables <- function(results, type = c("onlineSL_table", "forecasts")){
#   
#   window10 <- c("onlineSL_adapt_window10", "onlineSL_historical_window10", 
#                 "onlineSL_individual_window10")
#   window5 <- c("onlineSL_adapt_window5", "onlineSL_historical_window5", 
#                "onlineSL_individual_window5")
#   window1 <- c("onlineSL_adapt_window1", "onlineSL_historical_window1", 
#                "onlineSL_individual_window1")
#   window0 <- c("onlineSL_adapt", "onlineSL_historical", "onlineSL_individual")
#   
#   if(type == "onlineSL_table"){
#     window10 <- c(window10[1], paste0("MSE_", window10[1]),
#                   window10[2], paste0("MSE_", window10[2]),
#                   window10[3], paste0("MSE_", window10[3]))
#     window5 <- c(window5[1], paste0("MSE_", window5[1]),
#                   window5[2], paste0("MSE_", window5[2]),
#                   window5[3], paste0("MSE_", window5[3]))
#     window1 <- c(window1[1], paste0("MSE_", window1[1]),
#                  window1[2], paste0("MSE_", window1[2]),
#                  window1[3], paste0("MSE_", window1[3]))
#     window0 <- c(window0[1], paste0("MSE_", window0[1]),
#                  window0[2], paste0("MSE_", window0[2]),
#                  window0[3], paste0("MSE_", window0[3]))
#     tbl <- data.table(results$onlineSL_table)
#     m <- 1
#   } else if(type == "forecasts"){
#     locf <- ifelse(d[[args$outcome_locf]] == 1 | d[["row_locf"]] == 1, 1, 0)
#     locf <- locf[args$outcome_gap:nrow(d)]
#     locf <- c(locf, rep(NA, nrow(d)-length(locf)))
#     Y_tbl <- data.table(truth = d[[args$outcome]], LOCF = locf, 
#                         forecast_time = d[["min"]])
#     Y_tbl <- Y_tbl[forecast_time <= max(results$forecast_table$forecast_time),]
#     tbl <- data.table(merge(Y_tbl, results$forecast_table, all.x=F, all.y=T, 
#                             by="forecast_time"))
#     m <- 5
#   }
#   tbl_window0 <- tbl[, window0, with=F]
#   tbl_window0_1 <- tbl_window0[c(1:(1*m)), ]
#   colnames(tbl_window0_1) <- window1
#   tbl_window1 <- rbind(tbl_window0_1, tbl[-c(1:(1*m)), window1, with=F])
#   tbl_window0_5 <- tbl_window0[c(1:(5*m)), ]
#   colnames(tbl_window0_5) <- window5
#   tbl_window5 <- rbind(tbl_window0_5, tbl[-c(1:(5*m)), window5, with=F])
#   tbl_window0_10 <- tbl_window0[c(1:(10*m)), ]
#   colnames(tbl_window0_10) <- window10
#   tbl_window10 <- rbind(tbl_window0_10, tbl[-c(1:(10*m)), window10, with=F])
#   oSL <- data.table(cbind(tbl_window0, tbl_window1, tbl_window5, tbl_window10))
#   
#   if(type == "onlineSL_table"){
#     return(oSL)
#   } else if(type == "forecasts"){
#     tbl_other <- tbl[, -c(window0,window1,window5,window10), with=F]
#     return(cbind(oSL, tbl_other))
#   }
# }
# 
# onlineSL_table <- process_tables(results, "onlineSL_table")
# perf <- onlineSL_table[,grep("MSE", colnames(onlineSL_table)),with=F]
# perf <- sapply(perf, as.numeric)
# summary_mean <- colMeans(perf, na.rm=T)
# summary_median  <- apply(perf, 2, median, na.rm=T)
# n <- apply(perf, 2, function(x) sum(is.na(x)))
# summary_sum_standardized <- colSums(perf, na.rm=T)/n
# onlineCVrisk_summary <- rbind(summary_mean, summary_median, summary_sum_standardized)
# Summary_Type <- c("Mean", "Median", "Standardized Sum")
# onlineCVrisk_summary <- data.table(Summary_Type, onlineCVrisk_summary)
# 
# forecast_tbl <- process_tables(results, "forecasts")
# 

# 
# truth <- forecast_tbl[,grep("truth", colnames(forecast_tbl)), with=F]
# time <- forecast_tbl[,grep("forecast_time", colnames(forecast_tbl)), with=F]
# make_df <- function(pattern){
#   oSL <- forecast_tbl[,grep(pattern,colnames(forecast_tbl)),with=F]
#   df <- cbind.data.frame(time, truth, oSL)
#   df <- melt(df, id.vars = "forecast_time")
#   df$mysize <- rep(0.4, nrow(df))
#   df$mysize[df$variable=="truth"] <- 1
#   return(df)
# }
# make_plot <- function(df, title="Online SL Forecasts alongside True Mean BP"){
#   ggplot(data = df, aes(x=forecast_time, y=value, size=mysize, color=variable)) + 
#     geom_line() +
#     ylim(45,95) +
#     scale_size(range = c(0.5, 1), guide="none") +
#     labs(x="Time (min)", y="Mean BP", title=title, color="Type") +
#     scale_colour_brewer(palette="Dark2") 
# }
# 
# df <- make_df("onlineSL_adapt")
# pdf(file=paste0(dir, "/", "gap", args$outcome_gap, "min", "_onlineSLresults_adapt.pdf"), 
#     height=8, width=10)
# make_plot(df)
# dev.off()
# 
# df <- make_df("onlineSL_historical")
# pdf(file=paste0(dir, "/", "gap", args$outcome_gap, "min", "_onlineSLresults_historical.pdf"), 
#     height=8, width=10)
# make_plot(df)
# dev.off()
# 
# df <- make_df("onlineSL_individual")
# pdf(file=paste0(dir, "/", "gap", args$outcome_gap, "min", "_onlineSLresults_individual.pdf"), 
#     height=8, width=10)
# make_plot(df)
# dev.off()
# 
# truth <- forecast_tbl[["truth"]]
# time <- forecast_tbl[["forecast_time"]]
# not_lrnr <- c("forecast_time", "truth", "LOCF")
# lrnr_tbl <- forecast_tbl[,-not_lrnr, with=F]
# MSE_tbl <- apply(lrnr_tbl, 2, function(x) ((x-truth)^2))
# MSEs <- colSums(MSE_tbl, na.rm = T)
