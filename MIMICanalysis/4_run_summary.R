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

source(here::here("R", "v3", "utils_plotting.R"))

args <- list(abpmean_name = "abpmean", 
             individual_data = "individual30.Rdata", 
             individual_id = "791", 
             outcome = "Y15",
             outcome_locf = "abpmean_locf", 
             outcome_gap = 15)
################################## load data ###################################
data_path <- "~/Downloads/"

newfolder <- paste0("id", args$individual_id)
dir <- file.path(dirname(data_path), newfolder)
if(!dir.exists(dir)){
  dir.create(dir)
}

# subset to subject of interest
load(paste0(data_path, args$individual_data))
individual$id <- as.character(individual$id)
d <- individual[id == args$individual_id,]
rm(individual)
setorder(d, time_and_date)

pdf(file = paste0(dir, "/", "patient_chart.pdf"))
visualize(individual_data=d, smoothed_data=F)
dev.off()

results_path <- paste0(args$outcome, "_id", args$individual_id, ".Rdata")
load(file=paste0(data_path, results_path))


process_tables <- function(results, type = c("onlineSL_table", "forecasts")){
  
  window10 <- c("onlineSL_adapt_window10", "onlineSL_historical_window10", 
                "onlineSL_individual_window10")
  window5 <- c("onlineSL_adapt_window5", "onlineSL_historical_window5", 
               "onlineSL_individual_window5")
  window1 <- c("onlineSL_adapt_window1", "onlineSL_historical_window1", 
               "onlineSL_individual_window1")
  window0 <- c("onlineSL_adapt", "onlineSL_historical", "onlineSL_individual")
  
  if(type == "onlineSL_table"){
    window10 <- c(window10[1], paste0("MSE_", window10[1]),
                  window10[2], paste0("MSE_", window10[2]),
                  window10[3], paste0("MSE_", window10[3]))
    window5 <- c(window5[1], paste0("MSE_", window5[1]),
                  window5[2], paste0("MSE_", window5[2]),
                  window5[3], paste0("MSE_", window5[3]))
    window1 <- c(window1[1], paste0("MSE_", window1[1]),
                 window1[2], paste0("MSE_", window1[2]),
                 window1[3], paste0("MSE_", window1[3]))
    window0 <- c(window0[1], paste0("MSE_", window0[1]),
                 window0[2], paste0("MSE_", window0[2]),
                 window0[3], paste0("MSE_", window0[3]))
    tbl <- data.table(results$onlineSL_table)
    m <- 1
  } else if(type == "forecasts"){
    locf <- ifelse(d[[args$outcome_locf]] == 1 | d[["row_locf"]] == 1, 1, 0)
    locf <- locf[args$outcome_gap:nrow(d)]
    locf <- c(locf, rep(NA, nrow(d)-length(locf)))
    Y_tbl <- data.table(truth = d[[args$outcome]], LOCF = locf, 
                        forecast_time = d[["min"]])
    Y_tbl <- Y_tbl[forecast_time <= max(results$forecast_table$forecast_time),]
    tbl <- data.table(merge(Y_tbl, results$forecast_table, all.x=F, all.y=T, 
                            by="forecast_time"))
    m <- 5
  }
  tbl_window0 <- tbl[, window0, with=F]
  tbl_window0_1 <- tbl_window0[c(1:(1*m)), ]
  colnames(tbl_window0_1) <- window1
  tbl_window1 <- rbind(tbl_window0_1, tbl[-c(1:(1*m)), window1, with=F])
  tbl_window0_5 <- tbl_window0[c(1:(5*m)), ]
  colnames(tbl_window0_5) <- window5
  tbl_window5 <- rbind(tbl_window0_5, tbl[-c(1:(5*m)), window5, with=F])
  tbl_window0_10 <- tbl_window0[c(1:(10*m)), ]
  colnames(tbl_window0_10) <- window10
  tbl_window10 <- rbind(tbl_window0_10, tbl[-c(1:(10*m)), window10, with=F])
  oSL <- data.table(cbind(tbl_window0, tbl_window1, tbl_window5, tbl_window10))
  
  if(type == "onlineSL_table"){
    return(oSL)
  } else if(type == "forecasts"){
    tbl_other <- tbl[, -c(window0,window1,window5,window10), with=F]
    return(cbind(oSL, tbl_other))
  }
}

onlineSL_table <- process_tables(results, "onlineSL_table")
perf <- onlineSL_table[,grep("MSE", colnames(onlineSL_table)),with=F]
perf <- sapply(perf, as.numeric)
summary_mean <- colMeans(perf, na.rm=T)
summary_median  <- apply(perf, 2, median, na.rm=T)
n <- apply(perf, 2, function(x) sum(is.na(x)))
summary_sum_standardized <- colSums(perf, na.rm=T)/n
onlineCVrisk_summary <- rbind(summary_mean, summary_median, summary_sum_standardized)
Summary_Type <- c("Mean", "Median", "Standardized Sum")
onlineCVrisk_summary <- data.table(Summary_Type, onlineCVrisk_summary)

forecast_tbl <- process_tables(results, "forecasts")

path <- paste0(dir, "/", "gap", args$outcome_gap, "min", "_results_tables.xlsx")
write.xlsx(onlineCVrisk_summary, path, sheetName = "summary_onlineSL_performance", 
           row.names=F, showNA=F)
write.xlsx(onlineSL_table, path, sheetName = "onlineSL_performance", 
           row.names=F, showNA=F,append=T)
write.xlsx(forecast_tbl, path, sheetName="all_forecasts", append=T, 
           row.names=F, showNA=F)
write.xlsx(results$loss_table, path_tables, sheetName="onlineSL_loss", append=T, 
           row.names=F, showNA=F)

truth <- forecast_tbl[,grep("truth", colnames(forecast_tbl)), with=F]
time <- forecast_tbl[,grep("forecast_time", colnames(forecast_tbl)), with=F]
make_df <- function(pattern){
  oSL <- forecast_tbl[,grep(pattern,colnames(forecast_tbl)),with=F]
  df <- cbind.data.frame(time, truth, oSL)
  df <- melt(df, id.vars = "forecast_time")
  df$mysize <- rep(0.4, nrow(df))
  df$mysize[df$variable=="truth"] <- 1
  return(df)
}
make_plot <- function(df, title="Online SL Forecasts alongside True Mean BP"){
  ggplot(data = df, aes(x=forecast_time, y=value, size=mysize, color=variable)) + 
    geom_line() +
    ylim(45,95) +
    scale_size(range = c(0.5, 1), guide="none") +
    labs(x="Time (min)", y="Mean BP", title=title, color="Type") +
    scale_colour_brewer(palette="Dark2") 
}

df <- make_df("onlineSL_adapt")
pdf(file=paste0(dir, "/", "gap", args$outcome_gap, "min", "_onlineSLresults_adapt.pdf"), 
    height=8, width=10)
make_plot(df)
dev.off()

df <- make_df("onlineSL_historical")
pdf(file=paste0(dir, "/", "gap", args$outcome_gap, "min", "_onlineSLresults_historical.pdf"), 
    height=8, width=10)
make_plot(df)
dev.off()

df <- make_df("onlineSL_individual")
pdf(file=paste0(dir, "/", "gap", args$outcome_gap, "min", "_onlineSLresults_individual.pdf"), 
    height=8, width=10)
make_plot(df)
dev.off()
