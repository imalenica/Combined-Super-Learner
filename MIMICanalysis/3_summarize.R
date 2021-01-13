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
library(Hmisc)
library(gifski)

source(here::here("MIMICanalysis", "utils_plotting.R"))
source(here::here("MIMICanalysis", "utils_animation.R"))

################################################################################
summarize_results <- function(individual_data, ids, outcome,
                              summarize_full_fit = TRUE){
  # get forecast horizon from outcome name
  if(grepl("15", outcome)){
    horizon <- 15
  } else if(grepl("20", outcome)){
    horizon <- 20
  } else if(grepl("25", outcome)){
    horizon <- 25
  } else if(grepl("30", outcome)){
    horizon <- 30
  } else if(grepl("35", outcome)){
    horizon <- 35
  } else if(grepl("40", outcome)){
    horizon <- 40
  } else if(grepl("45", outcome)){
    horizon <- 45
  }

  if(grepl("mean", outcome)){
    smooth_type <- "mean"
    current_outcome <- "abpmean_lag5_mean"
  } else if(grepl("median", outcome)){
    smooth_type <- "median"
    current_outcome <- "abpmean_lag5_median"
  }

  for(i in 1:length(ids)){

    # subset data
    individual_id <- ids[i]
    cat("\nSummarizing results for id", individual_id, "with outcome", outcome, "\n")
    d <- individual_data[id == individual_id,]
    setorder(d, time_and_date)

    # retain true outcomes
    truth_tbl <- d[, c(current_outcome, outcome, "min"), with=F]
    colnames(truth_tbl) <- c("current_outcome","future_outcome","time")
    truth_tbl$future_outcome_time <- truth_tbl$time + horizon

    # make main folder and subfolders
    newfolder <- paste0("id", individual_id)
    dir <- paste0(save_path, newfolder)
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
    tables_folder <- file.path(dir, "tables")
    if(!dir.exists(tables_folder)){
      dir.create(tables_folder)
    }

    # patient chart
    chart_name <- paste0(smooth_type, "_smooth_chart.pdf")
    if(!file.exists(paste0(chart_folder, "/", chart_name))){
      chart <- make_patient_chart(individual_data = d, smooth_type = smooth_type)
      ggsave(filename = chart_name, plot = chart, device = "pdf",
             path = paste0(chart_folder, "/"), height = 8, width = 10)
    }

    # load results
    results_path <- paste0(outcome, "_", "id", individual_id,  ".Rdata")
    load(file = paste0(data_path, results_path))

    # movie
    forecast_tbl <- result$forecast_table
    movie_tbl <- suppressMessages(full_join(truth_tbl, forecast_tbl))
    osl_forecasts <- data.frame(time=movie_tbl$future_outcome_time,
                                forecast_time=movie_tbl$forecast_time_precise,
                                forecast=movie_tbl$onlineSL)
    truth <- data.frame(time=movie_tbl$future_outcome_time, truth=movie_tbl$future_outcome)
    if(summarize_full_fit){
      cols_list <- list()
      colnames(osl_forecasts)[3] <- "onlineSL"
      cols_list[[1]] <- list(forecast_time = colnames(osl_forecasts)[2],
                             forecast = colnames(osl_forecasts)[3])
      osl_forecasts$onlineSL_fullfit <- movie_tbl$onlineSL_full_fit
      cols_list[[2]] <- list(forecast_time = colnames(osl_forecasts)[2],
                             forecast = colnames(osl_forecasts)[4])
    }
    movie_name <- paste0("/", outcome, "_", "id", individual_id, ".mp4")
    movie_tbl <- suppressMessages(data.table(full_join(truth, osl_forecasts)))
    if(summarize_full_fit){
      make_multiforecast_movie_mp4(data.frame(movie_tbl), cols_list,
                                   name = paste0(movie_folder, movie_name))
    } else {
      make_forecast_movie(movie_tbl, name=paste0(movie_folder, movie_name))
    }

    # forecast plots
    forecast_tbl <- merge(truth_tbl, forecast_tbl, by="time", all.x=F, all.y=T)
    if(summarize_full_fit){
      learners <- c("onlineSL", "onlineSL_full_fit")
    } else {
      learners <- "onlineSL"
    }
    forecast_plot <- make_forecast_plot(forecast_tbl, learners = learners)
    plot_name <-  paste0(outcome, "_forecasts.pdf")
    ggsave(filename = plot_name, plot = forecast_plot, device = "pdf",
           path = paste0(plots_folder, "/"), height=8, width=10)

    # individual/historical learner proportion plots
    prop_tbl <- data.frame(Individual = result$sl_table$prop_individual,
                           Historical = result$sl_table$prop_historical,
                           time = result$sl_table$min_start)
    prop_plot <- make_prop_plot(prop_tbl)
    plot_name <- paste0(outcome, "_props.pdf")
    ggsave(filename = plot_name, plot = prop_plot, device = "pdf",
           path = paste0(plots_folder, "/"), height = 8, width = 10)

    # tables
    tbl_path <- paste0(tables_folder, "/", outcome, "_tables.xlsx")
    write.xlsx(result$sl_table, tbl_path, sheetName = "onlineSL", row.names=F, showNA=F)
    write.xlsx(forecast_tbl, tbl_path, sheetName="forecasts", append=T, row.names=F, showNA=F)
    MSE_tbl <- result$MSE_table
    MSE_tbl <- data.table(merge(truth_tbl, MSE_tbl, by="time", all.x=F, all.y=T))
    write.xlsx(MSE_tbl, tbl_path, sheetName = "MSE", row.names=F, showNA=F,append=T)
    write.xlsx(result$loss_table, tbl_path, sheetName="learner_losses", append=T,
               row.names=F, showNA=F)
  }
}

summarize_id <- function(individual_id, individual_mean_data,
                         individual_median_data, summarize_full_fit = TRUE){

  # make main folder and subfolders
  newfolder <- paste0("id", individual_id)
  dir <- paste0(save_path, newfolder)
  if(!dir.exists(dir)){
    dir.create(dir)
  }

  smooth_types <- c("mean", "median")

  if(summarize_full_fit){
    oSL_MSEsummary_smoothFF <- list()
  }
  oSL_MSEsummary_smooth <- list()
  meanMSE_smooth <- list()
  medianMSE_smooth <- list()
  NAsummary_smooth <- list()
  for (s in 1:length(smooth_types)){
    smooth_type <- smooth_types[s]
    horizons <- c(15, 20, 25, 30, 35, 40, 45)
    outcomes <- sapply(horizons, function(h) paste0("Y", h, "_lag5_", smooth_type))

    if(smooth_type == "mean"){
      current_outcome <- "abpmean_lag5_mean"
      d <- individual_mean_data[id == individual_id,]
    } else if(smooth_type == "median"){
      current_outcome <- "abpmean_lag5_median"
      d <- individual_median_data[id == individual_id,]
    }
    setorder(d, time_and_date)

    cols_list <- list()
    movie_tbls <- list()
    forecast_tbls <- list()
    meanMSE <- list()
    medianMSE <- list()
    NA_summary <- list()
    oSL_MSEsummary <- matrix(nrow=7, ncol=6)
    if(summarize_full_fit){
      cols_listFF <- list()
      movie_tblsFF <- list()
      forecast_tblsFF <- list()
      oSL_MSEsummaryFF <- matrix(nrow=7, ncol=6)
    }
    for(i in 1:length(outcomes)){
      outcome <- outcomes[i]
      results_path <- paste0(outcome, "_", "id", individual_id,  ".Rdata")
      load(file = paste0(data_path, results_path))

      # retain true outcomes
      truth_tbl <- d[, c(current_outcome, outcome, "min"), with=F]
      colnames(truth_tbl) <- c("current_outcome","future_outcome","time")
      truth_tbl$future_outcome_time <- truth_tbl$time + horizons[i]

      # movie
      forecast_tbl <- result$forecast_table
      movie_tbl <- suppressMessages(full_join(truth_tbl, forecast_tbl))
      osl_forecasts <- data.frame(time=movie_tbl$future_outcome_time,
                                  forecast_time=movie_tbl$forecast_time_precise,
                                  forecast=movie_tbl$onlineSL)
      colnames(osl_forecasts)[2] <- paste0("forecast_time_horizon", horizons[i])
      colnames(osl_forecasts)[3] <- paste0("horizon", horizons[i])
      cols_list[[i]] <- list(forecast_time = colnames(osl_forecasts)[2],
                             forecast = colnames(osl_forecasts)[3])
      truth <- data.frame(time=movie_tbl$future_outcome_time, truth=movie_tbl$future_outcome)
      movie_tbl <- suppressMessages(data.table(full_join(truth, osl_forecasts)))
      movie_tbls[[i]] <- movie_tbl

      # forecast plots
      forecast_tbl <- merge(truth_tbl, forecast_tbl, by="time", all.x=F, all.y=T)
      forecast_tbl <- data.table(future_outcome_time=forecast_tbl$future_outcome_time,
                                 future_outcome=forecast_tbl$future_outcome,
                                 forecast=forecast_tbl$onlineSL)
      colnames(forecast_tbl)[3] <- paste0("horizon", horizons[i])
      forecast_tbls[[i]] <- forecast_tbl

      # MSE summary
      time_idx <- grepl("time", colnames(result$MSE_table))
      MSE_tbl <- result$MSE_table[, -time_idx, with=F]
      meanMSE[[i]] <- data.table(t(round(colMeans(MSE_tbl, na.rm=T),4)))
      medianMSE[[i]] <- data.table(t(round(apply(MSE_tbl, 2, median, na.rm=T),4)))
      NA_summary[[i]] <- data.table(t(result$NA_summary))
      oSL_MSEsummary[i,] <- round(as.numeric(summary(MSE_tbl[["onlineSL"]], na.rm=T)),4)

      if(summarize_full_fit){
        # movie
        forecast_tbl <- result$forecast_table
        movie_tbl <- suppressMessages(full_join(truth_tbl, forecast_tbl))
        osl_forecasts <- data.frame(time=movie_tbl$future_outcome_time,
                                    forecast_time=movie_tbl$forecast_time_precise,
                                    forecast=movie_tbl$onlineSL_full_fit)
        colnames(osl_forecasts)[2] <- paste0("forecast_time_horizon", horizons[i])
        colnames(osl_forecasts)[3] <- paste0("horizon", horizons[i])
        cols_listFF[[i]] <- list(forecast_time = colnames(osl_forecasts)[2],
                                 forecast = colnames(osl_forecasts)[3])
        movie_tbl <- suppressMessages(data.table(full_join(truth, osl_forecasts)))
        movie_tblsFF[[i]] <- movie_tbl

        # forecast plots
        forecast_tbl <- merge(truth_tbl, forecast_tbl, by="time", all.x=F, all.y=T)
        forecast_tbl <- data.table(future_outcome_time=forecast_tbl$future_outcome_time,
                                   future_outcome=forecast_tbl$future_outcome,
                                   forecast=forecast_tbl$onlineSL_full_fit)
        colnames(forecast_tbl)[3] <- paste0("horizon", horizons[i])
        forecast_tblsFF[[i]] <- forecast_tbl

        # MSE summary
        oSL_MSEsummaryFF[i,] <- round(as.numeric(summary(MSE_tbl[["onlineSL_full_fit"]], na.rm=T)),4)
      }

    }
    base_dt <- data.table(id=rep(individual_id, 7), horizon=horizons,
                          smooth_type=rep(smooth_type, 7))
    oSL_MSEsummary <- data.table(base_dt, oSL_MSEsummary)
    colnames(oSL_MSEsummary)[4:9] <- c("Min","1stQ","Median","Mean","3rdQ","Max")
    oSL_MSEsummary_smooth[[s]] <- oSL_MSEsummary
    meanMSE_smooth[[s]] <- data.table(base_dt, rbindlist(meanMSE, fill=T))
    medianMSE_smooth[[s]] <- data.table(base_dt, rbindlist(medianMSE, fill=T))
    NAsummary_smooth[[s]] <- data.table(base_dt, rbindlist(NA_summary, fill=T))

    movie_tbl <- suppressMessages(data.frame(movie_tbls %>% reduce(full_join)))
    movie_name <- paste0("/", "movie_smooth_", smooth_type, "_id", individual_id, ".mp4")
    make_multiforecast_movie_mp4(movie_tbl, cols_list, name = paste0(dir, movie_name))

    forecast_tbl <- suppressMessages(data.table(forecast_tbls %>% reduce(full_join)))
    plot_name <-  paste0("forecasts_smooth_", smooth_type, "_id", individual_id, ".pdf")
    learners <- sapply(horizons, function(h) paste0("horizon", h))
    forecast_plot <- make_forecast_plot(forecast_tbl, learners = learners,
                                        truth_line_size = 0.3,
                                        forecast_line_size = 0.2)
    ggsave(filename = plot_name, plot = forecast_plot, device = "pdf",
           path = paste0(dir, "/"), height=8, width=10)

    if(summarize_full_fit){
      oSL_MSEsummary <- data.table(base_dt, oSL_MSEsummaryFF)
      colnames(oSL_MSEsummary)[4:9] <- c("Min","1stQ","Median","Mean","3rdQ","Max")
      oSL_MSEsummary_smoothFF[[s]] <- oSL_MSEsummary

      movie_tbl <- suppressMessages(data.frame(movie_tblsFF %>% reduce(full_join)))
      movie_name <- paste0("/", "movie_fullfit_smooth_", smooth_type, "_id", individual_id, ".mp4")
      make_multiforecast_movie_mp4(movie_tbl, cols_list, name = paste0(dir, movie_name))

      forecast_tbl <- suppressMessages(data.table(forecast_tblsFF %>% reduce(full_join)))
      plot_name <-  paste0("forecasts_fullfit_smooth_", smooth_type, "_id", individual_id, ".pdf")
      learners <- sapply(horizons, function(h) paste0("horizon", h))
      forecast_plot <- make_forecast_plot(forecast_tbl, learners = learners,
                                          truth_line_size = 0.3,
                                          forecast_line_size = 0.2)
      ggsave(filename = plot_name, plot = forecast_plot, device = "pdf",
             path = paste0(dir, "/"), height=8, width=10)
    }
  }
  summ <- rbindlist(oSL_MSEsummary_smooth, fill=T)
  setorder(summ, horizon)
  meanMSE <- rbindlist(meanMSE_smooth, fill=T)
  setorder(meanMSE, horizon)
  medianMSE <- rbindlist(medianMSE_smooth, fill=T)
  setorder(medianMSE, horizon)
  NAsummary <- rbindlist(NAsummary_smooth, fill=T)
  setorder(NAsummary, horizon)

  tbl_path <- paste0(dir, "/", "summarized_performance_id", individual_id, ".xlsx")
  write.xlsx(meanMSE, tbl_path, sheetName="mean_MSE", row.names=F, showNA=F)
  write.xlsx(medianMSE, tbl_path, sheetName="median_MSE", append=T, row.names=F, showNA=F)
  write.xlsx(summ, tbl_path, sheetName="onlineSL_MSE_summary", row.names=F, showNA=F, append=T)
  if(summarize_full_fit){
    summFF <- rbindlist(oSL_MSEsummary_smoothFF, fill=T)
    setorder(summFF, horizon)
    write.xlsx(summFF, tbl_path, sheetName="onlineSL_fullfit_MSE_summary",
               row.names=F, showNA=F, append=T)
  }
  write.xlsx(NAsummary, tbl_path, sheetName="NA_summary", row.names=F, showNA=F, append=T)
}
################################################################################
data_path <- "~/Box/symphony-data/"
save_path <- "~/Google Drive File Stream/My Drive/onlineSL results MIMIC/"

load(paste0(data_path, "individual_mean.Rdata"))
individual$id <- as.character(individual$id)
individual_mean <- individual

load(paste0(data_path, "individual_median.Rdata"))
individual$id <- as.character(individual$id)
individual_median <- individual

ids <- c("791","9275","3094", "26451")
for(i in 1:length(ids)){
  summarize_id(ids[i], individual_mean, individual_median)
}

horizons <- c(15, 20, 25, 30, 35, 40, 45)
lapply(horizons, function(horizon){
  summarize_results(individual_mean, ids = ids, outcome = paste0("Y", horizon, "_lag5_mean"))
  summarize_results(individual_median, ids = ids, outcome = paste0("Y", horizon, "_lag5_median"))
})
