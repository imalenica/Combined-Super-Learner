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
library(RColorBrewer)

source(here::here("R", "utils_plotting.R"))
source(here::here("R", "utils_animation.R"))

summarize_outcome <- function(individual_data, ids, outcome, data_path,
                              save_path, summarize_full_fit = TRUE){

  # get forecast horizon from outcome name
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

  if(any(grepl("abpmean_lag5_mean", colnames(individual_data)))){
    smooth_type <- "mean"
  } else if(any(grepl("abpmean_lag5_median", colnames(individual_data)))){
    smooth_type <- "median"
  }

  if(grepl("AHE", outcome)){
    current_outcome <- "AHE"
    outcome_clear <- paste0(outcome, "_", smooth_type)
    outcome_type <- "binomial"
  } else {
    outcome_clear <- outcome
    if(grepl("mean", outcome)){
      current_outcome <- "abpmean_lag5_mean"
    } else if(grepl("median", outcome)){
      current_outcome <- "abpmean_lag5_median"
    }
    outcome_type <- "continuous"
  }


  for(i in 1:length(ids)){

    # subset data
    individual_id <- ids[i]
    cat("\nSummarizing results for id", individual_id, "with outcome",
        outcome_clear, "\n")
    d <- individual_data[id == individual_id,]
    setorder(d, time_and_date)

    # retain true outcomes
    truth_tbl <- d[, c(current_outcome,"min", outcome), with=F]
    colnames(truth_tbl) <- c("currentY", "time", "futureY")
    truth_tbl$futureY_time <- truth_tbl$time + horizon

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
    results_path <- paste0(outcome_clear, "_id", individual_id, ".Rdata")
    load(file = paste0(data_path, results_path))

    # movie
    forecast_tbl <- result$forecast_table
    movie_tbl <- suppressMessages(full_join(truth_tbl, forecast_tbl))
    osl_forecasts <- data.frame(time=movie_tbl$futureY_time,
                                forecast_time=movie_tbl$forecast_time,
                                forecast=movie_tbl$onlineSL)
    truth <- data.frame(time=movie_tbl$futureY_time, truth=movie_tbl$futureY)
    if(summarize_full_fit){
      cols_list <- list()
      colnames(osl_forecasts)[3] <- "onlineSL"
      cols_list[[1]] <- list(forecast_time = "forecast_time",
                             forecast = "onlineSL")
      osl_forecasts$onlineSL_full_fit <- movie_tbl$onlineSL_full_fit
      cols_list[[2]] <- list(forecast_time = "forecast_time",
                             forecast = "onlineSL_full_fit")
    }
    movie_name <- paste0("/", outcome_clear, "_", "id", individual_id, ".mp4")
    movie_tbl <- suppressMessages(data.table(full_join(truth, osl_forecasts)))
    if(summarize_full_fit){
      make_multiforecast_movie_mp4(data.frame(movie_tbl), cols_list,
                                   name = paste0(movie_folder, movie_name),
                                   outcome_type = outcome_type)
    } else {
      make_forecast_movie(movie_tbl, name=paste0(movie_folder, movie_name),
                          outcome_type = outcome_type)
    }

    # forecast plots
    forecast_tbl <- merge(truth_tbl, forecast_tbl, by="time", all.x=F, all.y=T)
    if(summarize_full_fit){
      learners <- c("onlineSL", "onlineSL_full_fit")
    } else {
      learners <- "onlineSL"
    }
    forecast_plot <- make_forecast_plot(data.table(forecast_tbl),
                                        learners = learners,
                                        outcome_type = outcome_type)
    plot_name <-  paste0(outcome_clear, "_forecasts.pdf")
    ggsave(filename = plot_name, plot = forecast_plot, device = "pdf",
           path = paste0(plots_folder, "/"), height=8, width=10)

    # individual/historical learner proportion plots
    prop_tbl <- data.frame(Individual = result$sl_table$prop_individual,
                           Historical = result$sl_table$prop_historical,
                           time = result$sl_table$time)
    prop_plot <- make_prop_plot(prop_tbl)
    plot_name <- paste0(outcome_clear, "_props.pdf")
    ggsave(filename = plot_name, plot = prop_plot, device = "pdf",
           path = paste0(plots_folder, "/"), height = 8, width = 10)

    # tables
    tbl_path <- paste0(tables_folder, "/", outcome_clear, "_long_tables.xlsx")
    write.xlsx(result$sl_table, tbl_path, sheetName = "onlineSL", row.names=F,
               showNA=F)
    write.xlsx(forecast_tbl, tbl_path, sheetName="forecasts", append=T,
               row.names=F, showNA=F)

    MSE_tbl <- result$performance_table
    if(identical(MSE_tbl[["onlineSL"]], MSE_tbl[["onlineSL_full_fit"]])){
      print(paste0("Fixing losses for ", outcome_clear))
      forecast_tbl <- data.table(forecast_tbl)
      extra_cols <- c("time", "currentY", "futureY","futureY_time", "obs",
                      "outcome_time", "forecast_time", "forecast_time_precise")
      forecasts <- forecast_tbl[, -extra_cols, with=F]
      obs <- forecast_tbl[["obs"]]
      MSE_tbl <- data.table(apply(forecasts, 2, function(pred) (pred-obs)^2))
      MSE_tbl <- cbind(time = forecast_tbl[["time"]], MSE_tbl)
    }
    MSE_tbl <- data.table(merge(truth_tbl, MSE_tbl, by="time", all.x=F, all.y=T))

    write.xlsx(MSE_tbl, tbl_path, sheetName = "performance", row.names=F,
               showNA=F, append=T)
    write.xlsx(result$loss_table, tbl_path, sheetName="learner_losses",
               append=T, row.names=F, showNA=F)
    write.xlsx(result$risk_table, tbl_path, sheetName="learner_cv_risks",
               append=T, row.names=F, showNA=F)
    write.xlsx(result$timer_matrix, tbl_path, sheetName="timers",
               append=T, row.names=F, showNA=F)
    write.xlsx(result$missing_summary, tbl_path, sheetName="missing_summary",
               append=T, row.names=F, showNA=F)

    tbl_path <- paste0(tables_folder, "/", outcome_clear, "_performance_summary.csv")
    write.csv(result$performance_summary, tbl_path, row.names=F)

    if(current_outcome == "AHE"){
      tbl_path <- paste0(tables_folder, "/", outcome_clear, "_AUC.csv")
      write.csv(result$AUC_tbl, tbl_path, row.names=F)
    }
  }
}

summarize_id <- function(individual_id, individual_mean_data,
                         individual_median_data, horizons,
                         data_path, save_path, summarize_full_fit = TRUE){

  cat("\n************* Summarizing results for id", individual_id, "*************\n")

  # make main folder and subfolders
  newfolder <- paste0("id", individual_id)
  dir <- paste0(save_path, newfolder)
  if(!dir.exists(dir)){
    dir.create(dir)
  }

  outcome_types <- c("binomial", "continuous")
  mean_performance_Ytype <- list()
  median_performance_Ytype <- list()
  NAsummary_Ytype <- list()
  for (y in 1:length(outcome_types)){
    outcome_type <- outcome_types[y]
    smooth_types <- c("mean", "median")

    if(outcome_type == "binomial"){
      AUC_smooth <- list()
      AHEprev_smooth <- list()
    }
    mean_performance_smooth <- list()
    median_performance_smooth <- list()
    NAsummary_smooth <- list()
    for (s in 1:length(smooth_types)){
      smooth_type <- smooth_types[s]
      if(outcome_type == "binomial"){
        outcomes <- sapply(horizons, function(h) paste0("Y", h, "_AHE"))
      } else if (outcome_type == "continuous"){
        outcomes <- sapply(horizons, function(h) paste0("Y", h, "_lag5_", smooth_type))
      }

      if(smooth_type == "mean"){
        d <- individual_mean_data[id == individual_id,]
      } else if(smooth_type == "median"){
        d <- individual_median_data[id == individual_id,]
      }
      setorder(d, time_and_date)

      cols_list <- list()
      movie_tbls <- list()
      forecast_tbls <- list()
      mean_performance <- list()
      median_performance <- list()
      NA_summary <- list()
      if(summarize_full_fit){
        cols_listFF <- list()
        movie_tblsFF <- list()
        forecast_tblsFF <- list()
      }
      if(outcome_type == "binomial"){
        AUCsummary <- list()
        AHEprev <- list()
      }
      for(i in 1:length(outcomes)){
        outcome <- outcomes[i]
        if(grepl("AHE", outcome)){
          current_outcome <- "AHE"
          outcome_clear <- paste0(outcome, "_", smooth_type)
        } else {
          outcome_clear <- outcome
          if(grepl("mean", outcome)){
            current_outcome <- "abpmean_lag5_mean"
          } else if(grepl("median", outcome)){
            current_outcome <- "abpmean_lag5_median"
          }
        }
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

        results_path <- paste0(outcome_clear, "_id", individual_id, ".Rdata")
        load(file = paste0(data_path, results_path))

        # retain true outcomes
        truth_tbl <- d[, c(current_outcome, "min", outcome), with=F]
        colnames(truth_tbl) <- c("currentY","time", "futureY")
        truth_tbl$futureY_time <- truth_tbl$time + horizon

        # movie
        forecast_tbl <- result$forecast_table
        movie_tbl <- suppressMessages(full_join(truth_tbl, forecast_tbl))
        osl_forecasts <- data.frame(time=movie_tbl$futureY_time,
                                    forecast_time=movie_tbl$forecast_time,
                                    forecast=movie_tbl$onlineSL)
        colnames(osl_forecasts)[2] <- paste0("forecast_time_horizon", horizon)
        colnames(osl_forecasts)[3] <- paste0("horizon", horizon)
        cols_list[[i]] <- list(forecast_time = colnames(osl_forecasts)[2],
                               forecast = colnames(osl_forecasts)[3])
        truth <- data.frame(time=movie_tbl$futureY_time, truth=movie_tbl$futureY)
        movie_tbl <- suppressMessages(data.table(full_join(truth, osl_forecasts)))
        movie_tbls[[i]] <- movie_tbl

        # forecast plots
        forecast_tbl <- merge(truth_tbl, forecast_tbl, by="time", all.x=F, all.y=T)
        forecast_tbl <- data.table(futureY_time=forecast_tbl$futureY_time,
                                   futureY=forecast_tbl$futureY,
                                   forecast=forecast_tbl$onlineSL)
        colnames(forecast_tbl)[3] <- paste0("horizon", horizon)
        forecast_tbls[[i]] <- forecast_tbl

        # MSE summary
        perf_tbl <- result$performance_table
        if(identical(perf_tbl[["onlineSL"]], perf_tbl[["onlineSL_full_fit"]])){
          print(paste0("Fixing losses for outcome ", outcome_clear))
          forecast_tbl2 <- data.table(result$forecast_table)
          extra_cols <- c("time", "obs", "outcome_time", "forecast_time", "forecast_time_precise")
          forecasts <- forecast_tbl2[, -extra_cols, with=F]
          obs <- forecast_tbl2[["obs"]]
          perf_tbl <- data.table(apply(forecasts, 2, function(pred) (pred-obs)^2))
        } else {
          time_idx <- grepl("time", colnames(perf_tbl))
          perf_tbl <- perf_tbl[, -time_idx, with=F]
        }
        mean_performance[[i]] <- data.table(t(colMeans(perf_tbl, na.rm=T)))
        median_performance[[i]] <- data.table(t(apply(perf_tbl, 2, median, na.rm=T)))
        NA_summary[[i]] <- data.table(t(result$missing_summary))

        if(outcome_type == "binomial"){
          AUCsummary[[i]] <- result$AUC_tbl
          AHEprev[[i]] <- result$AHE_prevalence
        }

        if(summarize_full_fit){
          # movie
          forecast_tbl <- result$forecast_table
          movie_tbl <- suppressMessages(full_join(truth_tbl, forecast_tbl))
          osl_forecasts <- data.frame(time=movie_tbl$futureY_time,
                                      forecast_time=movie_tbl$forecast_time,
                                      forecast=movie_tbl$onlineSL_full_fit)
          colnames(osl_forecasts)[2] <- paste0("forecast_time_horizon", horizon)
          colnames(osl_forecasts)[3] <- paste0("horizon", horizon)
          cols_listFF[[i]] <- list(forecast_time = colnames(osl_forecasts)[2],
                                   forecast = colnames(osl_forecasts)[3])
          movie_tbl <- suppressMessages(data.table(full_join(truth, osl_forecasts)))
          movie_tblsFF[[i]] <- movie_tbl

          # forecast plots
          forecast_tbl <- merge(truth_tbl, forecast_tbl, by="time", all.x=F, all.y=T)
          forecast_tbl <- data.table(futureY_time=forecast_tbl$futureY_time,
                                     futureY=forecast_tbl$futureY,
                                     forecast=forecast_tbl$onlineSL_full_fit)
          colnames(forecast_tbl)[3] <- paste0("horizon", horizon)
          forecast_tblsFF[[i]] <- forecast_tbl
        }

      }
      base_dt <- data.table(
        id = rep(individual_id, length(horizons)),
        horizon = horizons, smooth_type=rep(smooth_type, length(horizons)),
        outcome_type = rep(outcome_type, length(horizons))
      )
      mean_performance_smooth[[s]] <- data.table(base_dt, rbindlist(mean_performance, fill=T))
      median_performance_smooth[[s]] <- data.table(base_dt, rbindlist(median_performance, fill=T))
      NAsummary_smooth[[s]] <- data.table(base_dt, rbindlist(NA_summary, fill=T))

      if(outcome_type == "binomial"){
        base_dt12 <- data.table(
          id=rep(individual_id, length(horizons)*3),
          horizon=rep(horizons, each=3),
          smooth_type=rep(smooth_type, length(horizons)*3),
          outcome_type = rep(outcome_type, length(horizons)*3)
        )
        AUC_smooth[[s]] <- data.table(base_dt12, rbindlist(AUCsummary, fill=T))
        AHEprev_smooth[[s]] <- data.table(base_dt, "AHEprev" = unlist(AHEprev))
      }

      movie_tbl <- suppressMessages(data.frame(movie_tbls %>% reduce(full_join)))
      movie_name <- paste0("/movie_Ytype_", outcome_type , "_smooth_",
                           smooth_type, "_id", individual_id, ".mp4")
      make_multiforecast_movie_mp4(movie_tbl, cols_list,
                                   name = paste0(dir, movie_name),
                                   outcome_type = outcome_type)

      forecast_tbl <- suppressMessages(data.table(forecast_tbls %>% reduce(full_join)))
      plot_name <-  paste0("forecasts_Ytype_", outcome_type , "_smooth_",
                           smooth_type, "_id", individual_id, ".pdf")
      learners <- sapply(horizons, function(h) paste0("horizon", h))
      forecast_plot <- make_forecast_plot(forecast_tbl, learners = learners,
                                          truth_line_size = 0.3,
                                          forecast_line_size = 0.2,
                                          outcome_type = outcome_type)
      ggsave(filename = plot_name, plot = forecast_plot, device = "pdf",
             path = paste0(dir, "/"), height=8, width=10)

      if(summarize_full_fit){
        movie_tbl <- suppressMessages(data.frame(movie_tblsFF %>% reduce(full_join)))
        movie_name <- paste0("/movieFF_Ytype_", outcome_type , "_smooth_",
                             smooth_type, "_id", individual_id, ".mp4")
        make_multiforecast_movie_mp4(movie_tbl, cols_list,
                                     name = paste0(dir, movie_name),
                                     outcome_type = outcome_type)

        forecast_tbl <- suppressMessages(data.table(forecast_tblsFF %>% reduce(full_join)))
        plot_name <- paste0("forecastsFF_Ytype_", outcome_type , "_smooth_",
                            smooth_type, "_id", individual_id, ".pdf")
        learners <- sapply(horizons, function(h) paste0("horizon", h))
        forecast_plot <- make_forecast_plot(forecast_tbl, learners = learners,
                                            truth_line_size = 0.3,
                                            forecast_line_size = 0.2,
                                            outcome_type = outcome_type)
        ggsave(filename = plot_name, plot = forecast_plot, device = "pdf",
               path = paste0(dir, "/"), height=8, width=10)
      }
    }
    mean_performance_Ytype[[y]] <- rbindlist(mean_performance_smooth, fill=T)
    median_performance_Ytype[[y]] <- rbindlist(median_performance_smooth, fill=T)
    NAsummary_Ytype[[y]] <- rbindlist(NAsummary_smooth, fill=T)

    if(outcome_type == "binomial"){
      tbl_path <- paste0(dir, "/", "summarized_performance_id",
                         individual_id, ".xlsx")
      AUCsummary <- rbindlist(AUC_smooth, fill=T)
      non_nums <- c("id", "horizon", "smooth_type", "outcome_type", "Metric")
      nums <- AUCsummary[,-non_nums,with=F]
      nums <- data.table(apply(nums, 2, function(x) as.numeric(x)))
      AUCsummary <- cbind(AUCsummary[,non_nums,with=F], nums)
      setorder(AUCsummary, smooth_type, horizon)
      AUC <- AUCsummary[Metric == "AUC",]
      AUCPR <- AUCsummary[Metric == "AUCPR",]
      nonNA_PREV <- AUCsummary[Metric == "nonNA_PREV",]
      AHEprev <- rbindlist(AHEprev_smooth, fill=T)
      write.xlsx(AHEprev, tbl_path, sheetName="AHEprev", row.names=F,
                 showNA=F)
      write.xlsx(AUC, tbl_path, sheetName="AUC", append=T, row.names=F,
                 showNA=F)
      write.xlsx(AUCPR, tbl_path, sheetName="AUCPR", append=T, row.names=F,
                 showNA=F)
      write.xlsx(nonNA_PREV, tbl_path, sheetName="nonNA_PREV", append=T,
                 row.names=F, showNA=F)
    }
  }
  mean_performance <- rbindlist(mean_performance_Ytype, fill=T)
  setorder(mean_performance, smooth_type, horizon)
  mean_performance_binY <- mean_performance[outcome_type == "binomial",]
  mean_performance_conY <- mean_performance[outcome_type == "continuous",]

  median_performance <- rbindlist(median_performance_Ytype, fill=T)
  setorder(median_performance, smooth_type, horizon)
  median_performance_binY <- median_performance[outcome_type == "binomial",]
  median_performance_conY <- median_performance[outcome_type == "continuous",]

  write.xlsx(mean_performance_binY, tbl_path, sheetName="mean_performance_binaryY",
             append = T, row.names=F, showNA=F)
  write.xlsx(median_performance_binY, tbl_path,
             sheetName="median_performance_binaryY", append=T, row.names=F,
             showNA=F)
  write.xlsx(mean_performance_conY, tbl_path,
             sheetName="mean_performance_continuousY",
             append = T, row.names=F, showNA=F)
  write.xlsx(median_performance_conY, tbl_path,
             sheetName="median_performance_continuousY", append=T, row.names=F,
             showNA=F)
  NAsummary <- rbindlist(NAsummary_Ytype, fill=T)
  setorder(NAsummary, outcome_type, smooth_type, horizon)
  write.xlsx(NAsummary, tbl_path, sheetName="NAsummary", row.names=F, showNA=F,
             append=T)
  print(paste0("Done summarizing id ", individual_id, "!!"))
}

summarize_all_ids <- function(results_path){
  # get ids
  ids <- list.files(results_path)[grep("^id", list.files(results_path))]
  ids <- gsub("id", "", ids)

  # result-type specific tables
  mean_performance_binaryY <- list()
  median_performance_binaryY <- list()
  mean_performance_continuousY <- list()
  median_performance_continuousY <- list()
  AUC <- list()
  AUCPR <- list()
  for(i in 1:length(ids)){
    id <- ids[i]
    path <- paste0(results_path, "id", ids[i], "/summarized_performance_id", ids[i], ".xlsx")
    wb <- loadWorkbook(path)
    wb <- readWorksheet(wb, sheet = getSheets(wb))
    AUC[[i]] <- wb$AUC
    AUCPR[[i]] <- wb$AUCPR
    mean_performance_binaryY[[i]] <- wb$mean_performance_binaryY
    median_performance_binaryY[[i]] <- wb$median_performance_binaryY
    mean_performance_continuousY[[i]] <- wb$mean_performance_continuousY
    median_performance_continuousY[[i]] <- wb$mean_performance_continuousY
  }
  summarize_list <- function(list_object, ignore = c("Metric", "id")){
    dt <- data.table(do.call(rbind, list_object))
    num_ids <- length(unique(dt[["id"]]))
    groups <- c("horizon", "smooth_type", "outcome_type")
    cols <- colnames(dt)[-which(colnames(dt) %in% c(groups, ignore))]
    dt <- dt[, lapply(.SD, mean), by=groups, .SDcols=cols]
    return(cbind(num_ids = rep(num_ids, nrow(dt)), dt))
  }
  AUC <- summarize_list(AUC)
  AUCPR <- summarize_list(AUCPR)
  mean_performance_binaryY <- summarize_list(mean_performance_binaryY, ignore = "id")
  median_performance_binaryY <- summarize_list(median_performance_binaryY, ignore = "id")
  mean_performance_continuousY <- summarize_list(mean_performance_continuousY, ignore = "id")
  median_performance_continuousY <- summarize_list(median_performance_continuousY, ignore = "id")

  sheets <- list(
    "AUC" = AUC, "AUCPR" = AUCPR,
    "mean_performance_binaryY" = mean_performance_binaryY,
    "median_performance_binaryY" = median_performance_binaryY,
    "mean_performance_continuousY" = mean_performance_continuousY,
    "median_performance_continuousY" = median_performance_continuousY
  )
  write_xlsx(sheets, paste0(results_path, "summary_all.xlsx"))
  print("Done!")
}

summarize_id_continuousY <- function(individual_id, individual_mean_data,
                                     individual_median_data, horizons,
                                     data_path, save_path,
                                     summarize_full_fit = TRUE){

  cat("\n************* Summarizing results for id", individual_id, "*************\n")

  # make main folder and subfolders
  newfolder <- paste0("id", individual_id)
  dir <- paste0(save_path, newfolder)
  if(!dir.exists(dir)){
    dir.create(dir)
  }

  outcome_types <- "continuous"
  mean_performance_Ytype <- list()
  median_performance_Ytype <- list()
  NAsummary_Ytype <- list()
  for (y in 1:length(outcome_types)){
    outcome_type <- outcome_types[y]
    smooth_types <- c("mean", "median")

    if(outcome_type == "binomial"){
      AUC_smooth <- list()
      AHEprev_smooth <- list()
    }
    mean_performance_smooth <- list()
    median_performance_smooth <- list()
    NAsummary_smooth <- list()
    for (s in 1:length(smooth_types)){
      smooth_type <- smooth_types[s]
      if(outcome_type == "binomial"){
        outcomes <- sapply(horizons, function(h) paste0("Y", h, "_AHE"))
      } else if (outcome_type == "continuous"){
        outcomes <- sapply(horizons, function(h) paste0("Y", h, "_lag5_", smooth_type))
      }

      if(smooth_type == "mean"){
        d <- individual_mean_data[id == individual_id,]
      } else if(smooth_type == "median"){
        d <- individual_median_data[id == individual_id,]
      }
      setorder(d, time_and_date)

      cols_list <- list()
      movie_tbls <- list()
      forecast_tbls <- list()
      mean_performance <- list()
      median_performance <- list()
      NA_summary <- list()
      if(summarize_full_fit){
        cols_listFF <- list()
        movie_tblsFF <- list()
        forecast_tblsFF <- list()
      }
      if(outcome_type == "binomial"){
        AUCsummary <- list()
        AHEprev <- list()
      }
      for(i in 1:length(outcomes)){
        outcome <- outcomes[i]
        if(grepl("AHE", outcome)){
          current_outcome <- "AHE"
          outcome_clear <- paste0(outcome, "_", smooth_type)
        } else {
          outcome_clear <- outcome
          if(grepl("mean", outcome)){
            current_outcome <- "abpmean_lag5_mean"
          } else if(grepl("median", outcome)){
            current_outcome <- "abpmean_lag5_median"
          }
        }
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

        results_path <- paste0(outcome_clear, "_id", individual_id, ".Rdata")
        load(file = paste0(data_path, results_path))

        # retain true outcomes
        truth_tbl <- d[, c(current_outcome, "min", outcome), with=F]
        colnames(truth_tbl) <- c("currentY","time", "futureY")
        truth_tbl$futureY_time <- truth_tbl$time + horizon

        # movie
        forecast_tbl <- result$forecast_table
        movie_tbl <- suppressMessages(full_join(truth_tbl, forecast_tbl))
        osl_forecasts <- data.frame(time=movie_tbl$futureY_time,
                                    forecast_time=movie_tbl$forecast_time,
                                    forecast=movie_tbl$onlineSL)
        colnames(osl_forecasts)[2] <- paste0("forecast_time_horizon", horizon)
        colnames(osl_forecasts)[3] <- paste0("horizon", horizon)
        cols_list[[i]] <- list(forecast_time = colnames(osl_forecasts)[2],
                               forecast = colnames(osl_forecasts)[3])
        truth <- data.frame(time=movie_tbl$futureY_time, truth=movie_tbl$futureY)
        movie_tbl <- suppressMessages(data.table(full_join(truth, osl_forecasts)))
        movie_tbls[[i]] <- movie_tbl

        # forecast plots
        forecast_tbl <- merge(truth_tbl, forecast_tbl, by="time", all.x=F, all.y=T)
        forecast_tbl <- data.table(futureY_time=forecast_tbl$futureY_time,
                                   futureY=forecast_tbl$futureY,
                                   forecast=forecast_tbl$onlineSL)
        colnames(forecast_tbl)[3] <- paste0("horizon", horizon)
        forecast_tbls[[i]] <- forecast_tbl

        # MSE summary
        perf_tbl <- result$performance_table
        if(identical(perf_tbl[["onlineSL"]], perf_tbl[["onlineSL_full_fit"]])){
          print(paste0("Fixing losses for outcome ", outcome_clear))
          forecast_tbl2 <- data.table(result$forecast_table)
          extra_cols <- c("time", "obs", "outcome_time", "forecast_time", "forecast_time_precise")
          forecasts <- forecast_tbl2[, -extra_cols, with=F]
          obs <- forecast_tbl2[["obs"]]
          perf_tbl <- data.table(apply(forecasts, 2, function(pred) (pred-obs)^2))
        } else {
          time_idx <- grepl("time", colnames(perf_tbl))
          perf_tbl <- perf_tbl[, -time_idx, with=F]
        }
        mean_performance[[i]] <- data.table(t(colMeans(perf_tbl, na.rm=T)))
        median_performance[[i]] <- data.table(t(apply(perf_tbl, 2, median, na.rm=T)))
        NA_summary[[i]] <- data.table(t(result$missing_summary))

        if(outcome_type == "binomial"){
          AUCsummary[[i]] <- result$AUC_tbl
          AHEprev[[i]] <- result$AHE_prevalence
        }

        if(summarize_full_fit){
          # movie
          forecast_tbl <- result$forecast_table
          movie_tbl <- suppressMessages(full_join(truth_tbl, forecast_tbl))
          osl_forecasts <- data.frame(time=movie_tbl$futureY_time,
                                      forecast_time=movie_tbl$forecast_time,
                                      forecast=movie_tbl$onlineSL_full_fit)
          colnames(osl_forecasts)[2] <- paste0("forecast_time_horizon", horizon)
          colnames(osl_forecasts)[3] <- paste0("horizon", horizon)
          cols_listFF[[i]] <- list(forecast_time = colnames(osl_forecasts)[2],
                                   forecast = colnames(osl_forecasts)[3])
          movie_tbl <- suppressMessages(data.table(full_join(truth, osl_forecasts)))
          movie_tblsFF[[i]] <- movie_tbl

          # forecast plots
          forecast_tbl <- merge(truth_tbl, forecast_tbl, by="time", all.x=F, all.y=T)
          forecast_tbl <- data.table(futureY_time=forecast_tbl$futureY_time,
                                     futureY=forecast_tbl$futureY,
                                     forecast=forecast_tbl$onlineSL_full_fit)
          colnames(forecast_tbl)[3] <- paste0("horizon", horizon)
          forecast_tblsFF[[i]] <- forecast_tbl
        }

      }
      base_dt <- data.table(
        id = rep(individual_id, length(horizons)),
        horizon = horizons, smooth_type=rep(smooth_type, length(horizons)),
        outcome_type = rep(outcome_type, length(horizons))
      )
      mean_performance_smooth[[s]] <- data.table(base_dt, rbindlist(mean_performance, fill=T))
      median_performance_smooth[[s]] <- data.table(base_dt, rbindlist(median_performance, fill=T))
      NAsummary_smooth[[s]] <- data.table(base_dt, rbindlist(NA_summary, fill=T))

      if(outcome_type == "binomial"){
        base_dt12 <- data.table(
          id=rep(individual_id, length(horizons)*3),
          horizon=rep(horizons, each=3),
          smooth_type=rep(smooth_type, length(horizons)*3),
          outcome_type = rep(outcome_type, length(horizons)*3)
        )
        AUC_smooth[[s]] <- data.table(base_dt12, rbindlist(AUCsummary, fill=T))
        AHEprev_smooth[[s]] <- data.table(base_dt, "AHEprev" = unlist(AHEprev))
      }

      movie_tbl <- suppressMessages(data.frame(movie_tbls %>% reduce(full_join)))
      movie_name <- paste0("/movie_Ytype_", outcome_type , "_smooth_",
                           smooth_type, "_id", individual_id, ".mp4")
      make_multiforecast_movie_mp4(movie_tbl, cols_list,
                                   name = paste0(dir, movie_name),
                                   outcome_type = outcome_type)

      forecast_tbl <- suppressMessages(data.table(forecast_tbls %>% reduce(full_join)))
      plot_name <-  paste0("forecasts_Ytype_", outcome_type , "_smooth_",
                           smooth_type, "_id", individual_id, ".pdf")
      learners <- sapply(horizons, function(h) paste0("horizon", h))
      forecast_plot <- make_forecast_plot(forecast_tbl, learners = learners,
                                          truth_line_size = 0.3,
                                          forecast_line_size = 0.2,
                                          outcome_type = outcome_type)
      ggsave(filename = plot_name, plot = forecast_plot, device = "pdf",
             path = paste0(dir, "/"), height=8, width=10)

      if(summarize_full_fit){
        movie_tbl <- suppressMessages(data.frame(movie_tblsFF %>% reduce(full_join)))
        movie_name <- paste0("/movieFF_Ytype_", outcome_type , "_smooth_",
                             smooth_type, "_id", individual_id, ".mp4")
        make_multiforecast_movie_mp4(movie_tbl, cols_list,
                                     name = paste0(dir, movie_name),
                                     outcome_type = outcome_type)

        forecast_tbl <- suppressMessages(data.table(forecast_tblsFF %>% reduce(full_join)))
        plot_name <- paste0("forecastsFF_Ytype_", outcome_type , "_smooth_",
                            smooth_type, "_id", individual_id, ".pdf")
        learners <- sapply(horizons, function(h) paste0("horizon", h))
        forecast_plot <- make_forecast_plot(forecast_tbl, learners = learners,
                                            truth_line_size = 0.3,
                                            forecast_line_size = 0.2,
                                            outcome_type = outcome_type)
        ggsave(filename = plot_name, plot = forecast_plot, device = "pdf",
               path = paste0(dir, "/"), height=8, width=10)
      }
    }
    mean_performance_Ytype[[y]] <- rbindlist(mean_performance_smooth, fill=T)
    median_performance_Ytype[[y]] <- rbindlist(median_performance_smooth, fill=T)
    NAsummary_Ytype[[y]] <- rbindlist(NAsummary_smooth, fill=T)
  }
  tbl_path <- paste0(dir, "/", "summarized_performance_id",
                     individual_id, ".xlsx")
  mean_performance <- rbindlist(mean_performance_Ytype, fill=T)
  setorder(mean_performance, smooth_type, horizon)
  write.xlsx(mean_performance, tbl_path,
             sheetName="mean_performance_continuousY",
             row.names=F, showNA=F)

  median_performance <- rbindlist(median_performance_Ytype, fill=T)
  setorder(median_performance, smooth_type, horizon)
  write.xlsx(median_performance, tbl_path,
             sheetName="median_performance_continuousY", append=T, row.names=F,
             showNA=F)
  NAsummary <- rbindlist(NAsummary_Ytype, fill=T)
  setorder(NAsummary, smooth_type, horizon)
  write.xlsx(NAsummary, tbl_path, sheetName="NAsummary", row.names=F, showNA=F,
             append=T)
  print(paste0("Done summarizing id ", individual_id, "!!"))
}

