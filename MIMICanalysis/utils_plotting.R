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

source(here::here("R", "v3", "utils_mimic.R"))

visualize <- function(individual_data, smoothed_data=F, smooth_type=NULL){
  
  if(smoothed_data){
    y_name_bp_plot <- paste0("5-Min Smoothed ", smooth_type, " Mean BP")
    y_name_vitals_plot <- paste0("5-Min Smoothed ", smooth_type, " Vitals")
    abpmean_name <- paste0("abpmean_lag5_", smooth_type)
    spo2_name <- paste0("spo2_lag5_", smooth_type)
    hr_name <- paste0("hr_lag5_", smooth_type)
    abpsys_name <- paste0("abpsys_lag5_", smooth_type)
    abpdias_name <- paste0("abpdias_lag5_", smooth_type)
  } else {
    y_name_bp_plot <- "Mean BP"
    y_name_vitals_plot <- "Vitals"
    abpmean_name <- "abpmean"
    abpsys_name <- "abpsys"
    abpdias_name <- "abpdias"
    spo2_name <- "spo2"
    hr_name <- "hr"
  }
  
  
  # make summary tables
  bp_tbl_summary <- make_bp_tbl_summary(individual_data, abpmean_name)
  missing <- make_missing(bp_tbl_summary)
  
  # make plots 
  p1 <- bp_plot(individual_data, abpmean_name, y_name=y_name_bp_plot)
  p2 <- trt_plot(individual_data) 
  p3 <- vitals_plot(individual_data, abpsys_name, abpdias_name, 
                    spo2_name, hr_name, y_name=y_name_vitals_plot)
  p4 <- locf_plot(bp_tbl_summary)
  col1 <- ggarrange(p1, p3, nrow = 2, ncol = 1)
  col2 <- ggarrange(p2, p4, nrow = 2, ncol = 1)
  plots <- ggarrange(col1, col2, nrow = 1, ncol = 2, widths=c(1.5,1))
  
  # make table of W for bottom of plots
  base <- c("id", "rank_icu", "sex", "age", "sapsi_first", "sofa_first", 
            "bmi", "care_unit", "admission_type_descr", "subject_id")
  tblW <- unique(individual_data[, base, with = FALSE])
  tbl <- cbind(tblW, missing)
  tbl$bmi <- round(tbl$bmi, 2)
  tbl <- ggtexttable(tbl[,-1], rows=NULL, 
                     theme=ttheme(base_size=7, base_colour="grey80"))
  fig <- ggarrange(plots, tbl, ncol=1, nrow=2, heights=c(6,1), widths=c(1.5,1))
  title <- paste0("Subject ID: ",tblW$subject_id," & Subject Stay ID: ",tblW$id)
  annotate_figure(fig, top=text_grob(title, size=10, face="bold"))
}

### hourly summaries 
make_bp_tbl_summary <- function(individual_data, abpmean_name){
  bp_tbl_summary <- individual_data %>% 
    dplyr::mutate(hour = as.integer(min/60) + 1) %>%
    dplyr::group_by(hour) %>%
    dplyr::summarize_at(.vars = c(abpmean_name), 
                        .funs = list(abpmean_mean = mean, abpmean_median = median, 
                                     abpmean_min = min, abpmean_max = max))
  bp_tbl_summary2 <- individual_data %>%  
    dplyr::mutate(abpmean_locf = as.numeric(as.character(abpmean_locf))) %>%
    dplyr::mutate(amine = as.numeric(as.character(amine))) %>%
    dplyr::mutate(sedation = as.numeric(as.character(sedation))) %>%
    dplyr::mutate(ventilation = as.numeric(as.character(ventilation))) %>%
    dplyr::mutate(row_locf = as.numeric(as.character(row_locf))) %>%
    dplyr::mutate(hour = as.integer(min/60) + 1) %>%
    dplyr::group_by(hour) %>%
    dplyr::summarize_at(.vars = c("abpmean_locf", "amine", "sedation", "row_locf",
                                  "ventilation", "hr", "spo2"), 
                        .funs = list(mean = mean))
  bp_tbl_summary <- merge(bp_tbl_summary, bp_tbl_summary2, by = c("hour"))
  bp_tbl_summary$hour <- as.numeric(as.factor(bp_tbl_summary$hour))
  bp_tbl_summary[order(bp_tbl_summary$hour),]
}

### gaps of data
make_missing <- function(bp_tbl_summary){
  max_hour <- bp_tbl_summary %>%
    summarize(hour_max = max(hour))
  len <- bp_tbl_summary %>%
    summarize(hour_count = n())
  missing <- cbind(max_hour, len)
  missing$diff <- missing$hour_max - missing$hour_count
  return(missing)
}

##################### plotting functions to visualize patterns #################
bp_plot <- function(individual_data, abpmean_name, y_name){
  sub_hypo <- individual_data %>%
    mutate(hour = as.numeric(min/60 + 1)) %>%
    select(y=all_of(abpmean_name), x=hour)
  ggplot(sub_hypo, aes(x=x, y=y)) +
    geom_line(size=.4) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    geom_hline(yintercept = 65, color = "red",size=.4) +
    labs(x = "Time (hours)", y = y_name) + 
    theme(text = element_text(size = 8))
}

locf_plot <- function(bp_tbl_summary){
  ggplot(bp_tbl_summary, aes(x = hour, y = row_locf_mean)) +
    geom_line() +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(x = "Time (hours)", y = "Hourly Prop of Row-wise LOCF") + 
    theme(text = element_text(size = 8))
}

vitals_plot <- function(individual_data, abpsys_name, abpdias_name, spo2_name,
                        hr_name, y_name){
  sub_hypo <- individual_data %>%
    mutate(hour = as.numeric(min/60 + 1)) %>%
    select(SBP=all_of(abpsys_name), DBP=all_of(abpdias_name), 
           SpO2=all_of(spo2_name), HR=all_of(hr_name), hour)
  sub_hypo <- melt(sub_hypo, id.vars = "hour", 
                   measure.vars = c("SBP", "DBP", "SpO2", "HR"))
  #hypo_line <- data.frame(variable = "abpmean", threshold = 65)
  ggplot(sub_hypo, aes(x = hour, y = value)) +
    geom_line(size=.4) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(x = "Time (hours)", y=y_name) +
    facet_grid(rows = vars(variable), scales = "free") +
    #geom_hline(data = hypo_line, aes(yintercept = threshold), color = "red") + 
    theme(text = element_text(size = 8))
}

trt_plot <- function(individual_data){
  sub_hypo <- individual_data %>%
    mutate(hour = as.numeric(min/60 + 1)) %>%
    dplyr::mutate(amine = as.numeric(as.character(amine))) %>%
    dplyr::mutate(sedation = as.numeric(as.character(sedation))) %>%
    dplyr::mutate(ventilation = as.numeric(as.character(ventilation))) %>%
    select(c("hour", "amine", "sedation", "ventilation"))
  sub_hypo <- melt(sub_hypo, id.vars = "hour", 
                   measure.vars = c("amine", "sedation", "ventilation"))
  ggplot(sub_hypo, aes(x = hour, y = value)) +
    geom_line() +
    scale_y_continuous(breaks = c(0,1), labels = c(0,1)) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(x = "Time (hours)", y = "Treatmeant") +
    facet_grid(rows = vars(variable)) + 
    theme(text = element_text(size = 8))
}