################################# forecast plot ################################
make_forecast_plot <- function(forecast_dt, learners = "onlineSL",
                               forecast_line_size = 0.4, truth_line_size=1,
                               outcome_type){
  cols <- c("#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#E6AB02",
            "#1F78B4", "#A6761D")

  learner_dt <- forecast_dt[, learners, with=F]
  df <- cbind.data.frame(time=forecast_dt$futureY_time,
                         truth=forecast_dt$futureY, learner_dt)

  df <- melt(df, id.vars = "time")
  df$mysize <- rep(forecast_line_size, nrow(df))
  df$mysize[df$variable == "truth"] <- truth_line_size

  if(outcome_type == "continuous"){
    plot <- ggplot(data = df, aes(x=time, y=value, size=mysize, color=variable)) +
      geom_line() +
      ylim(0,120) +
      scale_size(range = c(forecast_line_size, truth_line_size), guide="none") +
      labs(x="Time (minutes)", y="MAP",
           title="Online SL Forecasts Overlapping True MAP", color="Type") +
      scale_color_manual(values = cols) +
      theme(legend.title = element_blank())
  } else if(outcome_type == "binomial"){
    plot <- ggplot(data = df, aes(x=time, y=value, size=mysize, color=variable)) +
      geom_line() +
      ylim(0,1) +
      scale_size(range = c(forecast_line_size, truth_line_size), guide="none") +
      labs(x="Time (minutes)", y="AHE",
           title="Online SL Forecasts Overlapping True AHE", color="Type") +
      scale_color_manual(values = cols) +
      theme(legend.title = element_blank())
  }

  return(plot)
}

############################ learner proportion chart ##########################
make_prop_plot <- function(prop_tbl){
  df <- melt(prop_tbl, id.vars = "time")
  plot <- ggplot(data = df, aes(x=time, y=value, color=variable)) +
    geom_point(size=0.4) +
    labs(x = "Time (minutes)", y = "Contribution (%)",
         title = "Learner Type Contributions to Online SL", color = "Type") +
    scale_color_brewer(palette = "Dark2") +
    scale_y_continuous(breaks = c(0,.25,.5,.75,1), labels = c(0,25,50,75,100)) +
    theme(legend.position = "top")
  return(plot)
}
################################# patient chart ################################
make_patient_chart <- function(individual_data, smoothed_data=TRUE, smooth_type=NULL){

  if(smoothed_data){
    y_name_bp_plot <- paste0("5-Min Smoothed ", smooth_type, " MAP")
    y_name_vitals_plot <- paste0("5-Min Smoothed ", smooth_type, " Vitals")
    abpmean_name <- paste0("abpmean_lag5_", smooth_type)
    spo2_name <- paste0("spo2_lag5_", smooth_type)
    hr_name <- paste0("hr_lag5_", smooth_type)
    abpsys_name <- paste0("abpsys_lag5_", smooth_type)
    abpdias_name <- paste0("abpdias_lag5_", smooth_type)
  } else {
    y_name_bp_plot <- "MAP"
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
  base <- c("id", "sex", "age", "sapsi_first", "sofa_first",
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
  max_hour <- max(bp_tbl_summary$hour)
  hour_count <- nrow(bp_tbl_summary)
  missing <- data.frame(max_hour, hour_count)
  missing$diff <- missing$max_hour - missing$hour_count
  return(missing)
}

##################### plotting functions to visualize patterns #################
bp_plot <- function(individual_data, abpmean_name, y_name){
  sub_hypo <- individual_data %>%
    select(y=all_of(abpmean_name), x=min)
  ggplot(sub_hypo, aes(x=x, y=y)) +
    geom_line(size=.4) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    geom_hline(yintercept = 65, color = "red",size=.4) +
    labs(x = "Time (minutes)", y = y_name) +
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
    select(SBP=all_of(abpsys_name), DBP=all_of(abpdias_name),
           SpO2=all_of(spo2_name), HR=all_of(hr_name), min)
  sub_hypo <- melt(sub_hypo, id.vars = "min",
                   measure.vars = c("SBP", "DBP", "SpO2", "HR"))
  #hypo_line <- data.frame(variable = "abpmean", threshold = 65)
  ggplot(sub_hypo, aes(x = min, y = value)) +
    geom_line(size=.4) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(x = "Time (minutes)", y=y_name) +
    facet_grid(rows = vars(variable), scales = "free") +
    #geom_hline(data = hypo_line, aes(yintercept = threshold), color = "red") +
    theme(text = element_text(size = 8))
}

trt_plot <- function(individual_data){
  sub_hypo <- individual_data %>%
    dplyr::mutate(amine = as.numeric(as.character(amine))) %>%
    dplyr::mutate(sedation = as.numeric(as.character(sedation))) %>%
    dplyr::mutate(ventilation = as.numeric(as.character(ventilation))) %>%
    select(c("min", "amine", "sedation", "ventilation"))
  sub_hypo <- melt(sub_hypo, id.vars = "min",
                   measure.vars = c("amine", "sedation", "ventilation"))
  ggplot(sub_hypo, aes(x = min, y = value)) +
    geom_line() +
    scale_y_continuous(breaks = c(0,1), labels = c(0,1)) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(x = "Time (minutes)", y = "Treatmeant") +
    facet_grid(rows = vars(variable)) +
    theme(text = element_text(size = 8))
}
