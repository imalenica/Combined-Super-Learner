library(tidyverse)
library(data.table)
library(here)

load(here::here("Data","mimic.Rdata"))

# fill in gaps in min_elapsed by creating new row and carring last row forward 
dat <- data.table(mimic)
setkey(dat, subject_id, min_elapsed)
dat <- dat[CJ(unique(subject_id), seq(min(min_elapsed), max(min_elapsed))),
           roll = TRUE] 

# create `all_locf` column to indicate whether or not last row was carried forward
anti_dat <- anti_join(dat, mimic)
anti_dat$all_locf <- rep(1, nrow(anti_dat))
semi_dat <- semi_join(dat, mimic)
semi_dat$all_locf <- rep(0, nrow(semi_dat))
dat <- full_join(anti_dat, semi_dat)

dat <- dat[order(dat$subject_id, dat$time_and_date), ]

############################# add leaded outcome ###############################
data <- dat %>%
  group_by(subject_id) %>%
  mutate(Y_15 = dplyr::lead(abpmean, n = 15, default = NA)) %>%
  mutate(Y_20 = dplyr::lead(abpmean, n = 20, default = NA)) %>%
  mutate(Y_25 = dplyr::lead(abpmean, n = 25, default = NA)) %>%
  mutate(Y_30 = dplyr::lead(abpmean, n = 30, default = NA)) 

########### add time series summaries / history to current time data ###########

### for time-varying binary variables, we may be interested in on/off switch:
# (1) identify min when trt was switched, and 
# (2) if it was a switch on trt (diff=0-1=-1) or off trt (diff=1-0=1)

data$amine <- as.numeric(data$amine)-1
data$sedation <- as.numeric(data$sedation)-1
data$ventilation <- as.numeric(data$ventilation)-1

dat_trt <- data %>%
  group_by(subject_id) %>%
  mutate(a_last = dplyr::lag(amine, n = 1, default = NA)) %>% 
  mutate(amine_switch1 = ifelse(amine-a_last == 1, min_elapsed-1, NA)) %>%
  mutate(amine_switch0 = ifelse(amine-a_last == -1, min_elapsed-1, NA)) %>%
  mutate(s_last = dplyr::lag(sedation, n = 1, default = NA)) %>% 
  mutate(sedation_switch1 = ifelse(sedation-s_last == 1, min_elapsed-1, NA)) %>%
  mutate(sedation_switch0 = ifelse(sedation-s_last == -1, min_elapsed-1, NA)) %>%
  mutate(v_last = dplyr::lag(ventilation, n = 1, default = NA)) %>% 
  mutate(ventilation_switch1 = ifelse(ventilation-v_last == 1, min_elapsed-1, NA)) %>%
  mutate(ventilation_switch0 = ifelse(ventilation-v_last == -1, min_elapsed-1, NA)) %>%
  select(subject_id, min_elapsed, time_and_date,
         amine, amine_switch1, amine_switch0, 
         sedation, sedation_switch1, sedation_switch0,
         ventilation, ventilation_switch1, ventilation_switch0)

dat_full <- left_join(dat_trt, data)
dat <- dat_full[order(dat_full$subject_id, dat_full$time_and_date), ]

##### time-varying binary trt summaries
summarize_history_binary <- function(history_indicator_vector,
                                     history_switch1_time_vector,
                                     history_switch0_time_vector,
                                     current_time, name){
  
  # 1. proportion of time ON in history (for intercurrent event = 1)
  mean <- mean(history_indicator_vector)

  # 3. most recent switch ON time in history (in terms of min_elapsed)
  switch1 <- suppressWarnings(max(history_switch1_time_vector, na.rm = TRUE)) 
  switch1_time <- ifelse(is.infinite(switch1), 0, switch1)
  # 4. most recent switch ON time in history relative to current time 
  #    (i.e. min elapsed since switch)
  switch1_time_relative <- ifelse(is.infinite(switch1), 0, current_time-switch1)
  # 5. whether or not there was switch ON in the history
  switch1 <- ifelse(is.infinite(switch1), 0, 1)
  # 6. how many switches ON in the history
  switch1_count <- sum(ifelse(history_switch1_time_vector != 0, 1, 0), na.rm = TRUE)
  
  # 7. most recent switch OFF time in history (in terms of min_elapsed)
  switch0 <- suppressWarnings(max(history_switch0_time_vector, na.rm = TRUE)) 
  switch0_time <- ifelse(is.infinite(switch0), 0, switch0)
  # 8. most recent switch OFF time in history relative to current time 
  #    (i.e. min elapsed since switch)
  switch0_time_relative <- ifelse(is.infinite(switch0), 0, current_time-switch0)
  # 9. whether or not there was switch OFF in the history
  switch0 <- ifelse(is.infinite(switch0), 0, 1)
  # 10. how many switches OFF in the history
  switch0_count <- sum(ifelse(history_switch0_time_vector != 0, 1, 0), na.rm = TRUE)
  
  return_vec <- c(mean = mean, switch1 = switch1, switch1_time = switch1_time, 
                  switch1_time_relative = switch1_time_relative, 
                  switch1_count = switch1_count, switch0 = switch0, 
                  switch0_time = switch0_time, 
                  switch0_time_relative = switch0_time_relative, 
                  switch0_count = switch0_count)
  names(return_vec) <- paste0(name, names(return_vec))
  return(return_vec)
}

add_history_binary <- function(ind_data, min_history){
  dat_ordered <- ind_data[order(ind_data$time_and_date), ]
  
  current_history_binary <- list()
  for(i in 2:nrow(dat_ordered)){
    history_full <- dat_ordered[c(1:(i-1)),]
    # subset to history if i > min_history
    if(i <= min_history){
      history <- history_full
    } else {
      history <- history_full[c(
        (nrow(history_full)-min_history+1):(nrow(history_full))),]
    }
    current <- dat_ordered[i,]
    amine <- summarize_history_binary(history$amine, history$amine_switch1, 
                                      history$amine_switch0, current$min_elapsed,
                                      "amine_history_")
    sed <- summarize_history_binary(history$sedation, history$sedation_switch1, 
                                    history$sedation_switch0, current$min_elapsed,
                                    "sedation_history_")
    vent <- summarize_history_binary(history$ventilation, history$ventilation_switch1, 
                                     history$ventilation_switch0, current$min_elapsed,
                                     "ventilation_history_")
    current_history_binary[[i]] <- c(current, amine, sed, vent)
  }
  df <- do.call(rbind, current_history_binary)
  left_join(return_obj, dat_ordered[1,])
  return(return_obj)
}

N <- length(unique(dat$subject_id))
registerDoParallel(cores = (detectCores()-1)) # set up parallelization
getDoParWorkers() 
binary_history_list <- foreach(n = 1:N) %dopar% {
  id <- dat$subject_id[N]
  ind_dat <- dplyr::filter(dat, subject_id == id)
  
  # considering history as past 30 min and past 60 min
  history60 <- add_history_binary(ind_dat, 30)
  history30 <- add_history_binary(ind_dat, 60)
  
  return_list <- list(history60 = history60,
                      history30 = history30)
  return(return_list)
}

binary_history60_all <- do.call(rbind, lapply(binary_history_list, "[[", "history60"))
binary_history30_all <- do.call(rbind, lapply(binary_history_list, "[[", "history30"))