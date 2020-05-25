library(tidyverse)
library(data.table)
library(here)
library(lubridate)
library(tsibble)
library(vrtest)
library(feasts)
library(forecast)

source(here::here("R", "v3", "utils_mimic.R"))

load(here::here("Data", "all_history30.Rdata"))
load(here::here("Data", "all_history60.Rdata"))

new30 <- all_history30 %>% 
  group_by(subject_id) %>%
  mutate(init_split = min_elapsed - lag(min_elapsed, n = 1) > 600) %>%
  select(c("subject_id", "min_elapsed", "init_split"))

ids_mult_visits <- new30 %>% 
  group_by(subject_id) %>%
  summarize(mult_split = sum(init_split, na.rm = T)) %>%
  filter(mult_split > 0)

ids_mult <- ids_mult_visits$subject_id

list_multi_stay_data <- lapply(ids_mult, function(id){
  dat <- new30 %>% filter(subject_id == id) 
  split_points <- append((which(dat$init_split == TRUE)-1), nrow(dat))
  return_dat <- matrix(NA, ncol = 5, nrow = nrow(dat))
  for(i in 1:length(split_points)){
    split_i <- split_points[i]
    if(i == 1){
      split_init <- 1
    } else {
      split_init <- split_points[i-1]+1
    }
    dat_split_i <- dat[split_init:split_i,]
    dat_split_i$subject_stay_id <- rep(paste0(id, "_", i), nrow(dat_split_i))
    dat_split_i$min <- 1:(split_i-split_init+1)
    return_dat[(split_init:split_i),] <- as.matrix(dat_split_i)
  }
  return_dat <- data.frame(return_dat)[,-3]
  colnames(return_dat) <- c("subject_id", "min_elapsed", "subject_stay_id", "min")
  return_dat$min_elapsed <- as.numeric(as.character(return_dat$min_elapsed))
  return_dat$min <- as.numeric(as.character(return_dat$min))
  return(return_dat)
})
multi_stay_data <- do.call(rbind, list_multi_stay_data)

single_stay_data <- new30 %>% 
  filter(subject_id %nin% ids_mult) %>%
  group_by(subject_id) %>%
  mutate(subject_stay_id = subject_id) %>%
  mutate(min = min_elapsed) %>%
  select(c("subject_id", "min_elapsed", "subject_stay_id", "min"))

fixed_data <- rbind.data.frame(multi_stay_data, single_stay_data)

d30 <- data.table(merge(fixed_data, all_history30, all.y = TRUE,
                          by = c("subject_id", "min_elapsed")))
d60 <- data.table(merge(fixed_data, all_history60, all.y = TRUE,
                        by = c("subject_id", "min_elapsed")))

d30 <- setorder(d30, subject_id, min_elapsed)
dim(d30) 
# [1] 5035962     163
d60 <- setorder(d60, subject_id, min_elapsed)

# fuck -- just realized the history is all fucked for the multiple stays.. 
d30_firststay <- d30 %>% filter(!grepl("_2|_3|_4", subject_stay_id)) 
d60_firststay <- d60 %>% filter(!grepl("_2|_3|_4", subject_stay_id))
dim(d30_firststay) 
# [1] 4434999     163
# data loss ~ 12%

################## redo history only on stays after the first ##################
source(here::here("R", "preprocess", "history.R"))

d30_laterstay <- d30 %>% filter(grepl("_2|_3|_4", subject_stay_id)) 
nrow(d30_laterstay) # 600963
d <- d30_laterstay[,1:30]
d$subject_stay_id <- droplevels(d$subject_stay_id)

dat_trt <- d %>%
  group_by(subject_stay_id) %>%
  mutate(a_last = dplyr::lag(amine, n = 1, default = NA)) %>% 
  mutate(amine_switch1 = ifelse(amine-a_last == 1, min, NA)) %>%
  mutate(amine_switch0 = ifelse(amine-a_last == -1, min, NA)) %>%
  mutate(s_last = dplyr::lag(sedation, n = 1, default = NA)) %>% 
  mutate(sedation_switch1 = ifelse(sedation-s_last == 1, min, NA)) %>%
  mutate(sedation_switch0 = ifelse(sedation-s_last == -1, min, NA)) %>%
  mutate(v_last = dplyr::lag(ventilation, n = 1, default = NA)) %>% 
  mutate(ventilation_switch1 = ifelse(ventilation-v_last == 1, min, NA)) %>%
  mutate(ventilation_switch0 = ifelse(ventilation-v_last == -1, min, NA)) %>%
  select(subject_id, min_elapsed, time_and_date,
         amine, amine_switch1, amine_switch0, 
         sedation, sedation_switch1, sedation_switch0,
         ventilation, ventilation_switch1, ventilation_switch0)
dat_full <- data.table(left_join(d, dat_trt))
dat <- dat_full[order(dat_full$subject_id, dat_full$time_and_date), ]

# remove subject_stay_ids with 
dat_subset <- dat %>%
  group_by(subject_stay_id) %>%
  filter(n() >= 30)
nrow(dat_subset) #600775
dat_subset$subject_stay_id <- droplevels(dat_subset$subject_stay_id)

N <- length(unique(dat_subset$subject_stay_id))
registerDoParallel(cores = (detectCores()-14)) # set up parallelization
getDoParWorkers() 
binary_history_list <- foreach(n = 1:N) %dopar% {
  id <- levels(dat_subset$subject_stay_id)[n]
  ind_dat <- dplyr::filter(dat_subset, subject_stay_id == id)
  
  # considering history as past 30 min and past 60 min
  history30 <- add_history_binary(ind_dat, 30)
  history60 <- add_history_binary(ind_dat, 60)
  
  if(nrow(history30) != nrow(ind_dat)){
    print(paste0("nrow fuck up for subject_stay_id level ", n))
  }
  
  return_list <- list(history30 = history30,
                      history60 = history60)
  
  print(paste0("Finishing subject_id level ", n))
  
  return(return_list)
}

binary_history60 <- do.call(rbind, lapply(binary_history_list, "[[", "history60"))
binary_history30 <- do.call(rbind, lapply(binary_history_list, "[[", "history30"))

comparison30 <- compare(binary_history30[,c(1:36)], dat_subset, allowAll = TRUE)
comparison30$result
comparison60 <- compare(binary_history60[,c(1:36)], dat_subset, allowAll = TRUE)
comparison60$result

getDoParWorkers() 
continuous_history_list <- foreach(n = 1:N) %dopar% {
  id <- levels(binary_history30$subject_stay_id)[n]
  ind_dat_30 <- dplyr::filter(binary_history30, subject_stay_id == id)
  ind_dat_60 <- dplyr::filter(binary_history60, subject_stay_id == id)
  
  # considering history as past 30 min and past 60 min
  history30 <- add_history_continuous(ind_dat_30, 30)
  history60 <- add_history_continuous(ind_dat_60, 60)
  
  if(nrow(history30) != nrow(ind_dat_30)){
    print(paste0("nrow fucked up for subject_stay_id level ", n))
  }
  
  return_list <- list(history30 = history30,
                      history60 = history60)
  
  print(paste0("Finishing subject_id level ", n))
  
  return(return_list)
}

history30 <- lapply(continuous_history_list, '[[', 'history30')
history60 <- lapply(continuous_history_list, '[[', 'history60')

# check all column names equal
apply(do.call(rbind,lapply(history30,colnames)), 2, function(x) length(unique(x)) == 1)
apply(do.call(rbind,lapply(history60,colnames)), 2, function(x) length(unique(x)) == 1)

multi_history30 <- data.table(do.call(rbind, history30))
multi_history60 <- data.table(do.call(rbind, history60))

# comparison check
comparison30 <- compare(multi_history30[,c(1:63)], binary_history30, allowAll = TRUE)
comparison30$result
comparison60 <- compare(multi_history60[,c(1:63)], binary_history60, allowAll = TRUE)
comparison60$result

# ensure successful merge, there was an issue with the summary stats
colnames(d30_firststay)[which(colnames(d30_firststay) %nin% colnames(multi_history30))]
name <- c("abpsys_Min", "abpsys_1stQ", "abpsys_Median", "abpsys_Mean", "abpsys_3rdQ",
          "abpsys_Max", "abpdias_Min", "abpdias_1stQ", "abpdias_Median", "abpdias_Mean",   
          "abpdias_3rdQ", "abpdias_Max", "abpmean_Min", "abpmean_1stQ", "abpmean_Median",
          "abpmean_Mean", "abpmean_3rdQ", "abpmean_Max", "spo2_Min", "spo2_1stQ", 
          "spo2_Median", "spo2_Mean", "spo2_3rdQ", "spo2_Max")
colnames(multi_history30)[c(64:69,79:84,94:99,109:114)] <- name
colnames(multi_history60)[c(64:69,79:84,94:99,109:114)] <- name

multi_fix30 <- multi_history30[,lapply(.SD, function(x) gsub("Min.|Max.", "", x)), 
                               .SDcols = name]
multi_fix30 <- multi_fix30[,lapply(.SD, function(x) gsub("1st Qu.|3rd Qu.", "", x)), 
                           .SDcols = name]
multi_fix30 <- multi_fix30[,lapply(.SD, function(x) gsub("[^0-9.-]", "", x)), 
                           .SDcols = name]
multi_fix30 <- multi_fix30[, lapply(.SD, as.numeric), .SDcols = name]
fixed30 <- cbind(multi_fix30, multi_history30[,-name,with=FALSE])

multi_fix60 <- multi_history60[,lapply(.SD, function(x) gsub("Min.|Max.", "", x)), 
                               .SDcols = name]
multi_fix60 <- multi_fix60[,lapply(.SD, function(x) gsub("1st Qu.|3rd Qu.", "", x)), 
                           .SDcols = name]
multi_fix60 <- multi_fix60[,lapply(.SD, function(x) gsub("[^0-9.-]", "", x)), 
                           .SDcols = name]
multi_fix60 <- multi_fix60[, lapply(.SD, as.numeric), .SDcols = name]
fixed60 <- cbind(multi_fix60, multi_history60[,-name,with=FALSE])

d30 <- data.table(rbind(d30_firststay, fixed30))
d60 <- data.table(rbind(d60_firststay, fixed60))

mimic30 <- setorder(d30, "subject_id", "min_elapsed")
mimic60 <- setorder(d60, "subject_id", "min_elapsed")

# remove cols with features corresponding to abpsys
dups <- c("abpdias_trend_strength", "abpdias_spikiness", "abpdias_linearity", 
          "abpdias_curvature", "abpdias_stl_e_acf1", "abpdias_stl_e_acf10",
          "abpdias_spectral_entropy", "abpdias_n_flat_spots", 
          "abpdias_coef_hurst", "abpmean_trend_strength", "abpmean_spikiness",
          "abpmean_linearity", "abpmean_curvature", "abpmean_stl_e_acf1",
          "abpmean_stl_e_acf10", "abpmean_spectral_entropy", 
          "abpmean_n_flat_spots", "abpmean_coef_hurst", "spo2_trend_strength",
          "spo2_spikiness", "spo2_linearity", "spo2_curvature", 
          "spo2_stl_e_acf1", "spo2_stl_e_acf10", "spo2_spectral_entropy",
          "spo2_n_flat_spots", "spo2_coef_hurst")
mimic30 <- mimic30[,-dups, with=FALSE]
mimic60 <- mimic60[,-dups, with=FALSE]

# remove cols from binary history that are no longer relevant
irrelevant <- c("amine_switch1", "amine_switch0", "sedation_switch1", 
                "sedation_switch0", "ventilation_switch1", "ventilation_switch0")
mimic30 <- mimic30[,-irrelevant, with=FALSE]
mimic60 <- mimic60[,-irrelevant, with=FALSE]

# change leads to lags
mimic30 <- setorder(mimic30, "subject_id", "min_elapsed")
mimic60 <- setorder(mimic60, "subject_id", "min_elapsed")

vars <- c("abpsys", "abpdias", "abpmean", "spo2")
lag_names <- lapply(seq(10), function(x) paste0(names, "_lag_", x))
leads_to_lags <- function(dt){
  dt[, (lag_names[[1]]) := shift(.SD, n=1), by=subject_stay_id, .SDcols=vars]
  dt[, (lag_names[[2]]) := shift(.SD, n=2), by=subject_stay_id, .SDcols=vars]
  dt[, (lag_names[[3]]) := shift(.SD, n=3), by=subject_stay_id, .SDcols=vars]
  dt[, (lag_names[[4]]) := shift(.SD, n=4), by=subject_stay_id, .SDcols=vars]
  dt[, (lag_names[[5]]) := shift(.SD, n=5), by=subject_stay_id, .SDcols=vars]
  dt[, (lag_names[[6]]) := shift(.SD, n=6), by=subject_stay_id, .SDcols=vars]
  dt[, (lag_names[[7]]) := shift(.SD, n=7), by=subject_stay_id, .SDcols=vars]
  dt[, (lag_names[[8]]) := shift(.SD, n=8), by=subject_stay_id, .SDcols=vars]
  dt[, (lag_names[[9]]) := shift(.SD, n=9), by=subject_stay_id, .SDcols=vars]
  dt[, (lag_names[[10]]) := shift(.SD, n=10), by=subject_stay_id, .SDcols=vars]
  return(dt)
}

m30 <- leads_to_lags(mimic30)
m60 <- leads_to_lags(mimic60)

save(mimic30, file = here::here("Data","mimic30.Rdata"), compress = TRUE)
save(mimic60, file = here::here("Data","mimic60.Rdata"), compress = TRUE)

################################ history functions #############################
summarize_history_binary <- function(history_indicator_vector,
                                     history_switch1_time_vector,
                                     history_switch0_time_vector,
                                     current_time, name){
  
  # 1. proportion of time ON in history (for intercurrent event = 1)
  mean <- mean(history_indicator_vector)
  
  # 3. most recent switch ON time in history (in terms of min)
  switch1 <- suppressWarnings(max(history_switch1_time_vector, na.rm = TRUE)) 
  switch1_time <- ifelse(is.infinite(switch1), 0, switch1)
  # 4. most recent switch ON time in history relative to current time 
  #    (i.e. min elapsed since switch)
  switch1_time_relative <- ifelse(is.infinite(switch1), 0, current_time-switch1)
  # 5. whether or not there was switch ON in the history
  switch1 <- ifelse(is.infinite(switch1), 0, 1)
  # 6. how many switches ON in the history
  switch1_count <- sum(ifelse(history_switch1_time_vector != 0, 1, 0), na.rm = TRUE)
  
  # 7. most recent switch OFF time in history (in terms of min)
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
                                      history$amine_switch0, current$min,
                                      "amine_history_")
    sed <- summarize_history_binary(history$sedation, history$sedation_switch1, 
                                    history$sedation_switch0, current$min,
                                    "sedation_history_")
    vent <- summarize_history_binary(history$ventilation, history$ventilation_switch1,
                                     history$ventilation_switch0, current$min,
                                     "ventilation_history_")
    current_history_binary[[i]] <- data.frame(current, t(amine), t(sed), t(vent))
  }
  df <- do.call(rbind, current_history_binary)
  return_df <- suppressMessages(full_join(dat_ordered[1,], df))
  return(return_df)
}

min.true <- function(x){
  min(which(x==!TRUE))-1
}

summarize_history_continuous <- function(history, current, var="abpsys"){
  
  hist <- nrow(history)
  
  full_data <- rbind.data.frame(history, current)
  full_data <- full_data[order((hist+1):1),]
  
  #Basic statistics:
  sum_hist <- unclass(summary(history[,var]))
  sum_hist_names <- c("Min","1stQ","Median","Mean","3rdQ","Max")
  names(sum_hist) <- paste0(var, "_", sum_hist_names)

  #Time-series tools
  history_ts <- as_tsibble(history, index=min)
  full_data_ts<-as_tsibble(full_data, index=min)

  #Extract time-series features:
  #1) Features related to STL decomposition of the series (strength of trend, seasonality)
  #2) Spectral entropy of a time-series
  #3) Number of flat points
  #4) Hurst coefficient (level of fractional differencing of a series)
  
  all_names_features <- c("trend_strength","spikiness","linearity","curvature","stl_e_acf1",     
                          "stl_e_acf10","spectral_entropy","n_flat_spots","coef_hurst")
  all_names_features <- paste0(var, "_", all_names_features)
  
  if(nrow(history_ts) == 1){
    features_ts <- rep(NA, 9)
    features_ts <- data.frame(t(features_ts))
    names(features_ts) <- all_names_features
  }else{
    features_ts <- suppressWarnings(history_ts %>% 
                                      features(abpsys, list(feat_stl, feat_spectral, 
                                                            n_flat_spots, coef_hurst))) 
    names_features <- names(features_ts)  
    names(features_ts) <- paste0(var, "_", names_features)
    
  }
  
  #Check all worked
  if(length(features_ts) != 9){
    Missing <- setdiff(all_names_features, names(features_ts))
    features_ts[Missing] <- NA
    features_ts <- features_ts[,all_names_features]
  }
  
  return(cbind(t(sum_hist), features_ts))
}

add_history_continuous <- function(ind_data, min_history){
  dat_ordered <- ind_data[order(ind_data$time_and_date), ]
  
  current_history <- list()
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
    
    abpsys <- summarize_history_continuous(history=history, current=current, 
                                           var="abpsys")
    abpdias <- summarize_history_continuous(history=history, current=current, 
                                            var="abpdias")
    abpmean <- summarize_history_continuous(history=history, current=current, 
                                            var="abpmean")
    spo2 <- summarize_history_continuous(history=history, current=current, 
                                         var="spo2")
    
    current_history[[i]] <- data.frame(current, abpsys, abpdias, abpmean, spo2)
  }
  df <- do.call(rbind, current_history)
  return_df <- suppressMessages(full_join(dat_ordered[1,], df))
  
  #Add lags now:
  lags_sys <- data.frame(matrix(NA, nrow = nrow(return_df), 10))
  lags_dias <- data.frame(matrix(NA, nrow = nrow(return_df), 10))
  lags_mean <- data.frame(matrix(NA, nrow = nrow(return_df), 10))
  lags_spo2 <- data.frame(matrix(NA, nrow = nrow(return_df), 10))
  
  for(l in 1:10){
    lags_sys[,l]  <- lead(return_df$abpsys,l)
    lags_dias[,l] <- lead(return_df$abpdias,l)
    lags_mean[,l] <- lead(return_df$abpmean,l)
    lags_spo2[,l] <- lead(return_df$spo2,l)
  }
  
  names(lags_sys) <- paste0("abpsys_lag_", seq(10))
  names(lags_dias) <- paste0("abpdias_lag_", seq(10))
  names(lags_mean) <- paste0("abpmean_lag_", seq(10))
  names(lags_spo2) <- paste0("spo2_lag_", seq(10))
  
  return_df <- cbind.data.frame(return_df,lags_sys,lags_dias,lags_mean,lags_spo2)  
  
  return(return_df)
}