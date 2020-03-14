library(tidyverse)
library(data.table)
library(here)
library(doParallel)
library(compare)

load(here::here("Data","mimic.Rdata"))

################################################################################
# add leaded outcome
################################################################################

dat <- data.table(mimic)
dat_list <- split(dat, dat$subject_id)
registerDoParallel(cores = (detectCores()-1))
getDoParWorkers()

##############################################

added_outcome_list_all <- foreach(n = 1:length(dat_list)) %dopar% {
  x <- dat_list[[n]]
  added_outcome_list <- list()
  for(i in 1:nrow(x)){
    current <- x[i,]
    gaps <- c(15,20,25,30)
    future_outcomes_list <- list()
    for(j in 1:length(gaps)){
      # if the row for X-min ahead exists, add X-min ahead abpmean to current row
      future <- current$min_elapsed + gaps[j]
      if(length(which(x$min_elapsed == future)) > 0){
        future_outcomes_list[[j]] <- x[which(x$min_elapsed == future),]$abpmean
      } else {
        future_outcomes_list[[j]] <- NA
      }
    }
    future_outcomes <- unlist(future_outcomes_list)
    names(future_outcomes) <- c(paste0("Y_", gaps))
    added_outcome_list[[i]] <- data.frame(current, t(future_outcomes))
  }
  names(added_outcome_list) <- 1:nrow(x)
  added_outcome <- bind_rows(added_outcome_list)
  
  comparison <- compare::compare(added_outcome[,-c(25:28)], x, allowAll = TRUE)
  if(!(comparison$result)){
    print(paste0("Non-equal comparison for subject_id", x$subject_id))
  }
  
  if(length(unique(added_outcome$time_and_date)) != length(unique(x$time_and_date))){
    print(paste0("Length of time_and_date not equal for subject_id", x$subject_id))
  }
  
  if(length(unique(added_outcome$min_elapsed)) != length(unique(x$min_elapsed))){
    print(paste0("Length of min_elapsed not equal for subject_id", x$subject_id))
  }
  return(added_outcome[order(added_outcome$time_and_date), ])
}

table(unlist(lapply(added_outcome_list_all, ncol)))

mimic_outcome <- do.call(rbind, added_outcome_list_all)
save(mimic_outcome, file = here("Data", "mimic_outcome.Rdata"), compress = TRUE)

################################################################################
# add history of BINARY time-varying covariates 
################################################################################

load(here("Data", "mimic_outcome.Rdata"))
data <- mimic_outcome

### for time-varying binary variables, we may be interested in on/off switch:
# (1) identify min when trt was switched, and 
# (2) if it was a switch on trt (diff=1-0=1) or off trt (diff=0-1=-1)

data$amine <- as.numeric(data$amine)-1
data$sedation <- as.numeric(data$sedation)-1
data$ventilation <- as.numeric(data$ventilation)-1

dat_trt <- data %>%
  group_by(subject_id) %>%
  mutate(a_last = dplyr::lag(amine, n = 1, default = NA)) %>% 
  mutate(amine_switch1 = ifelse(amine-a_last == 1, min_elapsed, NA)) %>%
  mutate(amine_switch0 = ifelse(amine-a_last == -1, min_elapsed, NA)) %>%
  mutate(s_last = dplyr::lag(sedation, n = 1, default = NA)) %>% 
  mutate(sedation_switch1 = ifelse(sedation-s_last == 1, min_elapsed, NA)) %>%
  mutate(sedation_switch0 = ifelse(sedation-s_last == -1, min_elapsed, NA)) %>%
  mutate(v_last = dplyr::lag(ventilation, n = 1, default = NA)) %>% 
  mutate(ventilation_switch1 = ifelse(ventilation-v_last == 1, min_elapsed, NA)) %>%
  mutate(ventilation_switch0 = ifelse(ventilation-v_last == -1, min_elapsed, NA)) %>%
  select(subject_id, min_elapsed, time_and_date,
         amine, amine_switch1, amine_switch0, 
         sedation, sedation_switch1, sedation_switch0,
         ventilation, ventilation_switch1, ventilation_switch0)

dat_full <- data.table(left_join(data, dat_trt))
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
    current_history_binary[[i]] <- data.frame(current, t(amine), t(sed), t(vent))
  }
  df <- do.call(rbind, current_history_binary)
  return_df <- suppressMessages(full_join(dat_ordered[1,], df))
  return(return_df)
}

ptm <- proc.time()
N <- length(unique(dat$subject_id))
registerDoParallel(cores = (detectCores()-1)) # set up parallelization
getDoParWorkers() 
binary_history_list <- foreach(n = 1:N) %dopar% {
  id <- levels(dat$subject_id)[n]
  ind_dat <- dplyr::filter(dat, subject_id == id)
  
  # considering history as past 30 min and past 60 min
  history30 <- add_history_binary(ind_dat, 30)
  history60 <- add_history_binary(ind_dat, 60)
  
  print(paste0("Finishing subject_id level ", n))
  
  return_list <- list(history30 = history30,
                      history60 = history60)
  
  return(return_list)
}
ptm-proc.time()

binary_history60 <- do.call(rbind, lapply(binary_history_list, "[[", "history60"))
binary_history30 <- do.call(rbind, lapply(binary_history_list, "[[", "history30"))

comparison30 <- compare(binary_history30[,c(1:28)], dat, allowAll = TRUE)
comparison30$result
comparison60 <- compare(binary_history60[,c(1:28)], dat, allowAll = TRUE)
comparison60$result

save(binary_history60, file = here("Data", "binary_history60.Rdata"),
     compress = TRUE)
save(binary_history30, file = here("Data", "binary_history30.Rdata"),
     compress = TRUE)

################################################################################
# add history of CONTINUOUS time-varying covariates 
################################################################################
library(tsibble)
library(vrtest)
library(feasts)
library(forecast)

#Relevant functions
summarize_history_continuous <- function(history, current, var="abpsys", lags=10){
  
  hist <- nrow(history)
  
  full_data <- rbind.data.frame(history, current)
  full_data <- full_data[order((hist+1):1),]
  
  #Basic statistics:
  sum_hist <- unclass(summary(history[,var]))
  sum_hist_names <- c("Min","1stQ","Median","Mean","3rdQ","Max")
  names(sum_hist) <- paste0(var, "_", sum_hist_names)
  
  #Include full lags (of cardinality lags):
  lags_data <- history[1:lags,var]
  names(lags_data) <- paste0("lag_", seq(lags))
  
  #Is the history stationary (Variance Ratio Test)
  #stat <- vrtest::Auto.VR(history[,var])[[1]]
  #names(stat) <- paste0(var, "_stationary")

  #Time-series tools
  history_ts <- as_tsibble(history, index=min_elapsed)
  full_data_ts<-as_tsibble(full_data, index=min_elapsed)
  
  #Autocorrelation
  #if(nrow(full_data_ts) < 3){
  #  acfs <- rep(NA, 10)
  #}else{
  #  acfs <- full_data_ts %>% ACF(abpsys, lag_max=10)
  #  acfs <- acfs$acf
  #}
  
  #if(length(acfs)<10){
  #  acfs <- c(acfs, rep(NA, 10-length(acfs)))
  #}
  
  #names_acf <- paste0(var, "_acf_lag", seq(10))
  #names(acfs) <- names_acf
  
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
  
  return(cbind(t(sum_hist), t(lags_data), features_ts))
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
                                           var="abpsys", lags = 10)
    abpdias <- summarize_history_continuous(history=history, current=current, 
                                            var="abpdias", lags = 10)
    abpmean <- summarize_history_continuous(history=history, current=current, 
                                            var="abpmean", lags = 10)
    spo2 <- summarize_history_continuous(history=history, current=current, 
                                         var="spo2", lags = 10)
    
    current_history[[i]] <- data.frame(current, abpsys, abpdias, abpmean, spo2)
  }
  df <- do.call(rbind, current_history)
  return_df <- suppressMessages(full_join(dat_ordered[1,], df))
  return(return_df)
}

#############################################################################
#### Load data, with binary histories already added:
#############################################################################

load(file = here("Data", "binary_history30.Rdata"))
load(file = here("Data", "binary_history60.Rdata"))

#############################################################################
#### Check ACF, trend and seasonality for the whole time-series, per sample:
#############################################################################

min.true <- function(x){
  min(which(x==!TRUE))-1
}

N <- length(unique(binary_history30$subject_id))
registerDoParallel(cores = (detectCores()-1)) # set up parallelization
getDoParWorkers() 

check_ts_properties <- foreach(n = 1:N) %dopar% {
  id <- levels(binary_history30$subject_id)[n]
  ind_dat <- dplyr::filter(binary_history30, subject_id == id)
  
  ind_dat$time <- seq(1:nrow(ind_dat))
  
  ind_dat<-as_tsibble(ind_dat, index=time)
  
  #ACF indicates possibly some seasonality and trend in the series...
  #Allow for rapid changes:
  ind_dat_season <- ind_dat %>%
    model(STL(abpmean ~ season(window = 10))) %>% 
    components() #%>% autoplot()
  
  ind_dat_stl <- suppressWarnings(ind_dat %>% 
                                    features(abpsys, list(feat_stl))) 
  
  acfs_abpsys <- ind_dat %>% ACF(abpsys, lag_max=nrow(ind_dat))
  acfs_abpdias <- ind_dat %>% ACF(abpdias, lag_max=nrow(ind_dat))
  acfs_abpmean <- ind_dat %>% ACF(abpmean, lag_max=nrow(ind_dat))
  acfs_spo2 <- ind_dat %>% ACF(spo2, lag_max=nrow(ind_dat))
  acfs=cbind.data.frame(abpsys=acfs_abpsys, abpdias=acfs_abpdias,
                        abpmean=acfs_abpmean, spo2=acfs_spo2)
  
  sig_abpsys <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(ind_dat)))
  
  lags <- cbind.data.frame(lag_abpsys=min.true((acfs_abpsys$acf > sig_abpsys)),
                           lag_abpdias=min.true((acfs_abpdias$acf > sig_abpsys)),
                           lag_abpmean=min.true((acfs_abpmean$acf > sig_abpsys)),
                           lag_spo2=min.true((acfs_spo2$acf > sig_abpsys)))
  
  print(paste0("Finishing subject_id ", id, " number: ", n))
  
  return_list <- list(ind_dat_stl=ind_dat_stl,
                      ind_dat_season=ind_dat_season,
                      lags=lags,
                      acfs=acfs)
  
  return(return_list)
}  

#!!!
stls <- do.call(rbind, lapply(check_ts_properties, "[[", "ind_dat_stl"))
lags <- do.call(rbind, lapply(check_ts_properties, "[[", "lags"))        
colMeans(lags)

#############################################################################
#### Add 30 and 60 minute history, per sample:
#############################################################################

N <- length(unique(binary_history30$subject_id))
registerDoParallel(cores = (detectCores()-1)) # set up parallelization
getDoParWorkers() 

continuous_history_list <- foreach(n = 1:N) %dopar% {
  id <- levels(binary_history30$subject_id)[n]
  ind_dat_30 <- dplyr::filter(binary_history30, subject_id == id)
  ind_dat_60 <- dplyr::filter(binary_history60, subject_id == id)
  
  # considering history as past 30 min and past 60 min
  history30 <- add_history_continuous(ind_dat_30, 30)
  history60 <- add_history_continuous(ind_dat_60, 60)
  
  print(paste0("Finishing subject_id level ", n))
  
  return_list <- list(history30 = history30,
                      history60 = history60)
  
  return(return_list)
}


