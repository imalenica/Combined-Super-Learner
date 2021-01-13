library(tidyverse)
library(data.table)
library(here)
library(lubridate)
library(doParallel)
library(compare)

library(tsibble)
library(vrtest)
library(feasts)
library(forecast)
library(ForeCA)


##### time-varying binary trt summaries
### for time-varying binary variables, we may be interested in on/off switch:
# (1) identify min when trt was switched, and 
# (2) if it was a switch on trt (diff=1-0=1) or off trt (diff=0-1=-1)

init_history_binary <- function(data){
  data$amine <- as.numeric(data$amine)-1
  data$sedation <- as.numeric(data$sedation)-1
  data$ventilation <- as.numeric(data$ventilation)-1
  
  dat_trt <- data %>%
    group_by(id) %>%
    mutate(a_last = dplyr::lag(amine, n = 1, default = NA)) %>% 
    mutate(amine_switch1 = ifelse(amine-a_last == 1, min, NA)) %>%
    mutate(amine_switch0 = ifelse(amine-a_last == -1, min, NA)) %>%
    mutate(s_last = dplyr::lag(sedation, n = 1, default = NA)) %>% 
    mutate(sedation_switch1 = ifelse(sedation-s_last == 1, min, NA)) %>%
    mutate(sedation_switch0 = ifelse(sedation-s_last == -1, min, NA)) %>%
    mutate(v_last = dplyr::lag(ventilation, n = 1, default = NA)) %>% 
    mutate(ventilation_switch1 = ifelse(ventilation-v_last == 1, min, NA)) %>%
    mutate(ventilation_switch0 = ifelse(ventilation-v_last == -1, min, NA)) %>%
    select(id, min, time_and_date,
           amine, amine_switch1, amine_switch0, 
           sedation, sedation_switch1, sedation_switch0,
           ventilation, ventilation_switch1, ventilation_switch0)
  
  dat_full <- data.table(left_join(data, dat_trt))
  dat <- dat_full[order(dat_full$id, dat_full$time_and_date), ]
  return(dat)
}

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

add_all_history_binary <- function(dat, no.cores){
  
  N <- length(unique(dat$id))
  
  registerDoParallel(cores = no.cores) # set up parallelization
  getDoParWorkers() 
  binary_history_list <- foreach(n = 1:N) %dopar% {
    i <- levels(dat$id)[n]
    ind_dat <- dplyr::filter(dat, id == i)
    
    # considering history as past 30 min and past 60 min
    history30 <- add_history_binary(ind_dat, 30)
    history60 <- add_history_binary(ind_dat, 60)
    
    print(paste0("Finishing id level ", n))
    
    return_list <- list(history30 = history30,
                        history60 = history60)
    
    return(return_list)
  }
  
  binary_history60 <- do.call(rbind, lapply(binary_history_list, "[[", "history60"))
  binary_history30 <- do.call(rbind, lapply(binary_history_list, "[[", "history30"))
  return(list(binary_history30 = binary_history30, 
              binary_history60 = binary_history60))
}

################################################################################
# add history of BINARY time-varying covariates 
################################################################################
load(here::here("Data","mimic.Rdata"))
dat <- init_binary_history(mimic)
ptm <- proc.time()
binary_history_list <- add_all_history_binary(dat, 20)
proc.time() - ptm
# user    system   elapsed
# 76124.020  1088.452  5495.081
binary_history30 <- binary_history_list$binary_history30
binary_history60 <- binary_history_list$binary_history60
comparison30 <- compare(binary_history30[,c(1:38)], dat, allowAll = TRUE)
comparison30$result
comparison60 <- compare(binary_history60[,c(1:38)], dat, allowAll = TRUE)
comparison60$result
irrelevant <- c("amine_switch1", "amine_switch0", "sedation_switch1", 
                "sedation_switch0", "ventilation_switch1", "ventilation_switch0")
binary_history30 <- data.table(binary_history30)[,-irrelevant, with=FALSE]
binary_history60 <- data.table(binary_history60)[,-irrelevant, with=FALSE]
save(binary_history60, file = here::here("Data", "binary_history60.Rdata"),
     compress = TRUE)
save(binary_history30, file = here::here("Data", "binary_history30.Rdata"),
     compress = TRUE)

load(here::here("Data","mimic_smooth_mean.Rdata"))
dat_smooth_mean <- init_history_binary(mimic_smooth_mean)
ptm <- proc.time()
binary_history_list_smooth_mean <- add_all_history_binary(dat_smooth_mean, 20)
proc.time() - ptm
# user    system   elapsed
# 78022.315  1541.958  6105.254
binary_history30_smooth_mean <- binary_history_list_smooth_mean$binary_history30
binary_history60_smooth_mean <- binary_history_list_smooth_mean$binary_history60
comparison30 <- compare(binary_history30_smooth_mean[,c(1:38)], dat_smooth_mean, 
                        allowAll = TRUE)
comparison30$result
comparison60 <- compare(binary_history60_smooth_mean[,c(1:38)], dat_smooth_mean, 
                        allowAll = TRUE)
comparison60$result
binary_history60_smooth_mean <- data.table(binary_history30_smooth_mean)[,-irrelevant, with=FALSE]
binary_history30_smooth_mean <- data.table(binary_history60_smooth_mean)[,-irrelevant, with=FALSE]
save(binary_history30_smooth_mean, compress = TRUE,
     file = here::here("Data", "binary_history30_smooth_mean.Rdata"))
save(binary_history60_smooth_mean, compress = TRUE,
     file = here::here("Data", "binary_history60_smooth_mean.Rdata"))


load(here::here("Data","mimic_smooth_median.Rdata"))
dat_smooth_median <- init_history_binary(mimic_smooth_median)
ptm <- proc.time()
binary_history_list_smooth_median <- add_all_history_binary(dat_smooth_median, 20)
proc.time()-ptm

binary_history30_smooth_median <- binary_history_list_smooth_median$binary_history30
binary_history60_smooth_median <- binary_history_list_smooth_median$binary_history60
comparison30 <- compare(binary_history30_smooth_median[,c(1:38)], dat_smooth_median, 
                        allowAll = TRUE)
comparison30$result
comparison60 <- compare(binary_history60_smooth_median[,c(1:38)], dat_smooth_median, 
                        allowAll = TRUE)
comparison60$result
binary_history60_smooth_median <- data.table(binary_history30_smooth_median)[,-irrelevant, with=FALSE]
binary_history30_smooth_median <- data.table(binary_history60_smooth_median)[,-irrelevant, with=FALSE]
save(binary_history30_smooth_median, compress = TRUE,
     file = here::here("Data", "binary_history30_smooth_median.Rdata"))
save(binary_history60_smooth_median, compress = TRUE,
     file = here::here("Data", "binary_history60_smooth_median.Rdata"))

################################################################################
# add history of CONTINUOUS time-varying covariates 
################################################################################

#Relevant functions
summarize_history_continuous <- function(history, current, vars){
  
  #Basic statistics:
  summ_names <- c("Min","1stQ","Median","Mean","3rdQ","Max")
  summ_hist <- c(unname(summary(unlist(unname(dplyr::select(history, vars[1]))))),
                 unname(summary(unlist(unname(dplyr::select(history, vars[2]))))),
                 unname(summary(unlist(unname(dplyr::select(history, vars[3]))))),
                 unname(summary(unlist(unname(dplyr::select(history, vars[4]))))),
                 unname(summary(unlist(unname(dplyr::select(history, vars[5]))))))
  names(summ_hist) <- unlist(lapply(vars, function(x) paste0(x, "_", summ_names)))

  #Time-series tools
  history_ts <- as_tsibble(history, index = min)

  #Extract time-series features:
  #1) Features related to STL decomposition of the series (strength of trend, seasonality)
  #2) Spectral entropy of a time-series
  #3) Number of flat points
  #4) Hurst coefficient (level of fractional differencing of a series)
  
  names_features <- c("trend_strength","spikiness","linearity","curvature",
                      "stl_e_acf1", "stl_e_acf10","spectral_entropy",
                      "n_flat_spots","coef_hurst")
  all_names_features <- unlist(lapply(vars, function(x) paste0(x, "_", names_features)))
  
  if(nrow(history_ts) == 1){
    features_ts <- rep(NA, 45)
    features_ts <- data.frame(t(features_ts))
    names(features_ts) <- all_names_features
  }else{
    feature_idx <- which(colnames(history_ts) == vars[1])
    history_ts$feature <- unname(unlist(history_ts[, feature_idx]))
    feat_abpsys <- suppressWarnings(history_ts %>%                             
      features(feature, list(feat_stl, feat_spectral, n_flat_spots, coef_hurst)))
    names(feat_abpsys) <- paste(vars[1], names(feat_abpsys), sep = "_")
    
    feature_idx <- which(colnames(history_ts) == vars[2])
    history_ts$feature <- unname(unlist(history_ts[, feature_idx]))
    feat_abpdias <- suppressWarnings(history_ts %>%                             
      features(feature, list(feat_stl, feat_spectral, n_flat_spots, coef_hurst)))
    names(feat_abpdias) <- paste(vars[2], names(feat_abpdias), sep = "_")
    
    feature_idx <- which(colnames(history_ts) == vars[3])
    history_ts$feature <- unname(unlist(history_ts[, feature_idx]))
    feat_abpmean <- suppressWarnings(history_ts %>%                             
      features(feature, list(feat_stl, feat_spectral, n_flat_spots, coef_hurst)))
    names(feat_abpmean) <- paste(vars[3], names(feat_abpmean), sep = "_")
    
    feature_idx <- which(colnames(history_ts) == vars[4])
    history_ts$feature <- unname(unlist(history_ts[, feature_idx]))
    feat_spo2 <- suppressWarnings(history_ts %>%                             
      features(feature, list(feat_stl, feat_spectral, n_flat_spots, coef_hurst)))
    names(feat_spo2) <- paste(vars[4], names(feat_spo2), sep = "_")
    
    feature_idx <- which(colnames(history_ts) == vars[5])
    history_ts$feature <- unname(unlist(history_ts[, feature_idx]))
    feat_hr <- suppressWarnings(history_ts %>%                             
      features(feature, list(feat_stl, feat_spectral, n_flat_spots, coef_hurst)))
    names(feat_hr) <- paste(vars[5], names(feat_hr), sep = "_")
                              
    features_ts <- cbind(feat_abpsys, feat_abpdias, feat_abpmean, feat_spo2, 
                         feat_hr)
    
    if(length(features_ts) != 45){
      Missing <- setdiff(all_names_features, names(features_ts))
      features_ts[Missing] <- NA
      features_ts <- features_ts[, all_names_features]
    }
  }

  return(data.frame(t(summ_hist), features_ts))
}

add_history_continuous <- function(ind_data, min_history, vars){
  dat_ordered <- ind_data[order(ind_data$time_and_date), ]
  
  current_history <- list()
  for(i in 2:nrow(dat_ordered)){
    history_full <- data.frame(dat_ordered[c(1:(i-1)),])
    # subset to history if i > min_history
    if(i <= min_history){
      history <- history_full
    } else {
      history <- history_full[c(
        (nrow(history_full)-min_history+1):(nrow(history_full))),]
    }
    current <- dat_ordered[i,]
    history <- summarize_history_continuous(history, current, vars)
    current_history[[i]] <- data.frame(current, history)
  }
  df <- do.call(rbind, current_history)
  return_df <- suppressMessages(full_join(dat_ordered[1,], df))
  return(return_df)
}

dt_reorder <- function(dt){
  dt$subject_id <- as.numeric(as.character(dt$subject_id))
  dt$icustay_id <- as.numeric(as.character(dt$icustay_id))
  dt$min_elapsed <- as.numeric(as.character(dt$min_elapsed))
  dt <- data.table(dt)
  dt_ord <- setorder(dt, "subject_id", "icustay_id", "min_elapsed")
  dt_ord$subject_id <- as.factor(as.character(dt_ord$subject_id))
  dt_ord$icustay_id <- as.factor(as.character(dt_ord$icustay_id))
  dt_ord <- data.table(dt_ord)
  return(dt_ord)
}

add_all_history_continuous <- function(binary_history30, binary_history60, 
                                       no.cores, vars){
  N <- length(unique(binary_history30$id))
  registerDoParallel(cores = no.cores) # set up parallelization
  getDoParWorkers() 
  
  continuous_history_list <- foreach(n = 1:N) %dopar% {
    i <- levels(binary_history30$id)[n]
    ind_dat_30 <- dplyr::filter(binary_history30, id == i)
    ind_dat_60 <- dplyr::filter(binary_history60, id == i)
    if( any(duplicated(dplyr::select(ind_dat_30,min))) | 
        any(duplicated(dplyr::select(ind_dat_30,min))) ) {
      return(NULL)
    } else {
      # considering history as past 30 min and past 60 min
      add_hist <- function (d,t,v) {
        out <- tryCatch(add_history_continuous(d,t,v), error = function(e) NULL)
        return(out)
      }
      history30 <- add_hist(ind_dat_30, 30, vars)
      history60 <- add_hist(ind_dat_60, 60, vars)
      
      print(paste0("Finishing id level ", n))
      
      return_list <- list(history30 = history30,
                          history60 = history60)
      
      return(return_list)
    }
  }
  
  history_list <- continuous_history_list[!sapply(continuous_history_list, is.null)] 
  all_history30 <- data.table(do.call(rbind, lapply(history_list, '[[', 'history30')))
  all_history60 <- data.table(do.call(rbind, lapply(history_list, '[[', 'history60')))

  lag_names <- lapply(seq(10), function(x) paste0(vars, "_lag_", x))
  add_lags <- function(dt){
    dt[, (lag_names[[1]]) := shift(.SD, n=1), by=id, .SDcols=vars]
    dt[, (lag_names[[2]]) := shift(.SD, n=2), by=id, .SDcols=vars]
    dt[, (lag_names[[3]]) := shift(.SD, n=3), by=id, .SDcols=vars]
    dt[, (lag_names[[4]]) := shift(.SD, n=4), by=id, .SDcols=vars]
    dt[, (lag_names[[5]]) := shift(.SD, n=5), by=id, .SDcols=vars]
    dt[, (lag_names[[6]]) := shift(.SD, n=6), by=id, .SDcols=vars]
    dt[, (lag_names[[7]]) := shift(.SD, n=7), by=id, .SDcols=vars]
    dt[, (lag_names[[8]]) := shift(.SD, n=8), by=id, .SDcols=vars]
    dt[, (lag_names[[9]]) := shift(.SD, n=9), by=id, .SDcols=vars]
    dt[, (lag_names[[10]]) := shift(.SD, n=10), by=id, .SDcols=vars]
    return(dt)
  }
  history30 <- add_lags(all_history30)
  history60 <- add_lags(all_history60)
  history30 <- dt_reorder(history30)
  history60 <- dt_reorder(history60)
  
  return(list(history30 = history30, history60 = history60))
}


#############################################################################
#### Load data, with binary histories already added:
#############################################################################

load(file = here::here("Data", "binary_history30.Rdata"))
length(unique(binary_history30$id)) # 1162
load(file = here::here("Data", "binary_history60.Rdata"))
length(unique(binary_history60$id)) # 1162
forlag <- c("abpsys", "abpdias", "abpmean", "spo2", "hr")
ptm <- proc.time()
history_list <- add_all_history_continuous(binary_history30 = binary_history30, 
                                           binary_history60 = binary_history60, 
                                           no.cores = 20, vars = forlag)
# user     system    elapsed
# 934860.202   5937.411  60333.823
proc.time()-ptm

mimic30 <- history_list$history30
mimic60 <- history_list$history60

retained_ids <- unique(mimic30$id)
bin_hist30 <- binary_history30[id %in% retained_ids,]
retained_ids <- unique(mimic60$id)
bin_hist60 <- binary_history60[id %in% retained_ids,]
comparison30 <- compare(mimic30[,c(1:64)], bin_hist30, allowAll = TRUE)
comparison30$result
comparison60 <- compare(mimic60[,c(1:64)], bin_hist60, allowAll = TRUE)
comparison60$result
save(mimic30, file = here::here("Data", "mimic30.Rdata"), compress = T)
save(mimic60, file = here::here("Data", "mimic60.Rdata"), compress = T)

load(file = here::here("Data", "binary_history30_smooth_mean.Rdata"))
load(file = here::here("Data", "binary_history60_smooth_mean.Rdata"))
forlag <- c("abpsys_lag5_mean", "abpdias_lag5_mean", "abpmean_lag5_mean", "spo2_lag5_mean", 
            "hr_lag5_mean")
ptm <- proc.time()
history_list_smooth_mean <- add_all_history_continuous(
  binary_history30 = binary_history30_smooth_mean, 
  binary_history60 = binary_history60_smooth_mean, no.cores = 20, vars = forlag)
proc.time()-ptm
# user    system   elapsed
# 981022.85  13017.40  63186.77
mimic30_smooth_mean <- history_list_smooth_mean$history30
mimic60_smooth_mean <- history_list_smooth_mean$history60
retained_ids <- unique(mimic30_smooth_mean$id)
binary_history30_smooth_mean <- binary_history30_smooth_mean[id %in% retained_ids,]
retained_ids <- unique(mimic60_smooth_mean$id)
binary_history60_smooth_mean <- binary_history60_smooth_mean[id %in% retained_ids,]
comparison30 <- compare(mimic30_smooth_mean[,c(1:64)], binary_history30_smooth_mean, allowAll = TRUE)
comparison30$result
comparison60 <- compare(mimic60_smooth_mean[,c(1:64)], binary_history60_smooth_mean, allowAll = TRUE)
comparison60$result
save(mimic30_smooth_mean, file = here::here("Data", "mimic30_smooth_mean.Rdata"), 
     compress = TRUE)
save(mimic60_smooth_mean, file = here::here("Data", "mimic60_smooth_mean.Rdata"), 
     compress = TRUE)

load(file = here::here("Data", "binary_history30_smooth_median.Rdata"))
load(file = here::here("Data", "binary_history60_smooth_median.Rdata"))
forlag <- c("abpsys_lag5_median", "abpdias_lag5_median", "abpmean_lag5_median", 
            "spo2_lag5_median", "hr_lag5_median")
ptm <- proc.time()
history_list_smooth_median <- add_all_history_continuous(
  binary_history30 = binary_history30_smooth_median, 
  binary_history60 = binary_history60_smooth_median, no.cores = 20, vars = forlag)
proc.time()-ptm
# user    system   elapsed
# 995518.33  15965.76  64338.66
# 64338/60/60
# [1] 17.87167
mimic30_smooth_median <- history_list_smooth_median$history30
mimic60_smooth_median <- history_list_smooth_median$history60
retained_ids <- unique(mimic30_smooth_median$id)
binary_history30_smooth_median <- data.table(binary_history30_smooth_median)[id %in% retained_ids,]
retained_ids <- unique(mimic60_smooth_median$id)
binary_history60_smooth_median <- data.table(binary_history60_smooth_median)[id %in% retained_ids,]
comparison30 <- compare(mimic30_smooth_median[,c(1:64)], 
                        binary_history30_smooth_median, allowAll = TRUE)
comparison30$result
comparison60 <- compare(mimic60_smooth_median[,c(1:64)], 
                        binary_history60_smooth_median, allowAll = TRUE)
comparison60$result
save(mimic30_smooth_median, file = here::here("Data", "mimic30_smooth_median.Rdata"), 
     compress = TRUE)
save(mimic60_smooth_median, file = here::here("Data", "mimic60_smooth_median.Rdata"), 
     compress = TRUE)
#############################################################################
#### Check ACF, trend and seasonality for the whole time-series, per sample:
#############################################################################

# min.true <- function(x){
#   min(which(x==!TRUE))-1
# }

# check_ts_properties <- foreach(n = 1:N) %dopar% {
#   id <- levels(binary_history30$subject_id)[n]
#   ind_dat <- dplyr::filter(binary_history30, subject_id == id)
#   
#   ind_dat$time <- seq(1:nrow(ind_dat))
#   
#   ind_dat<-as_tsibble(ind_dat, index=time)
#   
#   #ACF indicates possibly some seasonality and trend in the series...
#   #Allow for rapid changes:
#   ind_dat_season <- ind_dat %>%
#     model(STL(abpmean ~ season(window = 10))) %>% 
#     components() #%>% autoplot()
#   
#   ind_dat_stl <- suppressWarnings(ind_dat %>% 
#                                     features(abpsys, list(feat_stl))) 
#   
#   acfs_abpsys <- ind_dat %>% ACF(abpsys, lag_max=nrow(ind_dat))
#   acfs_abpdias <- ind_dat %>% ACF(abpdias, lag_max=nrow(ind_dat))
#   acfs_abpmean <- ind_dat %>% ACF(abpmean, lag_max=nrow(ind_dat))
#   acfs_spo2 <- ind_dat %>% ACF(spo2, lag_max=nrow(ind_dat))
#   acfs=cbind.data.frame(abpsys=acfs_abpsys, abpdias=acfs_abpdias,
#                         abpmean=acfs_abpmean, spo2=acfs_spo2)
#   
#   sig_abpsys <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(ind_dat)))
#   
#   lags <- cbind.data.frame(lag_abpsys=min.true((acfs_abpsys$acf > sig_abpsys)),
#                            lag_abpdias=min.true((acfs_abpdias$acf > sig_abpsys)),
#                            lag_abpmean=min.true((acfs_abpmean$acf > sig_abpsys)),
#                            lag_spo2=min.true((acfs_spo2$acf > sig_abpsys)))
#   
#   print(paste0("Finishing id ", id, " number: ", n))
#   
#   return_list <- list(ind_dat_stl=ind_dat_stl,
#                       ind_dat_season=ind_dat_season,
#                       lags=lags,
#                       acfs=acfs)
#   
#   return(return_list)
# }  
# 
# #!!!
# stls <- do.call(rbind, lapply(check_ts_properties, "[[", "ind_dat_stl"))
# lags <- do.call(rbind, lapply(check_ts_properties, "[[", "lags"))        
# colMeans(lags)