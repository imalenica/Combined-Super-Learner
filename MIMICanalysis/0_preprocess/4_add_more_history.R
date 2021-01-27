library(data.table)
library(dplyr)
library(here)
library(doParallel)
library(foreach)
source(here("R", "v3", "impute_data.R"))
save_path <- "~/Box/Data_Shared/"

add_more_history <- function(d, vars, no.cores=22){

  setorder(d, subject_id, id, min)
  lags <- seq(1,180)

  # remove previous 6-number summaries & lags
  summary_names <- c("Min","1stQ","Median","Mean","3rdQ","Max")
  matches <- unique(grep(paste(c(summary_names,"lag_"), collapse="|"),
                         colnames(d), value=T))
  d <- d[, -matches, with=F]

  ##################################### lags ###################################
  # add new lags for generating new 6-number summaries
  lag_names <- lapply(lags, function(x) paste0(vars, "_lag_", x))
  for(i in 1:length(lag_names)){
    d[, (lag_names[[i]]) := shift(.SD, n=lags[i]), by=id, .SDcols=vars]
  }

  # go as far back with the lags as possible at the current time
  d_early <- d[d[, .I[1:(max(lags))], by=id]$V1, ] # select rows <= max lag
  d_late <- d[!d_early, on=.(id, min)]
  if(nrow(d_late) + nrow(d_early) != nrow(d)){
    stop("something went wrong with splitting data")
  }

  registerDoParallel(cores = no.cores) # set up parallelization
  getDoParWorkers()
  rows <- foreach(i = 1:length(lags)) %dopar% {
    if(i == 1){
      rowi <- d_early[d_early[, .I[lags[i]], by=id]$V1, ]
    } else {
      # set all lag cols to only lag to thier max possible
      row1toi <- d_early[d_early[, .I[1:lags[i]], by=id]$V1, ]
      for(j in 2:length(lag_names)){
        if(j < i){
          row1toi[,(lag_names[[j]]):=shift(.SD, n=lags[j]), by=id, .SDcols=vars]
        } else {
          row1toi[,(lag_names[[j]]):=shift(.SD, n=lags[(i-1)]), by=id, .SDcols=vars]
        }
      }
      rowi <- row1toi[row1toi[, .I[lags[i]], by=id]$V1, ]
    }
    return(rowi)
  }

  # a few checks to make sure we can build back original d
  d_early2 <- do.call(rbind, rows)
  setorder(d_early2, subject_id, id, min)
  all_names <- unlist(lag_names)
  if(!identical(d_early2[,-all_names,with=F], d_early[,-all_names,with=F])){
    stop("issue filling in NA lags")
  }
  d2 <- rbind(d_early2, d_late)
  setorder(d2, subject_id, id, min)
  if(!identical(d2[, -all_names, with=F], d[, -all_names, with=F])){
    stop("issue pairing filled NA lags with original lagged data")
  }
  d <- d2

  # remove excessively large objects
  rm(list=c("d2", "d_early", "d_late", "d_early2"))

  ##################################### AHE ####################################
  # add whether or not they had AHE, TODO (for binary Y): t since last episode
  Y_idx <- grep("abpmean", vars)
  Y_var <- vars[Y_idx]
  d[, "AHE0" := ifelse(.SD <= 65, 1, 0), by=id, .SDcols=Y_var]
  d[, "AHE1" := shift(AHE0, n=1, fill=0), by=id]
  d[, "AHE2" := shift(AHE0, n=2, fill=0), by=id]
  d[, "AHE3" := shift(AHE0, n=3, fill=0), by=id]
  d[, "AHE4" := shift(AHE0, n=4, fill=0), by=id]
  d[, "AHE5" := shift(AHE0, n=5, fill=0), by=id]
  d[, "AHEsum" := AHE0+as.numeric(AHE1)+as.numeric(AHE2)+as.numeric(AHE3)+
      as.numeric(AHE4)+as.numeric(AHE5)]
  d[, "AHE" := ifelse(AHEsum >= 5, 1, 0), by=id]
  d <- d[,-c("AHE0", "AHE1", "AHE2", "AHE3", "AHE4", "AHE5", "AHEsum"), with=F]

  ################################### summaries ################################
  # add hourly 6-number summaries (3 hours for Y, 1 hour for other vitals)
  registerDoParallel(cores = no.cores) # set up parallelization
  getDoParWorkers()
  summaries_list <- foreach(i = 1:length(vars)) %dopar% {
    var <- vars[[i]]
    lags <- list(1:30, 31:60)
    if(grepl("abpmean", var)){
      lags <- c(lags, list(61:90, 91:120, 121:150, 151:180))
    }
    return(data.table(do.call(cbind, lapply(seq_along(lags), function(i){
      lag_names <- paste0(var, "_lag_", lags[[i]])
      d_lag <- d[-1, lag_names, with=F]
      lag_summary <- data.table(do.call(rbind, apply(d_lag, 1, summary, na.rm=T)))
      if(any(colnames(lag_summary) == "NA's")){
        NA_idx <- which(colnames(lag_summary) == "NA's")
        lag_summary <- lag_summary[, -NA_idx, with=F]
      }
      min_lag <- min(lags[[i]])
      max_lag <- max(lags[[i]])
      name <- paste0(var, "_", summary_names, "_", min_lag, "_", max_lag)
      names(lag_summary) <- name
      row1 <- data.table(t(rep(NA, 6)))
      names(row1) <- name
      return(data.table(rbind(row1, lag_summary)))
    }))))
  }
  summaries <- data.table(do.call(cbind, summaries_list))

  # remove lags we no longer want & bind data with the summaries
  vars_lag60 <- vars[c(grep("abpmean", vars), grep("hr", vars))]
  vars_lag30 <- vars[-c(grep("abpmean", vars), grep("hr", vars))]
  to_rm <- c(
    unlist(lapply(vars_lag30, function(x) paste0(x, "_lag_", 31:180))),
    unlist(lapply(vars_lag60, function(x) paste0(x, "_lag_", 61:180)))
  )
  d <- data.table(cbind(d[, -to_rm, with=F], summaries))
  setorder(d, subject_id, id, min)
  return(d)
}

add_lost_rows <- function(df_mimic, dt_split){
  # add prev sliced out rows to dt_split, ie historical/individual data
  df_first <- df_mimic %>% group_by(id) %>% slice(1:10)
  df_last <- df_mimic %>% group_by(id) %>% slice((n()-45):n())
  dt_sliced <- rbind(data.table(df_first), data.table(df_last))
  dt <- rbind(dt_sliced, dt_split, fill=T)
  setorder(dt, subject_id, id, time_and_date)
  return(dt)
}

mini_process <- function(d){
  # remove duplicates
  if(any(duplicated(d[,c("id", "min"),with=F]))){
    dups <- which(duplicated(d))
    if(sum(duplicated(d[,c("id", "min"),with=F])) == length(dups)){
      print(paste0("Removing ", length(dups), " duplicates"))
      d <- d[-dups,]
    } else {
      stop("Check duplicates")
    }
  } else {
    print("No duplicates detected")
  }

  # remove treatment cols that we don't need
  rm <- c("amine_switch1", "amine_switch0", "sedation_switch1",
          "sedation_switch0", "ventilation_switch1", "ventilation_switch0",
          "amine_history_switch1_time", "amine_history_switch0_time",
          "sedation_history_switch1_time", "sedation_history_switch0_time",
          "ventilation_history_switch1_time", "ventilation_history_switch0_time")
  matches <- unique(grep(paste(rm, collapse="|"), colnames(d), value=T))
  d <- d[, -matches, with=F]

  return(d)
}

remove_low_support_ids <- function(d, max_lag=180){
  # remove subjects with too low observations
  id_length <- d[, .N, by=id]
  bad_ids <- id_length[N <= max_lag,]
  bad_ids <- bad_ids[["id"]]
  if(length(bad_ids) > 0){
    print(paste0("Removing ", length(bad_ids)," ids, not enough observations"))
    print(as.character(bad_ids))
    d <- d[!id %in% bad_ids, ]
  }
  return(d)
}

slice_rows <- function(dt){
  df <- data.frame(dt)
  data <- df %>% dplyr::group_by(id) %>% slice(11:n())
  return(data.table(data))
}

process_cols <- function(d){
  # remove columns with many NA
  manyNA <- c("trend_strength", "stl_e_acf1", "stl_e_acf10", "spectral_entropy")
  matches <- unique(grep(paste(manyNA, collapse="|"), colnames(d), value=T))
  d[, -matches, with=F]
}

vars <- c("abpsys", "abpdias", "abpmean", "spo2", "hr")
vars_mean <- paste0(vars, "_lag5_mean")
vars_median <- paste0(vars, "_lag5_median")

################################ mean smoothing ################################

# smoothed data after split
load(here("Data", "individual60_smooth_mean.Rdata"))
individual <- remove_low_support_ids(individual)
# [1] "Removing 57 ids, not enough observations"
# [1] "10670_1" "11420"   "11990"   "12004"   "1213"    "13990_2" "14081"   "14131"   "15245"
# [10] "15314"   "16434"   "17403"   "17956"   "18118"   "18645"   "19042"   "19166"   "19448"
# [19] "19841"   "19951"   "20154"   "21201"   "22156"   "24748"   "25075_1" "2563"    "25842"
# [28] "26261"   "27006"   "27075"   "28609"   "28948"   "29306"   "30311_1" "30544_2" "30710"
# [37] "30734"   "3076"    "31922"   "32332"   "3418_1"  "3460"    "3836"    "4117"    "4498_1"
# [46] "5563"    "573"     "5924"    "6026"    "6077"    "7105"    "7736"    "9029"    "9066_1"
# [55] "9066_2"  "9134"    "9841_1"
ind_ids <- unique(individual$id)
length(ind_ids) # 402

load(here("Data", "historical60_smooth_mean.Rdata"))
historical <- remove_low_support_ids(historical)
hist_ids <- unique(historical$id)
# "Removing 9 ids, not enough observations"
# [1] "12396_2" "13288_1" "17253"   "29505"   "30244"   "32504_1" "5963_2"  "7648_3"  "8981_2"
length(hist_ids) # 368

# BIG smoothed data before split
load(here("Data", "mimic60_smooth_mean.Rdata"))
df_mimic_mean_ind <- data.frame(mimic60_smooth_mean[id %in% ind_ids,])
mimic_mean_hist <- mimic60_smooth_mean[id %in% hist_ids,]
rm(mimic60_smooth_mean)

# add more history to data for individualized fits
individual <- add_lost_rows(df_mimic_mean_ind, individual)
individual <- mini_process(individual)
rm(df_mimic_mean_ind)
individual <- add_more_history(individual, vars_mean)
save(individual, file=here(save_path, "individual_mean_full.Rdata"), compress=T)

# add more history to data for historical fit
historical <- mini_process(mimic_mean_hist)
rm(mimic_mean_hist)
historical <- add_more_history(historical, vars_mean)
save(historical, file=here(save_path, "historical_mean_full.Rdata"), compress=T)

################################ median smoothing ##############################

# smoothed data after split
load(here("Data", "individual60_smooth_median.Rdata"))
individual <- remove_low_support_ids(individual)
# [1] "Removing 57 ids, not enough observations"
# [1] "10670_1" "11420"   "11990"   "12004"   "1213"    "13990_2" "14081"
# [8] "14131"   "15245"   "15314"   "16434"   "17403"   "17956"   "18118"
# [15] "18645"   "19042"   "19166"   "19448"   "19841"   "19951"   "20154"
# [22] "21201"   "22156"   "24748"   "25075_1" "2563"    "25842"   "26261"
# [29] "27006"   "27075"   "28609"   "28948"   "29306"   "30311_1" "30544_2"
# [36] "30710"   "30734"   "3076"    "31922"   "32332"   "3418_1"  "3460"
# [43] "3836"    "4117"    "4498_1"  "5563"    "573"     "5924"    "6026"
# [50] "6077"    "7105"    "7736"    "9029"    "9066_1"  "9066_2"  "9134"
# [57] "9841_1"
ind_ids <- unique(individual$id)
length(ind_ids) #402

load(here("Data", "historical60_smooth_median.Rdata"))
historical <- remove_low_support_ids(historical)
# [1] "Removing 9 ids, not enough observations"
# [1] "12396_2" "13288_1" "17253"   "29505"   "30244"   "32504_1" "5963_2"
# [8] "7648_3"  "8981_2"
hist_ids <- unique(historical$id)
length(hist_ids) #368

# BIG smoothed data before split
load(here("Data", "mimic60_smooth_median.Rdata"))
df_mimic_median_ind <- data.frame(mimic60_smooth_median[id %in% ind_ids,])
mimic_median_hist <- mimic60_smooth_median[id %in% hist_ids,]
rm(mimic60_smooth_median)

# add more history to data for individualized fits
individual <- add_lost_rows(df_mimic_median_ind, individual)
individual <- mini_process(individual)
rm(df_mimic_median_ind)
individual <- add_more_history(individual, vars_median)
save(individual, file=here(save_path, "individual_median_full.Rdata"), compress=T)

# add more history to data for historical fit
historical <- mini_process(mimic_median_hist)
rm(mimic_median_hist)
historical <- add_more_history(historical, vars_median)
save(historical, file=here(save_path, "historical_median_full.Rdata"), compress=T)

##################################### slice ####################################
load(here(save_path, "historical_mean_full.Rdata"))
historical <- slice_rows(historical)
historical <- process_cols(historical)
outcomes <- colnames(historical)[grepl("Y", colnames(historical))]
predictors <- historical[, -outcomes, with=F]
predictors_imputed <- impute_data(predictors, strata=c("admission_type_descr","sex"))
historical <- data.table(predictors_imputed, historical[, outcomes, with=F])
historical <- historical[, -"rank_icu", with=F]
save(historical, file=here(save_path, "historical_mean.Rdata"), compress=T)
rm(historical)

load(here(save_path, "individual_mean_full.Rdata"))
individual <- slice_rows(individual)
individual <- process_cols(individual)
individual <- individual[, -"rank_icu", with=F]
save(individual, file=here(save_path, "individual_mean.Rdata"), compress=T)
rm(individual)

load(here(save_path, "historical_median_full.Rdata"))
historical <- slice_rows(historical)
historical <- process_cols(historical)
outcomes <- colnames(historical)[grepl("Y", colnames(historical))]
predictors <- historical[, -outcomes, with=F]
predictors_imputed <- impute_data(predictors, strata=c("admission_type_descr","sex"))
historical <- data.table(predictors_imputed, historical[, outcomes, with=F])
historical <- historical[, -"rank_icu", with=F]
save(historical, file=here(save_path, "historical_median.Rdata"), compress=T)
rm(historical)

load(here(save_path, "individual_median_full.Rdata"))
individual <- slice_rows(individual)
individual <- process_cols(individual)
individual <- individual[, -"rank_icu", with=F]
save(individual, file=here(save_path, "individual_median.Rdata"), compress=T)
rm(individual)
