library(data.table)
library(dplyr)
library(here)

add_more_history <- function(d, vars, lags = seq(11,180)){
  
  # remove previous 6-number summaries
  summary_names <- c("Min","1stQ","Median","Mean","3rdQ","Max")
  matches <- unique(grep(paste(summary_names, collapse="|"), colnames(d), value=T))
  d <- d[, -matches, with=F]
  
  # add new lags for generating new 6-number summaries
  lag_names <- lapply(lags, function(x) paste0(vars, "_lag_", x))
  for(i in 1:length(lag_names)){
    d[, (lag_names[[i]]) := shift(.SD, n=lags[i]), by=id, .SDcols=vars]
  }
  
  # add hourly 6-number summaries
  lag_seq_list <- list(1:60, 1:120, 1:180)
  summaries <- data.table(do.call(cbind, lapply(vars, function(var){
    data.table(do.call(cbind, lapply(seq_along(lag_seq_list), function(i){
      lag_seq <- lag_seq_list[[i]]
      lag_names <- paste0(var, "_lag_", lag_seq)
      d_lag <- d[-1, lag_names, with=F]
      lag_summary <- data.table(do.call(rbind, apply(d_lag, 1, summary, na.rm=T)))
      if(any(colnames(lag_summary) == "NA's")){
        NA_idx <- which(colnames(lag_summary) == "NA's")
        lag_summary <- lag_summary[, -NA_idx, with=F]
      }
      names_lag_summary <- paste0(var, "_", summary_names, "_prev", i, "hrs")
      names(lag_summary) <- names_lag_summary
      row1 <- data.table(t(rep(NA, 6)))
      names(row1) <- names_lag_summary
      return(data.table(rbind(row1, lag_summary)))
    })))
  })))
  
  # remove lags we no longer want & bind with the summaries
  Y_idx <- grep("abpmean", vars)
  to_rm <- unlist(lapply(vars[-Y_idx], function(x) paste0(x, "_lag_", 11:180)))
  to_rm <- c(to_rm, paste0(vars[Y_idx], "_lag_", 31:180))
  return(data.table(cbind(d[, -to_rm, with=F], summaries)))
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

slice_rows <- function(dt){
  df <- data.frame(dt)
  data <- df %>% dplyr::group_by(id) %>% slice(31:(n()-46))
  return(data.table(data))
}

vars <- c("abpsys", "abpdias", "abpmean", "spo2", "hr")
vars_mean <- paste0(vars, "_lag5_mean")
vars_median <- paste0(vars, "_lag5_median")

################################ mean smoothing ################################

# smoothed data after split
load("~/Downloads/individual60_smooth_mean.Rdata")
ind_ids <- unique(individual$id)
load("~/Downloads/historical60_smooth_mean.Rdata")
hist_ids <- unique(historical$id)

# BIG smoothed data before split
load(here("Data", "mimic60_smooth_mean.Rdata"))
df_mimic_mean_ind <- data.frame(mimic60_smooth_mean[id %in% ind_ids,])
df_mimic_mean_hist <- data.frame(mimic60_smooth_mean[id %in% hist_ids,])
rm(mimic60_smooth_mean)

# add more history to data for individualized fits
start <- proc.time()
individual <- add_lost_rows(df_mimic_mean_ind, individual)
rm(df_mimic_mean_ind)
individual <- add_more_history(individual, vars_mean)
proc.time()-start
save(individual, file="~/Downloads/individual_mean_full.Rdata", compress=T)
individual <- slice_rows(individual)
save(individual, file="~/Downloads/individual_mean.Rdata", compress=T)
rm(individual)

# add more history to data for historical fit
start <- proc.time()
historical <- add_lost_rows(df_mimic_mean_hist, historical)
rm(df_mimic_mean_hist)
historical <- add_more_history(historical, vars_mean)
proc.time()-start
save(historical, file="~/Downloads/historical_mean_full.Rdata", compress=T)
historical <- slice_rows(historical)
save(historical, file="~/Downloads/historical_mean.Rdata", compress=T)
rm(historical)

################################ median smoothing ##############################

# smoothed data after split
load("~/Downloads/individual60_smooth_median.Rdata")
ind_ids <- unique(individual$id)
load("~/Downloads/historical60_smooth_median.Rdata")
hist_ids <- unique(historical$id)

# BIG smoothed data before split
load(here("Data", "mimic60_smooth_median.Rdata"))
df_mimic_median_ind <- data.frame(mimic60_smooth_median[id %in% ind_ids,])
df_mimic_median_hist <- data.frame(mimic60_smooth_median[id %in% hist_ids,])
rm(mimic60_smooth_median)

# add more history to data for individualized fits
start <- proc.time()
individual <- add_lost_rows(df_mimic_median_ind, individual)
rm(df_mimic_median_ind)
individual <- add_more_history(individual, vars_median)
proc.time()-start
save(individual, file="~/Downloads/individual_median_full.Rdata", compress=T)
individual <- slice_rows(individual)
save(individual, file="~/Downloads/individual_median.Rdata", compress=T)
rm(individual)

# add more history to data for historical fit
start <- proc.time()
historical <- add_lost_rows(df_mimic_median_hist, historical)
rm(df_mimic_median_hist)
historical <- add_more_history(historical, vars_median)
proc.time()-start
save(historical, file="~/Downloads/historical_median_full.Rdata", compress=T)
historical <- slice_rows(historical)
save(historical, file="~/Downloads/historical_median.Rdata", compress=T)
rm(historical)
