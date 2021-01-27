library(data.table)

# ==============================================================================
add_shorter_horizon <- function(d, smooth_type){

  # set up ordering
  l <- ncol(d)
  Ycols <- grep("Y", colnames(d))
  Ystart <- min(Ycols)

  setorder(d, subject_id, id, min)
  if(smooth_type == "mean"){
    d[, Y5_lag5_mean := shift(abpmean_lag5_mean, n=5, type="lead"), by=id]
    d[, Y10_lag5_mean := shift(abpmean_lag5_mean, n=10, type="lead"), by=id]
  } else if(smooth_type == "median"){
    d[, Y5_lag5_median := shift(abpmean_lag5_median, n=5, type="lead"), by=id]
    d[, Y10_lag5_median := shift(abpmean_lag5_median, n=10, type="lead"), by=id]
  }
  d[, c(1:(Ystart-1), l+1, l+2, Ystart:l), with=F]
}

add_longer_horizon <- function(d, smooth_type){

  # set up ordering
  l <- ncol(d)
  Ycols <- grep("Y", colnames(d))
  Yend <- max(Ycols)

  setorder(d, subject_id, id, min)
  if(smooth_type == "mean"){
    d[, Y50_lag5_mean := shift(abpmean_lag5_mean, n=50, type="lead"), by=id]
    d[, Y55_lag5_mean := shift(abpmean_lag5_mean, n=55, type="lead"), by=id]
    d[, Y60_lag5_mean := shift(abpmean_lag5_mean, n=60, type="lead"), by=id]
  } else if(smooth_type == "median"){
    d[, Y50_lag5_median := shift(abpmean_lag5_median, n=50, type="lead"), by=id]
    d[, Y55_lag5_median := shift(abpmean_lag5_median, n=55, type="lead"), by=id]
    d[, Y60_lag5_median := shift(abpmean_lag5_median, n=60, type="lead"), by=id]
  }
  d[, c(1:Yend, l+1, l+2, l+3, (Yend+1):l), with=F]
}

add_discrete_outcomes <- function(d, smooth_type){

  setorder(d, subject_id, id, min)

  # binary Y
  d[, Y5_AHE := shift(AHE, n=5, type="lead"), by=id]
  d[, Y10_AHE := shift(AHE, n=10, type="lead"), by=id]
  d[, Y15_AHE := shift(AHE, n=15, type="lead"), by=id]
  d[, Y20_AHE := shift(AHE, n=20, type="lead"), by=id]
  d[, Y25_AHE := shift(AHE, n=25, type="lead"), by=id]
  d[, Y30_AHE := shift(AHE, n=30, type="lead"), by=id]
  d[, Y35_AHE := shift(AHE, n=35, type="lead"), by=id]
  d[, Y40_AHE := shift(AHE, n=40, type="lead"), by=id]
  d[, Y45_AHE := shift(AHE, n=45, type="lead"), by=id]
  d[, Y50_AHE := shift(AHE, n=50, type="lead"), by=id]
  d[, Y55_AHE := shift(AHE, n=55, type="lead"), by=id]
  d[, Y60_AHE := shift(AHE, n=60, type="lead"), by=id]

  # categorical, ordinal Y
  if(smooth_type == "mean"){
    Y <- d[["abpmean_lag5_mean"]]
  } else {
    Y <- d[["abpmean_lag5_median"]]
  }
  AHEcat <- rep(NA, nrow(d))
  AHEcat[Y <= 60] <- "veryhigh"
  AHEcat[Y > 60 & Y <= 65] <- "high"
  AHEcat[Y > 65 & Y <= 70] <- "moderate"
  AHEcat[Y > 70] <- "low"
  AHEcat <- factor(AHEcat, levels = c("low", "moderate", "high", "veryhigh"))
  d$AHEcat <- AHEcat
  d[, Y5_AHEcat := shift(AHEcat, n=5, type="lead"), by=id]
  d[, Y10_AHEcat := shift(AHEcat, n=10, type="lead"), by=id]
  d[, Y15_AHEcat := shift(AHEcat, n=15, type="lead"), by=id]
  d[, Y20_AHEcat := shift(AHEcat, n=20, type="lead"), by=id]
  d[, Y25_AHEcat := shift(AHEcat, n=25, type="lead"), by=id]
  d[, Y30_AHEcat := shift(AHEcat, n=30, type="lead"), by=id]
  d[, Y35_AHEcat := shift(AHEcat, n=35, type="lead"), by=id]
  d[, Y40_AHEcat := shift(AHEcat, n=40, type="lead"), by=id]
  d[, Y45_AHEcat := shift(AHEcat, n=45, type="lead"), by=id]
  d[, Y50_AHEcat := shift(AHEcat, n=50, type="lead"), by=id]
  d[, Y55_AHEcat := shift(AHEcat, n=55, type="lead"), by=id]
  d[, Y60_AHEcat := shift(AHEcat, n=60, type="lead"), by=id]
}
# ==============================================================================
data_path <- "~/Box/Data_Shared/"
save_path <- "~/Box/symphony-data/"

load(paste0(data_path, "historical_median.Rdata"))
historical <- add_shorter_horizon(historical, "median")
historical <- add_longer_horizon(historical, "median")
historical <- add_discrete_outcomes(historical, "median")
save(historical, file=paste0(save_path, "historical_median.Rdata"), compress=T)
save(historical, file=paste0(data_path, "historical_median_update.Rdata"), compress=T)
rm(historical)

load(paste0(data_path, "individual_median.Rdata"))
individual <- add_shorter_horizon(individual, "median")
individual <- add_longer_horizon(individual, "median")
individual <- add_discrete_outcomes(individual, "median")
save(individual, file=paste0(save_path, "individual_median.Rdata"), compress=T)
save(individual, file=paste0(data_path, "individual_median_update.Rdata"), compress=T)
rm(individual)

load(paste0(data_path, "historical_mean.Rdata"))
historical <- add_shorter_horizon(historical, "mean")
historical <- add_longer_horizon(historical, "mean")
historical <- add_discrete_outcomes(historical, "mean")
save(historical, file=paste0(save_path, "historical_mean.Rdata"), compress=T)
save(historical, file=paste0(data_path, "historical_mean_update.Rdata"), compress=T)
rm(historical)

load(paste0(data_path, "individual_mean.Rdata"))
individual <- add_shorter_horizon(individual, "mean")
individual <- add_longer_horizon(individual, "mean")
individual <- add_discrete_outcomes(individual, "mean")
save(individual, file=paste0(save_path, "individual_mean.Rdata"), compress=T)
save(individual, file=paste0(data_path, "individual_mean_update.Rdata"), compress=T)
rm(individual)
