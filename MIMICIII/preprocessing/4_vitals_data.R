library(data.table)
library(dplyr)
library(arrow)
library(feasts)
library(tsibble)
library(DescTools)

box <- "/Users/Rachael/Library/CloudStorage/Box-Box/"
vitals_files <- list.files(paste0(box, "MIMIC-III/posl_numerics_trunc48hrs/"))
vitals_path <- paste0(box, "MIMIC-III/posl_numerics_trunc48hrs/")

files <- sub(".", "", vitals_files)
subject_ids <- gsub("_.*","", files)
subject_ids <- as.numeric(subject_ids)

patients_keep <- read.csv(paste0(box, "MIMIC-III/output/patients_keep.csv"))
setDT(patients_keep)

dopamine <- read.csv(paste0(box, "MIMIC-III/output/dopamine.csv"))
epinephrine <- read.csv(paste0(box, "MIMIC-III/output/epinephrine.csv"))
phenylephrine <- read.csv(paste0(box, "MIMIC-III/output/phenylephrine.csv"))
vasopressin <- read.csv(paste0(box, "MIMIC-III/output/vasopressin.csv"))
norepinephrine <- read.csv(paste0(box, "MIMIC-III/output/norepinephrine.csv"))
ventilation <- read.csv(paste0(box, "MIMIC-III/output/ventilation.csv"))
patient_data <- read.csv(paste0(box, "MIMIC-III/output/processed_patient_data.csv"))

sparse_vitals <- c(
  "N_ABP", "N_NBP SYS", "N_NBP DIAS", "N_NBP MEAN", "N_ART", "N_ART SYS",
  "N_ART DIAS", "N_ART MEAN", "N_PLSNBP", "N_IRREGULAR HR PERCENT",
  "N_HR VARIABILITY", "N_NBP", "N_CPP"
)
# ==============================================================================
# function for adding time-varying treatment data
join_tables <- function(tbl, patient_info, patient_vitals, newcol_name,
                        colvalue = NULL){
  added_col <- rep(0, nrow(patient_vitals))
  joined_tbl <- suppressMessages(data.table(inner_join(patient_info, tbl)))
  if(nrow(joined_tbl) > 0){
    for(j in 1:nrow(joined_tbl)){
      vitals_start <- min(patient_vitals$CHARTTIME)
      vitals_end <- min(patient_vitals$CHARTTIME)
      start_j <- joined_tbl[j,][["STARTTIME"]]
      end_j <- joined_tbl[j,][["ENDTIME"]]
      # don't include tbl info before/after vitals start/end
      if ( !((vitals_start > end_j) || (vitals_end < start_j)) ) {
        start_idx <- ifelse(
          vitals_start < start_j, which(patient_vitals$CHARTTIME == start_j), 1
        )
        end_idx <- ifelse(
          vitals_end > end_j, which(patient_vitals$CHARTTIME == end_j),
          nrow(patient_vitals)
        )
        if(is.null(colvalue)){
          added_col[start_idx:end_idx] <- 1
        } else {
          added_col[start_idx:end_idx] <- joined_tbl[j,][[colvalue]]
        }
      }
    }
  }
  patient_vitals$newcol <- added_col
  newcol_idx <- which(colnames(patient_vitals) == "newcol")
  colnames(patient_vitals)[newcol_idx] <- newcol_name
  return(patient_vitals)
}

# add new covariates based on summary measures of the history
add_historical_features <- function(vitals_history, feature, history_time = 1,
                                    basic_summary = TRUE){

  setDT(vitals_history)
  if(basic_summary) {
    feature_summary <- c(
      "median" = median(vitals_history[[feature]], na.rm = T),
      "mean" = mean(vitals_history[[feature]], na.rm = T),
      "var" = var(vitals_history[[feature]], na.rm = T),
      "NA" = mean(is.na(vitals_history[[feature]]), na.rm = T)
    )
    if (history_time != 1) {
      names(feature_summary) <- paste0(
        feature, "_history", history_time, "_", names(feature_summary)
      )
    } else {
      names(feature_summary) <- paste0(feature, "_", names(feature_summary))
    }

    return_dt <- data.table(t(unclass(feature_summary)))
  } else {
    feature_summary <- unname(summary(vitals_history[[feature]]))
    if(length(feature_summary) == 6) {
      feature_summary <- c(feature_summary[1:6], 0)
    } else if(length(feature_summary) == 7) {
      feature_summary <- c(feature_summary[1:6], feature_summary[7]/nrow(vitals_history))
    }
    feature_summary <- c(feature_summary, var(vitals_history[[feature]], na.rm = T))
    names(feature_summary) <- c("min", "Q1", "median", "mean", "Q3", "max", "NA", "var")
    if (history_time != 1) {
      names(feature_summary) <- paste0(
        feature, "_history", history_time, "_", names(feature_summary)
      )
    } else {
      names(feature_summary) <- paste0(feature, "_", names(feature_summary))
    }

    feature_ts <- as_tsibble(ts(vitals_history[[feature]]))
    feature_ts$feature <- unname(unlist(feature_ts[, "value"]))
    features_stl <- suppressWarnings(suppressMessages(tryCatch(
      feature_ts %>% features(feature, list(feat_stl)),
      error = function(err) NA,
      warning = function(warn) NA
    )))
    if((length(features_stl) == 1 && is.na(features_stl)) |
       length(features_stl) == 0){
      features_stl <- t(c(
        "trend_strength"=NA, "spikiness"=NA, "linearity"=NA, "curvature"=NA,
        "stl_e_acf1"=NA, "stl_e_acf10"=NA
      ))
    }
    spectral_entropy <- suppressWarnings(suppressMessages(tryCatch(
      feature_ts %>% features(feature, list(feat_spectral)),
      error = function(err) NA,
      warning = function(warn) NA
    )))
    longest_flat_spot <- suppressWarnings(suppressMessages(tryCatch(
      longest_flat_spot(x = ts(vitals_history[[feature]])),
      error = function(err) NA,
      warning = function(warn) NA
    )))
    coef_hurst <- suppressWarnings(suppressMessages(tryCatch(
      coef_hurst(x = ts(vitals_history[[feature]])),
      error = function(err) NA,
      warning = function(warn) NA
    )))
    feature_ts <- data.table(
      features_stl, spectral_entropy, longest_flat_spot, coef_hurst
    )
    if (history_time != 1) {
      colnames(feature_ts) <- paste0(
        feature, "_history", history_time, "_", colnames(feature_ts)
      )
    } else {
      colnames(feature_ts) <- paste0(feature, "_", colnames(feature_ts))
    }

    return_dt <- data.table(feature_ts, data.table(t(unclass(feature_summary))))
  }
  return(return_dt)
}
# ==============================================================================

col_ordering <- c(
  "SUBJECT_ID", "ICUSTAY_ID", "HEIGHT", "WEIGHT", "BMI", "AGE",
  "INSURANCE_PRIV", "SEX_M", "ETHNICITY_WHITE", "ETHNICITY_NONWHITE",
  "FIRST_CAREUNIT_TSICU", "FIRST_CAREUNIT_SICU", "FIRST_CAREUNIT_MICU",
  "FIRST_CAREUNIT_CSRU", "date_time_j", "minute", "MAP_trend_strength",
  "MAP_spikiness", "MAP_linearity", "MAP_curvature", "MAP_stl_e_acf1",
  "MAP_stl_e_acf10", "MAP_spectral_entropy", "MAP_longest_flat_spot",
  "MAP_coef_hurst", "MAP_min", "MAP_Q1", "MAP_median", "MAP_mean",
  "MAP_Q3", "MAP_max", "MAP_NA", "MAP_var", "HR_median", "HR_mean",
  "HR_var", "HR_NA", "PULSE_median", "PULSE_mean", "PULSE_var",
  "PULSE_NA", "SYS_median", "SYS_mean", "SYS_var", "SYS_NA", "DIAS_median",
  "DIAS_mean", "DIAS_var", "DIAS_NA", "RESP_median", "RESP_mean",
  "RESP_var", "RESP_NA", "SPO2_median", "SPO2_mean", "SPO2_var",
  "SPO2_NA", "DA", "NE", "AVP", "EPI", "PE", "VENT", "MAP_history5_trend_strength",
  "MAP_history5_spikiness", "MAP_history5_linearity", "MAP_history5_curvature",
  "MAP_history5_stl_e_acf1", "MAP_history5_stl_e_acf10", "MAP_history5_spectral_entropy",
  "MAP_history5_longest_flat_spot", "MAP_history5_coef_hurst",
  "MAP_history5_min", "MAP_history5_Q1", "MAP_history5_median",
  "MAP_history5_mean", "MAP_history5_Q3", "MAP_history5_max", "MAP_history5_NA",
  "MAP_history5_var", "HR_history5_median", "HR_history5_mean",
  "HR_history5_var", "HR_history5_NA", "PULSE_history5_median",
  "PULSE_history5_mean", "PULSE_history5_var", "PULSE_history5_NA",
  "SYS_history5_median", "SYS_history5_mean", "SYS_history5_var",
  "SYS_history5_NA", "DIAS_history5_median", "DIAS_history5_mean",
  "DIAS_history5_var", "DIAS_history5_NA", "RESP_history5_median",
  "RESP_history5_mean", "RESP_history5_var", "RESP_history5_NA",
  "SPO2_history5_median", "SPO2_history5_mean", "SPO2_history5_var",
  "SPO2_history5_NA", "DA_history5", "NE_history5", "AVP_history5",
  "EPI_history5", "PE_history5", "VENT_history5", "MAP_history15_trend_strength",
  "MAP_history15_spikiness", "MAP_history15_linearity", "MAP_history15_curvature",
  "MAP_history15_stl_e_acf1", "MAP_history15_stl_e_acf10", "MAP_history15_spectral_entropy",
  "MAP_history15_longest_flat_spot", "MAP_history15_coef_hurst",
  "MAP_history15_min", "MAP_history15_Q1", "MAP_history15_median",
  "MAP_history15_mean", "MAP_history15_Q3", "MAP_history15_max",
  "MAP_history15_NA", "MAP_history15_var", "HR_history15_median",
  "HR_history15_mean", "HR_history15_var", "HR_history15_NA", "PULSE_history15_median",
  "PULSE_history15_mean", "PULSE_history15_var", "PULSE_history15_NA",
  "SYS_history15_median", "SYS_history15_mean", "SYS_history15_var",
  "SYS_history15_NA", "DIAS_history15_median", "DIAS_history15_mean",
  "DIAS_history15_var", "DIAS_history15_NA", "RESP_history15_median",
  "RESP_history15_mean", "RESP_history15_var", "RESP_history15_NA",
  "SPO2_history15_median", "SPO2_history15_mean", "SPO2_history15_var",
  "SPO2_history15_NA", "DA_history15", "NE_history15", "AVP_history15",
  "EPI_history15", "PE_history15", "VENT_history15"
)

remaining_na_summary <- list()
cols <- list()

start_time <- proc.time()
for(i in 1:nrow(patients_keep)){

  patient_i <- patients_keep[i,]
  patient_i_data <- suppressMessages(left_join(patient_i, patient_data))

  idx <- which(subject_ids == as.numeric(patient_i$SUBJECT_ID))
  vitals_i <- arrow::read_parquet(paste0(vitals_path, vitals_files[[idx]]))
  if (any(colnames(vitals_i) == "N_%SPO2")) {
    colnames(vitals_i)[which(colnames(vitals_i) == "N_%SPO2")] <- "N_SPO2"
  }
  vitals_i <- vitals_i[,-which(colnames(vitals_i) %in% sparse_vitals), with=FALSE]
  vitals_i$ICUSTAY_ID <- patient_i$ICUSTAY_ID

  vitals_i <- join_tables(
    dopamine, patient_i, vitals_i, "dopamine", "dopamine_NEE"
  )
  vitals_i <- join_tables(
    norepinephrine, patient_i, vitals_i, "norepinephrine", "norepinephrine_NEE"
  )
  vitals_i <- join_tables(
    vasopressin, patient_i, vitals_i, "vasopressin", "vasopressin_NEE"
  )
  vitals_i <- join_tables(
    epinephrine, patient_i, vitals_i, "epinephrine", "epinephrine_NEE"
  )
  vitals_i <- join_tables(
    phenylephrine, patient_i, vitals_i, "phenylephrine", "phenylephrine_NEE"
  )
  vitals_i <- join_tables(
    ventilation, patient_i, vitals_i, "ventilation"
  )

  # format time with seconds removed so we can group by each minute
  vitals_i <- vitals_i[order(vitals_i$CHARTTIME), ]
  vitals_i$date_time <- format(vitals_i$CHARTTIME, format = '%Y%m%d %H:%M')
  setDT(vitals_i)
  summ_cols <- c(
    "N_HR", "N_PULSE", "N_ABP SYS", "N_ABP DIAS", "N_ABP MEAN", "N_RESP",
    "N_SPO2", "dopamine", "norepinephrine", "vasopressin", "epinephrine",
    "phenylephrine", "ventilation"
  )
  vitals_i <- vitals_i[, c(summ_cols,"date_time"), with=FALSE]
  colnames(vitals_i) <- c(
    "HR", "PULSE", "SYS", "DIAS", "MAP", "RESP", "SPO2", "DA", "NE", "AVP",
    "EPI", "PE", "VENT", "date_time"
  )

  data_i_list <- list()
  for(j in 1:length(unique(vitals_i$date_time))){
    print(paste0(
      "adding history for minute ", j, " of ", length(unique(vitals_i$date_time)),
      " for patient ", i, " of ", nrow(patients_keep)
    ))

    ########### summarize time-varying covariates over each minute #############
    date_time_j <- unique(vitals_i$date_time)[j]
    vitals_i_history_1min <- vitals_i[vitals_i$date_time == date_time_j,]
    vitals_i_j <- data.table(
      add_historical_features(vitals_i_history_1min, "MAP", basic_summary=F),
      add_historical_features(vitals_i_history_1min, "HR"),
      add_historical_features(vitals_i_history_1min, "PULSE"),
      add_historical_features(vitals_i_history_1min, "SYS"),
      add_historical_features(vitals_i_history_1min, "DIAS"),
      add_historical_features(vitals_i_history_1min, "RESP"),
      add_historical_features(vitals_i_history_1min, "SPO2"),
      "DA" = mean(vitals_i_history_1min[["DA"]]),
      "NE" = mean(vitals_i_history_1min[["NE"]]),
      "AVP" = mean(vitals_i_history_1min[["AVP"]]),
      "EPI" = mean(vitals_i_history_1min[["EPI"]]),
      "PE" = mean(vitals_i_history_1min[["PE"]]),
      "VENT" = mean(vitals_i_history_1min[["VENT"]])
    )

    if (j >= 5) {
      history_5min <- c(
        date_time_j, unique(vitals_i$date_time)[j-1],
        unique(vitals_i$date_time)[j-2], unique(vitals_i$date_time)[j-3],
        unique(vitals_i$date_time)[j-4]
      )
      vitals_i_history_5min <- vitals_i[vitals_i$date_time %in% history_5min,]
      vitals_i_j5 <- data.table(
        add_historical_features(vitals_i_history_5min, "MAP", 5, FALSE),
        add_historical_features(vitals_i_history_5min, "HR", 5),
        add_historical_features(vitals_i_history_5min, "PULSE", 5),
        add_historical_features(vitals_i_history_5min, "SYS", 5),
        add_historical_features(vitals_i_history_5min, "DIAS", 5),
        add_historical_features(vitals_i_history_5min, "RESP", 5),
        add_historical_features(vitals_i_history_5min, "SPO2", 5),
        "DA_history5" = mean(vitals_i_history_5min[["DA"]]),
        "NE_history5" = mean(vitals_i_history_5min[["NE"]]),
        "AVP_history5" = mean(vitals_i_history_5min[["AVP"]]),
        "EPI_history5" = mean(vitals_i_history_5min[["EPI"]]),
        "PE_history5" = mean(vitals_i_history_5min[["PE"]]),
        "VENT_history5" = mean(vitals_i_history_5min[["VENT"]])
      )
      vitals_i_j <- data.table(vitals_i_j, vitals_i_j5)
    }
    if (j >= 15) {
      history_15min <- c(
        history_5min, unique(vitals_i$date_time)[j-5],
        unique(vitals_i$date_time)[j-6], unique(vitals_i$date_time)[j-7],
        unique(vitals_i$date_time)[j-8], unique(vitals_i$date_time)[j-9],
        unique(vitals_i$date_time)[j-10], unique(vitals_i$date_time)[j-11],
        unique(vitals_i$date_time)[j-12], unique(vitals_i$date_time)[j-13],
        unique(vitals_i$date_time)[j-14]
      )
      vitals_i_history_15min <- vitals_i[vitals_i$date_time %in% history_15min,]
      vitals_i_j15 <- data.table(
        add_historical_features(vitals_i_history_15min, "MAP", 15, FALSE),
        add_historical_features(vitals_i_history_15min, "HR", 15),
        add_historical_features(vitals_i_history_15min, "PULSE", 15),
        add_historical_features(vitals_i_history_15min, "SYS", 15),
        add_historical_features(vitals_i_history_15min, "DIAS", 15),
        add_historical_features(vitals_i_history_15min, "RESP", 15),
        add_historical_features(vitals_i_history_15min, "SPO2", 15),
        "DA_history15" = mean(vitals_i_history_15min[["DA"]]),
        "NE_history15" = mean(vitals_i_history_15min[["NE"]]),
        "AVP_history15" = mean(vitals_i_history_15min[["AVP"]]),
        "EPI_history15" = mean(vitals_i_history_15min[["EPI"]]),
        "PE_history15" = mean(vitals_i_history_15min[["PE"]]),
        "VENT_history15" = mean(vitals_i_history_15min[["VENT"]])
      )
      vitals_i_j <- data.table(vitals_i_j, vitals_i_j15)
    }

    minute <- as.numeric(difftime(
      as.POSIXct(unique(vitals_i$date_time)[j], format = '%Y%m%d %H:%M'),
      as.POSIXct(unique(vitals_i$date_time)[1], format = '%Y%m%d %H:%M'),
      units = "mins"
    )) + 1
    data_i_list[[j]] <- data.table(
      patient_i_data, date_time_j, minute, vitals_i_j
    )
  }
  # bind all rows
  data_i <- rbindlist(data_i_list, fill = TRUE)
  setDT(data_i)

  # save cols, so we can later make sure they're identical for all subjects
  cols[[i]] <- colnames(data_i)

  # final check for NA, if so try LOCF
  remaining_na_summary[[i]] <- data.table(
    "SUBJECT_ID" = patient_i_data$SUBJECT_ID,
    "ICUSTAY_ID" = patient_i_data$ICUSTAY_ID, "na_remaining" = FALSE,
    "LOCF" = FALSE
  )
  na_by_col <- colSums(is.na(data_i[15:nrow(data_i),]))
  if (any(na_by_col > 0)) {
    # flag LOCF and summarize missingness after first 15 rows
    remaining_na_summary[[i]]$LOCF <- TRUE
    remaining_na_summary[[i]] <- data.table(
      remaining_na_summary[[i]], data.table(t(na_by_col[which(na_by_col > 0)]))
    )
    # seperate no missingness cols from those w missingness
    d_complete <- data_i[, which(na_by_col == 0), with = F]
    d_incomplete <- data_i[, which(na_by_col > 0), with = F]
    # try LOCF, flagging if it doesn't work
    d_incomplete_locf <- data.table(
      apply(d_incomplete, 2, function(col) DescTools::LOCF(col))
    )
    if (any(colSums(is.na(d_incomplete_locf[15:nrow(data_i),])) > 0)) {
      remaining_na_summary[[i]]$na_remaining <- TRUE
    }
    # put it all back together
    data_i <- data.table(d_incomplete_locf, d_complete)
  }
  data_i <- data_i[, col_ordering, with = FALSE]

  # save patient i's data
  save(
    data_i, compress = TRUE,
    file = paste0(box, "MIMIC-III/output/patients/", patient_i$SUBJECT_ID, "_",
                  patient_i$ICUSTAY_ID, ".Rdata")
  )
}
end_time <- proc.time() - start_time
end_time
# user    system   elapsed
# 35310.095   949.506 36930.538

names(cols) <- paste0(patients_keep$SUBJECT_ID, "_", patients_keep$ICUSTAY_ID)
save(cols, compress = TRUE, file = paste0(box, "MIMIC-III/output/cols_list.Rdata"))
length(unique(cols))


remaining_na_summary_tbl <- rbindlist(remaining_na_summary, fill = TRUE)
setDT(remaining_na_summary_tbl)
write.csv(
  remaining_na_summary_tbl, row.names = FALSE,
  file = paste0(box, "MIMIC-III/output/remaining_na_summary.csv")
)
dput(colnames(remaining_na_summary_tbl)[which(colSums(remaining_na_summary_tbl, na.rm = T) > 142700)])
to_rm <- c(
  "MAP_trend_strength", "MAP_stl_e_acf1", "MAP_stl_e_acf10",
  "MAP_spectral_entropy", "MAP_spikiness", "MAP_linearity",
  "MAP_curvature", "MAP_coef_hurst", "MAP_var", "HR_var", "PULSE_var",
  "SYS_var", "DIAS_var", "RESP_var", "SPO2_var", "MAP_history5_stl_e_acf10"
)
remaining_na_filtered <- remaining_na_summary_tbl[,-to_rm,with=F]

bad_id <- which(apply(remaining_na_filtered[,-c(1:4)], 1, function(row) any(row >= 1427)))
remaining_na_filtered[bad_id,]
# SUBJECT_ID ICUSTAY_ID
#     76568     290983
patients_keep_filtered <- patients_keep[-which(patients_keep$SUBJECT_ID == 76568)]

# make sure no NA with to_rm columns removed
na_by_col_list <- list()
any_na_list <- list()
for(i in 1:nrow(patients_keep_filtered)){
  patient_i <- patients_keep_filtered[i,]
  load(paste0(
    box, "MIMIC-III/output/patients/", patient_i$SUBJECT_ID, "_",
    patient_i$ICUSTAY_ID, ".Rdata"
  ))
  data_i <- data_i[, -to_rm, with=F]
  na_by_col_list[[i]] <- colSums(is.na(data_i[15:nrow(data_i),]))
  any_na_list[[i]] <- any(colSums(is.na(data_i[15:nrow(data_i),])) > 0)
}
ids_with_na <- which(any_na_list==T)
patients_keep_final <- patients_keep_filtered[-ids_with_na,]
write.csv(patients_keep_final,
          file = paste0(box, "MIMIC-III/output/patients_keep_final.csv"),
          row.names = FALSE)

covs <- col_ordering[-which(col_ordering %in% to_rm)]
covariates <- covs[-c(1:2)]
dput(covariates)
# c("HEIGHT", "WEIGHT", "BMI", "AGE", "INSURANCE_PRIV", "SEX_M",
#   "ETHNICITY_WHITE", "ETHNICITY_NONWHITE", "FIRST_CAREUNIT_TSICU",
#   "FIRST_CAREUNIT_SICU", "FIRST_CAREUNIT_MICU", "FIRST_CAREUNIT_CSRU",
#   "date_time_j", "minute", "MAP_longest_flat_spot", "MAP_min",
#   "MAP_Q1", "MAP_median", "MAP_mean", "MAP_Q3", "MAP_max", "MAP_NA",
#   "HR_median", "HR_mean", "HR_NA", "PULSE_median", "PULSE_mean",
#   "PULSE_NA", "SYS_median", "SYS_mean", "SYS_NA", "DIAS_median",
#   "DIAS_mean", "DIAS_NA", "RESP_median", "RESP_mean", "RESP_NA",
#   "SPO2_median", "SPO2_mean", "SPO2_NA", "DA", "NE", "AVP", "EPI",
#   "PE", "VENT", "MAP_history5_trend_strength", "MAP_history5_spikiness",
#   "MAP_history5_linearity", "MAP_history5_curvature", "MAP_history5_stl_e_acf1",
#   "MAP_history5_spectral_entropy", "MAP_history5_longest_flat_spot",
#   "MAP_history5_coef_hurst", "MAP_history5_min", "MAP_history5_Q1",
#   "MAP_history5_median", "MAP_history5_mean", "MAP_history5_Q3",
#   "MAP_history5_max", "MAP_history5_NA", "MAP_history5_var", "HR_history5_median",
#   "HR_history5_mean", "HR_history5_var", "HR_history5_NA", "PULSE_history5_median",
#   "PULSE_history5_mean", "PULSE_history5_var", "PULSE_history5_NA",
#   "SYS_history5_median", "SYS_history5_mean", "SYS_history5_var",
#   "SYS_history5_NA", "DIAS_history5_median", "DIAS_history5_mean",
#   "DIAS_history5_var", "DIAS_history5_NA", "RESP_history5_median",
#   "RESP_history5_mean", "RESP_history5_var", "RESP_history5_NA",
#   "SPO2_history5_median", "SPO2_history5_mean", "SPO2_history5_var",
#   "SPO2_history5_NA", "DA_history5", "NE_history5", "AVP_history5",
#   "EPI_history5", "PE_history5", "VENT_history5", "MAP_history15_trend_strength",
#   "MAP_history15_spikiness", "MAP_history15_linearity", "MAP_history15_curvature",
#   "MAP_history15_stl_e_acf1", "MAP_history15_stl_e_acf10", "MAP_history15_spectral_entropy",
#   "MAP_history15_longest_flat_spot", "MAP_history15_coef_hurst",
#   "MAP_history15_min", "MAP_history15_Q1", "MAP_history15_median",
#   "MAP_history15_mean", "MAP_history15_Q3", "MAP_history15_max",
#   "MAP_history15_NA", "MAP_history15_var", "HR_history15_median",
#   "HR_history15_mean", "HR_history15_var", "HR_history15_NA", "PULSE_history15_median",
#   "PULSE_history15_mean", "PULSE_history15_var", "PULSE_history15_NA",
#   "SYS_history15_median", "SYS_history15_mean", "SYS_history15_var",
#   "SYS_history15_NA", "DIAS_history15_median", "DIAS_history15_mean",
#   "DIAS_history15_var", "DIAS_history15_NA", "RESP_history15_median",
#   "RESP_history15_mean", "RESP_history15_var", "RESP_history15_NA",
#   "SPO2_history15_median", "SPO2_history15_mean", "SPO2_history15_var",
#   "SPO2_history15_NA", "DA_history15", "NE_history15", "AVP_history15",
#   "EPI_history15", "PE_history15", "VENT_history15")
