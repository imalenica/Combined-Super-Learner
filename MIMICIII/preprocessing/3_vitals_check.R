library(data.table)
library(dplyr)
library(arrow)

box <- "/Users/Rachael/Library/CloudStorage/Box-Box/"
vitals_files <- list.files(paste0(box, "MIMIC-III/posl_numerics_trunc48hrs/"))

files <- sub(".", "", vitals_files)
subject_ids <- gsub("_.*","", files)
subject_ids <- as.numeric(subject_ids)

icustay_ids <- substring(files, 8)
icustay_ids <- gsub("_.*","", icustay_ids)
icustay_ids <- as.numeric(icustay_ids)
patients_keep <- data.table("SUBJECT_ID" = subject_ids,
                            "ICUSTAY_ID" = icustay_ids)

dopamine <- read.csv(paste0(box, "MIMIC-III/output/dopamine.csv"))
epinephrine <- read.csv(paste0(box, "MIMIC-III/output/epinephrine.csv"))
phenylephrine <- read.csv(paste0(box, "MIMIC-III/output/phenylephrine.csv"))
vasopressin <- read.csv(paste0(box, "MIMIC-III/output/vasopressin.csv"))
norepinephrine <- read.csv(paste0(box, "MIMIC-III/output/norepinephrine.csv"))
ventilation <- read.csv(paste0(box, "MIMIC-III/output/ventilation.csv"))
patient_data <- read.csv(paste0(box, "MIMIC-III/output/processed_patient_data.csv"))

vitals_path <- paste0(box, "MIMIC-III/posl_numerics_trunc48hrs/")

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
# ==============================================================================

############################ 1. examine missinginess ###########################
patient_missingness_list <- list()
for(i in 1:nrow(patients_keep)){

  vitals_i <- arrow::read_parquet(paste0(vitals_path, vitals_files[[i]]))
  if (any(colnames(vitals_i) == "N_%SPO2")) {
    colnames(vitals_i)[which(colnames(vitals_i) == "N_%SPO2")] <- "N_SPO2"
  }
  patient_i <- patients_keep[i,]
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
  data_i <- vitals_i[ , lapply(.SD, mean, na.rm = TRUE), by = c("date_time")]
  patient_i <- suppressMessages(left_join(patient_i, patient_data))
  data_i <- data_i[,-which(colnames(data_i) %in% colnames(patient_i)),with=F]
  J <- unique(data_i$date_time)[length(unique(data_i$date_time))]
  max_minute <- as.numeric(difftime(
    as.POSIXct(J, format = '%Y%m%d %H:%M'),
    as.POSIXct(unique(data_i$date_time)[1], format = '%Y%m%d %H:%M'),
    units = "mins"
  ))
  data_i <- data.table(patient_i, data_i)
  missingness_tbl <- colSums(is.na(data_i[, -c(1:2), with=FALSE]))/nrow(data_i)
  patient_missingness_list[[i]] <- data.table(
    "SUBJECT_ID" = unique(data_i$SUBJECT_ID),
    "ICUSTAY_ID" = unique(data_i$ICUSTAY_ID),
    data.table(t(missingness_tbl)),
    "unique_min" = length(unique(vitals_i$date_time)),
    "total_min" = max_minute+1
  )
}
tbl <- rbindlist(patient_missingness_list, fill = TRUE)
tbl <- apply(tbl, 2, function(x) ifelse(is.na(x), 1, x))
tbl <- data.table(tbl)
write.csv(tbl, file = paste0(box, "MIMIC-III/output/patient_missingness.csv"),
          row.names = FALSE)

# all patients were supposed to have 24 hours of data
any(tbl$unique_min != tbl$total_min)
table(tbl$unique_min)
tbl <- tbl[which(tbl$unique_min == 1441), ]
tbl <- tbl[-which(AGE == 1 | HEIGHT == 1),]
patients_keep <- tbl[,c(1:2), with = FALSE]
write.csv(patients_keep, file = paste0(box, "MIMIC-III/output/patients_keep.csv"),
          row.names = FALSE)
dput(names(which(colMeans(tbl) > 0.1)))
sparse_vitals <- c(
  "N_ABP", "N_NBP SYS", "N_NBP DIAS", "N_NBP MEAN", "N_ART", "N_ART SYS",
  "N_ART DIAS", "N_ART MEAN", "N_PLSNBP", "N_IRREGULAR HR PERCENT",
  "N_HR VARIABILITY", "N_NBP", "N_CPP"
)
