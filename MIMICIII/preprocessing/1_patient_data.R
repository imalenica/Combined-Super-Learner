library(data.table)
library(dplyr)
library(feasts)
library(tsibble)

box <- "/Users/Rachael/Library/CloudStorage/Box-Box/"
vitals_files <- list.files(paste0(box, "MIMIC-III/posl_numerics_trunc48hrs/"))

files <- sub(".", "", vitals_files)
subject_ids <- gsub("_.*","", files)
subject_ids <- as.numeric(subject_ids)

icustay_ids <- substring(files, 8)
icustay_ids <- gsub("_.*","", icustay_ids)
icustay_ids <- as.numeric(icustay_ids)
patients_keep <- data.table("SUBJECT_ID" = subject_ids, "ICUSTAY_ID" = icustay_ids)

######################## time-invariant information ############################
# import time-invariant covariate data, removing ROW_ID column indicating row number
icustays <- read.csv(paste0(box, "MIMIC-III/ICUSTAYS.csv"))[,-1]
patients <- read.csv(paste0(box, "MIMIC-III/PATIENTS.csv"))[,-1]
admissions <- read.csv(paste0(box, "MIMIC-III/ADMISSIONS.csv"))[,-1]
baseline <- full_join(full_join(admissions, patients), icustays)
setDT(baseline)
# baseline <- baseline[DBSOURCE == "metavision", ]
baseline <- inner_join(baseline, patients_keep)

chartevents <- read.csv(paste0(box, "MIMIC-III/filtered_chartevents.csv"))
setDT(chartevents)
chartevents <- left_join(patients_keep, chartevents)
keep <- c("SUBJECT_ID", "HADM_ID", "ICUSTAY_ID", "VALUENUM")

height_cm <- chartevents[LABEL == "Height (cm)", keep, with = F]
height_cm[, "height_meters" := VALUENUM*0.01]
height_in <- chartevents[LABEL == "Height", keep, with = F]
height_in[, "height_meters" := VALUENUM*0.0254]
height <- rbind(height_in[, -"VALUENUM", with=F], height_cm[, -"VALUENUM", with=F])
height <- height[!duplicated(height[,-"height_meters",with=F]),]

weight_kg <- chartevents[LABEL == "Admission Weight (Kg)", keep, with = F]
weight_kg[, "weight_kg" := VALUENUM]
weight_lbs <- chartevents[LABEL == "Admission Weight (lbs.)", keep, with = F]
weight_lbs[, "weight_kg" := VALUENUM*0.45359237]
weight <- rbind(weight_kg[, -"VALUENUM", with=F], weight_lbs[, -"VALUENUM", with=F])
weight <- weight[!duplicated(weight[,-"weight_kg",with=F]),]
height_weight <- full_join(height, weight)
height_weight$BMI <- height_weight$weight_kg/(height_weight$height_meters^2)

baseline <- full_join(baseline, height_weight)
baseline <- baseline[DBSOURCE == "metavision", ]
baseline[, "age" := as.numeric(as.POSIXct(ADMITTIME)-as.POSIXct(DOB))/365.25]
baseline <- inner_join(patients_keep, baseline)
write.csv(baseline, file = paste0(box, "MIMIC-III/output/patient_data.csv"),
          row.names = F)

patient_data <- read.csv(paste0(box, "MIMIC-III/output/patient_data.csv"))
setDT(patient_data)
patient_data <- patient_data[,c(
  "SUBJECT_ID", "ICUSTAY_ID", "INSURANCE", "ETHNICITY",
  "GENDER", "FIRST_CAREUNIT", "height_meters", "weight_kg", "BMI", "age"
), with = FALSE]
patient_data[, "INSURANCE_PRIV" := ifelse(INSURANCE == "Private", 1, 0)]
patient_data[, "SEX_M" := ifelse(GENDER == "M", 1, 0)]
unk <- c("OTHER", "UNKNOWN/NOT SPECIFIED", "UNABLE TO OBTAIN",
         "PATIENT DECLINED TO ANSWER")
white <- c("WHITE - RUSSIAN", "WHITE")
patient_data[,"ETHNICITY" := ifelse(
  is.na(ETHNICITY) | ETHNICITY %in% unk, "UNKNOWN",
  ifelse(ETHNICITY %in% white, "WHITE", "NON-WHITE")
)]
patient_data[, "ETHNICITY_WHITE" := ifelse(ETHNICITY == "WHITE", 1, 0)]
patient_data[, "ETHNICITY_NONWHITE" := ifelse(ETHNICITY == "NON-WHITE", 1, 0)]
patient_data[, "FIRST_CAREUNIT_TSICU" := ifelse(FIRST_CAREUNIT == "TSICU", 1, 0)]
patient_data[, "FIRST_CAREUNIT_SICU" := ifelse(FIRST_CAREUNIT == "SICU", 1, 0)]
patient_data[, "FIRST_CAREUNIT_MICU" := ifelse(FIRST_CAREUNIT == "MICU", 1, 0)]
patient_data[, "FIRST_CAREUNIT_CSRU" := ifelse(FIRST_CAREUNIT == "CSRU", 1, 0)]
patient_data[,"AGE":= ifelse(age > 100, NA, age)]
patient_data <- patient_data[,c(
  "SUBJECT_ID", "ICUSTAY_ID", "height_meters", "weight_kg", "BMI", "AGE",
  "INSURANCE_PRIV", "SEX_M", "ETHNICITY_WHITE", "ETHNICITY_NONWHITE",
  "FIRST_CAREUNIT_TSICU", "FIRST_CAREUNIT_SICU", "FIRST_CAREUNIT_MICU",
  "FIRST_CAREUNIT_CSRU"
), with=FALSE]
colnames(patient_data)[3:4] <- c("HEIGHT", "WEIGHT")
write.csv(
  patient_data, file = paste0(box, "MIMIC-III/output/processed_patient_data.csv"),
  row.names = FALSE
)
patient_data_complete <- patient_data[complete.cases(patient_data),]
write.csv(
  patient_data_complete,
  file = paste0(box, "MIMIC-III/output/processed_complete_patient_data.csv"),
  row.names = FALSE
)
