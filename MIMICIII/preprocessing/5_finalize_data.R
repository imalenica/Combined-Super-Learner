library(data.table)
library(dplyr)

box <- "/Users/Rachael/Library/CloudStorage/Box-Box/"
vitals_files <- list.files(paste0(box, "MIMIC-III/posl_numerics_trunc48hrs/"))
vitals_path <- paste0(box, "MIMIC-III/posl_numerics_trunc48hrs/")

patients_keep_final <- read.csv(paste0(box, "MIMIC-III/output/patients_keep_final.csv"))
setDT(patients_keep_final)

covariates <- c(
  "HEIGHT", "WEIGHT", "BMI", "AGE", "INSURANCE_PRIV", "SEX_M",
  "ETHNICITY_WHITE", "ETHNICITY_NONWHITE", "FIRST_CAREUNIT_TSICU",
  "FIRST_CAREUNIT_SICU", "FIRST_CAREUNIT_MICU", "FIRST_CAREUNIT_CSRU",
  "MAP_longest_flat_spot", "MAP_min",
  "MAP_Q1", "MAP_median", "MAP_mean", "MAP_Q3", "MAP_max", "MAP_NA",
  "HR_median", "HR_mean", "HR_NA", "PULSE_median", "PULSE_mean",
  "PULSE_NA", "SYS_median", "SYS_mean", "SYS_NA", "DIAS_median",
  "DIAS_mean", "DIAS_NA", "RESP_median", "RESP_mean", "RESP_NA",
  "SPO2_median", "SPO2_mean", "SPO2_NA", "DA", "NE", "AVP", "EPI",
  "PE", "VENT", "MAP_history5_trend_strength", "MAP_history5_spikiness",
  "MAP_history5_linearity", "MAP_history5_curvature", "MAP_history5_stl_e_acf1",
  "MAP_history5_spectral_entropy", "MAP_history5_longest_flat_spot",
  "MAP_history5_coef_hurst", "MAP_history5_min", "MAP_history5_Q1",
  "MAP_history5_median", "MAP_history5_mean", "MAP_history5_Q3",
  "MAP_history5_max", "MAP_history5_NA", "MAP_history5_var", "HR_history5_median",
  "HR_history5_mean", "HR_history5_var", "HR_history5_NA", "PULSE_history5_median",
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

any_na_list <- list()
for(i in 1:nrow(patients_keep_final)){
  patient_i <- patients_keep_final[i,]
  load(paste0(
    box, "MIMIC-III/output/patients/", patient_i$SUBJECT_ID, "_",
    patient_i$ICUSTAY_ID, ".Rdata"
  ))
  data_i <- data_i[, c("SUBJECT_ID", "ICUSTAY_ID", "minute", covariates), with=F]
  data_i[, Y := shift(MAP_median, n = 5, type = "lag")]
  data_i <- data_i[-c(1:14), ]
  any_na_list[[i]] <- any(colSums(is.na(data_i)) > 0)
  save(
    data_i, compress = TRUE,
    file = paste0(box, "MIMIC-III/output/data/", patient_i$SUBJECT_ID, "_",
                  patient_i$ICUSTAY_ID, ".Rdata")
  )
}
any(any_na_list == TRUE)

# split subjects in half
set.seed(9528)
half1 <- sample(1:nrow(patients_keep_final), nrow(patients_keep_final)/2)
half2 <- seq(1:nrow(patients_keep_final))[-half1]
any(half1 %in% half2)
dput(as.numeric(half1))
c(50, 199, 48, 162, 66, 79, 93, 51, 184, 21, 196, 170, 23, 203,
  28, 152, 47, 1, 62, 99, 59, 60, 78, 149, 132, 14, 177, 68, 43,
  158, 84, 165, 88, 70, 151, 87, 166, 141, 56, 154, 16, 212, 115,
  91, 213, 181, 7, 37, 127, 26, 178, 33, 123, 32, 111, 126, 195,
  135, 117, 211, 94, 188, 133, 80, 194, 221, 143, 224, 5, 12, 98,
  138, 96, 57, 210, 3, 121, 171, 69, 164, 128, 217, 58, 142, 228,
  52, 40, 82, 49, 54, 53, 95, 144, 226, 153, 4, 22, 169, 76, 27,
  119, 6, 108, 218, 72, 175, 134, 145, 225, 61, 125, 187, 120,
  202, 85)
dput(as.numeric(half2))
c(2, 8, 9, 10, 11, 13, 15, 17, 18, 19, 20, 24, 25, 29, 30, 31,
  34, 35, 36, 38, 39, 41, 42, 44, 45, 46, 55, 63, 64, 65, 67, 71,
  73, 74, 75, 77, 81, 83, 86, 89, 90, 92, 97, 100, 101, 102, 103,
  104, 105, 106, 107, 109, 110, 112, 113, 114, 116, 118, 122, 124,
  129, 130, 131, 136, 137, 139, 140, 146, 147, 148, 150, 155, 156,
  157, 159, 160, 161, 163, 167, 168, 172, 173, 174, 176, 179, 180,
  182, 183, 185, 186, 189, 190, 191, 192, 193, 197, 198, 200, 201,
  204, 205, 206, 207, 208, 209, 214, 215, 216, 219, 220, 222, 223,
  227, 229, 230)
