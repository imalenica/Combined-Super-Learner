library(data.table)
library(origami)
library(sl3)
library(here)

### read args
half_idx <- 1

source("~/symphony/R/slstream.R")
source("~/symphony/R/run_slstream.R")
source("~/symphony/R/setup.R")
source("~/symphony/R/process_task.R")
source("~/symphony/R/get_weights.R")

cv_stack <- make_individual_stack(online = TRUE)

files_loc <-"/Users/Rachael/Library/CloudStorage/Box-Box/MIMIC-III/output/"
patients <- read.csv(paste0(files_loc, "patients_keep_final.csv"))

if(half_idx == 1){
  individual_patients <- patients[half1,]
  load(paste0(files_loc, "fits-savio/historical_half2_pooled.Rdata"))
} else if (half_idx == 2){
  individual_patients <- patients[half2,]
  load(paste0(files_loc, "fits-savio/historical_half1_pooled.Rdata"))
}

slstream_weights_control <- list(
  "wts.all1" = list(window = NULL, delay_decay = NULL, rate_decay = NULL),
  "wts.broad" = list(window = NULL, delay_decay = 90, rate_decay = 0.001),
  "wts.med" = list(window = 120, delay_decay = 30, rate_decay = 0.001),
  "wts.narrow" = list(window = 60, delay_decay = 10, rate_decay = 0.005),
  "wts.narrowest" = list(window = 30, delay_decay = 10, rate_decay = 0.01)
)

i <- 1
patient_i <- individual_patients[i,]
patient_i_loc <- paste0(patient_i$SUBJECT_ID, "_", patient_i$ICUSTAY_ID)
load(paste0(files_loc, "data/", patient_i_loc, ".Rdata"))
print(paste0("Running patient-specific individual fit for patient ", i))
start_time <- proc.time()
set.seed(715)
result_i <- run_slstream(
  individual_data = data_i, individual_id = patient_i$SUBJECT_ID,
  outcome = "Y", covariates = covariates, cv_stack = cv_stack,
  slstream_weights_control = slstream_weights_control,
  historical_fit = historical_fit_pooled, forecast_with_full_fit = TRUE,
  horizon = 5, outcome_bounds = c(10, 160), batch = 1,
  time = "minute", individual_burn_in = 60, max_stop_time = 1440
)
end_time <- proc.time()
end_time - start_time
