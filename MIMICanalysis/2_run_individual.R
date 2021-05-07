######################## Fit individualized SL on Savio ########################
.libPaths("/global/scratch/rachelvphillips/R")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

################################## read args ##################################
args <- R.utils::commandArgs(T)
print(args)
ind_id <- as.character(args[1])
smooth_type <- as.character(args[2])
outcome_type <- as.character(args[3])

############################ load libraries, data, source ######################
library(data.table)
library(origami)
library(sl3)
library(here)

# commented out source for Savio, since it's available in this repo
# source(here::here("R", "setup.R"))
# source(here::here("R", "run_slstream.R"))
# source(here::here("R", "slstream.R"))
# source(here::here("R", "process_task.R"))
# source(here::here("R", "get_weights.R"))
source(here("MIMICanalysis", "sl3_setup.R")) # get sl3 libs, covariates
source(here("R", "v3", "run_adapt_sl.R")) # create data streams & continual updates
source(here("R", "v3", "make_adapt_sl.R")) # main function
source(here::here("R", "v3", "process_task.R")) # dependency for main function
source(here::here("R", "v3", "get_weights.R")) # dependency for main function

file_path <- "/global/scratch/rachelvphillips/symphony-data/"

cv_stack <- make_individual_cv_stack()
weights_control <- list(window = 120, delay_decay = 10, rate_decay = 0.01)
outcomes <- c("Y5", "Y10", "Y15", "Y20", "Y30")
if(outcome_type == "both"){
  outcomes <- c(paste0(outcomes, "_lag5_", smooth_type), paste0(outcomes, "_AHE"))
} else if (outcome_type == "binary"){
  outcomes <- paste0(outcomes, "_AHE")
} else if (outcome_type == "continuous"){
  outcomes <- paste0(outcomes, "_lag5_", smooth_type)
}
print(outcomes)

load(paste0(file_path, "individual_", smooth_type, ".Rdata"))
individual$id <- as.character(individual$id)
individual$time <- individual$min
data <- individual[id == ind_id,]
rm("individual")
setorder(data, min)

covs <- get_covariates(data)

run_slstream_allY(
  outcomes = outcomes, individual_data = data, individual_id = ind_id, 
  covariates = covs, cv_stack = cv_stack, parallel = FALSE,
  slstream_weights_control = weights_control, forecast_with_full_fit = TRUE,
  file_path = file_path, print_timers = FALSE
)