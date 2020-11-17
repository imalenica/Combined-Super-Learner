##### Savio
# if (grepl("savio2", Sys.info()["nodename"])) {
#   .libPaths("/global/scratch/rachelvphillips/R")
#   Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
# }
# 
# ################################## read args ##################################
# args <- R.utils::commandArgs(
#   trailingOnly=T, asValues=T, 
#   defaults=list()
# )

##### locally
args <- list(outcome="Y15_lag5_mean", 
             ids=c("12477", "16621", "20934", "23986_1", "26451", "26558", 
                   "2698", "29585_2", "3094", "8042_2", "791", "9275", "10520"))

# these ids have mean AHE > 0.5 for the duration of their stay 

############################ load libraries, data, source ######################
library(data.table)
library(origami)
library(sl3)
library(here)

##### Savio
# source(here::here("R", "make_adapt_sl.R"))
# source(here::here("R", "sl3_setup.R"))
# source(here::here("R", "process_task.R"))
# source(here::here("R", "get_weights.R"))
# file_path <- "/global/scratch/rachelvphillips/symphony-data/"

##### locally 
source(here::here("R", "v3", "make_adapt_sl.R"))
source(here::here("R", "v3", "run_adapt_sl.R"))
source(here::here("MIMICanalysis", "sl3_setup.R"))
source(here::here("R", "v3", "process_task.R"))
source(here::here("R", "v3", "get_weights.R"))
file_path <- "~/Downloads/"

get_all_outcomes <- function(outcome_gaps = c(10,15,20,25,30,35,40), smooth){
  sapply(outcome_gaps, function(gap) paste0("Y", gap+5, "_lag5_", smooth))
}

if(grepl("mean", args$outcome)){
  smooth <- "mean"
} else if(grepl("median", args$outcome)){
  smooth <- "median"
}

outcomes <- get_all_outcomes(smooth = smooth)

load(paste0(file_path, "individual_", smooth, ".Rdata"))
individual$id <- as.character(individual$id)
individual <- individual[id %in% args$ids,]

covs <- get_covariates(individual)
covs

cv_stack <- get_cv_stack_individual()
weights_control <- list(window = 120, delay_decay = 10, rate_decay = 0.01)

outcome_specific_results <- run_slstream_allID(
  multi_individual_data = individual, ids = args$ids, covariates = covs, 
  outcome = args$outcome, file_path = file_path, cv_stack = cv_stack, 
  slstream_weights_control = weights_control, cores = args$cores
)
results_path <- paste0(args$outcome, ".Rdata")
save(outcome_specific_results, file=paste0(file_path, results_path), compress=T)
