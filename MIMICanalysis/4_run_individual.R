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
args <- list(outcome="Y15_lag5_mean", ids=TBD)

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

options(sl3.verbose = FALSE)


if(grepl("mean", args$outcome)){
  smooth <- "mean"
} else if(grepl("median", args$outcome)){
  smooth <- "median"
}

get_all_outcomes <- function(outcome_gaps = c(10,15,20,25,30,35,40), smooth){
  sapply(outcome_gaps, function(gap) paste0("Y", gap+5, "_lag5_", smooth))
}
outcomes <- get_all_outcomes(smooth = smooth)

load(paste0(file_path, "individual_", smooth, ".Rdata"))
individual$id <- as.character(individual$id)
individual <- individual[id %in% args$ids,]

covs <- get_covariates(individual)
cv_stack <- get_cv_stack_individual()
weights_control <- list(window = 120, delay_decay = 10, rate_decay = 0.01)

plot_decay <- function(weights_control, max_time = 240){
  
  weights <- rep(1, max_time)
  times <- seq(from = max_time, to = 1)
  lags <- max_time - times
    
  if (!is.null(weights_control$window)) {
    window <- max_time - weights_control$window
    weights <- weights * ifelse(times <= window, 0, 1)
  } else {
    weights_control$window <- "NULL"
  }
  
  if (!is.null(weights_control$rate_decay)) {
    lags <- max_time - times
    if (!is.null(weights_control$delay_decay)) {
      lags_delayed <- lags - weights_control$delay_decay
      lags <- ifelse(lags_delayed < 0, 0, lags_delayed)
    } else {
      weights_control$delay_decay <- "NULL"
    }
    weights <- weights * (1 - weights_control$rate_decay)^lags
  } else {
    weights_control$rate_decay <- "NULL"
  }
  title <- paste0("Decay in time wrt window=", weights_control$window, 
                  ", delay_decay=", weights_control$delay_decay,
                  ", rate_decay=", weights_control$rate_decay)
  plot(weights, xlab = "Distance from current observation (minutes)", 
       ylab = "Weight applied to learner losses", main = title, type = "l")
}


outcome_specific_results <- run_slstream_allID(
  multi_individual_data = individual, ids = args$ids, covariates = covs, 
  outcome = args$outcome, file_path = file_path, cv_stack = cv_stack, 
  slstream_weights_control = weights_control, cores = NULL
)
results_path <- paste0(args$outcome, ".Rdata")
save(outcome_specific_results, file=paste0(file_path, results_path), compress=T)
