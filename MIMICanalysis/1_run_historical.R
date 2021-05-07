######################## Fit historical stacks on Savio ########################
.libPaths("/global/scratch/rachelvphillips/R")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
options("sl3.verbose" = TRUE)
data_path <- "/global/scratch/rachelvphillips/symphony-data/"

### read args
args <- R.utils::commandArgs(T)
print(args)
outcome <- as.character(args[1])
smooth_type <- as.character(args[2])

outcome_type <- ifelse(grepl("AHE", outcome), "binomial", "continuous")

### load libraries, data, source
library(data.table)
library(origami)
library(sl3)

source(here::here("MIMICanalysis", "sl3_setup.R"))

historical_file <- paste0("historical_", smooth_type, ".Rdata")
load(paste0(data_path, historical_file))

### make learner stack 
stack <- make_historical_stack(outcome_type)

### make ML task
covs <- get_covariates(historical)
d <- data.table::copy(historical)
d[, weights := 1/.N, by = id]
d$normalized_weights <- d$weights/sum(d$weights)            
task <- make_sl3_Task(data=d, outcome=outcome, id="id", covariates=covs,
                      drop_missing_outcome=T, weights="normalized_weights")
rm(list=c("historical", "get_covariates", "d", "covs", "make_historical_stack"))

### train stack on task
set.seed(715)
t <- proc.time()
historical_fit <- stack$train(task)
historical_fit$is_trained
timer <- proc.time() - t
print(timer)

print(historical_fit$fit_object$learner_fits[[1]]$fit_object$selected)

### save
fit_path <- paste0(data_path, outcome, "_", historical_file)
save(historical_fit, file=fit_path, compress=T)
rm(list=c("task", "historical_fit"))
sessionInfo()