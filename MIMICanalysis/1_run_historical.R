# if (grepl("savio2", Sys.info()["nodename"])) {
#   .libPaths("/global/scratch/rachelvphillips/R")
#   Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
# }
# 
# ### read args
# args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE,
#                              defaults = list(file = "history30.Rdata",
#                                              outcome = "Y10"))

args <- list(file = "historical60_mean.Rdata", outcome = "Y10_lag5_mean")

### load libraries, data, source
library(data.table)
library(origami)
library(sl3)
source(here::here("MIMICanalysis","sl3_setup.R"))
# source(here::here("R", "setup.R"))
options(sl3.verbose = FALSE)

# prep
sl3_debug_mode()
stack <- make_historical_stack()
# d <- load_and_prep_historical_data(args$file)
load("~/Downloads/historical60_mean.Rdata")
d <- historical

covs <- get_covariates(d)
rm(get_covariates)
rm(make_historical_stack)

# task
task <- make_sl3_Task(data = d, outcome = args$outcome, id = "id", 
                      covariates = covs)
rm(d)
rm(covs)

# train
set.seed(715)
t <- proc.time()
fit <- suppressMessages(stack$train(task))
rm(task)
rm(stack)
fit$is_trained
fit$learner_fits[[1]]$learner_fits[[1]]$fit_object$selected
proc.time() - t

# save
fit_name <- paste0(args$outcome, "_", args$file)
fit_path <- paste0("/global/scratch/rachelvphillips/symphony-data/", fit_name)
save(fit, file = fit_path, compress = TRUE)