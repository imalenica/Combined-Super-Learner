################################### read args ##################################
args <- R.utils::commandArgs(trailingOnly = TRUE, 
                             asValues = TRUE,
                             defaults = list(
                               outcome = "Y15",
                               historical_data = "historical30.Rdata",
                               individual_data = "individual30.Rdata",
                               individual_id = "791",
                               outcome_locf = "abpmean_locf",
                               outcome_gap = 15
                               ))
args <- list(outcome = "Y15", historical_data = "historical30.Rdata",
             individual_data = "individual30.Rdata", individual_id = "791")

############################ load libraries, data, source ######################
library(data.table)
library(origami)
library(sl3)
library(earth)

# source(here::here("R", "setup.R"))
# source(here::here("R", "sl_funs.R"))
source(here::here("R", "v3", "2_sl3_setup.R"))
source(here::here("R", "v3", "make_adapt_and_online_sl.R"))
source(here::here("R", "v3", "utils_sl.R"))
source(here::here("R", "v3", "process_task.R"))
source(here::here("R", "v3", "get_weights.R"))
# data_path <- "/global/scratch/rachelvphillips/symphony-data/"
data_path <- "~/Downloads/"
# data_path <- "~/symphony/data/"
options(sl3.verbose = FALSE)

############################ load & subset individual data #####################
individual_data_path <- paste0(data_path, args$individual_data)
load(individual_data_path)

# subset to subject of interest
individual$id <- as.character(individual$id)
d <- individual[id == args$individual_id,]
rm(individual)
setorder(d, time_and_date)

############################## load historical fit #############################
historical_fit_path <- paste0(data_path, args$outcome, "_", args$historical_data)
load(historical_fit_path)
historical_fit <- fit # rename
rm(fit)

############################### retain covariates ##############################
W <- c("rank_icu", "sex", "age", "care_unit", "sapsi_first", "sofa_first", 
       "bmi", "admission_type_descr")
if(any(is.na(d[, W, with=FALSE]))){
  stop("NA baseline covariates")
}
covs <- get_covariates(d)
rm(get_covariates)

########################### set up pseudo-online data ##########################
max_stop_time <- 1440
initial_burn_in <- 30
batch <- 5

max_times <- c(max_stop_time, nrow(d))
max_time <- max_times[which.min(max_times)]
if(max_time < 100) {
  stop(paste0("Insufficient time for subject :", args$individual_id))
}

splits <- seq((initial_burn_in+batch), max_time, batch)
training_data_list <- lapply(splits, function(x) data.table(d[1:x, ]))
forecast_data_list <- lapply(splits, function(x) data.table(d[(x+1):(x+batch), ]))

# function to determine new cv stacks based on sample size
get_cv_stack_group <- function(n) ifelse(n<50, 1, ifelse(n>=50 && n<200, 2, 3))

# function to increase burn_in (i.e., initial training size) w complex cv stacks
get_burn_in <- function(n) ifelse(n<200, initial_burn_in, initial_burn_in*5)

################################## run #########################################

loss_tables <- list()
onlineSL_tables <- list()
onlineSL_forecasts <- list()
individual_fits <- list()
all_forecasts <- list()

sl3_debug_mode()
# N <- length(splits)-1
N <- 2
start_time <- proc.time()
for(i in 1:N){
  
  # quick clean of results we no longer need
  if(i > 2){
    loss_tables[[i-2]] <- NA
    onlineSL_tables[[i-2]] <- NA
    onlineSL_forecasts[[i-2]] <- NA
    individual_fits[[i-2]] <- NA
  }
  
  # initiate "online" training and forecast data
  training_data <- data.table(training_data_list[[i]])
  forecast_data <- data.table(forecast_data_list[[i]])
  burn_in <- get_burn_in(splits[[i]])

  if(i == 1){
    cv_stack <- make_individual_cv_stack(n = splits[[i]])
    fit <- invisible(make_adapt_and_online_sl(
      individual_training_data = training_data, outcome = args$outcome, 
      id = args$individual_id, covariates = covs, 
      historical_fit = historical_fit, batch = batch, first_window = burn_in, 
      individual_stack = cv_stack, individual_forecast_data = forecast_data, 
      first_fit = TRUE
    ))
  } else {
    # don't make new stack unless we're jumping to new cv_stack group
    previous_group <- get_cv_stack_group(splits[[i-1]])
    current_group <- get_cv_stack_group(splits[[i]])
    if(previous_group == current_group) {
      cv_stack <- NULL
    } else {
      cv_stack <- make_individual_cv_stack(n = splits[[i]])
    }
    fit <- invisible(make_adapt_and_online_sl(
      individual_training_data = training_data, outcome = args$outcome, 
      id = args$individual_id, covariates = covs, batch = batch, 
      first_window = burn_in, historical_fit = historical_fit, 
      individual_forecast_data = forecast_data, individual_stack = cv_stack, 
      first_fit = FALSE, individual_fit = individual_fits[[i-1]], 
      loss_table = loss_tables[[i-1]], onlineSL_table = onlineSL_tables[[i-1]],
      previous_onlineSL_forecasts = onlineSL_forecasts[[i-1]]
    ))
  }
  loss_tables[[i]] <- fit$loss_table
  onlineSL_tables[[i]] <- fit$onlineSL_table
  onlineSL_forecasts[[i]] <- fit$onlineSL_forecasts
  individual_fits[[i]] <- fit$individual_fit
  all_forecasts[[i]] <- fit$all_forecasts
}
end_time <- proc.time() - start_time

# quick summary of the individual
d <- training_data_list[[N]]
summary_table <- data.table(
  subject_id = unique(d$subject_id), icustay_id = unique(d$icustay_id), 
  id = unique(d$id), abpmean_mean = mean(d$abpmean), amine_mean = mean(d$amine), 
  sex = unique(d$sex), age = unique(d$age), bmi = unique(d$bmi), 
  sapsi_first = unique(d$sapsi_first), sofa_first = unique(d$sofa_first), 
  care_unit = unique(d$care_unit), 
  admission_type_descr = unique(d$admission_type_desc)
)

locf <- ifelse(d[[args$outcome_locf]] == 1|| d[["row_locf"]] == 1, 1, 0)
locf <- locf[args$outcome_gap:nrow(d)]
locf <- c(locf, rep(NA, nrow(d)-length(locf)))
Y_tbl <- data.table(truth = d[[args$outcome]], locf = locf, time = d[["min"]])
pred_tbl <- rbindlist(all_forecasts, fill=T)
forecast_tbl <- data.table(merge(Y_tbl, pred_tbl, all.x=F, all.y=T, by="time"))

results <- list(
  id_summary = summary_table,
  onlineSL_table = onlineSL_tables[[N]],
  loss_table = loss_tables[[N]],
  forecast_table = forecast_tbl
)

################################## save ########################################
path <- paste0(args$outcome, "_id", args$individual_id, ".Rdata")
save(results, file=paste0(data_path, path), compress=T)