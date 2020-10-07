get_covariates <- function(data){
  outcomes <- colnames(data)[grepl("Y", colnames(data))]
  not_covs <- c(outcomes, "time_and_date", "init_time_and_date", "subject_id", 
                "icustay_id", "id")
  colnames(data)[-which(colnames(data) %in% not_covs)]
}

args <- list(outcome = "Y15", individual_data = "id791_mimic.Rdata")

############################ load libraries, data, source ######################
library(data.table)
library(origami)
library(sl3)

# data_path <- "/global/scratch/rachelvphillips/symphony-data/"
data_path <- "~/Downloads/"
# data_path <- "~/symphony/data/"
options(sl3.verbose = FALSE)

############################ load & subset individual data #####################
individual_data_path <- paste0(data_path, args$individual_data)
load(individual_data_path)
setorder(d, time_and_date)
# save(d, file=paste0(data_path,"id791_mimic.Rdata"), compress=T)

covs <- get_covariates(d)

########################### set up pseudo-online data ##########################
max_stop_time <- 1440
initial_burn_in <- 30
batch <- 5

max_times <- c(max_stop_time, nrow(d))
max_time <- max_times[which.min(max_times)]
splits <- seq((initial_burn_in+batch), max_time, batch)
training_data_list <- lapply(splits, function(x) data.table(d[1:x, ]))
forecast_data_list <- lapply(splits, function(x) data.table(d[(x+1):(x+batch), ]))

# initiate "online" training and forecast data
i <- 1
training_data <- data.table(training_data_list[[i]])
forecast_data <- data.table(forecast_data_list[[i]])
burn_in <- get_burn_in(splits[[i]])

folds <- make_folds(training_data, fold_fun = folds_rolling_origin, gap = 0,
                    first_window = burn_in, validation_size = batch, 
                    batch = batch)
task <- make_sl3_Task(data = training_data, covariates = covs, time = "min",
                      folds = folds, outcome = args$outcome)

# some timeseries learners
arima_aicc <- Lrnr_arima$new(stepwise=FALSE, approximation=FALSE, ic="aicc")
arima_bic <- Lrnr_arima$new(stepwise=FALSE, approximation=FALSE, ic="bic")
arima_stack <- make_learner(Stack, arima_aicc, arima_bic)
lasso <- Lrnr_glmnet$new()
lasso_screener <- Lrnr_screener_coefs$new(lasso, threshold=0)
lasso_pipe <- make_learner(Pipeline, lasso_screener, arima_stack)

nlts_lstar <- Lrnr_tsDyn$new(learner="lstar")
nlts_setar <- Lrnr_tsDyn$new(learner="setar", model="MTAR")
lstm_relu <- Lrnr_lstm$new(activation="relu", early_stopping=T, patience=10)
lstm_selu <- Lrnr_lstm$new(activation="selu", early_stopping=T, patience=10)
ts_stack <- make_learner(Stack, lasso_pipe, nlts_lstar, nlts_setar, lstm_relu,
                         lstm_selu)

fit <- ts_stack$train(task)

                         