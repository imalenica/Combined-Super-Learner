#devtools::install_github("tlverse/sl3@timeseries-overhaul")
library(tidyverse)
library(data.table)
library(here)
library(ggplot2)
library(sl3)
library(origami)

source(here::here("R", "v3", "utils_mimic.R"))
source(here::here("R", "v3", "make_historical_sl.R"))

load(here::here("Data", "historical_data30.Rdata"))
load(here::here("Data", "historical_data60.Rdata"))

############################### set up learners ###############################
# regular learners
grid_params <- list(max_depth = c(3, 8), eta = c(0.05, 0.2), nrounds = c(50)) 
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
params_default <- list(nthread = getOption("sl.cores.learners", 1))
xgb_learners <- apply(grid, MARGIN = 1, function(params_tune) {
  do.call(Lrnr_xgboost$new, c(params_default, as.list(params_tune)))
})
names(xgb_learners) <- sprintf("xgboost_%d",seq_along(xgb_learners))

grid_params <- list(alpha = seq(0, 1, 0.1), nfolds = 5)
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
params_default <- list(nthread = getOption("sl.cores.learners", 1))
glmnet_learners <- apply(grid, MARGIN = 1, function(params_tune) {
  do.call(Lrnr_glmnet$new, c(params_default, as.list(params_tune)))
})
names(glmnet_learners) <- sprintf("glmnet_%d", seq_along(glmnet_learners))

lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- Lrnr_glm_fast$new()
lrnr_ranger <- make_learner(Lrnr_ranger)
lrnr_polspline <- make_learner(Lrnr_polspline)

# time series learners
lrnr_arima <- make_learner(Lrnr_arima)
lrnr_lstm10 <- make_learner(Lrnr_lstm, epochs = 500, batch_size = 10,
                            early_stop = TRUE)
lrnr_lstm1 <- make_learner(Lrnr_lstm, epochs = 500, batch_size = 1,
                           early_stop = TRUE)
lrnr_arima_strat <- Lrnr_multiple_ts$new(learner = lrnr_arima)
lrnr_lstm10_strat <- Lrnr_multiple_ts$new(learner = lrnr_lstm10)
lrnr_lstm1_strat <- Lrnr_multiple_ts$new(learner = lrnr_lstm1)

stack <- make_learner(Stack, unlist(list(xgb_learners, glmnet_learners, 
                                         lrnr_glm, lrnr_ranger, lrnr_polspline, 
                                         lrnr_arima_strat, lrnr_lstm10_strat, 
                                         lrnr_lstm1_strat), recursive = TRUE))

# screeners & pipelines
screener_cor <- make_learner(Lrnr_screener_corRank, rank = 20)
screener_rf <- make_learner(Lrnr_screener_randomForest, nVar = 20, ntree = 500)

screen_cor_pipe <- make_learner(Pipeline, screener_cor, stack)
screen_rf_pipe <- make_learner(Pipeline, screener_rf, stack)

# final stack
stack <- make_learner(Stack, screen_cor_pipe, screen_rf_pipe)

############################# fit historical data ##############################
notX <- c("subject_id","time_and_date","min_elapsed","Y_15","Y_20","Y_25","Y_30",
          "amine_switch1", "amine_switch0", "sedation_switch1", 
          "sedation_switch0", "ventilation_switch1", "ventilation_switch0")
X <- colnames(historical_dat)[-which(colnames(historical_dat) %in% notX)]

historical30_dat <- historical_data30 %>%
  group_by(subject_id) %>%
  slice(11:n())

historical30_fit <- make_historical_sl(
  historical_data = data.table(historical30_dat), 
  outcome = "Y", 
  covariates = X, 
  id = "subject_id", 
  historical_stack = learners,
  parallelize = TRUE, 
  cpus_logical = 23
)

historical60_dat <- historical_data60 %>%
  group_by(subject_id) %>%
  slice(11:n())

historical60_fit <- make_historical_sl(
  historical_data = data.table(historical60_dat), 
  outcome = "Y", 
  covariates = X, 
  id = "subject_id", 
  historical_stack = learners,
  parallelize = TRUE, 
  cpus_logical = 23
)
