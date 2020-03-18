# Test dependence with history included

library(tidyverse)
library(data.table)
library(here)
library(ggplot2)
library(sl3)
library(origami)

source(here::here("R", "utils_mimic.R"))
load(here::here("Data", "fin_history30.Rdata"))
load(here::here("Data", "fin_history60.Rdata"))

data_30<-fin_history30
data_60<-fin_history60

############################### set up learners ###############################

grid_params = list(max_depth = c(2,5,8),
                   eta = c(0.005, 0.1, 0.25))
grid = expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
params_default = list(nthread = getOption("sl.cores.learners", 1))
xgb_learners = apply(grid, MARGIN = 1, function(params_tune) {
  do.call(Lrnr_xgboost$new, c(params_default, as.list(params_tune)))
})
lrnr_lasso <- make_learner(Lrnr_glmnet, alpha = 1)
lrnr_glm <- make_learner(Lrnr_glm)
learners <- make_learner(Stack, unlist(list(xgb_learners, lrnr_glm, lrnr_lasso),
                                       recursive = TRUE))

metalearner <- make_learner(Lrnr_nnls)
sl <- Lrnr_sl$new(learners = learners, metalearner = metalearner)

lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm_simple <- make_learner(Lrnr_glm_fast)
individual_stack <- make_learner(Stack, lrnr_mean, lrnr_glm_simple)

#########################################################################
########        New Version of the Adaptive Super Learner        ########
#########################################################################

################# Learn the historical fit ###############################
set.seed(4197)

#Randomly sample
ids <- unique(data_30$subject_id)
subject_sample <- sample(ids, 40)
subject_sample_30 <- filter(data_30, subject_id %in% subject_sample)
subject_sample_60 <- filter(data_60, subject_id %in% subject_sample)
ind_sample <- sample(ids[!(ids %in% subject_sample)], 1)
ind_sample %in% subject_sample

all_covs = c("min_elapsed","abpsys","abpsys_locf",
             "abpdias","abpdias_locf","abpmean","abpmean_locf","spo2","spo2_locf","amine","sedation",
             "ventilation","rank_icu","sex","age","sapsi_first","sofa_first","bmi","care_unit",
             "admission_type_descr","amine_switch1","amine_switch0",
             "sedation_switch1","sedation_switch0","ventilation_switch1","ventilation_switch0",
             "amine_history_mean","amine_history_switch1_time",
             "sedation_history_mean","sedation_history_switch1_time",
             "ventilation_history_mean","ventilation_history_switch1_time",
             "abpsys_Min","abpsys_Mean","abpsys_Max",
             "abpsys_trend_strength","abpsys_spikiness","abpsys_linearity","abpsys_curvature",
             "abpsys_spectral_entropy","abpsys_n_flat_spots","abpsys_coef_hurst",
             "abpdias_Min","abpdias_Mean","abpdias_Max",
             "abpdias_trend_strength","abpdias_spikiness","abpdias_linearity","abpdias_curvature",
             "abpdias_spectral_entropy","abpdias_n_flat_spots","abpdias_coef_hurst",
             "abpmean_Min","abpmean_Mean","abpmean_Max","abpmean_trend_strength","abpmean_spikiness",
             "abpmean_linearity","abpmean_curvature","abpmean_spectral_entropy",
             "abpmean_n_flat_spots","abpmean_coef_hurst",
             "spo2_Min","spo2_Mean","spo2_Max","spo2_trend_strength","spo2_spikiness",
             "spo2_linearity","spo2_curvature","spo2_spectral_entropy","spo2_n_flat_spots","spo2_coef_hurst")

#Create historical fit (here, use xgboost)
historical_fit <- make_historical_fit(
  historical_data = data.table(subject_sample_30), 
  outcome = "Y_15", 
  covariates = all_covs,
  id = "subject_id", 
  historical_stack = learners
)

################# Learn the individual fit ###############################

#Chose the one sample we were interested in
ind_data <- filter(data_30, subject_id == ind_sample)
#Split data into 5 row chuncks
splits <- seq(30,300,60)
split_data <- lapply(splits, function(x) data.table(ind_data[1:x,]))

#Save all results
result_list <- list()

#Fit initial iteration of the individual learner
result_list[[1]] <- make_adapt_sl(
  individual_training_data = split_data[[1]], 
  indiviual_forecast_data = ind_data[c((nrow(split_data[[1]])+1):(nrow(split_data[[1]])+5)),],
  outcome = "Y_15", 
  covariates = all_covs, 
  subject_id = ind_sample,
  historical_fit = historical_fit,
  individual_stack = individual_stack
)

result_list[[2]] <- make_adapt_sl(
  individual_training_data = split_data[[2]], 
  indiviual_forecast_data = ind_data[c((nrow(split_data[[1]])+1):(nrow(split_data[[1]])+5)),],
  outcome = "Y_15", 
  covariates = all_covs, 
  subject_id = ind_sample,
  historical_fit = historical_fit,
  individual_stack = individual_stack
)



