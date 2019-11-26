library(here)
library(tidyverse)
library(data.table)
library(sl3)
library(origami)
library(SuperLearner)

source(here::here("R", "utils_mimic.R"))
load(here::here("Data","mimic_all.Rdata"))

################################################################################
# only include subjects with 5 hours of data at most 1 min time gap between
# two consecutive measurements 
# having no gaps is necessary for pooled cross validation scheme, so each 
# subject has same number of rows of observations
################################################################################

# how many subjects with 5 hours of data have no more than 1 min gap between 
# consecutive observations? 
eval_missingness(min = 1, dataset = mimic_all, total_hrs = 8)[["num"]] 
# 682
dat <- eval_missingness(min = 1, dataset = mimic_all, total_hrs = 8)[["dat"]]

################################################################################
# prep for combined SL
################################################################################

# adding time column
dat <- dat %>%
  dplyr::group_by(subject_id) %>%
  dplyr::mutate(time = seq(1,n(), 1))

# outcome needs to be numeric for correlation screener to work
dat$hypo_event <- as.numeric(dat$hypo_event) - 1

outcome <- "hypo_event"

covars_baseline <- c("gender","age","care_unit", "admission_type_descr", 
                     "sapsi_first", "sofa_first", "bmi", "rank_icu",
                     "imputed_age", "imputed_bmi", "imputed_sofa", 
                     "imputed_sapsi")

covars_timevarying <- c("amine", "sedation", "ventilation", "spo2", "hr", 
                        "abpmean", "imputed_abpmean")

grid_params = list(max_depth=c(4,8),
                   eta = c(0.001, 0.01, 0.1, 0.2))
grid = expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
params_default = list(nthread = getOption("sl.cores.learners", 1))
xgb_learners = apply(grid, MARGIN = 1, function(params_tune) {
  do.call(Lrnr_xgboost$new, c(params_default, as.list(params_tune)))
})

stack <- make_learner(
  Stack, xgb_learners[[1]], xgb_learners[[2]], xgb_learners[[3]],
  xgb_learners[[4]], xgb_learners[[5]], xgb_learners[[6]], xgb_learners[[7]], 
  xgb_learners[[8]])

screen_cor <- Lrnr_pkg_SuperLearner_screener$new("screen.corP")
cor_pipeline <- make_learner(Pipeline, screen_cor, stack)
stack_screen <- make_learner(Stack, cor_pipeline, stack)

metalearner <- make_learner(Lrnr_nnls)
sl <- Lrnr_sl$new(learners = stack, metalearner = metalearner)

################################################################################
# run combined SL with varying cv schemes and training times
################################################################################

# we need at least gap + h time in test set, so we do not train on full 8 hrs
# we need at least gap + test_size time in training set, so we do not include 
# 30 min as a training time
times <- seq(30, 420, by = 30)

ptw <- proc.time()
comSL_cvrw1 <- lapply(times[1:7], function(x){
  combine_SL2(train_all = dat, outcome, t = x, stack_pool = stack, 
             stack_individual = stack, stack_screen = stack_screen, sl = sl,
             covars = covars_timevarying, covars_baseline = covars_baseline,
             cv = "folds_rolling_window", gap = 30, h = 15, 
             test_size = 15, mini_batch = 5, window_size = 15)
  })
names(comSL_cvrw1) <- paste0(times[1:7], " min training time")
ptw - proc.time()
save(comSL_cvrw1, file = here::here("Results", "comSL_cvrw1.Rdata"), 
     compress = TRUE)

ptw <- proc.time()
comSL_cvrw2 <- lapply(times[7:11], function(x){
  combine_SL2(train_all = dat, outcome, t = x, stack_pool = stack, 
              stack_individual = stack, stack_screen = stack_screen, sl = sl,
              covars = covars_timevarying, covars_baseline = covars_baseline,
              cv = "folds_rolling_window", gap = 30, h = 15, 
              test_size = 15, mini_batch = 5, window_size = 15)})
names(comSL_cvrw2) <- paste0(times[7:11], " min training time")
save(comSL_cvrw2, file = here::here("Results", "comSL_cvrw2.Rdata"), 
     compress = TRUE)
comSL_cvrw3 <- lapply(times[12:14], function(x){
  combine_SL2(train_all = dat, outcome, t = x, stack_pool = stack, 
              stack_individual = stack, stack_screen = stack_screen, sl = sl,
              covars = covars_timevarying, covars_baseline = covars_baseline,
              cv = "folds_rolling_window", gap = 30, h = 15, 
              test_size = 15, mini_batch = 5, window_size = 15)})
names(comSL_cvrw3) <- paste0(times[12:14], " min training time")
save(comSL_cvrw3, file = here::here("Results", "comSL_cvrw3.Rdata"), 
     compress = TRUE)
proc.time() - ptw
comSL_cvro2 <- lapply(times[7:14], function(x){
  combine_SL2(train_all = dat, outcome, t = x, stack_pool = stack, 
              stack_individual = stack, stack_screen = stack_screen, sl = sl,
              covars = covars_timevarying, covars_baseline = covars_baseline,
              cv = "folds_rolling_origin", gap = 30, h = 15, 
              test_size = 15, mini_batch = 5, window_size = 15)})
names(comSL_cvro2) <- paste0(times, " min training time")
save(comSL_cvro2, file = here::here("Results", "comSL_cvro2.Rdata"), 
     compress = TRUE)
proc.time() - ptw

comSL_cvro <- lapply(times, function(x){
  combine_SL(train_all = dat, outcome, t = x, stack_pool = stack, 
             stack_individual = stack, stack_screen = stack_screen, sl = sl,
             covars = covars_timevarying, covars_baseline = covars_baseline,
             cv = "folds_rolling_origin", gap = 30, h = 15, 
             test_size = 15, mini_batch = 5, window_size = 15)
})
names(comSL_cvro) <- paste0(times, " min training time")
save(comSL_cvro, file = here::here("Results", "comSL_cvro.Rdata"), 
     compress = TRUE)