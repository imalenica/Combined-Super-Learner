library(here)
library(tidyverse)
library(data.table)
library(sl3)
library(origami)
library(SuperLearner)

source(here::here("R", "utils_mimic.R"))
load(here::here("Data","mimic_all.Rdata"))

################################################################################
# only include subjects with 8 hours of data at most 3 min time gap between
# two consecutive measurements 
################################################################################

# how many subjects with 8 hours of data have no more than 3 min gap between 
# consecutive observations? 
eval_missingness(min = 3, dataset = mimic_all, total_hrs = 8)[["num"]] 
# 786
dat <- eval_missingness(min = 3, dataset = mimic_all, total_hrs = 8)[["dat"]] 

# out of these 786 subjects how many times do they having 2 minute or greater 
# gaps of missingness?  
eval <- dat
eval$tdiff <- unlist(tapply(eval$time_and_date, 
                            INDEX = eval$subject_id,
                            FUN = function(x) c(0, diff(as.numeric(x)))))
missing <- eval %>% dplyr::filter(tdiff > 60)
gaps_by_id <- missing %>%
  dplyr::select("subject_id") %>%
  dplyr::group_by(subject_id) %>%
  tally()
gaps_by_id$subject_id <- as.numeric(as.character(gaps_by_id$subject_id))
subset(gaps_by_id, rowSums(gaps_by_id) - gaps_by_id$subject_id > 2)
# not often, 3 times is the highest

################################################################################
# prep for combined SL
################################################################################

# need "time" column 
colnames(dat)[3] <- "time"

id <- "subject_id"

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
# run combined SL with varying cv schemes and gaps
################################################################################

# we need at least gap + h time in test set, so we do not train on full 8 hrs
# we need at least gap + test_size time in training set, so we do not include 30 min as a training time
times <- seq(30, 420, by = 30)[-1]

comSL_cvrw <- lapply(times, function(x){
  combine_SL(train_all = dat, outcome, t = x, stack_pool = stack, 
             stack_individual = stack, stack_screen = stack_screen, sl = sl,
             covars = covars_timevarying, covars_baseline = covars_baseline,
             id = id, cv = "folds_rolling_window", gap = 30, h = 15, 
             test_size = 15, mini_batch = 5, window_size = 15)
  })
names(comSL_cvrw) <- paste0(times, " min training time")
save(comSL_cvrw, file = here::here("Results", "comSL_cvrw.Rdata"), 
     compress = TRUE)

comSL_cvro <- lapply(times, function(x){
  combine_SL(train_all = dat, outcome, t = x, stack_pool = stack, 
             stack_individual = stack, stack_screen = stack_screen, sl = sl,
             covars = covars_timevarying, covars_baseline = covars_baseline,
             id = id, cv = "folds_rolling_origin", gap = 30, h = 15, 
             test_size = 15, mini_batch = 5, window_size = 15)
})
names(comSL_cvro) <- paste0(times, " min training time")
save(comSL_cvro, file = here::here("Results", "comSL_cvro.Rdata"), 
     compress = TRUE)