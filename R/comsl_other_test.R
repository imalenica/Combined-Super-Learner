# how much time is saved

#devtools::install_github("tlverse/sl3@devel")
library(tidyverse)
library(data.table)
library(here)
library(ggplot2)
library(sl3)
library(origami)

source(here::here("R", "utils_mimic.R"))
load(here::here("Data", "data_summaryY.Rdata"))

dat_filtered <- filter(data, 
                       imputed_abpmean == 0 & imputed_abpsys_abpdias == 0)
dat_filtered$hour <- as.integer(dat_filtered$min_elapsed/60) + 1
length(unique(dat_filtered$subject_id)) # 678

############################ identify guinea pig subject #######################

bp_tbl_summary <- dat_filtered %>% 
  group_by(subject_id, hour) %>%
  summarize(med_abpmean = median(abpmean), 
            mean_abpmean = mean(abpmean))
# potential hypo subjects: 10250, 10564, 10061, 11638
sub_hypo <- filter(dat_filtered, subject_id == 11638)
ggplot(sub_hypo, aes(x = min_elapsed, y = abpmean)) +
  geom_line() +
  geom_hline(yintercept = 65, color = "red") +
  labs(x = "Minutes Elapsed", y = "Mean Arterial Blood Pressure", 
       main = "Subject ID 11638") 
  
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

# create historical data and create historical fit
set.seed(4197)
subject_sample <- sample(unique(dat_filtered$subject_id), 100)
subject_sample <- filter(dat_filtered, subject_id %in% subject_sample)

historical_fit <- make_historical_fit(
  historical_data = data.table(subject_sample), 
  outcome = "Y", 
  covariates = c("gender","age","care_unit", "admission_type_descr",
                 "sapsi_first", "sofa_first", "bmi", "rank_icu", "imputed_age",
                 "imputed_bmi", "imputed_sofa", "imputed_sapsi",
                 "amine", "sedation", "ventilation", "spo2", "hr",
                 "abpmean", "abpsys", "abpdias", "min_elapsed"), 
  id = "subject_id", 
  historical_stack = learners
  )

# create individual data, fit, and combine with historical fit
ind_data <- filter(dat_filtered, subject_id == 11638)
ind_data_init <- data.table(filter(ind_data, min_elapsed < 62))
individual_stack <- make_learner(Stack, xgb_learners)
fit_initial <- make_combined_sl(
  individual_data = ind_data_init, 
  outcome = "Y", 
  covariates = c("gender","age","care_unit", "admission_type_descr",
                 "sapsi_first", "sofa_first", "bmi", "rank_icu", "imputed_age",
                 "imputed_bmi", "imputed_sofa", "imputed_sapsi",
                 "amine", "sedation", "ventilation", "spo2", "hr",
                 "abpmean", "abpsys", "abpdias", "min_elapsed"), 
  historical_fit = historical_fit,
  subject_id = "11638",
  individual_stack = individual_stack, 
  cv_type = "folds_rolling_origin"
  )

ind_data2 <- data.table(filter(ind_data, min_elapsed < 72))
fit2 <- make_combined_sl(
  individual_data = ind_data2, 
  outcome = "Y", 
  covariates = c("gender","age","care_unit", "admission_type_descr",
                 "sapsi_first", "sofa_first", "bmi", "rank_icu", "imputed_age",
                 "imputed_bmi", "imputed_sofa", "imputed_sapsi",
                 "amine", "sedation", "ventilation", "spo2", "hr",
                 "abpmean", "abpsys", "abpdias", "min_elapsed"), 
  cv_type = "folds_rolling_origin", 
  historical_fit = historical_fit,
  subject_id = "11638", 
  past_individual_fit = fit_initial$individual_fit,
  past_cv_type = fit_initial$cv_type,
  past_sl_weights = fit_initial$sl_weights
)

ind_data3 <- data.table(filter(ind_data, min_elapsed < 82))
fit3 <- make_combined_sl(
  individual_data = ind_data3, 
  outcome = "Y", 
  covariates = c("gender","age","care_unit", "admission_type_descr",
                 "sapsi_first", "sofa_first", "bmi", "rank_icu", "imputed_age",
                 "imputed_bmi", "imputed_sofa", "imputed_sapsi",
                 "amine", "sedation", "ventilation", "spo2", "hr",
                 "abpmean", "abpsys", "abpdias", "min_elapsed"), 
  cv_type = "folds_rolling_origin", 
  historical_fit = historical_fit,
  subject_id = "11638", 
  past_individual_fit = fit2$individual_fit,
  past_cv_type = fit2$cv_type,
  past_sl_weights = fit2$sl_weights
)

fit2$sl_loss + fit3$sl_loss