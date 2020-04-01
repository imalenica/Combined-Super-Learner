# how much time is saved

#devtools::install_github("tlverse/sl3@devel")
library(tidyverse)
library(data.table)
library(here)
library(ggplot2)
library(sl3)
library(origami)

#source(here::here("R", "utils_mimic.R"))
#load(here::here("Data", "data_summaryY.Rdata"))

dat_filtered <- filter(data, 
                       imputed_abpmean == 0 & imputed_abpsys_abpdias == 0)
dat_filtered$hour <- as.integer(dat_filtered$min_elapsed/60) + 1
length(unique(dat_filtered$subject_id)) # 678
dat_filtered$time <- min_elapsed+30
############################ identify guinea pig subject #######################

bp_tbl_summary <- dat_filtered %>% 
  group_by(subject_id, hour) %>%
  summarize(med_abpmean = median(abpmean), 
            mean_abpmean = mean(abpmean))
# potential hypo subjects: 10250, 10564, 10061, 11638
sub_hypo <- filter(dat_filtered, subject_id == 11638)
ggplot(sub_hypo, aes(x = min_elapsed, y = abpmean)) +
  geom_line() +
  #geom_hline(yintercept = 65, color = "red") +
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
grid_params <- list(covar1 = c("amine", "sedation", "ventilation", "spo2", "hr",
                               "abpmean", "abpsys", "abpdias", "min_elapsed"),
                    covar2 = c("amine", "sedation", "ventilation", "spo2", "hr",
                               "abpmean", "abpsys", "abpdias", "min_elapsed"))
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
grid$dup <- ifelse(grid$covar1 == grid$covar2, 1, 0)
grid_unique <- filter(grid, dup != 1)[,-3]

formulas_2terms <- list()
for(i in 1:nrow(grid_unique)){
  formulas_2terms[[i]] <- paste("Y ~", 
                                paste(c(as.character(grid_unique[i,1]),
                                        as.character(grid_unique[i,2])), 
                                      collapse= "+"))
}

term <- c("amine", "sedation", "ventilation", "spo2", "hr", "abpmean", 
          "abpsys", "abpdias")
formulas_1term <- list()
for(i in 1:length(term)){
  formulas_1term[[i]] <- paste("Y~", term[i])
}

formula <- c(formulas_1term, formulas_2terms)

glm_learners <- lapply(formula, function(params_tune) {
  listpars <- list(as.character(params_tune))
  names(listpars) <- "formula"
  do.call(Lrnr_glm$new, c(params_default, listpars))
})

individual_stack <- make_learner(Stack, glm_learners)

ind_data <- filter(dat_filtered, subject_id == 11638)
splits <- seq(10,400,5)
split_data <- lapply(splits, function(x) data.table(ind_data[1:x,]))

indiviual_forecast_data$Y <- rep(NA,5)
incoming_data <- rbind(individual_training_data,indiviual_forecast_data)
incoming_data$Y_time <- incoming_data$min_elapsed + 30
incoming_data <- incoming_data[,-c(28:29,31:35)]
all_covs <- c("gender","age","care_unit", "admission_type_descr",
              "sapsi_first", "sofa_first", "bmi", "rank_icu", "imputed_age",
              "imputed_bmi", "imputed_sofa", "imputed_sapsi", "amine", 
              "sedation", "ventilation", "spo2", "hr", "abpmean", "abpsys",
              "abpdias", "min_elapsed")

lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm_simple <- make_learner(Lrnr_glm, formula = "Y ~ abpmean")
lrnr_glm_simple2 <- make_learner(Lrnr_glm, formula = "Y ~ abpmean + hr")
lrnr_glm_simple3 <- make_learner(Lrnr_glm, formula = "Y ~ abpmean + hr + spo2")
lrnr_lasso <- make_learner(Lrnr_glmnet, alpha = 1)
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_xgboost <- make_learner(Lrnr_xgboost)
lrnr_ranger <- make_learner(Lrnr_ranger)

result_list <- list()
i <- 1
while(i < 50){
  if(i == 1){
    individual_stack <- make_learner(Stack, lrnr_mean, lrnr_glm_simple)
    result_list[[i]] <- make_adapt_sl(
      individual_training_data = split_data[[i]], 
      indiviual_forecast_data = ind_data[
        c((nrow(split_data[[i]])+1):(nrow(split_data[[i]])+5)),],
      outcome = "Y", 
      covariates = all_covs, 
      subject_id = "11638",
      historical_fit = historical_fit,
      individual_stack = individual_stack
    )
  } else if(i %in% 2:10) {
    if(i == 2){
      individual_stack <- make_learner(Stack, lrnr_mean, lrnr_glm_simple, 
                                       lrnr_glm_simple2, lrnr_glm_simple3, 
                                       lrnr_lasso, lrnr_glm)
    } else {
      individual_stack <- NULL
    }
    result_list[[i]] <- make_adapt_sl(
      individual_training_data = split_data[[i]], 
      indiviual_forecast_data = ind_data[
        c((nrow(split_data[[i]])+1):(nrow(split_data[[i]])+5)),],
      outcome = "Y", 
      covariates = all_covs, 
      historical_fit = historical_fit,
      subject_id = "11638", 
      past_individual_fit = result_list[[i-1]]$individual_fit,
      individual_stack = individual_stack
    )
  } else if(i %in% 11:30) {
    if(i == 11){
    individual_stack <- make_learner(Stack, lrnr_mean, lrnr_glm_simple, 
                                     lrnr_glm_simple2, lrnr_glm_simple3, 
                                     lrnr_lasso, lrnr_glm, lrnr_xgboost, 
                                     lrnr_ranger)
    } else {
      individual_stack <- NULL
    }
    result_list[[i]] <- make_adapt_sl(
      individual_training_data = split_data[[i]], 
      indiviual_forecast_data = ind_data[
        c((nrow(split_data[[i]])+1):(nrow(split_data[[i]])+5)),],
      outcome = "Y", 
      covariates = all_covs, 
      historical_fit = historical_fit,
      subject_id = "11638", 
      past_individual_fit = result_list[[i-1]]$individual_fit,
      individual_stack = individual_stack
    )
  } else if(i >= 31) {
    if(i == 31){
      individual_stack <- make_learner(
        Stack, unlist(list(xgb_learners, lrnr_glm, lrnr_lasso), recursive = T)
        )
    } else { 
      individual_stack <- NULL
    }
    result_list[[i]] <- make_adapt_sl(
      individual_training_data = split_data[[i]], 
      indiviual_forecast_data = ind_data[
        c((nrow(split_data[[i]])+1):(nrow(split_data[[i]])+5)),],
      outcome = "Y", 
      covariates = all_covs, 
      historical_fit = historical_fit,
      subject_id = "11638", 
      past_individual_fit = result_list[[i-1]]$individual_fit,
      individual_stack = individual_stack
    )
  } 
  i <- i + 1
}

fit2$sl_loss + fit3$sl_loss