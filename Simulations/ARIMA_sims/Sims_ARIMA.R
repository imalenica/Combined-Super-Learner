### ARIMA simulations

options(warn=-1)
suppressMessages(library(xtable))
suppressMessages(library(here))
suppressMessages(library(pROC))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(forecast))
suppressMessages(library(kableExtra))
suppressMessages(library(sl3))
suppressMessages(library(origami))
suppressMessages(library(data.table))
options(xtable.comment = FALSE)

### Load all the files 

#Load Combined SL
file.sources = list.files(path=here("R/v3/"),pattern="*.R")
sapply(paste(here("R/v3"), file.sources, sep="/"),source,.GlobalEnv)
#Load simulation scripts
file.sources = list.files(path=here("Simulations/ARIMA_sims/R/"),pattern="*.R")
sapply(paste(here("Simulations/ARIMA_sims/R"), file.sources, sep="/"),source,.GlobalEnv)
#Load the subset data
load(here::here("Data/fin_history60_subset.Rdata"))

### set up learners 
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

learners_ind <- make_learner(Stack, lrnr_mean, lrnr_glm, lrnr_lasso)

#### Covariates to be included
W <-         c("min_elapsed","hr","hr_locf","rank_icu","sex","age","sapsi_first","sofa_first",
               "bmi","care_unit","admission_type_descr")

W_with_tv <- c("min_elapsed","hr","hr_locf","spo2","spo2_locf","amine","sedation","ventilation",
               "rank_icu","sex","age","sapsi_first","sofa_first","bmi","care_unit","admission_type_descr",
               "amine_switch1","amine_switch0","sedation_switch1","sedation_switch0",
               "ventilation_switch1","ventilation_switch0","amine_history_mean",
               "sedation_history_mean","ventilation_history_mean","spo2_Mean")

######################### Process data #########################
data <- fin_history60_subset

#Limit data to the first 5 hours
dat <- data %>%
  dplyr::group_by(subject_id) %>%
  dplyr::mutate(time=1:n()) %>%
  dplyr::filter(min_elapsed <= 300)

#Exclude samples with time gaps (can't model with ARIMA)
dat_no_gap <- dat %>% 
  dplyr::mutate(time_count = n()) %>%
  dplyr::filter(time_count == 300)

#Number of samples in the subset with at least 5 hours of data: (263)
no_samples <- length(unique(dat_no_gap$subject_id)) 

######################### Sample and fit ARIMA models #########################

fit_v1 <- sample_and_fit(dat=dat_no_gap, pool_size=no_samples, ind_size=no_samples)

#Save the fits
save.image(file=here("Simulations/ARIMA_sims/Results/fit_ARIMA_v1.Rdata"))

#Examine accuracy of one-step ahead forecasts
acc_all <- rbind.data.frame(fit_v1$fit_pooled$accuracy_pooled, fit_v1$fit_pooled$accuracy)
acc_all$subject_id<-as.numeric(levels(acc_all$subject_id))[acc_all$subject_id]
acc_all[1,1]<-"pooled"
acc_ind <- t(rbind_list(lapply(fit_v1$fit_individuals, function(x){x$accuracy})))
acc_ind <- cbind.data.frame(row.names(acc_ind), acc_ind)
colnames(acc_ind) <- names(acc_all)
row.names(acc_ind) <- NULL
acc_comb <- merge(x=acc_all[,c("subject_id", "MAE")], 
                  acc_ind[,c("subject_id", "MAE")],  by="subject_id")
names(acc_comb)[2:3] <- c("MAE pooled", "MAE individual")

######################### Simulate from ARIMA models #########################

#Example patient (very different fits then pooled)
id <- 13569
arima_simulate <- run_simulation(id=id, fit=fit_v1, df=dat_no_gap, W=W)

#Save the simulation
save(arima_simulate, file = here::here("Simulations/ARIMA_sims/Results/t=0", 
                                   paste0("arima_simulate", id, ".Rdata")), compress = TRUE)

############### Combined SL for ARIMA simulated time-series with t=0 ###############

### set up W,Y,datasets
data_ind <- arima_simulate$data_ind
data_hist <- arima_simulate$data_hist
covs <- W
outcome <- "Y_20"

### Learn the historical fit
historical_fit <- make_historical_fit(
  historical_data = data_hist,
  outcome = outcome, 
  covariates = covs, 
  id = "subject_id", 
  historical_stack = learners
)

### Learn the individual fit 
#Split data into time chuncks
splits <- seq(10,300,20)
split_data <- lapply(splits, function(x) data.table(data_ind[1:x,]))

#Save all results
result_list <- list()
result_weights <- list()
result_forecasts <- list()

for(i in 1:length(splits)){
  result_list[[i]] <- make_adapt_sl(
    individual_training_data = split_data[[i]], 
    indiviual_forecast_data = data_ind[c((nrow(split_data[[i]])+1):(nrow(split_data[[i]])+5)),],
    outcome = outcome, 
    covariates = covs, 
    subject_id = id,
    historical_fit = historical_fit,
    individual_stack = learners_ind
  )
  result_weights[[i]] <- t(result_list[[i]]$sl_weights)
  result_forecasts[[i]] <- t(result_list[[i]]$sl_forecasts)
}

save.image(file = here::here("Simulations/ARIMA_sims/Results", 
                             paste0("fit_combSL_", id, ".Rdata")), compress = TRUE)

############## Examine closer match between pool and individual ############## 

#Example patient (very different fits then pooled)
id <- 13569
arima_simulate <- run_simulation(id=id, fit=fit_v1, df=dat_no_gap, W=W)

#Save the simulation
save(arima_simulate, file = here::here("Simulations/ARIMA_sims/Results/t=0", 
                                       paste0("arima_simulate", id, ".Rdata")), compress = TRUE)

### set up W,Y,datasets
data_ind <- arima_simulate$data_ind
data_hist <- arima_simulate$data_hist
covs <- W
outcome <- "Y_20"

### Learn the historical fit
historical_fit <- make_historical_fit(
  historical_data = data_hist,
  outcome = outcome, 
  covariates = covs, 
  id = "subject_id", 
  historical_stack = learners
)

### Learn the individual fit 
#Split data into time chuncks
splits <- seq(10,300,20)
split_data <- lapply(splits, function(x) data.table(data_ind[1:x,]))

#Save all results
result_list <- list()
result_weights <- list()
result_forecasts <- list()

for(i in 1:length(splits)){
  result_list[[i]] <- make_adapt_sl(
    individual_training_data = split_data[[i]], 
    indiviual_forecast_data = data_ind[c((nrow(split_data[[i]])+1):(nrow(split_data[[i]])+5)),],
    outcome = outcome, 
    covariates = covs, 
    subject_id = id,
    historical_fit = historical_fit,
    individual_stack = learners_ind
  )
  result_weights[[i]] <- t(result_list[[i]]$sl_weights)
  result_forecasts[[i]] <- t(result_list[[i]]$sl_forecasts)
}

save.image(file = here::here("Simulations/ARIMA_sims/Results", 
                             paste0("fit_combSL_", id, ".Rdata")), compress = TRUE)
















