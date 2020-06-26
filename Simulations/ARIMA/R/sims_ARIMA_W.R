### ARIMA simulations with common W

options(warn=-1)
#remotes::install_github("tlverse/sl3@timeseries-overhaul")
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
require(smooth)

#Load Combined SL
file.sources = list.files(path=here("R/v3/"),pattern="*.R")
sapply(paste(here("R/v3"), file.sources, sep="/"),source,.GlobalEnv)
source(here("Simulations/ARIMA/R/sims_functions.R"))
source(here("Simulations/ARIMA/R/sims_data_gen.R"))

### Create a compatible arimax process:
#load(here::here("Data/fin_history60_subset.Rdata"))
#data <- fin_history60_subset

#Limit data to the first 5 hours
#dat <- data %>%
#  dplyr::group_by(subject_id) %>%
#  dplyr::mutate(time=1:n()) %>%
#  dplyr::filter(min_elapsed <= 300)

#Exclude samples with time gaps (can't model with ARIMA)
#dat_no_gap <- dat %>% 
#  dplyr::mutate(time_count = n()) %>%
#  dplyr::filter(time_count == 300)

#data <- dat_no_gap[,c("subject_id", "sex", "age", "care_unit", "abpmean")]
#data$sex <- ifelse(data$sex=="M", 1, 0)
#data$care_unit <- ifelse(data$care_unit=="CCU", 0, ifelse(data$care_unit=="CSRU", 1, 2))
#fit_arimax <- auto.arima(y = data$abpmean, xreg = as.matrix(data[,2:4]))

##########################################################################################
#                                 Run simulations 
##########################################################################################

set.seed(11)

### set up learners 
grid_params = list(max_depth = c(5,8),
                   eta = c(0.05, 0.1, 0.25), nrounds = c(50))
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
params_default <- list(nthread = getOption("sl.cores.learners", 1))
xgb_learners <- apply(grid, MARGIN = 1, function(params_tune) {
  do.call(Lrnr_xgboost$new, c(params_default, as.list(params_tune)))
})

#Tune neural net
grid_params = list(max_depth = c(5,8),
                   eta = c(0.05, 0.1, 0.25))
grid = expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
params_default = list(nthread = getOption("sl.cores.learners", 1))
xgb_learners = apply(grid, MARGIN = 1, function(params_tune) {
  do.call(Lrnr_xgboost$new, c(params_default, as.list(params_tune)))
})
lrnr_lasso <- make_learner(Lrnr_glmnet, alpha = 0.5)
lrnr_glm <- make_learner(Lrnr_glm_fast)
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_arima <- make_learner(Lrnr_arima)
lrnr_arima_strat <- Lrnr_multiple_ts$new(learner = lrnr_arima)
learners <- make_learner(Stack, unlist(list(xgb_learners, lrnr_glm),
                                       recursive = TRUE))

### set up parameters
t=60*9
n=30
MC=50

### set up variables
covs <- c("sex","age","care_unit","lag1","lag2","lag3","lag4","lag5")
outcome <- "Y"

ts_data <- list()
sl_discrete <- list()
sl_nnls_convex <- list()
sl_nnls <- list()
loss <- list()

for(m in 40:MC){
  paste0("Iteration: ", m)
  
  ### Simulate data:
  data <- data_gen_v1(n=n,t=t)
  splits <- seq(10,520,20)
  
  res <- run_posl(data=data, covs=covs, outcome=outcome, 
                  learners=learners, splits=splits)
  
  ts_data[[m]]        <- data
  sl_discrete[[m]]    <- res$res_discrete
  sl_nnls_convex[[m]] <- res$res_nnls_convex
  sl_nnls[[m]]        <- res$res_nnls 
  loss[[m]]           <- res$loss
}

save.image(file = here::here("Simulations/ARIMA/Results", 
                             paste0("fit_W_offset_stationary_v4.Rdata")), compress = TRUE)
