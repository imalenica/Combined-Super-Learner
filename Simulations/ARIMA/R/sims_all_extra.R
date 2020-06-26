### All DGP: get predictions

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

##########################################################################################
#                                 Run simulations 
##########################################################################################

set.seed(1111)

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

### set up variables
covs <- c("sex","age","care_unit","lag1","lag2","lag3","lag4","lag5")
outcome <- "Y"
splits <- seq(10,520,20)

### Simulation 0:
data <- data_gen_v0(n=n,t=t)
res_sim0 <- run_posl(data=data, covs=covs, outcome=outcome, 
                learners=learners, splits=splits, return_preds = TRUE)
preds_sim0 <- res_sim0$cv_preds

### Simulation 1:
data <- data_gen_v1(n=n,t=t)
res_sim1 <- run_posl(data=data, covs=covs, outcome=outcome, 
                     learners=learners, splits=splits, return_preds = TRUE)
preds_sim1 <- res_sim1$cv_preds

### Simulation 2:
data <- data_gen_v2(n=n,t=t)
res_sim2 <- run_posl(data=data, covs=covs, outcome=outcome, 
                     learners=learners, splits=splits, return_preds = TRUE)
preds_sim2 <- res_sim2$cv_preds

### Simulation 3:
data <- data_gen_v3(n=n,t=t)
res_sim3 <- run_posl(data=data, covs=covs, outcome=outcome, 
                     learners=learners, splits=splits, return_preds = TRUE)
preds_sim3 <- res_sim3$cv_preds

save.image(file = here::here("Simulations/ARIMA/Results", 
                             paste0("fit_preds_extra_v4.Rdata")), compress = TRUE)
