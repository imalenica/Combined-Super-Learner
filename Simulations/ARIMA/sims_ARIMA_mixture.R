### mixture ARIMA simulations with common W

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
source(here("Simulations/ARIMA/sims_functions.R"))

data_gen_v2 <- function(n,t){
  sim_historical <- list() 
  ts_offset <- function(W){
    start <- 0.5*W[,"care_unit"] + 0.02*W[,"age"] + 0.5*W[,"sex"]
    return(start)
  }
  
  #Construct historical time-series
  for(i in 1:n){
    W <- cbind.data.frame(id       = rep(i, t),
                          sex      = rbinom(n = t, size = 1, prob = 0.5),
                          age      = round(runif(n = t, min = 19, max=90)),
                          care_unit = round(runif(n = t, min = 0, max=2)))
    
    sim_ts <- as.numeric(75 + ts_offset(W) + 
                           arima.sim(model=list(ar=c(0.6,0.4,0.1,-0.1,-0.05)),n=t))
    sim_historical[[i]] <- cbind.data.frame(series   = sim_ts,
                                            time     = seq(1:t),
                                            lag1     = lead(sim_ts, n = 1),
                                            lag2     = lead(sim_ts, n = 2),
                                            lag3     = lead(sim_ts, n = 3),
                                            lag4     = lead(sim_ts, n = 4),
                                            lag5     = lead(sim_ts, n = 5),
                                            W,
                                            Y        = lead(sim_ts, n = 15))
    
  }
  #Remove NAs: issues with sl3
  sim_historical <- do.call(rbind,sim_historical)
  sim_hist_cc <- complete.cases(sim_historical)
  sim_historical <- sim_historical[sim_hist_cc,]
  
  #Construct individual time-series
  W <- cbind.data.frame(id       = rep(n+1, 1),
                        sex      = rbinom(n = 1, size = 1, prob = 0.5),
                        age      = round(runif(n = 1, min = 19, max=90)),
                        care_unit = round(runif(n = 1, min = 0, max=2)))
  
  sim_ts <- as.numeric(75 + ts_offset(W) + 
                         c(arima.sim(model=list(ma=c(-0.5,-0.3,-0.1,0.2,0.5)),n=(2*t/3)),
                           arima.sim(model=list(ar=c(0.6,0.4,0.1,-0.1,-0.05)),n=(1*t/3))))
  sim_individual <- cbind.data.frame(series   = sim_ts,
                                     time     = seq(1:t),
                                     lag1     = lead(sim_ts, n = 1),
                                     lag2     = lead(sim_ts, n = 2),
                                     lag3     = lead(sim_ts, n = 3),
                                     lag4     = lead(sim_ts, n = 4),
                                     lag5     = lead(sim_ts, n = 5),
                                     W,
                                     Y        = lead(sim_ts, n = 15))
  
  #Construct full data set
  sim_all <- rbind.data.frame(sim_historical, sim_individual)
  sim_all_cc <- complete.cases(sim_all)
  sim_all <- sim_all[sim_all_cc,]
  
  return(list(full=sim_all, 
              historical=sim_historical, 
              individual=sim_individual))
}

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

for(m in 1:MC){
  paste0("Iteration: ", m)
  
  ### Simulate data:
  data <- data_gen_v2(n=n,t=t)
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
                             paste0("fit_mixture_stationary_v3.Rdata")), compress = TRUE)
