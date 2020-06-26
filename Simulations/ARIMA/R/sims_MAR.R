### Mixture autoregressive models (MAR)
# consist of multiple stationary or non-stationary autoregressive components
# MAR(K; p1, p2, .... pk): finite mixture of K Gaussian AR models

#https://cran.r-project.org/web/packages/tsfeatures/vignettes/tsfeatures.html
#Entropy: measures the forecastability of the series
#Stability: based on non-overlapping windows; variance of the means
#max_level_shift: largest mean shift between two consecutive windows (along with time)
#max_var_shift: largest variance shift between two consecutive windows (along with time)
#max_kl_shift: largest shift in Kulback-Leibler divergence between two consecutive windows (along with time)
#crossing_points: number of times a ts crosses the mean line
#nperiods: number of seasonal periods (=1 for non-seasonal data)

#### STL decomposition:
# x_t = f_t + s_{1,t} + s_{M,t + e_t}
#trend: smoothed trend component (f_t)
#spike: spikiness of a time-series (variance of the leave-one-out variances of 
#                                   the remainder component (e_t))

#### heterogeneity:
#arch_acf: sum of squares of the first 12 autocorrelations of x_t^2
#garch_acf: sum of squares of the first 12 autocorrelations of z_t^2
#arch_r2: R^2 value of an AR model applied to x_t^2
#garch_r2: R^2 value of an AR model applied to z_t^2

#nonlinearity: coefficient computed using a modification of the statistic used in 
#              Terasvirta's nonlinearity test; large values = nonlinear

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
suppressMessages(library(gratis))
suppressMessages(library(tsfeatures))

#Load Combined SL
file.sources = list.files(path=here("R/v3/"),pattern="*.R")
sapply(paste(here("R/v3"), file.sources, sep="/"),source,.GlobalEnv)
source(here("Simulations/ARIMA/R/sims_functions.R"))
source(here("Simulations/ARIMA/R/sims_data_gen.R"))

###############################################################################################
#load(here::here("Data/fin_history60_subset.Rdata"))
#data <- fin_history60_subset

#Limit data to the first 7 hours
#dat <- data %>%
#  dplyr::group_by(subject_id) %>%
#  dplyr::mutate(time=1:n()) %>%
#  dplyr::filter(min_elapsed <= 420)

#Exclude samples with time gaps 
#dat_no_gap <- dat %>% 
#  dplyr::mutate(time_count = n()) %>%
#  dplyr::filter(time_count == 420)

#ts_abp <- dat_no_gap$abpmean
#ts_abp_features <- tsfeatures(ts_abp, features=c("entropy", "stability", "max_level_shift", 
#                                                 "max_var_shift","max_kl_shift", 
#                                                 "crossing_points", "stl_features",
#                                                 "heterogeneity","nonlinearity"))
###############################################################################################

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
MC=5

### set up variables
covs <- c("sex","age","care_unit","lag1","lag2","lag3","lag4","lag5")
outcome <- "Y"

ts_data <- list()
sl_discrete <- list()
sl_nnls_convex <- list()
sl_nnls <- list()
loss <- list()

for(m in 6:10){
  paste0("Iteration: ", m)
  
  ### Simulate data:
  data <- data_gen_v3(n=n,t=t)
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
                             paste0("fit_MAR_v4.Rdata")), compress = TRUE)


