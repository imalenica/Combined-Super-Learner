args <- list(outcome="Y15_lag5_mean", id_type="AHE_high", parallel=TRUE)

if(args$id_type == "AHE_high"){
  # these ids have AHE for > 50% of their stay 
  ids <- c("791", "9275", "26451", "3094", "20934", "10520", "26558", "16621", 
           "23986_1", "29585_2", "2698", "8042_2", "12477")
} else if (args$id_type == "AHE_medium"){
  # these ids have AHE for 25-50% of their stay 
  ids <- c("33044", "32798", "2738_2", "22690", "2254", "25496", "26038", 
           "2790", "23543", "3078", "17420", "6075", "6464", "22587", "5544", 
           "29585_1", "24605_2", "12878", "2759", "9313_1", "28711", "9429", 
           "26520", "29334", "17435", "14822", "27887", "1830", "2904", 
           "2974", "2380", "28695", "8042_1")
} else if (args$id_type == "AHE_none"){
  # these ids have no AHE during their stay 
  ids <- c("21707_1", "24832", "11548", "4362", "7826", "31262", "7934_1", 
           "18414", "4341", "18517", "11071", "15121", "32346", "18494", 
           "25296", "26325", "5058", "1749", "27580", "30357", "19349", 
           "4435", "30198", "14463", "20820", "25258", "4140", "9841_2")
}

source(here::here("R", "v3", "make_adapt_sl.R"))
source(here::here("R", "v3", "run_adapt_sl.R"))
source(here::here("MIMICanalysis", "sl3_setup.R"))
source(here::here("R", "v3", "process_task.R"))
source(here::here("R", "v3", "get_weights.R"))
file_path <- "~/Downloads/"

library(data.table)
library(origami)
library(sl3)
library(here)

if(grepl("mean", args$outcome)){
  smooth_type <- "mean"
} else if(grepl("median", args$outcome)){
  smooth_type <- "median"
}

load(paste0(file_path, "individual_", smooth_type, ".Rdata"))
individual$id <- as.character(individual$id)
individual <- individual[id %in% ids,]

covs <- get_covariates(individual)

cv_stack <- make_individual_cv_stack()
weights_control <- list(window = 120, delay_decay = 10, rate_decay = 0.01)

run_slstream_allID(
  multi_individual_data = individual, ids = ids, covariates = covs, 
  outcome = args$outcome, file_path = file_path, cv_stack = cv_stack, 
  slstream_weights_control = weights_control, parallel = args$parallel
)