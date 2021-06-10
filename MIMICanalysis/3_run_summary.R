############################## Run summary on Savio ############################
.libPaths("/global/scratch/rachelvphillips/R")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
args <- R.utils::commandArgs(T)
print(args)
ind_id <- as.character(args[1])

# data_path <- "~/Box/symphony-data/"
data_path <- "/global/scratch/rachelvphillips/symphony-data/"
load(paste0(data_path, "individual_mean.Rdata"))
individual$id <- as.character(individual$id)
d_mean <- individual
load(paste0(data_path, "individual_median.Rdata"))
individual$id <- as.character(individual$id)
d_med <- individual
rm("individual")

# source(here::here("MIMICanalysis", "utils_summary.R"))
# save_path <- "~/Google Drive/My Drive/onlineSL results MIMIC/"
source(here::here("R", "utils_summary.R"))
save_path <- "/global/scratch/rachelvphillips/symphony-results/"

horizons <- c(5, 10, 15, 20, 30)
summarize_id_continuousY(ind_id, d_mean, d_med, horizons, data_path, save_path)

for(i in 1:length(horizons)){
  print(paste0("Summarizing outcomes for horizon ", horizons[i]))
  # summarize_outcome(d_mean, ind_id, paste0("Y", horizons[i], "_AHE"),
  #                   data_path, save_path)
  # summarize_outcome(d_med, ind_id, paste0("Y", horizons[i], "_AHE"),
  #                   data_path, save_path)
  summarize_outcome(d_mean, ind_id, paste0("Y", horizons[i], "_lag5_mean"),
                    data_path, save_path)
  summarize_outcome(d_med, ind_id, paste0("Y", horizons[i], "_lag5_median"),
                    data_path, save_path)
}

