############################## Run summary on Savio ############################
.libPaths("/global/scratch/rachelvphillips/R")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
args <- R.utils::commandArgs(T)
print(args)
ind_id <- as.character(args[1])

#source(here::here("R", "utils_summary.R"))
source(here::here("MIMICanalysis", "utils_summary.R"))

data_path <- "/global/scratch/rachelvphillips/symphony-data/"
# save_path <- "~/Google Drive/My Drive/onlineSL results MIMIC/"
save_path <- "/global/scratch/rachelvphillips/symphony-results/"

load(paste0(data_path, "individual_mean.Rdata"))
individual$id <- as.character(individual$id)
d_mean <- individual

load(paste0(data_path, "individual_median.Rdata"))
individual$id <- as.character(individual$id)
d_median <- individual
rm("individual")

ids <- c("9275", "33044")
horizons <- c(5, 10, 15, 20, 25, 30)

for(i in 1:length(ids)){
  summarize_id(ids[i], d_mean, d_median, horizons, data_path, save_path)
}

for(i in 1:length(horizons)){
  summarize_outcome(d_mean, ids, paste0("Y", horizon[i], "_AHE"),
                    data_path, save_path)
  summarize_outcome(d_median, ids, paste0("Y", horizon[i], "_AHE"),
                    data_path, save_path)
  summarize_outcome(d_mean, ids, paste0("Y", horizon[i], "_lag5_mean"),
                    data_path, save_path)
  summarize_outcome(d_median, ids, paste0("Y", horizon[i], "_lag5_median"),
                    data_path, save_path)
}

library(here)
library(tidyverse)
library(data.table)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(xlsx)
library(animation)
library(Hmisc)
library(RColorBrewer)

pkgs <- c("ggplot2", "reshape2", "scales", "grid", "gridExtra", "ggpubr", 
          "cowplot", "xlsx", "RColorBrewer", "Hmisc", "animation", "gifski")
options(repos = structure(c(CRAN = "https://cran.rstudio.com/")))
lib <- "/global/scratch/rachelvphillips/R"
opts <- "--no-lock"
install.packages("magick", lib = lib, INSTALL_opts = opts)
module load imagemagick
imagemagick/7.0.8-29
