library(data.table)
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)
library(gridExtra)
library(ggpubr)
library(cowplot)

box <- "/Users/Rachael/Library/CloudStorage/Box-Box/"
patients_keep_final <- read.csv(paste0(box, "MIMIC-III/output/patients_keep_final.csv"))
setDT(patients_keep_final)
source("~/symphony/R/utils_plotting.R")

for(i in 1:nrow(patients_keep_final)){
  print(paste0("Starting patient ", i, " chart"))
  patient_i <- patients_keep_final[i,]
  load(paste0(
    box, "MIMIC-III/output/data/", patient_i$SUBJECT_ID, "_",
    patient_i$ICUSTAY_ID, ".Rdata"
  ))
  chart <- make_patient_chart(data_i, patient_i$SUBJECT_ID, patient_i$ICUSTAY_ID)
  ggsave(filename = paste0(patient_i$SUBJECT_ID, "_", patient_i$ICUSTAY_ID,".pdf"),
         plot = chart, device = "pdf", height = 12, width = 10,
         path = paste0(box, "MIMIC-III/output/patient_charts"))
}