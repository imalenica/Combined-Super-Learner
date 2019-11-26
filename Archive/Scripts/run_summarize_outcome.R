library(dplyr)
library(here)
library(data.table)

load(here("Data", "comSLargs.Rdata"))
source(here("R", "summarize_outcome.R"))
data <- comSLargs[[1]]

data1 <- add_summary_outcome(data, summary_measure = "median")
colnames(data1)[31] <- "Y_median"
data2 <- add_summary_outcome(data1, summary_measure = "mean")
colnames(data2)[32] <- "Y_mean"
data3 <- add_summary_outcome(data2, summary_measure = "min")
colnames(data3)[33] <- "Y_min"
data4 <- add_summary_outcome(data3, summary_measure = "max")
colnames(data4)[34] <- "Y_max"

data <- data.table(data4)
save(data, file = here("Data", "data_summaryY.Rdata"), compress = TRUE)