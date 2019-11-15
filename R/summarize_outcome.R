summary_function <- function(vector, function_name) { 
  FUN <- match.fun(function_name) 
  FUN(vector) 
}

add_summary_outcomes <- function(data, gap = 30, summary_window = 15, 
                                summary_measure = "median"){
  dat <- data %>%
  group_by(subject_id) %>%
  mutate(Y = dplyr::lead(abpmean, n = gap, default = NA))
  
  dat$Y_summary <- rep(NA, nrow(dat))
  
  for(i in unique(dat$subject_id)){
    sub <- dat[dat$subject_id == i,]
    init <- na.omit(sub[,"Y"])
    Y_summary <- rep(0, nrow(init)-summary_window)
    for(j in 1:(nrow(init)-summary_window)){
      vec <- data.frame(init[(j):(j+summary_window-1),])[,1]
      Y_summary[j] <- summary_function(vec, summary_measure)
    }
    out <- c(Y_summary, rep(NA, (nrow(sub)-length(Y_summary))))
    dat[dat$subject_id == i,"Y_summary"] <- out
  }
  return(dat)
}

library(dplyr)
library(here)
library(data.table)

load(here("Data", "comSLargs.Rdata"))
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
