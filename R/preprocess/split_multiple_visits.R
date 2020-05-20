library(tidyverse)
library(data.table)
library(here)
library(lubridate)

source(here::here("R", "v3", "utils_mimic.R"))

load(here::here("Data", "all_history30.Rdata"))
load(here::here("Data", "all_history60.Rdata"))

new30 <- all_history30 %>% 
  group_by(subject_id) %>%
  mutate(init_split = min_elapsed - lag(min_elapsed, n = 1) > 600) %>%
  select(c("subject_id", "min_elapsed", "init_split"))

ids_mult_visits <- new30 %>% 
  group_by(subject_id) %>%
  summarize(mult_split = sum(init_split, na.rm = T)) %>%
  filter(mult_split > 0)

ids_mult <- ids_mult_visits$subject_id

list_multi_stay_data <- lapply(ids_mult, function(id){
  dat <- new30 %>% filter(subject_id == id) 
  split_points <- append((which(dat$init_split == TRUE)-1), nrow(dat))
  return_dat <- matrix(NA, ncol = 5, nrow = nrow(dat))
  for(i in 1:length(split_points)){
    split_i <- split_points[i]
    if(i == 1){
      split_init <- 1
    } else {
      split_init <- split_points[i-1]+1
    }
    dat_split_i <- dat[split_init:split_i,]
    dat_split_i$subject_stay_id <- rep(paste0(id, "_", i), nrow(dat_split_i))
    dat_split_i$min <- 1:(split_i-split_init+1)
    return_dat[(split_init:split_i),] <- as.matrix(dat_split_i)
  }
  return_dat <- data.frame(return_dat)[,-3]
  colnames(return_dat) <- c("subject_id", "min_elapsed", "subject_stay_id", "min")
  return_dat$min_elapsed <- as.numeric(as.character(return_dat$min_elapsed))
  return_dat$min <- as.numeric(as.character(return_dat$min))
  return(return_dat)
})
multi_stay_data <- do.call(rbind, list_multi_stay_data)

single_stay_data <- new30 %>% 
  filter(subject_id %nin% ids_mult) %>%
  group_by(subject_id) %>%
  mutate(subject_stay_id = subject_id) %>%
  mutate(min = min_elapsed) %>%
  select(c("subject_id", "min_elapsed", "subject_stay_id", "min"))

fixed_data <- rbind.data.frame(multi_stay_data, single_stay_data)

d30 <- data.table(merge(fixed_data, all_history30, all.y = TRUE,
                          by = c("subject_id", "min_elapsed")))
d60 <- data.table(merge(fixed_data, all_history60, all.y = TRUE,
                        by = c("subject_id", "min_elapsed")))

mimic_history30 <- setorder(d30, subject_id, min_elapsed)
dim(mimic_history30) 
# [1] 5035962     163
mimic_history60 <- setorder(d60, subject_id, min_elapsed)

# fuck -- just realized the history is all fucked for the multiple stays.. 
# we would have to redo the history after dealing with the multi stay.. 
# and that takes forever so I am just going to keep the first stay, woohoo..
mimic_history30 <- mimic_history30 %>% filter(!grepl("_2|_3|_4", subject_stay_id)) 
mimic_history60 <- mimic_history60 %>% filter(!grepl("_2|_3|_4", subject_stay_id))
dim(mimic_history30) 
# [1] 4434999     163
# data loss ~ 12%

save(mimic_history30, file = here::here("Data","mimic_history30.Rdata"), 
     compress = TRUE)
save(mimic_history60, file = here::here("Data","mimic_history60.Rdata"), 
     compress = TRUE)