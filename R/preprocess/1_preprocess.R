library(here)
library(tidyverse)
library(parallel)
library(doParallel)
library(data.table)
source(here::here("R", "v3", "utils_mimic.R"))

################################################################################
# outcome data : fix repeated measure data and identify/impute outliers
################################################################################

load(here::here("Data", "infos_admission.RData"))
length(unique(infos_admission$subject_id)) # 32425
numerics <- readRDS(here::here("Data", "numerics.rds"))
length(unique(numerics$id)) #1353

################ account for subjects with multiple hospital stays #############
attribution_icustay_id <- function(df1, df2){
  temp = all_temp = NULL
  for(i in 1:length(unique(df1$icustay_id))){
    temp <- subset(df2, time_and_date >= df1$admission_hos[i] & 
                     time_and_date <= df1$discharge_hop[i])
    if(nrow(temp)){
      temp$icustay_id <- df1$icustay_id[i]
      all_temp <- rbind(all_temp, temp)
    }
  }
  return(all_temp)
}
no_cores <- 1
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)
v_id <- unique(numerics$id)
new_numerics <- foreach(sub = v_id,.combine = rbind) %dopar% {
  attribution_icustay_id(df1 = subset(infos_admission, subject_id == sub),
                         df2 = subset(numerics, id == sub))
}
colnames(new_numerics)[1] <- "subject_id"
length(unique(new_numerics$subject_id)) #1322

# there are 145 subjects whose icustay_id varyies in the same hospital stay
# these are patients who were moved in and out of the icu so their icustay_id 
# changes along with other covariates (typically sofa, sapsi, bmi)
# we can retain some of these 145 subjects because outcome/numerics data is not 
# available for all icustay_id's 

df <- merge(infos_admission[,c(1:2,4:11)], new_numerics,
            by = c("subject_id", "icustay_id"), all.y = TRUE)
levs_by_id <- df %>%
  dplyr::select(c("subject_id","gender", "age", "bmi","care_unit",
                  "admission_type_descr","sapsi_first", "sofa_first")) %>%
  dplyr::group_by(subject_id) %>%
  summarise_each(funs(n_distinct))
oddities <- subset(levs_by_id, rowSums(levs_by_id) - levs_by_id$subject_id != 7)
new_numerics1 <- new_numerics[!(new_numerics$subject_id %in% 
                                  oddities$subject_id),]

attribution_icustay_id2 <- function(df1, df2){
  temp = all_temp = NULL
  for( i in 1:length(unique(df1$icustay_id))){
    temp <- subset(df2, time_and_date >= df1$admission_icu[i] & 
                     time_and_date <= df1$discharge_icu[i])
    if(nrow(temp)){
      temp$icustay_id <- df1$icustay_id[i]
      all_temp <- rbind(all_temp, temp)
    }
  }
  return(all_temp)
}
v_id2 <- unique(oddities$subject_id)
new_numerics2 <- foreach(sub = v_id2,.combine = rbind) %dopar% {
  attribution_icustay_id2(df1 = subset(infos_admission, subject_id == sub),
                          df2 = subset(numerics, id == sub))
}
colnames(new_numerics2)[1] <- "subject_id"

new_numerics <- rbind(new_numerics1, new_numerics2)
new_numerics <- new_numerics[order(new_numerics$subject_id, 
                                   new_numerics$time_and_date), ]
levs_by_id <- new_numerics %>%
  dplyr::select(c("subject_id","icustay_id")) %>%
  dplyr::group_by(subject_id) %>%
  summarise_each(funs(n_distinct))
oddities2 <- subset(levs_by_id, rowSums(levs_by_id) - 
                      levs_by_id$subject_id != 1)
# oddities2 are 38 subjects whose icustay_id varyies in the same hospital stay
# and outcome/numerics data IS available for multiple icustay_id's
# we retain the data for both

num <- new_numerics %>%
  group_by(icustay_id) %>%
  mutate(init_time_and_date = min(time_and_date)) %>%
  mutate(min_elapsed = as.integer((time_and_date - init_time_and_date) / 60) + 1)
num <- num[order(num$icustay_id, num$min_elapsed), ]

new <- num %>% 
  group_by(icustay_id) %>%
  mutate(init_split = min_elapsed - lag(min_elapsed, n = 1) > 960) %>%
  select(c("icustay_id", "min_elapsed", "init_split"))

ids_mult_visits <- new %>% 
  group_by(icustay_id) %>%
  summarize(mult_split = sum(init_split, na.rm = T)) %>%
  filter(mult_split > 0)
ids_mult <- ids_mult_visits$icustay_id

list_multi_stay_data <- lapply(ids_mult, function(id){
  dat <- new %>% filter(icustay_id == id) 
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
    dat_split_i$id <- rep(paste0(id, "_", i), nrow(dat_split_i))
    dat_split_i$min <- 1:(split_i-split_init+1)
    return_dat[(split_init:split_i),] <- as.matrix(dat_split_i)
  }
  return_dat <- data.frame(return_dat)[,-3]
  colnames(return_dat) <- c("icustay_id", "min_elapsed", "id", "min")
  return_dat$min_elapsed <- as.numeric(as.character(return_dat$min_elapsed))
  return_dat$min <- as.numeric(as.character(return_dat$min))
  return(return_dat)
})
multi_stay_data <- do.call(rbind, list_multi_stay_data)

single_stay_data <- new %>% 
  filter(icustay_id %nin% ids_mult) %>%
  group_by(icustay_id) %>%
  mutate(id = icustay_id) %>%
  mutate(min = min_elapsed) %>%
  select(c("icustay_id", "min_elapsed", "id", "min"))
single_stay_data$icustay_id <- as.factor(as.character(single_stay_data$icustay_id))
single_stay_data$id <- as.factor(as.character(single_stay_data$id))

fixed_data <- rbind.data.frame(multi_stay_data, single_stay_data)
d <- data.table(merge(fixed_data, num, all.y = TRUE, by = c("icustay_id", "min_elapsed")))
d <- setorder(d, icustay_id, min_elapsed)
length(unique(d$subject_id)) #1320
length(unique(d$icustay_id)) #1361
length(unique(d$id)) #1509

d_new <- list()
for(i in 1:length(unique(d$id))){
  s <- unique(d$id)[i]
  dd <- dplyr::filter(d, id == s)
  n_occur <- data.frame(table(dd$min_elapsed))
  if(any(n_occur$Freq > 1)){
    index <- min(as.numeric(as.character(n_occur[which(n_occur$Freq > 1),]$Var1)))
    dsub <- dplyr::filter(dd, min_elapsed < index)
    d_new[[i]] <- dsub
    print(paste0(as.character(s), " subset ", min(dsub$min_elapsed), " to ", 
                 index, " min_elapsed"))
  } else {
    d_new[[i]] <- dd
  }
}

# [1] "4970_2 subset 2767 to 6618 min_elapsed"
# [1] "7934_2 subset 6303 to 14459 min_elapsed"
# [1] "8279_1 subset 1 to 2024 min_elapsed"
# [1] "10668_1 subset 1 to 6370 min_elapsed"
# [1] "11642_2 subset 22944 to 31101 min_elapsed"
# [1] "17234_1 subset 1 to 4439 min_elapsed"
# [1] "19927_1 subset 1 to 6618 min_elapsed"
# [1] "20039_2 subset 5130 to 6618 min_elapsed"
# [1] "11071 subset 1 to 828 min_elapsed"
# [1] "11705 subset 1 to 861 min_elapsed"
# [1] "1275 subset 1 to 1692 min_elapsed"
# [1] "13858 subset 1 to 611 min_elapsed"
# [1] "14083 subset 1 to 6617 min_elapsed"
# [1] "15669 subset 1 to 4583 min_elapsed"
# [1] "17611 subset 1 to 2744 min_elapsed"
# [1] "19212 subset 1 to 11534 min_elapsed"
# [1] "19564 subset 1 to 12326 min_elapsed"
# [1] "20328 subset 1 to 2831 min_elapsed"
# [1] "21704 subset 1 to 4449 min_elapsed"
# [1] "21891 subset 1 to 1701 min_elapsed"
# [1] "22364 subset 1 to 6005 min_elapsed"
# [1] "22702 subset 1 to 1734 min_elapsed"
# [1] "22925 subset 1 to 557 min_elapsed"
# [1] "23584 subset 1 to 3671 min_elapsed"
# [1] "24514 subset 1 to 2496 min_elapsed"
# [1] "2483 subset 1 to 2214 min_elapsed"
# [1] "26196 subset 1 to 4029 min_elapsed"
# [1] "28566 subset 1 to 3419 min_elapsed"
# [1] "30627 subset 1 to 272 min_elapsed"
# [1] "3163 subset 1 to 6623 min_elapsed"
# [1] "32876 subset 1 to 10873 min_elapsed"
# [1] "32908 subset 1 to 869 min_elapsed"
# [1] "4372 subset 1 to 2961 min_elapsed"
# [1] "4562 subset 1 to 108 min_elapsed"
# [1] "5402 subset 1 to 1151 min_elapsed"
# [1] "7874 subset 1 to 3702 min_elapsed"
# [1] "9510 subset 1 to 2158 min_elapsed"
# [1] "9578 subset 1 to 3520 min_elapsed"

# 38 ids with repeated min data, subset to not include repeats
d_new <- do.call(rbind, d_new)

d_new$icustay_id <- as.numeric(as.character(d_new$icustay_id))
d_new <- d_new[order(d_new$icustay_id, d_new$time_and_date), ]
d_new$icustay_id <- as.factor(as.character(d_new$icustay_id))

length(unique(d_new$subject_id)) # 1320
length(unique(d_new$icustay_id)) # 1361
length(unique(d_new$id)) # 1509

########################## outcome measurement error ###########################
# SBP > 250
# DBP > 130
# MAP > 150
# HR > 180
# SpO2 > 100

d_new$abpsys <- ifelse(d_new$abpsys <= 0 | d_new$abpsys > 250, NA, d_new$abpsys)
d_new$abpdias <- ifelse(d_new$abpdias <= 0 | d_new$abpdias > 130, NA, d_new$abpdias)
d_new$abpmean <- ifelse(d_new$abpmean <= 0 | d_new$abpmean > 150, NA, d_new$abpmean)
d_new$spo2 <- ifelse(d_new$spo2 <= 0 | d_new$spo2 > 180, NA, d_new$spo2)
d_new$hr <- ifelse(d_new$hr <= 0 | d_new$hr > 100, NA, d_new$hr)

d_new$abpdias_locf <- ifelse(is.na(d_new$abpdias), 1, 0)
d_new$abpsys_locf <- ifelse(is.na(d_new$abpsys), 1, 0)
d_new$abpmean_locf <- ifelse(is.na(d_new$abpmean), 1, 0)
d_new$spo2_locf <- ifelse(is.na(d_new$spo2), 1, 0)
d_new$hr_locf <- ifelse(is.na(d_new$hr), 1, 0)

full_num_list <- lapply(unique(d_new$id), function(i){
  ind <- dplyr::filter(d_new, id == i)
  sorted <- ind[order(ind$time_and_date), ]
  if(nrow(sorted) < 60){
    print(paste0("Insufficient data for subject_id = ", i))
    locf_dat <- NULL
  } else if(any(colSums(sorted[,c(13:17)]) == nrow(sorted))){
    name <- c(names(which(colSums(sorted[,c(13:17)]) == nrow(sorted))))
    print(paste0(name, " all NA for subject_id ", i))
    locf_dat <- NULL
  } else if(any(colSums(sorted[,c(13:17)]) >= .5*nrow(sorted))){
    name <- c(names(which(colSums(sorted[,c(13:17)]) >= .5*nrow(sorted))))
    print(paste0(name, " atleast 50% NA for subject ", i))
    locf_dat <- NULL
  } else {
    sorted$abpmean <- imputeTS::na_locf(sorted$abpmean, na_remaining = "keep")
    sorted$abpsys <- imputeTS::na_locf(sorted$abpsys, na_remaining = "keep")
    sorted$abpdias <- imputeTS::na_locf(sorted$abpdias, na_remaining = "keep")
    sorted$hr <- imputeTS::na_locf(sorted$hr, na_remaining = "keep")
    sorted$spo2 <- imputeTS::na_locf(sorted$spo2, na_remaining = "keep")
    locf_dat <- na.omit(sorted)
  }
  return(locf_dat)
})
num_list <- full_num_list[!sapply(full_num_list, is.null)] 
num_dat <- do.call(rbind, num_list)
length(unique(num_dat$subject_id)) #1058
length(unique(num_dat$icustay_id)) #1086
length(unique(num_dat$id)) #1163

num_dat <- num_dat[order(num_dat$subject_id, num_dat$icustay_id, 
                         num_dat$time_and_date), ]

################################################################################
# treatment data : incorporate to reflect minute of administration (not hour)
################################################################################

load(here::here("Data", "sedation.RData"))
sedation <- dplyr::select(sedation, -c("itemid", "duration", "startrealtime"))
sedation <- subset(sedation, sedation$icustay_id %in% num_dat$icustay_id)
sedation <- sedation[complete.cases(sedation),]
sedation <- sedation %>% distinct()
sedation <- sedation[order(c(sedation$icustay_id, sedation$starttime)), ]
sedation <- sedation[!duplicated(sedation[,c("icustay_id", "endtime")]),]
sedation <- sedation[!duplicated(sedation[,c("icustay_id", "starttime")], 
                                 fromLast = TRUE),]

load(here::here("Data", "amine.RData"))
amine <- dplyr::select(amine, -c("itemid", "duration", "startrealtime"))
amine <- subset(amine, amine$icustay_id %in% num_dat$icustay_id)
amine <- amine[complete.cases(amine),]
amine <- amine %>% distinct()
amine <- amine[order(c(amine$icustay_id, amine$starttime)), ]
amine <- amine[!duplicated(amine[,c("icustay_id", "endtime")]),]
amine <- amine[!duplicated(amine[,c("icustay_id", "starttime")], 
                           fromLast = TRUE),]

load(here::here("Data", "ventilation.RData"))
ventilation <- dplyr::select(ventilation, -c("extubated", "selfextubated"))
ventilation <- ventilation[complete.cases(ventilation),]
ventilation <- ventilation %>% distinct()
ventilation <- ventilation[order(c(ventilation$icustay_id, ventilation$realtime)), ]
ventilation <- ventilation[!duplicated(ventilation[,c("icustay_id", "realtime")]),]

############################# amine and sedation ###############################
attribution_treatment <- function(sub, trt, dat){
  df2 <- subset(dat, icustay_id == sub)
  
  if(trt %in% c("amine", "sedation")){
    if(trt == "amine"){
      df1 <- subset(amine, icustay_id == sub)
      }
    if(trt == "sedation"){
      df1 <- subset(sedation, icustay_id == sub)
      }
    if(nrow(df1) > 0){
      temp1 <- list()
      for(i in 1:nrow(df1)){
        temp1[[i]] <- subset(df2, time_and_date >= df1$starttime[i] & 
                                     time_and_date <= df1$endtime[i])
        if(nrow(temp1[[i]])){
          temp1[[i]]$txt <- 1
        }
      }
      temp1_df <- rbindlist(temp1, fill = TRUE)
      temp0 <- subset(df2, !(df2$time_and_date %in% temp1_df$time_and_date))
      if(nrow(temp0)){
        temp0$txt <- 0
      }
      all <- rbind.data.frame(temp1_df, temp0, fill = TRUE)
      }
    if(nrow(df1) == 0){
      df2$txt <- 0
      all <- df2
  }
  }
  if(trt == "ventilation"){
    
    df1 <- subset(ventilation, icustay_id == sub)
    df1 <- subset(df1, realtime >= min(df2$time_and_date) &
                    realtime <= max(df2$time_and_date))
    df1 <- df1[order(c(df1$realtime)), ]
    if(nrow(df1) > 1){
      temp10 <- list()
      temp01 <- list()
      temp11 <- list()
      temp00 <- list()
      for(i in 1:nrow(df1)){
        if(i == nrow(df1)){
          if(df1$mechvent[nrow(df1)] == 1){
            temp_end <- subset(df2, time_and_date >= df1$realtime[nrow(df1)])
            if(nrow(temp_end)){
              temp_end$txt <- 1
            }
          }
          if(df1$mechvent[nrow(df1)] == 0){
            temp_end <- subset(df2, time_and_date >= df1$realtime[nrow(df1)])
            if(nrow(temp_end)){
              temp_end$txt <- 0
            }
          }
        }
        if(i != nrow(df1)){
          if(i == 1){
            temp_start <- subset(df2, time_and_date < df1$realtime[1])
            if(nrow(temp_start)){
                temp_start$txt <- 0
              }
            }
          if(df1$mechvent[i] == 1 & df1$mechvent[i+1] == 0){
            temp10[[i]] <- subset(df2, time_and_date >= df1$realtime[i] & 
                                  time_and_date < df1$realtime[i+1])
            if(nrow(temp10[[i]])){
              temp10[[i]]$txt <- 1
            }
          }
          if(df1$mechvent[i] == 1 & df1$mechvent[i+1] == 1){
            temp11[[i]] <- subset(df2, time_and_date >= df1$realtime[i] & 
                                  time_and_date < df1$realtime[i+1])
            if(nrow(temp11[[i]])){
              temp11[[i]]$txt <- 1
            }
          }
          if(df1$mechvent[i] == 0 & df1$mechvent[i+1] == 0){
            temp00[[i]] <- subset(df2, time_and_date >= df1$realtime[i] & 
                                  time_and_date < df1$realtime[i+1])
            if(nrow(temp00[[i]])){
              temp00[[i]]$txt <- 0
            }
          }
          if(df1$mechvent[i] == 0 & df1$mechvent[i+1] == 1){
            temp01[[i]] <- subset(df2, time_and_date >= df1$realtime[i] & 
                                  time_and_date < df1$realtime[i+1])
            if(nrow(temp01[[i]])){
              temp01[[i]]$txt <- 0
            }
          }
        }
      }
      temp10_df <- rbindlist(temp10, fill = TRUE)
      temp11_df <- rbindlist(temp11, fill = TRUE)
      temp00_df <- rbindlist(temp00, fill = TRUE)
      temp01_df <- rbindlist(temp01, fill = TRUE)
      all <- rbind(temp_start, temp_end, temp10_df, temp11_df, temp00_df, 
                   temp01_df, fill = TRUE)
    }
    if(nrow(df1) == 0){
      df2$txt <- 0
      all <- df2
    }
    if(nrow(df1) == 1){
      if(df1$mechvent == 1){
        start <- subset(df2, time_and_date < df1$realtime)
        start$txt <- 0
        end <- subset(df2, time_and_date >= df1$realtime)
        end$txt <- 1
        all <- rbind(start, end)
      }
      if(df1$mechvent == 0){
        df2$txt <- 0
        all <- df2
      }
    }
  }
  all <- all[order(all$time_and_date), ]
  all <- all %>% distinct()
  if(nrow(all) > nrow(df2)){
    print(paste0("nrow(all) = ", nrow(all), " > nrow(df2) = ", nrow(df2),
                 " for icustay_id ",sub))
  }
  if(nrow(all) < nrow(df2)){
    print(paste0("nrow(all) = ", nrow(all), " < nrow(df2) = ", nrow(df2),
                 " for icustay_id ",sub))
  }
  return(all)
}

# amine
numerics_amine_list <- lapply(unique(num_dat$icustay_id), 
                              attribution_treatment, trt = "amine",
                              dat = num_dat)
numerics_amine <- rbindlist(numerics_amine_list)
colnames(numerics_amine)[18] <- "amine"

# sedation
numerics_sedation_list <- lapply(unique(num_dat$icustay_id), 
                                 attribution_treatment, trt = "sedation",
                                 dat = numerics_amine)
numerics_sedation <- rbindlist(numerics_sedation_list)
colnames(numerics_sedation)[19] <- "sedation"

# ventilation 
numerics_ventilation_list <- lapply(unique(num_dat$icustay_id), 
                                    attribution_treatment, trt = "ventilation",
                                    dat = numerics_sedation)
numerics_ventilation <- rbindlist(numerics_ventilation_list)
colnames(numerics_ventilation)[20] <- "ventilation"

dim(num_dat)
dim(numerics_ventilation)
common <- intersect(numerics_ventilation[,-c(18:20)], num_dat)
dim(common)

num_trt <- numerics_ventilation
save(num_trt, file = here::here("Data","num_trt.Rdata"), compress = TRUE)

################################################################################
# covariate data : incorporate relevant covariates
################################################################################
infos_admission$icustay_id <- as.factor(as.character(infos_admission$icustay_id))
df <- merge(num_trt, infos_admission[,c(1:2, 4:11)], 
            by = c("subject_id", "icustay_id"), all.x = TRUE)

df <- df[!(df$age < 18),]
df$age <- ifelse(df$age > 115, NA, df$age)
df$bmi <- ifelse(df$bmi > 60 | df$bmi < 10, NA, df$bmi)
colnames(df)[22] <- "sex"

colSums(is.na(df))

df$care_unit <- ifelse(df$care_unit == "1", "MICU",
                ifelse(df$care_unit == "2", "SICU",
                ifelse(df$care_unit == "4", "CSRU",
                ifelse(df$care_unit == "5", "NICU",
                ifelse(df$care_unit == "6", "CCU", NA)))))

############################## format column classes ##########################
cols_fac <- c("sex", "care_unit", "admission_type_descr", "amine", 
              "sedation", "ventilation", "rank_icu", "subject_id", "icustay_id",
              "id","abpdias_locf", "abpsys_locf", "abpmean_locf", 
              "spo2_locf", "hr_locf")
cols_num <- c("age", "sapsi_first","sofa_first", "hr", "spo2", "abpsys", 
              "abpdias","abpmean", "bmi", "min_elapsed", "min")
dat <- data.frame(df)
cols_all <- c(colnames(dat)[-c(6,12)])
dat[cols_all] <- lapply(dat[cols_all], as.character)
dat[cols_num] <- lapply(dat[cols_num], as.numeric)
dat[cols_fac] <- lapply(dat[cols_fac], as.factor)

dat$subject_id <- as.numeric(as.character(dat$subject_id))
dat$icustay_id <- as.numeric(as.character(dat$icustay_id))
mimic <- dat[order(dat$subject_id, dat$icustay_id, dat$time_and_date), ]
mimic$subject_id <- as.factor(as.character(mimic$subject_id))
mimic$icustay_id <- as.factor(as.character(mimic$icustay_id))

length(unique(mimic$subject_id)) # 1057
length(unique(mimic$icustay_id)) # 1085
length(unique(mimic$id)) # 1162
mimic <- droplevels(mimic)
save(mimic, file = here::here("Data","mimic.Rdata"), compress = TRUE)

################################################################################
# smooth vitals data with lagged 5-minute medians
################################################################################
load(file = here::here("Data","mimic.Rdata"))

# carry last obs forward for gaps in min
locf1 <- mimic %>%
  group_by(id) %>%
  select(-c("min", "min_elapsed")) %>%
  complete(time_and_date = seq(min(time_and_date), max(time_and_date), by = "min")) %>%
  mutate(row_locf = ifelse(is.na(subject_id), 1, 0)) %>%
  fill(c("subject_id", "icustay_id", "hr", "abpsys", "abpdias", "abpmean", "spo2", 
         "init_time_and_date", "abpdias_locf", "abpsys_locf", "abpmean_locf", 
         "spo2_locf", "hr_locf", "amine", "sedation", "ventilation", "rank_icu", 
         "sex", "age", "sapsi_first", "sofa_first", "bmi", "care_unit", 
         "admission_type_descr"))

locf2 <- locf1 %>%
  group_by(icustay_id) %>%
  mutate(min_elapsed = as.integer((time_and_date - min(time_and_date))/60) + 1)

locf3 <- locf2 %>%
  group_by(id) %>%
  mutate(min = as.integer((time_and_date - min(time_and_date))/60) + 1) 


d <- data.table(locf3)
d <- setorder(d, subject_id, icustay_id, id, min_elapsed)

vars <- c("abpsys", "abpdias", "abpmean", "spo2", "hr")
Ynames <- lapply(seq(), function(x) paste0(vars, "_lag_", x))
d[, Y10 := shift(abpmean, n=10, type="lead"), by=id]
d[, Y15 := shift(abpmean, n=15, type="lead"), by=id]
d[, Y20 := shift(abpmean, n=20, type="lead"), by=id]
d[, Y25 := shift(abpmean, n=25, type="lead"), by=id]
d[, Y30 := shift(abpmean, n=30, type="lead"), by=id]
d[, Y35 := shift(abpmean, n=35, type="lead"), by=id]
d[, Y40 := shift(abpmean, n=40, type="lead"), by=id]
d[, Y45 := shift(abpmean, n=45, type="lead"), by=id]

############################# initiate lags for all vitals #####################
vars <- c("abpsys", "abpdias", "abpmean", "spo2", "hr")
lag_names <- lapply(seq(10), function(x) paste0(vars, "_lag_", x))
d[, (lag_names[[1]]) := shift(.SD, n=1), by=id, .SDcols=vars]
d[, (lag_names[[2]]) := shift(.SD, n=2), by=id, .SDcols=vars]
d[, (lag_names[[3]]) := shift(.SD, n=3), by=id, .SDcols=vars]
d[, (lag_names[[4]]) := shift(.SD, n=4), by=id, .SDcols=vars]
d[, (lag_names[[5]]) := shift(.SD, n=5), by=id, .SDcols=vars]

### try 5-min median lags 
d$abpsys_lag5_median <- apply(d[,c("abpsys", "abpsys_lag_1", "abpsys_lag_2",
                                "abpsys_lag_3", "abpsys_lag_4", "abpsys_lag_5"), 
                             with=F], 1, median, na.rm = T)
d$abpdias_lag5_median <- apply(d[,c("abpdias", "abpdias_lag_1", "abpdias_lag_2", 
                                 "abpdias_lag_3", "abpdias_lag_4", 
                                 "abpdias_lag_5"), with=F], 1, median, na.rm = T)
d$abpmean_lag5_median <- apply(d[,c("abpmean", "abpmean_lag_1", "abpmean_lag_2", 
                                 "abpmean_lag_3", "abpmean_lag_4", "abpmean_lag_5"), 
                              with=F], 1, median, na.rm = T)
d$spo2_lag5_median <- apply(d[,c("spo2", "spo2_lag_1", "spo2_lag_2", "spo2_lag_3", 
                              "spo2_lag_4", "spo2_lag_5"), with=F], 1, median, 
                         na.rm = T)
d$hr_lag5_median <- apply(d[,c("hr", "hr_lag_1", "hr_lag_2", "hr_lag_3", "hr_lag_4",
                            "hr_lag_5"), with=F], 1, median, na.rm = T)

d[, Y15_lag5_median := shift(abpmean_lag5_median, n=15, type="lead"), by=id]
d[, Y20_lag5_median := shift(abpmean_lag5_median, n=20, type="lead"), by=id]
d[, Y25_lag5_median := shift(abpmean_lag5_median, n=25, type="lead"), by=id]
d[, Y30_lag5_median := shift(abpmean_lag5_median, n=30, type="lead"), by=id]
d[, Y35_lag5_median := shift(abpmean_lag5_median, n=35, type="lead"), by=id]
d[, Y40_lag5_median := shift(abpmean_lag5_median, n=40, type="lead"), by=id]
d[, Y45_lag5_median := shift(abpmean_lag5_median, n=45, type="lead"), by=id]

### try 5-min mean lags 
d$abpsys_lag5_mean <- apply(d[,c("abpsys", "abpsys_lag_1", "abpsys_lag_2", 
                                 "abpsys_lag_3", "abpsys_lag_4", "abpsys_lag_5"), 
                                  with=F], 1, mean, na.rm = T)
d$abpdias_lag5_mean <- apply(d[,c("abpdias", "abpdias_lag_1", "abpdias_lag_2", 
                                  "abpdias_lag_3", "abpdias_lag_4", "abpdias_lag_5"), 
                                   with=F], 1, mean, na.rm = T)
d$abpmean_lag5_mean <- apply(d[,c("abpmean", "abpmean_lag_1", "abpmean_lag_2", 
                                  "abpmean_lag_3", "abpmean_lag_4","abpmean_lag_5"), 
                                   with=F], 1, mean, na.rm = T)
d$spo2_lag5_mean <- apply(d[,c("spo2", "spo2_lag_1", "spo2_lag_2", "spo2_lag_3", 
                               "spo2_lag_4", "spo2_lag_5"), 
                                with=F], 1, mean, na.rm = T)
d$hr_lag5_mean <- apply(d[,c("hr", "hr_lag_1", "hr_lag_2", "hr_lag_3", "hr_lag_4",
                             "hr_lag_5"), with=F], 1, mean, na.rm = T)

d[, Y15_lag5_mean := shift(abpmean_lag5_mean, n=15, type="lead"), by=id]
d[, Y20_lag5_mean := shift(abpmean_lag5_mean, n=20, type="lead"), by=id]
d[, Y25_lag5_mean := shift(abpmean_lag5_mean, n=25, type="lead"), by=id]
d[, Y30_lag5_mean := shift(abpmean_lag5_mean, n=30, type="lead"), by=id]
d[, Y35_lag5_mean := shift(abpmean_lag5_mean, n=35, type="lead"), by=id]
d[, Y40_lag5_mean := shift(abpmean_lag5_mean, n=40, type="lead"), by=id]
d[, Y45_lag5_mean := shift(abpmean_lag5_mean, n=45, type="lead"), by=id]

################################### final save #################################

# remove lagged covariates, and set order
rm <- colnames(d)[grepl("lag_", colnames(d))]
d <- d[,-rm,with = F]
d$subject_id <- as.numeric(as.character(d$subject_id))
d$icustay_id <- as.numeric(as.character(d$icustay_id))
d$min_elapsed <- as.numeric(as.character(d$icustay_id))
d <- data.table(d)
d <- setorder(d, subject_id, icustay_id, min_elapsed)
d$subject_id <- as.factor(as.character(d$subject_id))
d$icustay_id <- as.factor(as.character(d$icustay_id))

### actual Y (no smoothing)
# only retain non-smoothed outcomes
rm <- colnames(d)[grepl("lag5_", colnames(d))]
mimic <- d[,-rm,with = F]
save(mimic, file=here("Data", "mimic.Rdata"), compress=T)

### smoothed Y 
lagged <- d[,-c("Y10","Y15","Y20","Y25","Y30","Y35","Y40","Y45"), with = F]
# only retain smoothed mean outcomes
rm <- colnames(lagged)[grepl("lag5_median", colnames(lagged))]
mimic_smooth_mean <- lagged[,-rm,with = F]
save(mimic_smooth_mean, file=here("Data", "mimic_smooth_mean.Rdata"), compress=T)
# only retain smoothed median outcomes
rm <- colnames(lagged)[grepl("lag5_mean", colnames(lagged))]
mimic_smooth_median <- lagged[,-rm,with = F]
save(mimic_smooth_median, file=here("Data", "mimic_smooth_median.Rdata"), compress=T)
