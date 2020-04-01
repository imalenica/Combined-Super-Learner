library(here)
library(tidyverse)
library(parallel)
library(doParallel)
library(data.table)
source(here::here("R", "utils_mimic.R"))

################################################################################
# outcome data : fix repeated measure data and identify/impute outliers
################################################################################

load(here::here("Data", "infos_admission.RData"))
numerics <- readRDS(here::here("Data", "numerics.rds"))

################ account for subjects with multiple hospital stays #############
attribution_icustay_id <- function(df1, df2){
  temp = all_temp = NULL
  for( i in 1:length(unique(df1$icustay_id))){
    temp <- subset(df2, time_and_date >= df1$admission_hos[i] & 
                     time_and_date <= df1$discharge_hop[i])
    if(nrow(temp)){
      temp$icustay_id <- df1$icustay_id[i]
      all_temp <- rbind(all_temp, temp)
    }
  }
  return(all_temp)
}
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)
v_id <- unique(numerics$id)
new_numerics <- foreach(sub = v_id,.combine = rbind) %dopar% {
  attribution_icustay_id(df1 = subset(infos_admission, subject_id == sub),
                         df2 = subset(numerics, id == sub))
}
colnames(new_numerics)[1] <- "subject_id"

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
v_id <- unique(oddities$subject_id)
new_numerics2 <- foreach(sub = v_id,.combine = rbind) %dopar% {
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
# we retain the data for the first icustay_id for these subjects
odd_numerics <- new_numerics[(new_numerics$subject_id %in% 
                                oddities2$subject_id),]
names <- (odd_numerics %>% group_by(subject_id) %>% 
            summarize(min(icustay_id)))[,2]
odd_numerics <- odd_numerics[(odd_numerics$icustay_id %in% names),]

new_numerics <- new_numerics[!(new_numerics$subject_id %in% 
                                 oddities2$subject_id),]
new_numerics <- rbind(new_numerics, odd_numerics)
full_num <- new_numerics[order(new_numerics$subject_id,
                                   new_numerics$time_and_date), ]
levs_by_id <- full_num %>%
  dplyr::select(c("subject_id","icustay_id")) %>%
  dplyr::group_by(subject_id) %>%
  summarise_each(funs(n_distinct))
subset(levs_by_id, rowSums(levs_by_id) - levs_by_id$subject_id != 1)

# SUMMARY - we retain all 1,282 subjects with "baseline" covariate information
# and numerics information

########################## outcome measurement error ###########################

abpsys <- ifelse(full_num$abpsys == full_num$abpdias | full_num$abpsys < 0 |
                   full_num$abpdias < 0, NA, full_num$abpsys)
abpdias <- ifelse(full_num$abpsys == full_num$abpdias | full_num$abpdias < 0 | 
                    full_num$abpsys < 0, NA, full_num$abpdias)
full_num$abpdias <- abpdias
full_num$abpsys <- abpsys

is.na.rle <- rle(is.na(full_num$abpdias))
is.na.rle$values <- is.na.rle$values & is.na.rle$lengths >= 3 
full_num$abpmean <- ifelse(inverse.rle(is.na.rle), NA, full_num$abpmean)

full_num$spo2 <- ifelse(full_num$spo2 == 0, NA, full_num$spo2)
full_num$hr <- ifelse(full_num$hr == 0, NA, full_num$hr)

full_num$abpdias_locf <- ifelse(is.na(full_num$abpdias), 1, 0)
full_num$abpsys_locf <- ifelse(is.na(full_num$abpsys), 1, 0)
full_num$abpmean_locf <- ifelse(is.na(full_num$abpmean), 1, 0)
full_num$spo2_locf <- ifelse(is.na(full_num$spo2), 1, 0)
full_num$hr_locf <- ifelse(is.na(full_num$hr), 1, 0)

full_num_list <- lapply(unique(full_num$subject_id), function(i){
  ind <- dplyr::filter(full_num, subject_id == i)
  sorted <- ind[order(ind$time_and_date), ]
  if(nrow(sorted) < 60){
    print(paste0("Insufficient data for subject_id = ", i))
    locf_dat <- NULL
  } else if(any(colSums(sorted[,c(9:13)]) == nrow(sorted))){
    name <- c(names(which(colSums(sorted[,c(9:13)]) == nrow(sorted))))
    print(paste0(name, " all NA for subject_id ", i))
    locf_dat <- NULL
  } else if(any(colSums(sorted[,c(9:13)]) >= .5*nrow(sorted))){
    name <- c(names(which(colSums(sorted[,c(9:13)]) >= .5*nrow(sorted))))
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

################### add column representing minutes elapsed ####################
num_clean <- num_dat %>%
  group_by(subject_id) %>%
  mutate(init_time_and_date = min(time_and_date)) %>%
  mutate(min_elapsed = as.integer((time_and_date - init_time_and_date) / 60) + 1)

numerics_clean <- num_clean[order(num_clean$subject_id, num_clean$time_and_date), ]
length(unique(numerics_clean$subject_id)) # 1154

################################################################################
# treatment data : incorporate to reflect minute of administration (not hour)
################################################################################

load(here::here("Data", "sedation.RData"))
sedation <- dplyr::select(sedation, -c("itemid", "duration", "startrealtime", 
                                       "icustay_id"))
sedation <- subset(sedation, sedation$subject_id %in% numerics_clean$subject_id)
sedation <- sedation[complete.cases(sedation),]
sedation <- sedation %>% distinct()
sedation <- sedation[order(c(sedation$subject_id, sedation$starttime)), ]
sedation <- sedation[!duplicated(sedation[,c("subject_id", "endtime")]),]
sedation <- sedation[!duplicated(sedation[,c("subject_id", "starttime")], 
                                 fromLast = TRUE),]

load(here::here("Data", "amine.RData"))
amine <- dplyr::select(amine, -c("itemid", "duration", "startrealtime", 
                                 "icustay_id"))
amine <- subset(amine, amine$subject_id %in% numerics_clean$subject_id)
amine <- amine[complete.cases(amine),]
amine <- amine %>% distinct()
amine <- amine[order(c(amine$subject_id, amine$starttime)), ]
amine <- amine[!duplicated(amine[,c("subject_id", "endtime")]),]
amine <- amine[!duplicated(amine[,c("subject_id", "starttime")], 
                           fromLast = TRUE),]

load(here::here("Data", "ventilation.RData"))
ventilation <- dplyr::select(ventilation, -c("extubated", "selfextubated", 
                                             "icustay_id"))
ventilation <- ventilation[complete.cases(ventilation),]
ventilation <- ventilation %>% distinct()
ventilation <- ventilation[order(c(ventilation$subject_id, 
                             ventilation$realtime)), ]
ventilation <- ventilation[!duplicated(ventilation[,c("subject_id", "realtime")]),]

############################# amine and sedation ###############################
attribution_treatment <- function(sub, trt, dat){
  df2 <- subset(dat, subject_id == sub)
  
  if(trt %in% c("amine", "sedation")){
    if(trt == "amine"){
      df1 <- subset(amine, subject_id == sub)
      }
    if(trt == "sedation"){
      df1 <- subset(sedation, subject_id == sub)
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
    
    df1 <- subset(ventilation, subject_id == sub)
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
                 " for subject_id ",sub))
  }
  if(nrow(all) < nrow(df2)){
    print(paste0("nrow(all) = ", nrow(all), " < nrow(df2) = ", nrow(df2),
                 " for subject_id ",sub))
  }
  return(all)
}

# amine
numerics_amine_list <- lapply(unique(numerics_clean$subject_id), 
                              attribution_treatment, trt = "amine",
                              dat = numerics_clean)
numerics_amine <- rbindlist(numerics_amine_list)
colnames(numerics_amine)[16] <- "amine"

# sedation
numerics_sedation_list <- lapply(unique(numerics_amine$subject_id), 
                                 attribution_treatment, trt = "sedation",
                                 dat = numerics_amine)
numerics_sedation <- rbindlist(numerics_sedation_list)
colnames(numerics_sedation)[17] <- "sedation"

# ventilation 
numerics_ventilation_list <- lapply(unique(numerics_sedation$subject_id), 
                                    attribution_treatment, trt = "ventilation",
                                    dat = numerics_sedation)
numerics_ventilation <- rbindlist(numerics_ventilation_list)
colnames(numerics_ventilation)[18] <- "ventilation"

dim(numerics_clean)
dim(numerics_ventilation)
common <- intersect(numerics_ventilation[,-c(16:18)], numerics_clean)
dim(common)

num_trt <- numerics_ventilation
save(num_trt, file = here::here("Data","num_trt.Rdata"), compress = TRUE)

################################################################################
# covariate data : incorporate relevant covariates
################################################################################

df <- merge(num_trt, infos_admission[,c(1:2, 4:11)], 
            by = c("subject_id", "icustay_id"), all.x = TRUE)

df <- df[!(df$age < 18),]
df$age <- ifelse(df$age > 115, NA, df$age)
df$bmi <- ifelse(df$bmi > 60 | df$bmi < 10, NA, df$bmi)

colSums(is.na(df))

df$care_unit <- ifelse(df$care_unit == "1", "MICU",
                       ifelse(df$care_unit == "2", "SICU",
                              ifelse(df$care_unit == "4", "CSRU",
                                     ifelse(df$care_unit == "5", "NICU",
                                            ifelse(df$care_unit == "6", "CCU", 
                                                   NA)))))
df_update <- df %>%
  group_by(subject_id) %>%
  mutate(init_time_and_date = min(time_and_date)) %>%
  mutate(min_elapsed = as.integer((time_and_date - init_time_and_date) / 60) + 1)

############################## format column classes ##########################
cols_fac <- c("gender", "care_unit", "admission_type_descr", "amine", 
              "sedation", "ventilation", "rank_icu", "subject_id", 
              "abpdias_locf", "abpsys_locf", "abpmean_locf", "spo2_locf", 
              "hr_locf")
cols_num <- c("age", "sapsi_first","sofa_first", "hr", "spo2", "abpsys", 
              "abpdias","abpmean", "bmi", "min_elapsed")
dat <- data.frame(df_update)
cols_all <- c(colnames(dat)[-3])
dat[cols_all] <- lapply(dat[cols_all], as.character)
dat[cols_num] <- lapply(dat[cols_num], as.numeric)
dat[cols_fac] <- lapply(dat[cols_fac], as.factor)

############## for now -- remove subjects with duplicated min elapsed ##########
drop_for_now <- list()
for(i in 1:length(unique(dat$subject_id))){
  s <- unique(dat$subject_id)[i]
  d <- dplyr::filter(dat, subject_id == s)
  n_occur <- data.frame(table(d$min))
  if(any(n_occur$Freq > 1)){
    # min_occur <- which(n_occur$Freq > 1)
    # if(any(diff(min_occur)) > 1){
      drop_for_now[[i]] <- as.character(s)
      print(as.character(s))
  #   } else {
  #     last <- d[(min(min_occur)-1), "sec"] # keep duplicate w same second term
  #     ahead <- d[(max(min_occur)+1), "sec"]
  #     for(j in 1:length(min_occur)){
  #
  # }
  }
}

bad_subjects_for_now <- unlist(drop_for_now)
mimic <- dat[!(dat$subject_id %in% bad_subjects_for_now),]
length(bad_subjects_for_now) #36

mimic <- mimic[order(mimic$subject_id, mimic$time_and_date), ]
mimic <- mimic[,c(1,3,15,4,13,5,10,6,9,7,11,8,12,16:26)]
colnames(mimic)[18] <- "sex"

length(unique(mimic$subject_id)) # 1117
mimic <- droplevels(mimic)
save(mimic, file = here::here("Data","mimic.Rdata"), compress = TRUE)