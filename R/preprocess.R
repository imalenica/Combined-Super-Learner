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

########################## insufficient numerics data ##########################

dat <- full_num %>%
  dplyr::group_by(subject_id) %>%
  mutate(init_time_and_date = min(time_and_date)) %>%
  mutate(min_elapsed = as.integer((time_and_date -
                                     init_time_and_date) / 60) + 1) %>%
  dplyr::filter(min_elapsed >= 90)
full_num <- full_num[(full_num$subject_id %in% dat$subject_id),]

########################## outcome measurement error ###########################

sum(is.na(new_numerics$abpmean))

detect_outliers <- function(id) {
  odd <- dplyr::filter(df_outcome_odd, subject_id == id)
  all <- dplyr::filter(full_num, subject_id == id)
  x <- all[!(all$time_and_date %in% odd$time_and_date),]
  qnt <- quantile(x$abpmean, probs=c(.25, .75), na.rm = TRUE)
  H <- 1.5 * IQR(x$abpmean, na.rm = TRUE)
  odd_out <- dplyr::mutate(odd, outlier_abpmean =
                             ifelse(abpmean < (qnt[1] - H) | 
                                      abpmean > (qnt[2] + H), 1, 0))
  odd_out$oddity_abp <- rep(1, nrow(odd_out))
  return(odd_out)
}

impute_oddities <- function(id) {
  
  oddities <- outcome_odd_df %>%
    dplyr::filter(subject_id == id) %>%
    dplyr::mutate(abpmean = ifelse(outlier_abpmean == 1, NA, abpmean))
  oddities$abpdias <- as.numeric(rep(NA, nrow(oddities)))
  oddities$abpsys <- as.numeric(rep(NA, nrow(oddities)))
  
  all <- dplyr::filter(full_num, subject_id == id)
  not_oddities <- all[!(all$time_and_date %in% oddities$time_and_date),]
  not_oddities$outlier_abpmean <- rep(0, nrow(not_oddities))
  not_oddities$oddity_abp <- rep(0, nrow(not_oddities))
  
  full <- full_join(not_oddities, oddities)
  full_sorted <- full[order(full$subject_id, full$time_and_date), ]
  full_sorted$abpmean <- as.numeric(full_sorted$abpmean)
  
  # last observation carried forward
  full_sorted$abpmean <- imputeTS::na_locf(full_sorted$abpmean)
  full_sorted$abpsys <- imputeTS::na_locf(full_sorted$abpsys)
  full_sorted$abpdias <- imputeTS::na_locf(full_sorted$abpdias)
  
  return(full_sorted)
}

df_outcome_odd <- full_num[(full_num$abpsys == full_num$abpdias),]
df_outcome_odd <- df_outcome_odd[!(df_outcome_odd$subject_id %in% 
                                     c(25373,26209)),]
outcome_odd_list <- lapply(unique(df_outcome_odd$subject_id), detect_outliers)
outcome_odd_df <- bind_rows(outcome_odd_list)
dat_clean_list <- lapply(unique(outcome_odd_df$subject_id), impute_oddities)
dat_clean2 <- bind_rows(dat_clean_list)

dat_new <- full_num[!(full_num$subject_id %in% c(25373,26209)),]
ids <- setdiff(dat_new$subject_id, dat_clean2$subject_id)
dat_new <- dat_new[(dat_new$subject_id %in% ids),]
dat_new$outlier_abpmean <- rep(0, nrow(dat_new))
dat_new$oddity_abp <- rep(0, nrow(dat_new))

dat_clean2 <- rbind(dat_clean2, dat_new)
numerics_clean <- dat_clean2[order(dat_clean2$subject_id, 
                                   dat_clean2$time_and_date), ]

# clearly indicate which columns were imputed 
colnames(numerics_clean)[9] <- "imputed_abpmean"
colnames(numerics_clean)[10] <- "imputed_abpsys_abpdias"

########################## classifying a hypotensive event #####################

numerics_clean <- new_Y_sol1(numerics_clean, cutoff = 65)
colnames(numerics_clean)[11] <- "hypo_event"

save(numerics_clean, file = here::here("Data","numerics_clean.Rdata"), compress = TRUE)

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
  print(paste0("subject_id ",sub))
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
      all <- rbind(temp1_df, temp0, fill = TRUE)
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
numerics_amine <- numerics_amine %>% distinct()
colnames(numerics_amine)[12] <- "amine"

# sedation
numerics_sedation_list <- lapply(unique(numerics_amine$subject_id), 
                              attribution_treatment, trt = "sedation",
                              dat = numerics_amine)
numerics_sedation <- rbindlist(numerics_sedation_list)
numerics_sedation <- numerics_sedation %>% distinct()
colnames(numerics_sedation)[13] <- "sedation"

# ventilation 
numerics_ventilation_list <- lapply(unique(numerics_sedation$subject_id), 
                                 attribution_treatment, trt = "ventilation",
                                 dat = numerics_sedation)
numerics_ventilation <- rbindlist(numerics_ventilation_list)
numerics_ventilation <- numerics_ventilation %>% distinct()
dim(numerics_clean)
dim(numerics_ventilation)
common <- intersect(numerics_ventilation[,-c(12:14)], numerics_clean)
dim(common)
colnames(numerics_ventilation)[14] <- "ventilation"

num_txt <- numerics_ventilation
save(num_txt, file = here::here("Data","num_txt.Rdata"), compress = TRUE)

################################################################################
# covariate data : incorporate relevant covariates
################################################################################

df <- merge(num_txt, infos_admission[,c(1:2, 4:11)], 
            by = c("subject_id", "icustay_id"), all.x = TRUE)

df <- df[!(df$age == 0),]
df$age <- ifelse(df$age == 200, NA, df$age)
df$bmi <- ifelse(df$bmi > 60 | df$bmi < 10, NA, df$bmi)

################### impute missing baseline characteristics ####################
colSums(is.na(df))

# how many subjects require imputing? 44
s1 <- df %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_sofa = sum(is.na(sofa_first)))
s2 <- df %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_sapsi = sum(is.na(sapsi_first)))
s3 <- df %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_age = sum(is.na(age)))
s4 <- df %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_age = sum(is.na(bmi)))
res <- Reduce(merge, list(s1, s2, s3, s4))
res <- res[apply(res[,-1], 1, function(x) !all(x==0)),]
s5 <- df %>% dplyr::group_by(subject_id) %>% dplyr::summarize(n = n())
res <- left_join(res, s5)
impute <- subset(res, rowSums(res[,c(2:5)]) != 0)

# add column indicating whether or not imputation was performed 
df$imputed_age <- ifelse(is.na(df$age) == TRUE, 1, 0)
df$imputed_bmi <- ifelse(is.na(df$bmi) == TRUE, 1, 0)
df$imputed_sofa <- ifelse(is.na(df$sofa_first) == TRUE, 1, 0)
df$imputed_sapsi <- ifelse(is.na(df$sapsi_first) == TRUE, 1, 0)

# impute with median
df$age <- ifelse(is.na(df$age) == TRUE, median(df$age, na.rm = TRUE), df$age)
df$bmi <- ifelse(is.na(df$bmi) == TRUE, median(df$bmi, na.rm = TRUE), df$bmi)
df$sofa_first <- ifelse(is.na(df$sofa_first) == TRUE, 
                        median(df$sofa_first, na.rm = TRUE), df$sofa_first)
df$sapsi_first <- ifelse(is.na(df$sapsi_first) == TRUE, 
                         median(df$sapsi_first, na.rm = TRUE), df$sapsi_first)

################### add column representing minutes elapsed ####################
df <- df %>%
  dplyr::group_by(subject_id) %>%
  mutate(min_elapsed = as.integer((time_and_date - min(time_and_date)) / 60)+1)

############################## final touches and save ##########################
df <- df[,c(1,3,27,4:26)]

df$care_unit <- ifelse(df$care_unit == "1", "MICU",
                       ifelse(df$care_unit == "2", "SICU",
                              ifelse(df$care_unit == "4", "CSRU",
                                     ifelse(df$care_unit == "5", "NICU",
                                            ifelse(df$care_unit == "6", "CCU", 
                                                   NA)))))

mimic <- run_class(data.frame(df),
                   cols_fac = c("gender", "care_unit", "admission_type_descr",
                                "amine", "sedation", "ventilation", "hypo_event",
                                "rank_icu", "imputed_abpmean", "subject_id",
                                "imputed_abpsys_abpdias", "imputed_age", 
                                "imputed_sofa", "imputed_sapsi", "imputed_bmi"),
                   cols_num = c("age", "sapsi_first","sofa_first", "hr", "spo2", 
                                "abpsys", "abpdias","abpmean", "bmi"))

mimic_ids <- mimic %>%
  dplyr::group_by(subject_id) %>%
  dplyr::select(abpmean) %>%
  summarise_each(funs(n_distinct))
mimic_ids_bad <- dplyr::filter(mimic_ids, abpmean == 1)
mimic_distinct <- mimic[!(mimic$subject_id %in% mimic_ids_bad$subject_id),]

mimic_all <- mimic_distinct[order(mimic_distinct$subject_id, 
                                  mimic_distinct$time_and_date), ]
save(mimic_all, file = here::here("Data","mimic_all.Rdata"), compress = TRUE)
