library(here)
library(tidyverse)
source(here::here("R", "utils_mimic.R"))

################################################################################
# Previously preprocessed data
################################################################################

load(here::here("Data", "data.Rdata"))
df <- df[order(df$subject_id, df$time_and_date), ]

df <- df[!with(df, is.na(abpmean)),]

df$care_unit <- ifelse(df$care_unit == 1, "MICU",
                  ifelse(df$care_unit == 2, "SICU",
                    ifelse(df$care_unit == 4, "CSRU",
                      ifelse(df$care_unit == 5, "NICU",
                        ifelse(df$care_unit == 6, "CCU", NA)))))

df <- run_class(df,
                   cols_fac = c("gender", "care_unit", "admission_type_descr",
                                "amine", "sedation", "ventilation", "event"),
                   cols_num = c("periode", "time", "age", "sapsi_first",
                                "sofa_first", "bmi", "los_icu", "los_hospital",
                                "hr", "spo2", "abpsys", "abpdias", "abpmean"))

#################### Subject IDs representing multiple patients ################

levs_by_id <- df %>%
  dplyr::select(c("subject_id", "gender", "age", "bmi", "care_unit",
                  "admission_type_descr", "los_icu", "los_hospital",
                  "sapsi_first", "sofa_first")) %>%
  dplyr::group_by(subject_id) %>%
  summarise_each(funs(n_distinct))

oddities <- subset(levs_by_id, rowSums(levs_by_id) - levs_by_id$subject_id != 9)
oddities
dat_new <- df[!(df$subject_id %in% oddities$subject_id),]

########################## Outcome Measurement Error ###########################

df_outcome_odd <- dat_new[(dat_new$abpsys == dat_new$abpdias),]

detect_outliers <- function(id) {
  odd <- dplyr::filter(df_outcome_odd, subject_id == id)
  all <- dplyr::filter(dat_new, subject_id == id)
  x <- all[!(all$time_and_date %in% odd$time_and_date),]
  qnt <- quantile(x$abpmean, probs=c(.25, .75), na.rm = TRUE)
  H <- 1.5 * IQR(x$abpmean, na.rm = TRUE)
  odd_out <- dplyr::mutate(odd, outlier =
  ifelse(abpmean < (qnt[1] - H) | abpmean > (qnt[2] + H), 1, 0))
  return(odd_out)
}
df_outcome_odd <- df_outcome_odd[!(df_outcome_odd$subject_id %in% c(25373,26209)),]
outcome_odd_list <- lapply(unique(df_outcome_odd$subject_id), detect_outliers)
outcome_odd_df <- bind_rows(outcome_odd_list)

remove_outliers <- function(id) {
  outlier <- outcome_odd_df %>%
    dplyr::filter(subject_id == id) %>%
    dplyr::filter(outlier == 1)
  all <- dplyr::filter(dat_new, subject_id == id)
  x <- all[!(all$time_and_date %in% outlier$time_and_date),]
  return(x)
}

dat_clean_list <- lapply(unique(outcome_odd_df$subject_id), remove_outliers)
dat_clean <- bind_rows(dat_clean_list)
dat_new <- dat_new[!(dat_new$subject_id %in% c(25373,26209)),]
ids <- setdiff(dat_new$subject_id, dat_clean$subject_id)
dat_clean <- rbind(dat_clean, dat_new[(dat_new$subject_id %in% ids),])

###################### Missing baseline characteristics ########################

colSums(is.na(dat_clean))

s1 <- dat_clean %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_bmi = sum(is.na(bmi)))
s2 <- dat_clean %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_sofa = sum(is.na(sofa_first)))
s3 <- dat_clean %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_sapsi = sum(is.na(sapsi_first)))
s4 <- dat_clean %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_ad = sum(is.na(admission_type_descr)))
s5 <- dat_clean %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_hospital = sum(is.na(los_hospital)))
res <- Reduce(merge, list(s1, s2, s3, s4, s5))
res <- res[apply(res[,-1], 1, function(x) !all(x==0)),]
s6 <- dat_clean %>% dplyr::group_by(subject_id) %>% dplyr::summarize(n = n())
res <- left_join(res, s6)
remove <- subset(res, rowSums(res[,c(3:6)]) != 0)

new_dat <- dat_clean[!(dat_clean$subject_id %in% remove$subject_id),]
mimic <- new_dat

colSums(is.na(mimic))

########################## Insufficient patient data ###########################

dat <- mimic %>%
          dplyr::group_by(subject_id) %>%
            mutate(init_time_and_date = min(time_and_date)) %>%
              mutate(min_elapsed = as.integer((time_and_date -
                init_time_and_date) / 60) + 1) %>%
                  dplyr::filter(min_elapsed >= 90)
mimic <- mimic[(mimic$subject_id %in% dat$subject_id),]

########################## Classifying a hypotensive event #####################

mimic <- new_Y_sol1(mimic, cutoff = 65)


mimic <- mimic[order(mimic$subject_id, mimic$time_and_date), ]
save(mimic, file = here::here("Data","mimic.Rdata"), compress = TRUE)

################################################################################
# Data numerics
################################################################################

########################### Omitting 30 minute gap periods #####################

load(here::here("Data", "data_numerics.Rdata"))
num_dat <- numerics[(numerics$subject_id %in% mimic$subject_id),]
num_outcome_odd <- num_dat[(num_dat$abpsys == num_dat$abpdias),]

detect_outliers <- function(id) {
  odd <- dplyr::filter(num_outcome_odd, subject_id == id)
  all <- dplyr::filter(num_dat, subject_id == id)
  x <- all[!(all$time_and_date %in% odd$time_and_date),]
  qnt <- quantile(x$abpmean, probs=c(.25, .75), na.rm = TRUE)
  H <- 1.5 * IQR(x$abpmean, na.rm = TRUE)
  odd_out <- dplyr::mutate(odd, outlier =
  ifelse(abpmean < (qnt[1] - H) | abpmean > (qnt[2] + H), 1, 0))
  return(odd_out)
}
num_odd_list <- lapply(unique(num_outcome_odd$subject_id), detect_outliers)
num_odd_df <- bind_rows(num_odd_list)


remove_outliers <- function(id) {
  outlier <- num_odd_df %>%
    dplyr::filter(subject_id == id) %>%
    dplyr::filter(outlier == 1)
  all <- dplyr::filter(num_dat, subject_id == id)
  x <- all[!(all$time_and_date %in% outlier$time_and_date),]
  return(x)
}

num_clean_list <- lapply(unique(num_outcome_odd$subject_id), remove_outliers)
num_clean <- bind_rows(num_clean_list)
ids <- setdiff(num_dat$subject_id, num_clean$subject_id)
num_clean <- rbind(num_clean, num_dat[(num_dat$subject_id %in% ids),])

merge_dat <- left_join(num_clean, mimic)

############################### Fill in missing values #########################

merged_dat <- merge_dat %>% dplyr::group_by(subject_id) %>%
             tidyr::fill(gender, age, sapsi_first, sofa_first, bmi,
                  care_unit, admission_type_descr, los_icu, los_hospital)

s1 <- merged_dat %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_age = sum(is.na(age)))
s2 <- merged_dat %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_gender = sum(is.na(gender)))
s3 <- merged_dat %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_sofa = sum(is.na(sofa_first)))
s4 <- merged_dat %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_sapsi = sum(is.na(sapsi_first)))
s5 <- merged_dat %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_ad = sum(is.na(admission_type_descr)))
s6 <- merged_dat %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_hospital = sum(is.na(los_hospital)))
s7 <- merged_dat %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_careunit = sum(is.na(care_unit)))
s8 <- merged_dat %>% dplyr::group_by(subject_id) %>%
  dplyr::summarize(na_icu = sum(is.na(los_icu)))
res <- Reduce(merge, list(s1, s2, s3, s4, s5, s6, s7, s8))
res <- res[apply(res[,-1], 1, function(x) !all(x==0)),]
s9 <- merged_dat %>% dplyr::group_by(subject_id) %>% dplyr::summarize(n = n())
res <- left_join(res, s9)

merged_dat <- merged_dat %>% dplyr::group_by(subject_id) %>%
             tidyr::fill(gender, age, sapsi_first, sofa_first, bmi,
                  care_unit, admission_type_descr, los_icu, los_hospital,
                  .direction = "up")

###################### Classifying a hypotensive event #########################

mimic_nogap <- new_Y_sol1(merged_dat, cutoff = 65)


mimic_nogap <- mimic_nogap[order(mimic_nogap$subject_id, mimic_nogap$time_and_date), ]
save(mimic_nogap, file = here::here("Data","mimic_nogap.Rdata"),compress = TRUE)


mimic_gap <- mimic_nogap[!is.na(mimic_nogap$event),]
save(mimic_gap, file = here::here("Data", "mimic_gap.Rdata"), compress = TRUE)
