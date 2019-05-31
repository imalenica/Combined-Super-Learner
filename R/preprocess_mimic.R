library(here)
library(tidyverse)
source(here::here("R", "utils_mimic.R"))

load(here::here("Data", "data.Rdata"))

df$care_unit <- ifelse(df$care_unit == 1, "MICU",
                  ifelse(df$care_unit == 2, "SICU",
                    ifelse(df$care_unit == 4, "CSRU",
                      ifelse(df$care_unit == 5, "NICU",
                        ifelse(df$care_unit == 6, "CCU", NA)))))

# PROBLEM: some samples are repeated?
# first remove the rows without our outcome of interest
new_df <- df[!with(df, is.na(abpdias) & is.na(abpsys)),]
rs <- rowSums(is.na(new_df))
df_ordered <- new_df[order(new_df$subject_id, new_df$time_and_date, rs), ]
dat <- df_ordered[!duplicated(df_ordered[c('subject_id', 'time_and_date')]), ]

# set time to be continuous
dat <- dat %>%
          group_by(subject_id) %>%
            mutate(init_time_and_date = min(time_and_date)) %>%
              mutate(min_elapsed = as.integer((time_and_date -
                init_time_and_date) / 60) + 1)

# PROBLEM: How to deal with missing values?
colSums(is.na(dat))

s1 <- dat %>% dplyr::group_by(subject_id) %>% dplyr::summarize(na_bmi = sum(is.na(bmi)))
s2 <- dat %>% dplyr::group_by(subject_id) %>% dplyr::summarize(na_sofa = sum(is.na(sofa_first)))
s3 <- dat %>% dplyr::group_by(subject_id) %>% dplyr::summarize(na_sapsi = sum(is.na(sapsi_first)))
s4 <- dat %>% dplyr::group_by(subject_id) %>% dplyr::summarize(na_ad= sum(is.na( admission_type_descr)))
res <- Reduce(merge, list(s1, s2, s3, s4))
res <- res[apply(res[,-1], 1, function(x) !all(x==0)),]
s5 <- dat %>% dplyr::group_by(subject_id) %>% dplyr::summarize(n = n())
res <- left_join(res, s5)
# res[!(res$na_bmi == res$n),]
# oddities <- c("3372","4685","9882","13715","16122","19620","21484","22130")
# odd_dat <- filter(dat, subject_id %in% oddities)

new_dat <- dat[!(dat$subject_id %in% res$subject_id),]
colSums(is.na(new_dat))
length(unique(new_dat$subject_id))

mimic <- run_class(new_dat,
                   cols_fac = c("subject_id", "gender", "care_unit",
                                "admission_type_descr", "amine", "sedation",
                                "ventilation", "event"),
                   cols_num = c("periode", "time", "age", "sapsi_first",
                                "sofa_first", "bmi", "los_icu", "los_hospital",
                                "hr", "spo2", "abpsys", "abpdias", "abpmean"))

save(mimic, file = here::here("Data", "mimic.Rdata"), compress = TRUE)
