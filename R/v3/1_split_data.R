library(here)
library(data.table)
library(tidyverse)
library(origami)
source(here("R", "v3", "process_data.R"))
#################### make sure data variations have same ids  ##################
load(here::here("data","mimic30_smooth_median.Rdata"))
load(here::here("data","mimic60_smooth_median.Rdata"))
mimic30_smooth_median_id <- unique(mimic30_smooth_median$id)
mimic60_smooth_median_id <- unique(mimic60_smooth_median$id)

load(here::here("data","mimic30_smooth_mean.Rdata"))
load(here::here("data","mimic60_smooth_mean.Rdata"))
mimic30_smooth_mean_id <- unique(mimic30_smooth_mean$id)
mimic60_smooth_mean_id <- unique(mimic60_smooth_mean$id)

load(here::here("data","mimic30.Rdata"))
load(here::here("data","mimic60.Rdata"))
mimic30_id <- unique(mimic30$id)
mimic60_id <- unique(mimic60$id)

all.equal(mimic30_id, mimic60_id)
all.equal(mimic30_id, mimic30_smooth_mean_id)
all.equal(mimic30_id, mimic60_smooth_mean_id)
all.equal(mimic30_id, mimic30_smooth_median_id)
all.equal(mimic30_id, mimic60_smooth_median_id)

#################### how long are people in the hospital?  #####################
max_hr_by_id <- mimic30 %>% group_by(id) %>% summarize(max_hr = max(min)/60)
hist(max_hr_by_id$max_hr)
summary(max_hr_by_id$max_hr)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.8833  20.9083  44.3500  68.4011  92.1958 441.8667 
subject_data <- distinct(dplyr::select(mimic30, c("subject_id", "id")))
maxtime_by_subject <- left_join(max_hr_by_id, subject_data)

potential_historical <- maxtime_by_subject %>% 
  group_by(subject_id) %>%
  mutate(subject_max = max(max_hr)) %>%
  mutate(candidate = ifelse(subject_max >= 3 & subject_max <= 96, 1, 0)) %>%
  filter(candidate == 1)
potential_historical_subject_ids <- potential_historical$subject_id

mimic30_potential <- mimic30 %>%
  filter(subject_id %in% potential_historical_subject_ids)
########## subset data based on hypotensive and vasopressor strata #############

# create hypotensive strata
data_hypo <- mimic30_potential %>%
  mutate(hypo = as.factor(ifelse(abpmean < 65, "1", "0"))) %>%
  group_by(subject_id, hypo) %>%
  summarise(n = n()) %>%
  mutate(freq = round(n / sum(n), 2)) %>%
  filter(hypo == 1) %>%
  mutate(hypo_strata = ifelse(freq < .5, 1, 2)) %>%
  select(c(subject_id, hypo_strata))

never_hypo <- potential_historical_subject_ids[
  -which(potential_historical_subject_ids %in% data_hypo$subject_id)]
data_never_hypo <- data.frame(subject_id = never_hypo, 
                              hypo_strata = rep(0, length(never_hypo)))
data_never_hypo$subject_id <- as.character(data_never_hypo$subject_id)
data_hypo <- data.frame(data_hypo)
data_hypo$subject_id <- as.character(data_hypo$subject_id)

AHE_summary <- rbind(data_hypo, data_never_hypo)
AHE_summary$hypo_strata <- as.factor(as.character(AHE_summary$hypo_strata))

data <- left_join(mimic30_potential, AHE_summary)

# create vasopressor strata
data_amine <- mimic30_potential %>%
  group_by(subject_id, amine) %>%
  summarise(n = n()) %>%
  mutate(freq = round(n / sum(n), 2)) %>%
  filter(amine == 1) %>%
  mutate(amine_strata = ifelse(freq < .5, 1,
                        ifelse(freq >= .5 & freq < 1, 2, 3))) %>%
  select(c(subject_id, amine_strata))

never_amine <- potential_historical_subject_ids[
  -which(potential_historical_subject_ids %in% data_amine$subject_id)]
data_never_amine <- data.frame(subject_id = never_amine, 
                               amine_strata = rep(0, length(never_amine)))
data_never_amine$subject_id <- as.character(data_never_amine$subject_id)
data_amine <- data.frame(data_amine)
data_amine$subject_id <- as.character(data_amine$subject_id)
amine_summary <- rbind(data_amine, data_never_amine)
amine_summary$amine_strata <- as.factor(as.character(amine_summary$amine_strata))

sub <- merge(amine_summary, AHE_summary, by = "subject_id", all = TRUE)

# take 50% of stratified subsets 
historical_subjects <- list()
set.seed(491)
for(i in 0:2){
  idx <- sub[which(sub$hypo_strata == i & sub$amine_strata == 0),]$subject_id
  if(length(idx) > 2){
    historical_subjects0 <- as.character(sample(idx, as.integer(0.5*length(idx))))
  } else {
    historical_subjects0 <- NULL
  }
  
  idx <- sub[which(sub$hypo_strata == i & sub$amine_strata == 1),]$subject_id
  if(length(idx) > 2){
    historical_subjects1 <- as.character(sample(idx, as.integer(0.5*length(idx))))
  } else {
    historical_subjects1 <- NULL
  }
  
  idx <- sub[which(sub$hypo_strata == i & sub$amine_strata == 2),]$subject_id
  if(length(idx) > 2){
    historical_subjects2 <- as.character(sample(idx, as.integer(0.5*length(idx))))
  } else {
    historical_subjects2 <- NULL
  }
  
  idx <- sub[which(sub$hypo_strata == i & sub$amine_strata == 3),]$subject_id
  if(length(idx) > 2){
    historical_subjects3 <- as.character(sample(idx, as.integer(0.5*length(idx))))
  } else {
    historical_subjects3 <- NULL
  }
  historical_subjects[[(i+1)]] <- list(historical_subjects0, 
                                       historical_subjects1,
                                       historical_subjects2,
                                       historical_subjects3)
  names(historical_subjects[[(i+1)]]) <- c(paste0("hypo", i, "_amine0"),
                                           paste0("hypo", i, "_amine1"),
                                           paste0("hypo", i, "_amine2"),
                                           paste0("hypo", i, "_amine3"))
}

historical_subjects <- unlist(historical_subjects, recursive = FALSE)
historical_subjects
# $hypo0_amine0
# [1] "5871"  "12267" "6145"  "26693" "1396"  "6145"  "16798" "21570" "26087" "20129"
# [11] "6145"  "19430" "3133"  "10995" "12231" "6145"  "17712" "16677" "17696" "19649"
# [21] "14187" "19330" "1650" 
# 
# $hypo0_amine1
# NULL
# 
# $hypo0_amine2
# NULL
# 
# $hypo0_amine3
# [1] "23782" "1818" 
# 
# $hypo1_amine0
# [1] "4136"  "21419" "10939" "5354"  "14167" "24559" "11372" "15464" "7632"  "14551"
# [11] "7487"  "18489" "10655" "14070" "8363"  "13720" "19618" "23318" "26069" "10651"
# [21] "18595" "24924" "14622" "20766" "4317"  "23907" "907"   "18982" "24431" "21481"
# [31] "9950"  "6233"  "13536" "26018" "9783"  "736"   "79"    "24822" "26156" "18377"
# [41] "11096" "5175"  "7849"  "17691" "17785" "19029" "22354" "19364" "24446" "19975"
# [51] "12187" "17748" "4565"  "9882"  "9268"  "6553"  "5496"  "22322" "15531" "13146"
# [61] "5198"  "23603" "8347"  "4270"  "15026" "4944"  "14486" "22859" "3821"  "13422"
# [71] "15270" "377"   "124"   "12461" "9595"  "408"   "22809" "24357" "328"   "24556"
# [81] "2265"  "1613"  "12171" "23321" "8272"  "6522"  "26604" "17865" "13837" "9426" 
# [91] "15208" "17092" "15052" "11907" "13099" "22766" "24556" "15302" "4805"  "16343"
# [101] "9001"  "18982" "24514" "8557"  "13096" "10651" "23922" "14057" "16337" "15817"
# [111] "19815" "17054" "8422"  "9426"  "263"   "8996"  "2834"  "22281" "8154"  "21305"
# [121] "20856" "15509" "22937" "2148"  "25271" "20181" "11787" "24556" "13868" "23869"
# [131] "3748"  "12663" "11684" "11641" "9389"  "17674" "14692" "16804"
# 
# $hypo1_amine1
# [1] "24666" "25328" "26233" "5818"  "18998" "7786"  "8723"  "10785" "23180" "20389"
# [11] "23270" "23627" "24577" "10419" "11698" "19125" "10384" "11945" "20860" "24562"
# [21] "16607" "22687" "14532" "1144"  "9112"  "7612"  "23459" "3884"  "16715" "11840"
# [31] "15514" "9430"  "20984" "7410"  "7886"  "6335"  "14123" "4059"  "5686"  "20"   
# [41] "4286" 
# 
# $hypo1_amine2
# [1] "2700"  "15333" "5237"  "15982" "21328" "6892"  "9685"  "3290"  "11763" "15631"
# [11] "14054" "7225"  "14561" "23292" "25627" "10320" "18681" "18269" "22306" "9330" 
# [21] "11244" "1586"  "9678"  "8890"  "16117" "10423" "4847"  "8138"  "4784"  "13002"
# [31] "3626"  "25107" "214"   "5879"  "18498" "6407"  "16881" "5343"  "25178" "7172" 
# [41] "25255" "16391" "11591" "10564" "11347" "3302"  "20448" "24152" "3218"  "9335" 
# [51] "25939" "23780" "9425"  "22603" "19513" "21317" "23749" "14922" "7542"  "11200"
# [61] "4852"  "18487" "16565" "26381" "1502"  "10710" "17847" "3995"  "26459" "9016" 
# 
# $hypo1_amine3
# [1] "19380" "5114"  "17690" "3914"  "4915"  "7183"  "25111" "17457" "4909"  "1449" 
# [11] "6455"  "20929" "1941"  "23015" "22138" "6673"  "17913" "16639" "7842"  "3863" 
# [21] "17803" "20246" "17216" "177"   "7705"  "18875" "2229"  "12586" "14899" "26054"
# [31] "22242" "21071" "22017" "19213" "22448" "24792" "22642" "14702" "9575"  "2361" 
# [41] "25721" "19620" "1950"  "3889"  "3992"  "7860"  "15279" "12807" "21504" "17722"
# [51] "8493"  "14584" "2514"  "3652"  "18952" "10925" "2598"  "13485" "12920" "17372"
# [61] "3171"  "3358"  "24837"
# 
# $hypo2_amine0
# [1] "3886"  "6254"  "25725" "20095" "19604" "17231" "18688"
# 
# $hypo2_amine1
# NULL
# 
# $hypo2_amine2
# [1] "3533"  "703"   "5784"  "10241"
# 
# $hypo2_amine3
# [1] "3365"  "16691" "6323"  "8013"  "17472" "14391" "19603" "25140" "23620"
historical_subjects <- unname(unlist(historical_subjects))
length(historical_subjects) # 357

################################## create folds ################################
make_data_by_type <- function(data, historical = TRUE){
  data <- data.table(data)
  if(historical) {
    data <- data[subject_id %in% historical_subjects, ]
  } else {
    data <- data[!(subject_id %in% historical_subjects), ]
  }
  data$subject_id <- as.numeric(as.character(data$subject_id))
  data <- setorder(data, subject_id, time_and_date)
  data$subject_id <- as.factor(as.character(data$subject_id))
  data <- data %>% group_by(id) %>% slice(11:(n()-46))
  return(data.table(data))
}
d <- make_data_by_type(data, historical = TRUE)
length(unique(d$id)) # 377
dim(d) # 885748    190
# 0.380309 reduction from 2329022 rows

set.seed(592)
folds <- make_folds(d, cluster_ids = d$subject_id, strata_ids = d$hypo_strata, 
                    fold_fun = folds_vfold)
save(folds, file = here::here("data", "folds.Rdata"), compress = TRUE)

########################## process historical data ################################
match_cols <- c("id", "time_and_date", "subject_id", "icustay_id", "hr", 
                "abpsys", "abpdias", "abpmean", "spo2", "init_time_and_date",
                "abpdias_locf", "abpsys_locf", "abpmean_locf", "spo2_locf", 
                "hr_locf", "amine", "sedation", "ventilation", "rank_icu", 
                "sex", "age", "sapsi_first", "sofa_first", "bmi", "care_unit", 
                "admission_type_descr", "row_locf", "min_elapsed", "min")
d <- data.table(d)[, match_cols, with = FALSE]

process_all_data <- function(data, outfile_individual, outfile_historical){
  
  # initialize historical and individual data
  data <- data.table(data)
  individual <- make_data_by_type(data, historical = FALSE)
  save(individual, file = here::here("Data", outfile_individual), compress = T)
  rm(individual)
  
  historical <- make_data_by_type(data, historical = TRUE)
  rm(data)
  
  comparison <- compare(historical[, match_cols, with = FALSE], d)
  print(paste0("compare result: ", comparison$result))
  rm(comparison)
  
  historical <- process_data(historical, strata=c("admission_type_descr","sex"))
  print(paste0("sum NA in imputed data: ", sum(is.na(historical))))
  save(historical, file = here::here("Data", outfile_historical), compress = T)
  rm(historical)
}

process_all_data(mimic30, "individual30.Rdata", "history30.Rdata")


load(here::here("data","mimic60.Rdata"))
process_historical_data(mimic60, "history60.Rdata")
rm(mimic60)

process_historical_data(mimic30_smooth_mean, "history30_smooth_mean.Rdata")
process_historical_data(mimic60_smooth_mean, "history60_smooth_mean.Rdata")
process_historical_data(mimic30_smooth_median, "history30_smooth_median.Rdata")
rm(mimic30_smooth_median)

load(here::here("data","mimic60_smooth_median.Rdata"))
process_historical_data(mimic60_smooth_median, "history60_smooth_median.Rdata")
rm(mimic60_smooth_median)



