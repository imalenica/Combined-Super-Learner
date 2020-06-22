library(here)
library(data.table)
library(tidyverse)
library(origami)

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

########## subset data based on hypotensive and vasopressor strata #############

# create hypotensive strata
data_hypo <- mimic30 %>%
  mutate(hypo = as.factor(ifelse(abpmean < 65, "1", "0"))) %>%
  group_by(subject_id, hypo) %>%
  summarise(n = n()) %>%
  mutate(freq = round(n / sum(n), 2)) %>%
  filter(hypo == 1) %>%
  mutate(hypo_strata = ifelse(freq < .5, 1, 2)) %>%
  select(c(subject_id, hypo_strata))
non_hypo <- unique(mimic30[-which(mimic30$subject_id %in% data_hypo$subject_id),]$subject_id)
remaining <- data.frame(subject_id = non_hypo, 
                        hypo_strata = rep(0, length(non_hypo)))
data_hypo <- rbind(data_hypo, remaining)
data_hypo$hypo_strata <- as.factor(as.character(data_hypo$hypo_strata))
data <- left_join(mimic30, data_hypo)

# create vasopressor strata
data_amine <- mimic30 %>%
  group_by(subject_id, amine) %>%
  summarise(n = n()) %>%
  mutate(freq = round(n / sum(n), 2)) %>%
  filter(amine == 1) %>%
  mutate(amine_strata = ifelse(freq < .5, 1,
                        ifelse(freq >= .5 & freq < 1, 2, 3))) %>%
  select(c(subject_id, amine_strata))
non_amine <- unique(mimic30[-which(mimic30$subject_id %in% data_amine$subject_id),]$subject_id)
remaining <- data.frame(subject_id = non_amine, 
                        amine_strata = rep(0, length(non_amine)))
data_amine <- rbind(data_amine, remaining)
data_amine$amine_strata <- as.factor(as.character(data_amine$amine_strata))

sub <- merge(data_amine, data_hypo, by = "subject_id", all = TRUE)

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
# [1] "22606" "11323" "3133"  "20303" "4771"  "12231" "5871"  "26693" "14822"
# [10] "19430" "20269" "7681"  "17696" "6145"  "16798" "2049"  "23590" "26087"
# [19] "10995" "11279" "3498"  "1650"  "14448" "15652" "17125" "16677" "14187"
# [28] "12267" "21570"
# 
# $hypo0_amine1
# NULL
# 
# $hypo0_amine2
# NULL
# 
# $hypo0_amine3
# [1] "23782" "1818"  "17929"
# 
# $hypo1_amine0
# [1] "787"   "24357" "18685" "11086" "2561"  "1378"  "22491" "11291" "14626"
# [10] "7894"  "2834"  "4778"  "14070" "26639" "16633" "10939" "4805"  "13716"
# [19] "4270"  "6204"  "3570"  "13308" "17748" "20986" "22980" "10926" "2264" 
# [28] "1408"  "18166" "21484" "9001"  "18846" "4406"  "26069" "6411"  "22859"
# [37] "7360"  "7683"  "79"    "22508" "4448"  "25271" "17285" "16161" "26133"
# [46] "11995" "17865" "9991"  "736"   "3138"  "14824" "12673" "6667"  "21305"
# [55] "8272"  "24431" "6553"  "20238" "23922" "13096" "5714"  "1222"  "12251"
# [64] "9732"  "22322" "23468" "5107"  "21138" "719"   "16337" "13536" "4481" 
# [73] "14551" "11801" "12663" "8832"  "17182" "4655"  "23869" "4870"  "408"  
# [82] "21517" "21873" "24559" "3266"  "14884" "17798" "9426"  "8269"  "8879" 
# [91] "12565" "5198"  "24822" "19364" "26024" "4565"  "15150" "18676" "18377"
# [100] "14692" "6464"  "6522"  "24514" "18696" "21419" "15817" "26097" "17948"
# [109] "11703" "22181" "14622" "24177" "1357"  "8347"  "4113"  "23344" "14579"
# [118] "12371" "24924" "19618" "11679" "22281" "24547" "11604" "23907" "9389" 
# [127] "15974" "14167" "15509" "13759" "7849"  "22809" "15903" "328"   "20856"
# [136] "23675" "26018" "10651" "21797" "9268"  "21817" "4862"  "13720" "10186"
# [145] "25073" "5175"  "3174"  "11096" "16455" "11907" "13099" "26576" "14057"
# [154] "21974" "14174" "907"   "6256"  "19891" "13033" "16343" "18489" "7632" 
# [163] "15302" "10973" "8422"  "15924" "25987" "17092" "5786"  "1414"  "18982"
# [172] "8557"  "19220"
# 
# $hypo1_amine1
# [1] "17629" "16038" "9178"  "5818"  "12113" "18088" "15701" "24746" "16607"
# [10] "25988" "4059"  "1046"  "20984" "20389" "24573" "9233"  "11945" "18695"
# [19] "7519"  "18988" "10552" "8368"  "14019" "369"   "23459" "7262"  "20546"
# [28] "6335"  "19208" "6889"  "3860"  "7410"  "3884"  "18998" "12112" "10785"
# [37] "12104" "4338"  "22389" "26519" "4401"  "9430"  "26233" "2332"  "8723" 
# [46] "21088" "23270" "16715" "7655"  "11840" "6478"  "20474" "19012" "26105"
# [55] "17516" "11318" "19977" "5686"  "2664"  "8989"  "20"   
# 
# $hypo1_amine2
# [1] "26459" "608"   "22241" "138"   "25699" "11380" "894"   "7910"  "16853"
# [10] "14131" "4847"  "24711" "23020" "22077" "14054" "15864" "5451"  "6407" 
# [19] "23591" "3675"  "14897" "17152" "10342" "11827" "23292" "21161" "16961"
# [28] "16122" "16565" "9425"  "20124" "22603" "25178" "3995"  "20922" "6892" 
# [37] "23617" "9615"  "12141" "3272"  "16511" "3798"  "4520"  "23401" "26575"
# [46] "24417" "21328" "5369"  "24597" "9335"  "11138" "3290"  "18219" "11591"
# [55] "23552" "7212"  "5879"  "13002" "7542"  "20354" "3302"  "6637"  "10611"
# [64] "15687" "21561" "925"   "4852"  "12215" "23584" "10205" "15631" "11877"
# [73] "18687" "1586"  "8249"  "8138"  "17810" "15333" "6535"  "21775" "20448"
# [82] "4266"  "10069" "217"   "565"   "22657" "22401" "13852" "22393" "24152"
# [91] "12400" "7172"  "9987"  "4656"  "16391"
# 
# $hypo1_amine3
# [1] "24792" "18358" "25222" "26054" "4252"  "19620" "17959" "10315" "2213" 
# [10] "1932"  "18875" "20658" "26688" "318"   "26380" "3171"  "21845" "1941" 
# [19] "2361"  "19734" "8532"  "15023" "3466"  "625"   "24837" "4451"  "5933" 
# [28] "23641" "8790"  "18952" "10209" "3358"  "17913" "14584" "6202"  "22017"
# [37] "12807" "12000" "6636"  "515"   "12920" "6349"  "439"   "13485" "5199" 
# [46] "17026" "668"   "24142" "3652"  "6042"  "10525" "2224"  "13195" "5960" 
# [55] "2229"  "17457" "9971"  "17803" "5114"  "24064" "3889"  "10925" "9575" 
# [64] "17372" "25016" "24457" "1244"  "13646" "15021" "549"   "177"   "22466"
# [73] "15124" "25851" "4685"  "4803"  "2075"  "1950"  "21710" "3279"  "17690"
# [82] "6944"  "22448" "26467" "3052"  "2619"  "11464" "5748"  "8461"  "18248"
# [91] "7860"  "3992"  "906"   "15569" "13101" "2395" 
# 
# $hypo2_amine0
# [1] "13993" "20095" "1944"  "19936" "3886"  "19604" "7651"  "2787"  "9043" 
# [10] "24922"
# 
# $hypo2_amine1
# [1] "871"  "9642"
# 
# $hypo2_amine2
# [1] "703"   "13353" "10241" "10485" "3533" 
# 
# $hypo2_amine3
# [1] "16691" "17028" "25140" "3365"  "1224"  "16827" "18738" "14391" "18401"
# [10] "19634" "507"   "2480" 
historical_subjects <- unname(unlist(historical_subjects))
length(historical_subjects) # 486

################################## create folds ################################
make_historical_data <- function(data){
  data <- setDT(data)[subject_id %in% historical_subjects, ]
  data$subject_id <- droplevels(data$subject_id)
  data$subject_id <- as.numeric(as.character(data$subject_id))
  data <- setorder(data, subject_id, time_and_date)
  data <- data %>% group_by(id) %>% slice(11:(n()-46))
  setDT(data)
}
d <- make_historical_data(data)
length(unique(d$id)) # 541
dim(d) # 2329022     190

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
d <- setDT(d)[,match_cols,with=FALSE]

process_historical_data <- function(history, outfile){
  history <- make_historical_data(history)
  
  comparison <- compare(history[, match_cols, with = FALSE], d)
  print(paste0("compare: ", comparison$result))
  
  history <- process_data(history, strata = c("admission_type_descr", "sex"))
  print(paste0("sum NA: ", sum(is.na(history))))
  save(history, file = here::here("data", outfile), compress = TRUE)
}

process_historical_data(mimic30, "history30.Rdata")
process_historical_data(mimic60, "history60.Rdata")
process_historical_data(mimic30_smooth_mean, "history30_smooth_mean.Rdata")
process_historical_data(mimic60_smooth_mean, "history60_smooth_mean.Rdata")
process_historical_data(mimic30_smooth_median, "history30_smooth_median.Rdata")
rm(mimic30_smooth_median)

load(here::here("data","mimic60_smooth_median.Rdata"))
process_historical_data(mimic60_smooth_median, "history60_smooth_median.Rdata")
rm(mimic60_smooth_median)

