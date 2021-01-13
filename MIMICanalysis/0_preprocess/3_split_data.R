library(here)
library(data.table)
library(tidyverse)
library(origami)
source(here("R", "v3", "impute_data.R"))

#################### make sure data variations have same ids  ##################
load(here::here("Data","mimic30_smooth_median.Rdata"))
load(here::here("Data","mimic60_smooth_median.Rdata"))
mimic30_smooth_median_id <- unique(mimic30_smooth_median$id)
rm(mimic30_smooth_median)
mimic60_smooth_median_id <- unique(mimic60_smooth_median$id)
rm(mimic60_smooth_median)

load(here::here("Data","mimic30_smooth_mean.Rdata"))
load(here::here("Data","mimic60_smooth_mean.Rdata"))
mimic30_smooth_mean_id <- unique(mimic30_smooth_mean$id)
rm(mimic30_smooth_mean)
mimic60_smooth_mean_id <- unique(mimic60_smooth_mean$id)
rm(mimic60_smooth_mean)

load(here::here("Data","mimic30.Rdata"))
load(here::here("Data","mimic60.Rdata"))
mimic30_id <- unique(mimic30$id)
mimic60_id <- unique(mimic60$id)
rm(mimic60)

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
token_subjects <- dplyr::filter(sub, amine_strata != 0 & hypo_strata != 0)$subject_id

# take 50% of stratified subsets
historical_subjects_list <- list()
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
  historical_subjects_list[[(i+1)]] <- list(historical_subjects0,
                                            historical_subjects1,
                                            historical_subjects2,
                                            historical_subjects3)
  names(historical_subjects_list[[(i+1)]]) <- c(paste0("hypo", i, "_amine0"),
                                                paste0("hypo", i, "_amine1"),
                                                paste0("hypo", i, "_amine2"),
                                                paste0("hypo", i, "_amine3"))
}

historical_subjects_list <- unlist(historical_subjects_list, recursive = FALSE)
historical_subjects_list
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
historical_subjects <- unique(unname(unlist(historical_subjects_list)))
length(historical_subjects) # 349
dput(historical_subjects)
# c("5871", "12267", "6145", "26693", "1396", "16798", "21570",
#   "26087", "20129", "19430", "3133", "10995", "12231", "17712",
#   "16677", "17696", "19649", "14187", "19330", "1650", "23782",
#   "1818", "4136", "21419", "10939", "5354", "14167", "24559", "11372",
#   "15464", "7632", "14551", "7487", "18489", "10655", "14070",
#   "8363", "13720", "19618", "23318", "26069", "10651", "18595",
#   "24924", "14622", "20766", "4317", "23907", "907", "18982", "24431",
#   "21481", "9950", "6233", "13536", "26018", "9783", "736", "79",
#   "24822", "26156", "18377", "11096", "5175", "7849", "17691",
#   "17785", "19029", "22354", "19364", "24446", "19975", "12187",
#   "17748", "4565", "9882", "9268", "6553", "5496", "22322", "15531",
#   "13146", "5198", "23603", "8347", "4270", "15026", "4944", "14486",
#   "22859", "3821", "13422", "15270", "377", "124", "12461", "9595",
#   "408", "22809", "24357", "328", "24556", "2265", "1613", "12171",
#   "23321", "8272", "6522", "26604", "17865", "13837", "9426", "15208",
#   "17092", "15052", "11907", "13099", "22766", "15302", "4805",
#   "16343", "9001", "24514", "8557", "13096", "23922", "14057",
#   "16337", "15817", "19815", "17054", "8422", "263", "8996", "2834",
#   "22281", "8154", "21305", "20856", "15509", "22937", "2148",
#   "25271", "20181", "11787", "13868", "23869", "3748", "12663",
#   "11684", "11641", "9389", "17674", "14692", "16804", "24666",
#   "25328", "26233", "5818", "18998", "7786", "8723", "10785", "23180",
#   "20389", "23270", "23627", "24577", "10419", "11698", "19125",
#   "10384", "11945", "20860", "24562", "16607", "22687", "14532",
#   "1144", "9112", "7612", "23459", "3884", "16715", "11840", "15514",
#   "9430", "20984", "7410", "7886", "6335", "14123", "4059", "5686",
#   "20", "4286", "2700", "15333", "5237", "15982", "21328", "6892",
#   "9685", "3290", "11763", "15631", "14054", "7225", "14561", "23292",
#   "25627", "10320", "18681", "18269", "22306", "9330", "11244",
#   "1586", "9678", "8890", "16117", "10423", "4847", "8138", "4784",
#   "13002", "3626", "25107", "214", "5879", "18498", "6407", "16881",
#   "5343", "25178", "7172", "25255", "16391", "11591", "10564",
#   "11347", "3302", "20448", "24152", "3218", "9335", "25939", "23780",
#   "9425", "22603", "19513", "21317", "23749", "14922", "7542",
#   "11200", "4852", "18487", "16565", "26381", "1502", "10710",
#   "17847", "3995", "26459", "9016", "19380", "5114", "17690", "3914",
#   "4915", "7183", "25111", "17457", "4909", "1449", "6455", "20929",
#   "1941", "23015", "22138", "6673", "17913", "16639", "7842", "3863",
#   "17803", "20246", "17216", "177", "7705", "18875", "2229", "12586",
#   "14899", "26054", "22242", "21071", "22017", "19213", "22448",
#   "24792", "22642", "14702", "9575", "2361", "25721", "19620",
#   "1950", "3889", "3992", "7860", "15279", "12807", "21504", "17722",
#   "8493", "14584", "2514", "3652", "18952", "10925", "2598", "13485",
#   "12920", "17372", "3171", "3358", "24837", "3886", "6254", "25725",
#   "20095", "19604", "17231", "18688", "3533", "703", "5784", "10241",
#   "3365", "16691", "6323", "8013", "17472", "14391", "19603", "25140",
#   "23620")

token_individuals <- token_subjects[-which(token_subjects %in% historical_subjects)]
length(token_individuals) # 191
token_individuals_tbl <- data.table(sub[sub$subject_id %in% token_individuals,])
setorder(token_individuals_tbl, -hypo_strata, -amine_strata)
token_individuals_tbl <- data.frame(left_join(token_individuals_tbl,
                                              maxtime_by_subject))
token_individuals_tbl$max_hr <- round(token_individuals_tbl$max_hr, 2)
token_individuals_tbl
# subject_id        amine_strata hypo_strata  id max_hr
# 1         1485            3           2    1830  64.07
# 2        16827            3           2   20934  18.88
# 3        18401            3           2   22910  76.28
# 4        18738            3           2   23322  41.27
# 5        21270            3           2   26451   9.87
# 6         2480            3           2    3078  29.38
# 7         6358            3           2    7900  77.78
# 8          650            3           2     791   6.75
# 9         8749            3           2   10885  26.72
# 10       13353            2           2   16621  40.58
# 11       18925            2           2   23543  28.80
# 12       23811            2           2 29585_1  40.45
# 13       23811            2           2 29585_2  44.35
# 14        5205            2           2    6464  31.43
# 15        8467            2           2   10520  20.80
# 16       10013            1           2   12477  51.45
# 17        2185            1           2    2698  44.87
# 18       10525            3           1   13112  31.22
# 19       11764            3           1   14660  84.08
# 20       11855            3           1 14775_1  66.97
# 21        1279            3           1    1577   8.18
# 22        1313            3           1    1618  77.20
# 23       13149            3           1   16377  46.77
# 24       13646            3           1   16990  69.30
# 25       14036            3           1   17456  16.78
# 26       14059            3           1   17482  67.13
# 27       15021            3           1   18654  31.90
# 28       15569            3           1   19351  68.05
# 29       16127            3           1 20045_2  16.67
# 30       16196            3           1   20131   9.40
# 31       16873            3           1   20995  15.52
# 32       16876            3           1   20998  74.78
# 33       17026            3           1   21199  87.53
# 34       17072            3           1   21260  30.12
# 35       17959            3           1   22364  88.10
# 36       18229            3           1   22702  28.88
# 37       18248            3           1   22723  33.45
# 38       18358            3           1   22858  29.48
# 39        1932            3           1    2380  68.07
# 40       19734            3           1 24605_1  24.23
# 41       19734            3           1 24605_2  40.92
# 42       20689            3           1   25759  83.33
# 43       21483            3           1   26736  55.42
# 44        2224            3           1    2759  44.57
# 45       22429            3           1   27887  60.93
# 46       23048            3           1   28647   4.35
# 47       23105            3           1   28711  47.75
# 48        2340            3           1    2912  87.93
# 49       23934            3           1   29734  11.47
# 50        2395            3           1    2974  64.40
# 51       24124            3           1   29961  22.95
# 52       24457            3           1   30357  23.70
# 53       25222            3           1   31323  24.68
# 54       25621            3           1   31833  41.65
# 55       25851            3           1   32125  90.15
# 56        2619            3           1    3252  29.20
# 57       26380            3           1   32798   9.58
# 58       26435            3           1   32865  42.20
# 59       26467            3           1 32901_1  64.20
# 60       26467            3           1 32901_2  27.03
# 61        3052            3           1    3787  28.43
# 62         318            3           1     384  44.47
# 63        3279            3           1    4068  17.65
# 64        3883            3           1    4825  44.03
# 65        4076            3           1    5077  33.03
# 66        4252            3           1    5291  25.65
# 67        4329            3           1    5380   4.67
# 68        4362            3           1    5419  40.57
# 69        4451            3           1    5527  11.85
# 70        4893            3           1    6075  30.17
# 71        5199            3           1    6456  29.37
# 72         549            3           1     674  69.47
# 73        5933            3           1    7394  76.52
# 74        6202            3           1    7717  77.50
# 75        6636            3           1    8239  44.43
# 76         695            3           1     848  33.73
# 77         710            3           1   864_1  63.55
# 78        7213            3           1    8966  27.47
# 79        7897            3           1  9841_1   2.43
# 80        7897            3           1  9841_2  72.63
# 81        8231            3           1   10240  21.80
# 82        8790            3           1   10930  43.38
# 83        9139            3           1   11377  73.30
# 84        9593            3           1   11952  49.95
# 85        1006            2           1    1246  70.55
# 86       10205            2           1   12696  20.98
# 87       10342            2           1   12878  43.07
# 88       10611            2           1   13209  88.88
# 89       10766            2           1   13439  19.98
# 90       11138            2           1   13910  20.52
# 91       11827            2           1   14739  60.55
# 92       11877            2           1   14822  59.90
# 93       12141            2           1   15150  43.62
# 94       12294            2           1   15333  92.92
# 95       12400            2           1   15448  42.32
# 96       12679            2           1   15771  25.83
# 97       13438            2           1   16725  80.83
# 98         138            2           1     168  18.70
# 99       13970            2           1 17376_1  39.07
# 100      13970            2           1 17376_2  75.20
# 101      13970            2           1   17377  60.25
# 102      14005            2           1   17420  30.63
# 103      14205            2           1   17659  17.28
# 104      15821            2           1   19657   4.08
# 105      16853            2           1   20969  29.95
# 106      16961            2           1   21112  69.22
# 107      17810            2           1 22172_1  27.48
# 108      18219            2           1   22690  20.47
# 109      18597            2           1   23151  40.93
# 110      18687            2           1   23263   6.28
# 111       1973            2           1    2432  66.22
# 112      20268            2           1   25258  70.22
# 113      20922            2           1   26038  24.60
# 114      21011            2           1   26151  22.42
# 115      21050            2           1   26196  65.80
# 116      21161            2           1   26325  20.22
# 117      21561            2           1   26834  41.55
# 118        217            2           1     268  19.78
# 119      21809            2           1   27127  75.47
# 120      21857            2           1   27190  80.68
# 121      22122            2           1   27513  91.77
# 122      22221            2           1   27637  52.97
# 123      22241            2           1   27660  14.28
# 124      22401            2           1   27851  23.20
# 125       2251            2           1    2790  25.83
# 126      22879            2           1   28430  23.92
# 127      23092            2           1   28695  69.27
# 128      23584            2           1   29299  20.68
# 129      23617            2           1   29334  49.15
# 130      24417            2           1 30311_1   3.23
# 131      24417            2           1 30311_2  37.55
# 132      24532            2           1   30457  38.45
# 133      24711            2           1   30684  36.60
# 134      26506            2           1   32959  56.98
# 135       3192            2           1    3963  20.18
# 136       3214            2           1    3985  53.92
# 137       3272            2           1    4060  22.50
# 138       3462            2           1    4297  68.75
# 139       3473            2           1    4307  13.98
# 140       3619            2           1  4498_1   0.97
# 141       3619            2           1  4498_2  28.70
# 142       3623            2           1    4502  22.60
# 143       3675            2           1    4565  66.77
# 144       4194            2           1    5228  26.95
# 145       4393            2           1    5452  24.85
# 146       5369            2           1  6671_1  92.18
# 147       5451            2           1    6777  14.02
# 148        565            2           1   690_1  26.92
# 149        565            2           1   690_2  36.52
# 150       7297            2           1    9074  17.63
# 151       7567            2           1    9429  47.08
# 152       7910            2           1    9860  11.38
# 153       7977            2           1    9936  27.90
# 154       8040            2           1   10008  30.13
# 155        894            2           1    1100  40.03
# 156        925            2           1    1145  66.87
# 157       9393            2           1   11706  51.13
# 158       9615            2           1   11974  93.72
# 159       9987            2           1   12448  28.13
# 160       1046            1           1    1300  38.53
# 161      10534            1           1 13124_2  39.28
# 162      11161            1           1   13935  91.75
# 163      11318            1           1   39274  36.85
# 164      12104            1           1   15104  80.60
# 165      12112            1           1   15114  21.93
# 166      13171            1           1   16400  44.38
# 167      13333            1           1   16599  91.62
# 168      13435            1           1   16719  23.90
# 169      13570            1           1   16891  18.07
# 170      14019            1           1   17435  53.62
# 171      15144            1           1   18803  17.00
# 172      15701            1           1   19509  27.83
# 173      17516            1           1   21791  76.10
# 174      18088            1           1   22529  17.00
# 175      18695            1           1   23271  86.22
# 176       1885            1           1    2320  39.05
# 177      19012            1           1 23661_1   6.43
# 178      19012            1           1 23661_2  47.18
# 179      19918            1           1   24832   4.70
# 180      19977            1           1   24900  27.48
# 181      20474            1           1   25496  24.33
# 182      22285            1           1   27711  41.93
# 183      23060            1           1   28660  21.53
# 184       2332            1           1    2904  64.47
# 185      23580            1           1 29294_2  67.45
# 186      23890            1           1   29687  33.02
# 187      24573            1           1   30507  95.48
# 188      25603            1           1   31814  60.33
# 189      25988            1           1   32290  27.93
# 190      26043            1           1   32358  43.67
# 191       2754            1           1    3417  21.37
# 192        369            1           1     449  20.28
# 193       4338            1           1    5392   4.57
# 194       4401            1           1    5465  59.62
# 195       5023            1           1    6229  20.20
# 196       5163            1           1    6411  74.30
# 197       5847            1           1    7280  45.78
# 198       6478            1           1    8053  20.68
# 199        772            1           1   941_2  45.70
# 200       8207            1           1 10206_1  34.42
# 201       9372            1           1   11679  51.25
selected <- c(token_individuals, historical_subjects)
other_individuals <- maxtime_by_subject %>%
  filter(!(subject_id %in% selected)) %>%
  filter(max_hr <= 24*3) # only consider other subjects with at most 3 days data
other_individuals <- as.character(other_individuals$subject_id)
individual_subjects <- c(other_individuals, token_individuals)
dput(individual_subjects)
# c("8249", "8269", "8368", "865", "8568", "871", "8879", "8905",
#   "9043", "9176", "9269", "9630", "9642", "981", "9962", "10186",
#   "1028", "10250", "10926", "11205", "11279", "11291", "11323",
#   "11604", "11658", "11679", "11703", "11850", "123", "12116",
#   "12217", "12251", "12277", "12371", "12536", "12565", "12581",
#   "12673", "13033", "13195", "13308", "13354", "1357", "1378",
#   "13716", "13759", "1408", "13993", "1414", "1418", "14131", "14174",
#   "14448", "14579", "14626", "14822", "14884", "14900", "14938",
#   "15013", "15150", "15329", "15426", "15567", "15652", "15831",
#   "15864", "15924", "15965", "16055", "16071", "16122", "16122",
#   "16210", "16399", "16561", "16581", "16633", "16740", "16915",
#   "17028", "17443", "17443", "17589", "1778", "17629", "17702",
#   "17798", "17920", "17948", "1824", "18126", "1854", "1861", "18676",
#   "18685", "18846", "19102", "1944", "19220", "19246", "19655",
#   "19848", "20062", "20124", "20269", "20303", "2049", "20354",
#   "20459", "208", "2075", "20763", "20986", "21108", "21138", "21147",
#   "21156", "21265", "21321", "21349", "21484", "2172", "21710",
#   "21766", "21797", "21805", "21817", "21974", "2213", "22181",
#   "2264", "22606", "22657", "22669", "22983", "23001", "23020",
#   "23047", "23299", "23344", "23351", "23590", "23619", "24133",
#   "24314", "24460", "24597", "24730", "24748", "2479", "2492",
#   "25168", "25207", "2561", "25659", "25699", "26024", "26037",
#   "26105", "26382", "26472", "26511", "26519", "26576", "26639",
#   "26688", "2755", "2787", "3097", "3138", "3174", "17125", "3266",
#   "3321", "3340", "3498", "3513", "3570", "3798", "3986", "4064",
#   "4113", "4248", "4465", "4481", "4520", "462", "4655", "4771",
#   "4808", "4853", "4862", "4894", "5107", "5506", "5642", "5714",
#   "5908", "5995", "6107", "6214", "6256", "6294", "6382", "6411",
#   "6470", "6667", "682", "6944", "6983", "7136", "7262", "7289",
#   "7289", "7347", "7360", "7448", "7478", "7651", "7681", "7683",
#   "10013", "1006", "10205", "10342", "1046", "10525", "10534",
#   "10611", "10766", "11138", "11161", "11318", "11764", "11827",
#   "11855", "11877", "12104", "12112", "12141", "12294", "12400",
#   "12679", "1279", "1313", "13149", "13171", "13333", "13353",
#   "13435", "13438", "13570", "13646", "138", "13970", "14005",
#   "14019", "14036", "14059", "14205", "1485", "15021", "15144",
#   "15569", "15701", "15821", "16127", "16196", "16827", "16853",
#   "16873", "16876", "16961", "17026", "17072", "17516", "17810",
#   "17959", "18088", "18219", "18229", "18248", "18358", "18401",
#   "18597", "18687", "18695", "18738", "1885", "18925", "19012",
#   "1932", "1973", "19734", "19918", "19977", "20268", "20474",
#   "20689", "20922", "21011", "21050", "21161", "21270", "21483",
#   "21561", "217", "21809", "2185", "21857", "22122", "22221", "2224",
#   "22241", "22285", "22401", "22429", "2251", "22879", "23048",
#   "23060", "23092", "23105", "2332", "2340", "23580", "23584",
#   "23617", "23811", "23890", "23934", "2395", "24124", "24417",
#   "24457", "24532", "24573", "24711", "2480", "25222", "25603",
#   "25621", "25851", "25988", "26043", "2619", "26380", "26435",
#   "26467", "26506", "2754", "3052", "318", "3192", "3214", "3272",
#   "3279", "3462", "3473", "3619", "3623", "3675", "369", "3883",
#   "4076", "4194", "4252", "4329", "4338", "4362", "4393", "4401",
#   "4451", "4893", "5023", "5163", "5199", "5205", "5369", "5451",
#   "549", "565", "5847", "5933", "6202", "6358", "6478", "650",
#   "6636", "695", "710", "7213", "7297", "7567", "772", "7897",
#   "7910", "7977", "8040", "8207", "8231", "8467", "8749", "8790",
#   "894", "9139", "925", "9372", "9393", "9593", "9615", "9987")

################################## create folds ################################
make_data_by_type <- function(data, historical = TRUE){
  data$subject_id <- as.character(data$subject_id)
  data <- data.table(data)
  if(historical) {
    data <- data[data$subject_id %in% historical_subjects,]
  } else {
    data <- data[data$subject_id %in% individual_subjects,]
  }
  data$subject_id <- as.numeric(data$subject_id)
  setorder(data, subject_id, time_and_date)
  data$subject_id <- as.factor(as.character(data$subject_id))
  data$id <- as.factor(as.character(data$id))
  data <- data %>% group_by(id) %>% slice(11:(n()-46))
  return(data.table(data))
}
d <- make_data_by_type(data, historical = TRUE)
length(unique(d$id)) # 377
dim(d) # 885748    190
# 0.380309 reduction from 2329022 rows

set.seed(58)
folds <- make_folds(d, cluster_ids = d$subject_id, strata_ids = d$hypo_strata,
                    fold_fun = folds_vfold)
save(folds, file = here::here("Data", "folds.Rdata"), compress = TRUE)

# no. AHE in 1st ~6 hour (rf varimp + W)
# motivating Q: can we assign trt / assess status at baseline?
# prop time <65
# people who enter ICU on vasopressors
# Y = prop time < 65 -- dplyr groupby, mutate

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

  historical <- impute_data(historical, strata=c("admission_type_descr","sex"))
  print(paste0("sum NA in imputed data: ", sum(is.na(historical))))
  save(historical, file = here::here("Data", outfile_historical), compress = T)
  rm(historical)
}

process_all_data(mimic30, "individual30.Rdata", "historical30.Rdata")
rm(mimic30)

load(here::here("data","mimic60.Rdata"))
process_all_data(mimic60,  "individual30.Rdata", "historical60.Rdata")
rm(mimic60)

load(here::here("data","mimic30_smooth_mean.Rdata"))
process_all_data(mimic30_smooth_mean,
                 "individual30_smooth_mean.Rdata",
                 "historical30_smooth_mean.Rdata")
rm(mimic30_smooth_mean)

load(here::here("data","mimic60_smooth_mean.Rdata"))
process_all_data(mimic60_smooth_mean,
                 "individual60_smooth_mean.Rdata",
                 "historical60_smooth_mean.Rdata")
rm(mimic60_smooth_mean)

load(here::here("data","mimic30_smooth_median.Rdata"))
process_all_data(mimic30_smooth_median,
                 "individual30_smooth_median.Rdata",
                 "historical30_smooth_median.Rdata")
rm(mimic30_smooth_median)

load(here::here("data","mimic60_smooth_median.Rdata"))
process_all_data(mimic60_smooth_median,
                 "individual60_smooth_median.Rdata",
                 "historical60_smooth_median.Rdata")
rm(mimic60_smooth_median)
