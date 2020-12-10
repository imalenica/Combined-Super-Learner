library(data.table)
library(dplyr)
source(here("MIMICanalysis", "utils_plotting.R"))
load("~/Downloads/individual_mean.Rdata")

individual$id <- as.character(individual$id)
setorder(individual, subject_id, id, time_and_date)
ids <- unique(individual$id)

id_summaries <- lapply(ids, function(i){
  d <- individual[id == i,]
  W <- c("rank_icu", "sex", "age", "care_unit", "sapsi_first", "sofa_first", 
         "bmi", "admission_type_descr")
  if(any(is.na(d[, W, with=FALSE]))){
    return(NULL)
  } else {
    bp_tbl_summary <- make_bp_tbl_summary(d, abpmean_name = "abpmean_lag5_mean")
    AHE_tbl <- suppressMessages(d %>% 
      mutate(hour = as.integer(min/60) + 1) %>%
      group_by(hour) %>% 
      summarize(mean_AHE = mean(AHE)))
    hourly_summary_tbl <- merge(AHE_tbl, bp_tbl_summary, by = "hour")
    
    base <- c("id", "subject_id", W)
    tblW <- unique(d[, base, with = FALSE])
    missing_hours <- make_missing(bp_tbl_summary)
    missing_rows <- (max(d$min)-min(d$min)+1) - nrow(d)
    mean_AHE <- mean(hourly_summary_tbl$mean_AHE)
    tbl <- cbind(tblW, missing_hours, missing_rows, mean_AHE)
    tbl$bmi <- round(tbl$bmi, 2)
    return(list(hourly_summary_tbl = hourly_summary_tbl, W_tbl = tbl))
  }
})
names(id_summaries) <- ids
id_summaries_valid <- Filter(Negate(is.null), id_summaries)
hourly_summaries <- lapply(id_summaries_valid, '[[', 'hourly_summary_tbl')
W_summary_tbl <- do.call(rbind, lapply(id_summaries_valid, '[[', 'W_tbl'))

# let's only consider at most 212 subjects with at most 3 day stays
token <- dplyr::filter(W_summary_tbl, hour_count <= 72)
nrow(token) 

# ids with AHE (0.50, 1] of the time
AHE_high <- token %>% filter(mean_AHE > 0.5) %>% arrange(hour_count)
dput(AHE_high$id)
c("791", "9275", "26451", "3094", "20934", "10520", "26558", 
  "16621", "23986_1", "29585_2", "2698", "8042_2", "12477")

# ids with AHE (0.25, 0.5] of the time
AHE_med <- token %>% filter(mean_AHE <= 0.5 & mean_AHE > 0.25)  %>% arrange(hour_count)
dput(AHE_med$id)
c("33044", "32798", "2738_2", "22690", "2254", "25496", "26038", 
  "2790", "23543", "3078", "17420", "6075", "6464", "22587", "5544", 
  "29585_1", "24605_2", "12878", "2759", "9313_1", "28711", "9429", 
  "26520", "29334", "17435", "14822", "27887", "1830", "2904", 
  "2974", "2380", "28695", "8042_1")

# ids with AHE (0, 0.25] of the time
AHE_low <- token %>% filter(mean_AHE <= 0.25 & mean_AHE > 0)  %>% arrange(hour_count)
dput(AHE_low$id)
c("19657", "28647", "5392", "4718_1", "23245_2", "23661_1", "1577", 
  "20131", "29734", "5527", "27660", "32908", "4307", "20995", 
  "24999", "5609_1", "6035", "17456", "17659", "18803", "32433_1", 
  "4068", "16891", "168", "12410", "7361", "13439", "13910", "268", 
  "29299", "3963", "449", "6229", "8053", "12696", "15765", "19667", 
  "28660", "3417", "14095", "15114", "26296", "4060", "10240", 
  "14578", "17619", "27851", "29961", "16719", "21707_2", "24605_1", 
  "28430", "31323", "14550", "5452", "15771", "22172_1", "24900", 
  "32901_2", "5228", "690_1", "8626_1", "8966", "1065", "19509", 
  "2807", "12448", "21942_1", "22702", "3252", "6456", "20969", 
  "18654", "19971", "22723", "29687", "5077", "8279_1", "848", 
  "690_2", "10403_2", "30311_2", "1300", "15419", "30457", "17376_1", 
  "19794", "2320", "1100", "23151", "26834", "29977", "31833", 
  "15448", "27711", "32865", "15150", "15639", "18183", "32358", 
  "16400", "17570", "4825", "829", "941_2", "2194", "7280", "16377", 
  "23661_2", "19711", "32975", "11952", "11706", "7463", "27637", 
  "3985", "26736", "14526", "28566", "32959", "5465", "14739", 
  "31814", "6344", "864_1", "32901_1", "2432", "14775_1", "17482", 
  "29294_2", "4565", "1145", "15662", "19351", "4297", "21112", 
  "25353", "3181_2")

# ids with no AHE 
AHE_none <- token %>% filter(mean_AHE == 0)  %>% arrange(hour_count)
dput(AHE_none$id)
c("21707_1", "24832", "11548", "4362", "7826", "31262", "7934_1", 
  "18414", "4341", "18517", "11071", "15121", "32346", "18494", 
  "25296", "26325", "5058", "1749", "27580", "30357", "19349", 
  "4435", "30198", "14463", "20820", "25258", "4140", "9841_2")