load(here::here("data","mimic30.Rdata"))
# do we want to do this only for historical patients? no
# L1? just amine?
# AHE = 5 consecutive min <65? fine
########## subset data based on hypotensive and vasopressor strata #############

W_cols <- c("subject_id","rank_icu", "sex", "age", "sapsi_first", 
            "sofa_first", "bmi", "care_unit", "admission_type_descr")
subjects <- unique(d$subject_id)

make_cov_tbl <- function(L1_type = c("none", "amine", "all_min1")){
  min1_list <- lapply(subjects, function(x){
    min1 <- d[which(d$subject_id == x & d$min == 1),]
    if(L1_type == "amine"){
      init_amine <- ifelse(min1$amine == 1, 1, 0)
      final <- min1[, c(W_cols,"init_amine"), with=F]
    } else if(L1_type == "all_min1"){
      init_amine <- ifelse(min1$amine == 1, 1, 0)
      init_sedation <- ifelse(min1$sedation == 1, 1, 0)
      init_ventilation <- ifelse(min1$ventilation == 1, 1, 0)
      init_abpmean <- min1$abpmean
      init_abpdias <- min1$abpdias 
      init_abpsys <- min1$abpsys
      init_hr <- min1$hr
      init_spo2 <- min1$spo2
      L1 <- c("init_amine", "init_sedation", "init_ventilation", 
              "init_abpmean", "init_abpdias", "init_abpsys", "init_hr",
              "init_spo2")
      final <- min1[, c(W_cols, L1), with=F]
    } else if(L1_type == "none"){
      final <- min1[, W_cols, with=F]
    }
    return(final)
  })
  covs_tbl <- do.call(rbind, min1_list)
  return(covs_tbl)
}

make_data <- function(hour = 1, covs_tbl){
  d_reduced <- d %>%
    filter(min <= 60*hour & min > 1) %>%
    mutate(hypo = as.numeric(ifelse(abpmean < 65, 1, 0))) %>%
    mutate(hypo1 = lag(hypo,1,0)) %>%
    mutate(hypo2 = lag(hypo,2,0)) %>%
    mutate(hypo3 = lag(hypo,3,0)) %>%
    mutate(hypo4 = lag(hypo,4,0)) %>%
    mutate(hypo5 = lag(hypo,5,0)) %>%
    mutate(hypo_lag_sum = hypo1 + hypo2 + hypo3 + hypo4 + hypo5) %>%
    mutate(AHE = as.factor(ifelse(hypo_lag_sum == 5, "1", "0"))) %>%
    group_by(subject_id, AHE) %>%
    summarise(n = n()) %>%
    mutate(freq = round(n / sum(n), 2)) %>%
    select(c(subject_id, AHE))
  return(merge(covs_tbl, d_reduced, by = "subject_id"))
}

process_data(d_reduced, strata=c("admission_type_descr","sex"))

