library(data.table)
source(here("MIMICanalysis", "utils_plotting.R"))
load("~/Downloads/individual_mean-2.Rdata")

individual$id <- as.character(individual$id)
setorder(individual, subject_id, id, time_and_date)

# pick token individuals 
ids <- unique(individual$id)
id_summaries <- lapply(ids, function(i){
  d <- individual[id == i,]
  W <- c("rank_icu", "sex", "age", "care_unit", "sapsi_first", "sofa_first", 
         "bmi", "admission_type_descr")
  if(any(is.na(d[, W, with=FALSE]))){
    return(NULL)
  } else {
    bp_tbl_summary <- make_bp_tbl_summary(d, abpmean_name = "Y45_lag5_mean")
    AHE <- suppressMessages(d %>% 
      dplyr::mutate(hour = as.integer(min/60) + 1) %>%
      dplyr::group_by(hour) %>% 
      dplyr::summarize(mean_AHE = mean(AHE)))
    hourly_summary_tbl <- merge(AHE, bp_tbl_summary, by = "hour")
    
    base <- c("id", "subject_id", W)
    tblW <- unique(d[, base, with = FALSE])
    missing <- make_missing(bp_tbl_summary)
    AHE <- d %>% dplyr::summarize(mean_AHE = mean(AHE))
    tbl <- cbind(tblW, missing, AHE)
    tbl$bmi <- round(tbl$bmi, 2)
    return(list(hourly_summary_tbl = hourly_summary_tbl, W_tbl = tbl))
  }
})
names(id_summaries) <- ids

id_summaries_valid <- Filter(Negate(is.null), id_summaries)
hourly_summaries <- lapply(id_summaries_valid, '[[', 'hourly_summary_tbl')
W_summary_tbl <- do.call(rbind, lapply(id_summaries_valid, '[[', 'W_tbl'))
ids_AHE_tbl <- W_summary_tbl %>% dplyr::filter(mean_AHE > 0.5)
dput(ids_AHE_tbl$id)
c("12477", "16621", "20934", "23986_1", "26451", "26558", "2698", 
  "29585_2", "3094", "8042_2", "791", "9275", "10520")