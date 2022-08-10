library(data.table)
library(dplyr)
library(feasts)
library(tsibble)

box <- "/Users/Rachael/Library/CloudStorage/Box-Box/"
vitals_files <- list.files(paste0(box, "MIMIC-III/posl_numerics_trunc48hrs/"))

files <- sub(".", "", vitals_files)
subject_ids <- gsub("_.*","", files)
subject_ids <- as.numeric(subject_ids)

icustay_ids <- substring(files, 8)
icustay_ids <- gsub("_.*","", icustay_ids)
icustay_ids <- as.numeric(icustay_ids)
patients_keep <- data.table("SUBJECT_ID" = subject_ids, "ICUSTAY_ID" = icustay_ids)

########################## time-varying information ############################
items_dictionary <- read.csv(paste0(box, "MIMIC-III/D_ITEMS.csv"))[,-1]
setDT(items_dictionary)
items_dictionary <- items_dictionary[DBSOURCE != "carevue",]

##### ventilation
items_ventilation <- items_dictionary[CATEGORY == "2-Ventilation",]
procedures <- read.csv(paste0(box,"MIMIC-III/PROCEDUREEVENTS_MV.csv"))[,-1]
setDT(procedures)
keep_cols <- c("SUBJECT_ID", "HADM_ID", "ICUSTAY_ID", "STARTTIME", "ENDTIME")
procedures <- procedures[ITEMID %in% items_ventilation[["ITEMID"]], keep_cols, with=F]
procedures <- inner_join(patients_keep, procedures)
write.csv(procedures, file = paste0(box, "MIMIC-III/output/ventilation.csv"),
          row.names = F)
##### vasopressors:
trt <- c("Vasopressin", "Phenylephrine", "Norepinephrine", "Epinephrine", "Dopamine")
keep_cols <- c("SUBJECT_ID", "HADM_ID", "ICUSTAY_ID", "STARTTIME", "ENDTIME", "RATE")
items_trt <- items_dictionary[LABEL %in% trt,]
inputevents <- read.csv(paste0(box,"MIMIC-III/INPUTEVENTS_MV.csv"))[,-1]
setDT(inputevents)
inputevents <- inputevents[ITEMID %in% items_trt[["ITEMID"]],]
inputevents <- inner_join(patients_keep, inputevents)
##### norepinephrine
ne <- inputevents[ITEMID == items_trt[LABEL == "Norepinephrine",][["ITEMID"]],
                  keep_cols, with=F]
ne$norepinephrine_NEE <- ne$RATE
ne <- ne[,-"RATE",with=F]
ne <- inner_join(patients_keep, ne)
write.csv(ne, file = paste0(box, "MIMIC-III/output/norepinephrine.csv"),
          row.names = F)
##### vasopressin
vas <- inputevents[ITEMID == items_trt[LABEL == "Vasopressin",][["ITEMID"]],
                   c(keep_cols, "RATEUOM"), with=F]
vas$vasopressin_NEE <- ifelse(
  vas$RATEUOM == "units/min", (vas$RATE*0.01)/0.04, ((vas$RATE/60)*0.01)/0.04
)
vas <- vas[,-c("RATE","RATEUOM"), with=F]
vas <- inner_join(patients_keep, vas)
write.csv(vas, file = paste0(box, "MIMIC-III/output/vasopressin.csv"),
          row.names = F)
##### phenylephrine
phen <- inputevents[ITEMID == items_trt[LABEL == "Phenylephrine",][["ITEMID"]],
                    keep_cols, with=F]
phen$phenylephrine_NEE <- (phen$RATE*0.01)/1
phen <- phen[,-"RATE", with=F]
phen <- inner_join(patients_keep, phen)
write.csv(phen, file = paste0(box, "MIMIC-III/output/phenylephrine.csv"),
          row.names = F)
##### epinephrine
epi <- inputevents[ITEMID == items_trt[LABEL == "Epinephrine",][["ITEMID"]],
                   keep_cols, with=F]
epi$epinephrine_NEE <- epi$RATE
epi <- epi[,-"RATE", with=F]
epi <- inner_join(patients_keep, epi)
write.csv(epi, file = paste0(box, "MIMIC-III/output/epinephrine.csv"),
          row.names = F)
##### dopamine
dop <- inputevents[ITEMID == items_trt[LABEL == "Dopamine",][["ITEMID"]],
                   keep_cols, with=F]
dop$dopamine_NEE <- (dop$RATE*0.01)/15
dop <- dop[,-"RATE", with=F]
dop <- inner_join(patients_keep, dop)
write.csv(dop, file = paste0(box, "MIMIC-III/output/dopamine.csv"),
          row.names = F)
