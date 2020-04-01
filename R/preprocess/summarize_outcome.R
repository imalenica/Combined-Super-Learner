summary_function <- function(vector, function_name) { 
  FUN <- match.fun(function_name) 
  FUN(vector) 
}

## Solution 1 - hypotensive episode for time t is defined as either:
# 1. abpmean at time t < 65 mmHg and there is 5-minute window around time t 
#    (i.e., 10 time-points) with at least 5 time-points abpmean < 65.
#    Once the episodes are identified for the entire period, we summarize. 

Y_sol1 = function(vec, window=5, cutoff=62){
  
  init <- ifelse(vec < cutoff,1,0)
  
  #Sum the lags around each possible event.
  lags<-rowSums(as.data.frame(lapply(seq(1:window), function(x) lag(init,x))), na.rm = TRUE)
  leads<-rowSums(as.data.frame(lapply(seq(1:window), function(x) lead(init,x))), na.rm = TRUE)
  mids<-lags+leads
  
  event<-ifelse((init==1 & mids >= 5),1,
                ifelse((init==0 & mids >= 8),1,0))
  
  #But now summarize over the entire period:
  event<-ifelse(sum(event)>0,1,0)
  return(event)
}

create_summary <- function(init, summary_window, summary_measure){
  cap <- nrow(init)-summary_window

  Y_summary <- unlist(lapply(seq(cap), function(j){
    vec <- data.frame(init[(j):(j+summary_window-1),])[,1]
    summary_function(vec, summary_measure)
  }))
  
  return(Y_summary)
}

#summary_measure: one of the standard stats functions, or custom function (Y_sol1)
add_summary_outcome <- function(data, gap = 30, summary_window = 15, 
                                summary_measure = "median"){
  dat <- data %>%
  group_by(subject_id) %>%
  mutate(Y = dplyr::lead(abpmean, n = gap, default = NA)) %>%
  mutate(Y_summary = NA)
  
  for(i in unique(dat$subject_id)){
    init <- na.omit(dat[dat$subject_id == i, "Y"])
    Y_summary <- create_summary(init=init, summary_window=summary_window, 
                                summary_measure=summary_measure)
    #Pad the rest with NAs
    cap <- (nrow(dat[dat$subject_id == i,])-length(Y_summary))
    Y_summary <- c(Y_summary, rep(NA, cap))
    dat[dat$subject_id == i,"Y_summary"] <- Y_summary
  }

  return(dat)
}