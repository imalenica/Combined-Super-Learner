summary_function <- function(vector, function_name) { 
  FUN <- match.fun(function_name) 
  FUN(vector) 
}

create_summary <- function(init, summary_window, summary_measure){
  cap <- nrow(init)-summary_window

  Y_summary <- unlist(lapply(seq(cap), function(j){
    vec <- data.frame(init[(j):(j+summary_window-1),])[,1]
    summary_function(vec, summary_measure)
  }))
  
  return(Y_summary)
}

add_summary_outcome <- function(data, gap = 30, summary_window = 15, 
                                summary_measure = "median"){
  dat <- data %>%
  group_by(subject_id) %>%
  mutate(Y = dplyr::lead(abpmean, n = gap, default = NA)) %>%
  mutate(Y_summary = NA)

  lapply(unique(dat$subject_id), function(i){
    init <- na.omit(dat[dat$subject_id == i, "Y"])
    Y_summary <- create_summary(init=init, summary_window=summary_window, 
                                summary_measure=summary_measure)
    #Pad the rest with NAs
    cap <- (nrow(dat[dat$subject_id == i,])-length(Y_summary))
    Y_summary <- c(Y_summary, rep(NA, cap))
    dat[dat$subject_id == i,"Y_summary"] <- out
  })

  return(dat)
}