## Sample n individuals with t periods
sample_n_t = function(df,n,t){

  #Pick t periods and find samples that have data for those periods:
  train_all_p <- lapply(1:t, function(x) df[df$periode==x,])
  samples_intersect <- Reduce(intersect,lapply(1:t, function(x) train_all_p[[x]]$subject_id))

  #Subset to n samples:
  samples<-sample(samples_intersect, n)
  train_all_c <- lapply(1:t, function(x) train_all_p[[x]][train_all_p[[x]]$subject_id %in% samples,])
  train_all <- do.call("rbind", train_all_c)

  train_all <- train_all %>% dplyr::group_by(subject_id) %>%
                dplyr::arrange(time,.by_group = TRUE) %>%
                  dplyr::arrange(periode,.by_group = TRUE)

  #PROBLEM: some samples are repeated?
  train_all <- train_all %>% distinct(subject_id, time, periode, .keep_all = TRUE)

  #Set time to be continuous
  train_all <- train_all %>% dplyr::group_by(subject_id) %>% dplyr::mutate(time = 1:(t*60))

  return(train_all)
}

### Create a new outcome

## Solution 1

new_Y_sol1 = function(train_all, window=5, cutoff=62){

  train_all$Y1<-rep(0,nrow(train_all))

  for(i in unique(train_all$subject_id)){

    sub<-train_all[train_all$subject_id==i,]

    init <- ifelse(sub[,"abpmean"] < cutoff,1,0)

    #Sum the lags around each possible event.
    lags<-rowSums(as.data.frame(lapply(seq(1:window), function(x) lag(init,x))), na.rm = TRUE)
    leads<-rowSums(as.data.frame(lapply(seq(1:window), function(x) lead(init,x))), na.rm = TRUE)
    mids<-lags+leads

    event<-ifelse((init==1 & mids >= 5),1,0)
    #test<-cbind.data.frame(init,lags,leads,mids,event)
    train_all[train_all$subject_id==i,"Y1"]<-event
  }

  return(train_all)
}

## Solution 2

# Acute Hypothensive Episode (AHE) is a period of WINDOW minutes
# during which at least 80% of the MAP measurements are
# no greater than CUTOFF.

new_Y_sol2 = function(train_all, window=5, cutoff=62){

  train_all$Y2<-rep(0,nrow(train_all))

  for(i in unique(train_all$subject_id)){

    sub<-train_all[train_all$subject_id==i,]
    t<-dim(sub)[1]
    windows <- seq(from = 1, to=t-window, by=window+1)

    for(j in windows){

      #Look at the "window"" minute interval:
      init<-ifelse(sub[(j:(j+window)),"abpmean"] < cutoff,1,0)

      if(sum(init)>(0.90*window)){
        sub[(j:(j+window)),"Y2"]<-data.frame(rep(1,length(init)))
      }else{
        sub[(j:(j+window)),"Y2"]<-data.frame(rep(0,length(init)))
      }
    }

    #sub[is.na(sub$Y),"Y"]<-ifelse(sub[is.na(sub$Y),"abpmean"]<cutoff,1,0)
    #Too short to be an episode- even if < cutoff
    sub[is.na(sub$Y2),"Y2"]<-0
    train_all[train_all$subject_id==i,]<-sub
  }

  return(train_all)
}

run_class <- function(df, cols_fac, cols_num) {
  cols_all <- c(cols_fac, cols_num)
  df[cols_all] <- lapply(df[cols_all], as.character)
  df[cols_num] <- lapply(df[cols_num], as.numeric)
  df[cols_fac] <- lapply(df[cols_fac], as.factor)
  return(df)
}

## Use proper CV structure (to be included in origami) 
folds_rolling_origin_pooled <- function(n, t, first_window, validation_size,
                                        gap = 0, batch = 1) {
  
  message(paste("Processing", n/t, "samples with", t, "time points."))
  
  #Index the observations
  dat <- cbind.data.frame(index=seq(n),time=rep(seq(t),n/t),id=rep(seq(n/t), each=t))
  ids <- unique(dat$id)
  
  # establish rolling origin forecast for time-series cross-validation
  rolling_origin_skeleton <- folds_rolling_origin(t, first_window,
                                                  validation_size, gap, batch)
  
  folds_rolling_origin <- lapply(rolling_origin_skeleton, function(h){
    train_indices <- lapply(ids, function(i){
      train <- dat[dat$id == i, ]
      train[h$training_set, ]$index
    })
    val_indices <- lapply(ids, function(j){
      val <- dat[dat$id == j, ]
      val[h$validation_set, ]$index
    })
    make_fold(v=h$v, training_set=unlist(train_indices), 
              validation_set=unlist(val_indices))
  })
  return(folds_rolling_origin)
}

folds_rolling_window_pooled <- function(n, t, window_size, validation_size,
                                        gap = 0, batch = 1) {
  
  message(paste("Processing", n/t, "samples with", t, "time points."))
  
  #Index the observations
  dat <- cbind.data.frame(index=seq(n),time=rep(seq(t),n/t),id=rep(seq(n/t), each=t))
  ids <- unique(dat$id)
  
  # establish rolling window forecast for time-series cross-validation
  rolling_window_skeleton <- folds_rolling_window(t, window_size,
                                                  validation_size, gap, batch)
  
  folds_rolling_window <- lapply(rolling_window_skeleton, function(h){
    train_indices <- lapply(ids, function(i){
      train <- dat[dat$id == i, ]
      train[h$training_set, ]$index
    })
    val_indices <- lapply(ids, function(j){
      val <- dat[dat$id == j, ]
      val[h$validation_set, ]$index
    })
    make_fold(v=h$v, training_set=unlist(train_indices), 
              validation_set=unlist(val_indices))
  })
  return(folds_rolling_window)
}
