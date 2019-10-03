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

## Solution 1 - hypotensive episode for time t is defined as either:
# 1. abpmean at time t < 65 mmHg and there is 5-minute window around time t 
#    (i.e., 10 time-points) with at least 5 time-points abpmean < 65.
# 2. there is 5-minute window around time t (i.e., 10 time-points) where at 
#    least 8 time-points abpmean < 65.
new_Y_sol1 = function(train_all, window=5, cutoff=62){

  train_all$Y1<-rep(0,nrow(train_all))

  for(i in unique(train_all$subject_id)){

    sub<-train_all[train_all$subject_id==i,]

    init <- ifelse(sub[,"abpmean"] < cutoff,1,0)

    #Sum the lags around each possible event.
    lags<-rowSums(as.data.frame(lapply(seq(1:window), function(x) lag(init,x))), na.rm = TRUE)
    leads<-rowSums(as.data.frame(lapply(seq(1:window), function(x) lead(init,x))), na.rm = TRUE)
    mids<-lags+leads

    event<-ifelse((init==1 & mids >= 5),1,
                  ifelse((init==0 & mids >= 8), 1,0))
    
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

eval_missingness <- function(min, dataset, total_hrs = 5) {
  sec <- min*60
  total_min <- total_hrs*60
  dataset <- dataset %>%
    dplyr::group_by(subject_id) %>%
    mutate(init_time_and_date = min(time_and_date)) %>%
    mutate(min_elapsed = as.integer((time_and_date - init_time_and_date) / 60) + 1)
  d <- dataset %>%
    dplyr::group_by(subject_id) %>%
    dplyr::filter(min_elapsed <= total_min)
  df_full <- d %>%
    dplyr::group_by(subject_id) %>%
    dplyr::filter(min_elapsed == total_min)
  d <- d[(d$subject_id %in% df_full$subject_id),]
  dd <- d[order(d$subject_id, d$time_and_date),]
  dd$tdiff <- unlist(tapply(dd$time_and_date, INDEX = dd$subject_id,
                            FUN = function(x) c(0, diff(as.numeric(x)))))
  df_bad <- dd %>% dplyr::filter(tdiff > sec)
  d_complete <- d[!(d$subject_id %in% df_bad$subject_id),]
  list(num = length(unique(d_complete$subject_id)), dat = d_complete)
}


# fxn to summarize combine_SL fxn results
calculations <- function(res){
  
  # get all predictions in a data-frame:
  preds_fin <- as.data.frame(res$preds_fin)
  preds_regSL <- as.data.frame(res$pred_regular_SL)
  truth <- as.data.frame(res$truth)
  names(preds_regSL) <- names(preds_fin)
  names(truth) <- names(preds_fin)
  #preds_indSL <- as.data.frame(res$pred_regular_individual_SL)
  #names(preds_indSL) <- names(preds_fin)
  
  # AUC over all validation time-points
  pred_fin <- data.frame(pred = unlist(preds_fin))
  pred_regSL <- data.frame(pred = unlist(preds_regSL))
  truth <- data.frame(truth = unlist(truth))
  #pred_indSL <- data.frame(pred = unlist(preds_indSL))
  
  # averaged over time and all samples
  avg_loss <- colMeans(res$losses_all)
  names(avg_loss) <- names(res$losses_all)
  
  # look at the weights as an average over all samples:
  fit_coef <- lapply(res$fit_coef, function(x){
    rownames(x) <- x$learners
    return(x[,2:3])
  })
  weights <- data.frame(bind_cols(fit_coef))
  toDelete_convexcoefs <- seq(2, ncol(weights), 2)
  weights <- weights[,-toDelete_convexcoefs]
  weight <- data.frame(rowMeans(weights, na.rm = TRUE))
  names(weight) <- "Coefficient"
  row.names(weight) <- res$fit_coef[[1]]$learners
  weight <- weight/sum(weight)
  #weight_lrn<-cbind.data.frame(Learner=res$fit_coef[[1]]$learners,Coefficients=weight)
  #weight<-weight_lrn[order(weight, decreasing = TRUE),1:2]  
  
  # make a more general category:
  #1. picks the best algorithm from each category? 
  #2. averages over the algorithms?
  weight$SL_type <- NA
  weight[grepl("GlobalSL_screenbaseline", row.names(weight), fixed = TRUE),"SL_type"] <- "GlobalSL_Screen_Baseline"
  weight[grepl("GlobalSL_baseline", row.names(weight), fixed = TRUE),"SL_type"] <- "GlobalSL_Baseline"
  weight[grepl("GlobalSL_screen_", row.names(weight), fixed = TRUE),"SL_type"] <- "GlobalSL_Screen"
  weight[grepl("IndividualSL_", row.names(weight), fixed = TRUE),"SL_type"] <- "IndividualSL"
  weight[is.na(weight$SL_type),"SL_type"] <- "GlobalSL"
  
  #1. Best algorithm from each category
  max_SL_type <- weight %>% 
                 dplyr::group_by(SL_type) %>% 
                 dplyr::summarise(max=max(Coefficient))
  #2. Average of algorithms for each category
  ave_SL_type <- weight %>% 
                 dplyr::group_by(SL_type) %>% 
                 dplyr::summarise(ave=mean(Coefficient))
  
  return(list(pred_fin = pred_fin, pred_regSL = pred_regSL, truth = truth,
              loss = loss, weight = weight, max_SL_type = max_SL_type,
              ave_SL_type = ave_SL_type))
}

