'%nin%' <- Negate('%in%')

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
  names(avg_loss) <- c("loss_combined_SL", 
                       "loss_regular_SL", 
                       "loss_individual_SL")
  
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
  avg_weight <- weight/sum(weight)
  #weight_lrn<-cbind.data.frame(Learner=res$fit_coef[[1]]$learners,Coefficients=weight)
  #weight<-weight_lrn[order(weight, decreasing = TRUE),1:2]  
  
  # make a more general category:
  #1. picks the best algorithm from each category? 
  #2. averages over the algorithms?
  avg_weight$SL_type <- NA
  avg_weight[grepl("GlobalSL_screenbaseline", row.names(avg_weight), fixed = TRUE),"SL_type"] <- "GlobalSL_Screen_Baseline"
  avg_weight[grepl("GlobalSL_baseline", row.names(avg_weight), fixed = TRUE),"SL_type"] <- "GlobalSL_Baseline"
  avg_weight[grepl("GlobalSL_screen_", row.names(avg_weight), fixed = TRUE),"SL_type"] <- "GlobalSL_Screen"
  avg_weight[grepl("IndividualSL_", row.names(avg_weight), fixed = TRUE),"SL_type"] <- "IndividualSL"
  avg_weight[is.na(avg_weight$SL_type),"SL_type"] <- "GlobalSL"
  
  #1. Best algorithm from each category
  max_SL_type <- avg_weight %>% 
                 dplyr::group_by(SL_type) %>% 
                 dplyr::summarise(max=max(Coefficient))
  #2. Average of algorithms for each category
  ave_SL_type <- avg_weight %>% 
                 dplyr::group_by(SL_type) %>% 
                 dplyr::summarise(ave=mean(Coefficient))
  
  comSL_summary <- list(pred_fin = pred_fin, pred_regSL = pred_regSL, 
                        truth = truth, avg_loss = avg_loss, 
                        avg_weight = avg_weight, max_SL_type = max_SL_type,
                        ave_SL_type = ave_SL_type)
  return(comSL_summary)
}

calculate_ind_AUCs <- function(res_all){
  
  ind_res <- list()
  ind_AUCs <- list()
  
  for(i in 1:length(res_all[[1]]$preds_fin)){
    truths <- lapply(res_all, function(res){
      factor(res$truth[[i]]$hypo_event, levels = c(0,1))})
    preds_com <- lapply(res_all, function(res) res$preds_fin[[i]][,1])
    preds_reg <- lapply(res_all, function(res) res$pred_regular_SL[[i]][,1])
    preds_ind <- lapply(res_all, function(res) res$pred_individual_SL[[i]][,1])
    id <- substring(colnames(res_all[[1]]$preds_fin[[i]]), 7)
    
    df <- list()
    for(j in 1:length(truths)){
      t <- rep(as.numeric(names(truths)[[j]]), length(truths[[j]]))
      df[[j]] <- data.frame(training_time = t, truth = truths[[j]],
                            pred_com = preds_com[[j]], 
                            pred_reg = preds_reg[[j]], 
                            pred_ind = preds_ind[[j]])
    }
    df <- bind_rows(df)
    subject_id <- rep(id, nrow(df))
    ind_res[[i]] <- data.table(subject_id, df)

    if(length(unique(ind_res[[i]]$truth)) == 1){
      AUC_com <- AUC_reg <- AUC_ind <- NA
      ind_AUCs[[i]] <- data.frame(subject_id = id, AUC_com, AUC_reg, AUC_ind)
    }
    if(length(unique(ind_res[[i]]$truth)) > 1){
      roc_com <- roc(ind_res[[i]]$truth, ind_res[[i]]$pred_com, ci = TRUE)
      roc_reg <- roc(ind_res[[i]]$truth, ind_res[[i]]$pred_reg, ci = TRUE)
      roc_ind <- roc(ind_res[[i]]$truth, ind_res[[i]]$pred_ind, ci = TRUE)
      AUC_com <- as.numeric(ci.auc(auc(roc_com)))[2]
      AUC_reg <- as.numeric(ci.auc(auc(roc_reg)))[2]
      AUC_ind <- as.numeric(ci.auc(auc(roc_ind)))[2]
      ind_AUCs[[i]] <- data.frame(subject_id = id, AUC_com, AUC_reg, AUC_ind)
    }
  }
  ind_AUCs_df <- bind_rows(ind_AUCs) %>%
    arrange(AUC_com)
  return(list(ind_AUCs_df = ind_AUCs_df,
              ind_res = ind_res))
}

plot_coefvtime <- function(comSL_summary_list, weight_grouping = c("max", "ave"), 
                           cv_type){
  if(weight_grouping == "max"){
    avg_weight_all <- t(cbind.data.frame(lapply(comSL_summary_list, function(z){
      z$max_SL_type[2]})))
    SL_types <- t(comSL_summary_list[[1]]$max_SL_type[1])
  }
  if(weight_grouping == "ave"){
    avg_weight_all <- t(cbind.data.frame(lapply(comSL_summary_list, function(z){
      z$ave_SL_type[2]})))
    SL_types <- t(comSL_summary_list[[1]]$ave_SL_type[1])
  }
  
  row.names(avg_weight_all) <- names(comSL_summary_list)
  colnames(avg_weight_all) <- t(comSL_summary_list[[1]]$max_SL_type[1])
  avg_weight_all <- cbind.data.frame(Time = row.names(avg_weight_all),
                                     avg_weight_all)
  avg_weight_all <- melt(avg_weight_all, id = "Time")
  avg_weight_all$Time <- as.numeric(levels(avg_weight_all$Time))[avg_weight_all$Time]
  obj <- ggplot(avg_weight_all, aes(x = Time, y = value, colour = variable)) +
    geom_line(aes(group = variable), size = 0.4) + 
    geom_point(shape = 1) + 
    ggtitle(paste0("SL Weights over Time for ", cv_type, " CV")) +
    labs(x = "Training Time (Minutes)", 
         y = "Average SL Coefficient",
         color = "SL Type") 
  return(obj)
}
  
calculate_AUCs <- function(comSL_summary_list, cv_type){
  
  truths <- lapply(comSL_summary_list, function(z) z$truth$truth)
  preds <- lapply(comSL_summary_list, function(z) z$pred_fin$pred)
  
  AUCs <- list()
  for(i in 1:length(comSL_summary_list)){
    roc_obj <- roc(truths[[i]], preds[[i]], ci = TRUE)
    AUCs[[i]] <- as.numeric(ci.auc(auc(roc_obj)))[2]
  }
  names(AUCs) <- names(comSL_summary_list)
  AUCs_df <- data.frame(unlist(AUCs))
  colnames(AUCs_df) <- "AUC"
  return(AUCs_df)
}


