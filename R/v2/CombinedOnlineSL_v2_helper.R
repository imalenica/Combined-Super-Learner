########### evaluate the loss:
# binary outcome: Binary Cross-Entropy / Log Loss
eval_loss <- function(ps, y){
  
  #Per row operation, so across all validation time-points.
  loss<-apply(ps, 1, function(p){
    #Avoids predictions out of range
    if(sum(as.numeric(p>1))>0){p[p>1]<-1-abs(1-p[p>1])}
    if(sum(as.numeric(p<0))>0){p[p<0]<-abs(p[p<0])}
    #Avoids inf loss
    p[(p==0)]<-0.000001; p[(p==1)]<-0.9999999
    #This is now sum over all h time-points
    estloss<-sum(-(y*log(p) + (1-y)*log(1-p)), na.rm=TRUE)
    if(is.na(estloss) & (p==1 || p==0)){estloss=0}
    estloss})
  
  return(loss=loss)
}

# continuous outcome: MSE
eval_loss_cont <- function(ps, y){
  
  #Per row operation, so across all validation time-points.
  loss<-apply(ps, 1, function(p){
    #This is now sum over all h time-points
    estloss<-sum((p-y)^2, na.rm = TRUE)
    estloss})
  
  return(loss=loss)
}

########### get the weights:
get_weights = function(ps,y,l=NULL,convex,discrete,interact=NULL){
  
  #Note: both interact and ps should be continuous, so should be ok
  #      to just multiply
  if(!is.null(interact)){
    ps <- ps * interact
  }
  
  if(discrete == TRUE & !is.null(l)){
    fit_coef<-as.numeric(l)
    fit_coef[which.min(l)]<-1
    fit_coef[-(which.min(l))]<-0
  }else if(discrete == FALSE & convex == TRUE){
    fit_coef <- lsei::pnnls(t(as.matrix(ps)), as.matrix(y), sum = 1)
    fit_coef <- fit_coef$x
  }else if(discrete == FALSE & convex == FALSE){
    #Gives a coefficient based on ALL the predictions (NOT time specific!)
    fit_coef <- nnls::nnls(t(as.matrix(ps)), as.matrix(y))
    fit_coef <- fit_coef$x
  }
  if(sum(fit_coef)==0){
    warning("All algorithms have zero weight", call. = FALSE)
  }
  
  return(fit_coef)
  
  #PROBLEM: If no weights given, distribute weights based on the loss?
  #if(sum(fit_coef)==0){
  #This "rewards" learnes with the smallest loss, but barely.
  #fit_coef<-as.numeric(l)
  #fit_coef<-sum(fit_coef)-fit_coef
  #fit_coef<-fit_coef/sum(fit_coef)
  
  #Give all weight to the discrete learner:
  #fit_coef[which.min(l)]<-1
  #}
  
  #TO DO: But why is it necessary for the coefficients to sum to 1?
  #if(sum(fit_coef)>0){
  #  fit_coef<-fit_coef/sum(fit_coef)
  #}else if(sum(fit_coef)==0){
  #  warning("All algorithms have zero weight", call. = FALSE)
  #Give all weight to the discrete learner:
  #fit_coef[which.min(l)]<-1
  ##NOTE: SuperLearner does not do this!  
  #Might cause issue with regular SL comparison?
  #}else if(sum(fit_coef)<1 && sum(fit_coef)!=0){
  #Give all weight to the discrete learner:
  #fit_coef<-rep(0,length(l))
  #fit_coef[which.min(l)]<-1 
  #  fit_coef<-fit_coef
  #}
}

#Function to return weights based on adapt_covs (baseline covariates):
get_weights_adaptive = function(train_all, adapt_covs, preds, truths, truths_training, outcome){
  
  #Create dataset with requested covariates
  sample_baseline_cov <- train_all %>% 
    group_by(subject_id) %>% 
    slice(1) %>%
    select(subject_id,adapt_covs)
  sample_baseline_cov <- sample_baseline_cov[,-1]
  
  #Create lists with repeated baseline covariates to match the number of time-points
  res<-split(sample_baseline_cov, f=row.names(sample_baseline_cov))
  res<-lapply(res, function(x){ splitstackshape::expandRows(x, count=nrow(truths_training[[1]]), 
                                                            count.is.col=FALSE) })
  res_val<-lapply(res, function(x){ splitstackshape::expandRows(x, count=nrow(truths[[1]]), 
                                                                count.is.col=FALSE) })
  
  #Create a dataset for glm (trained on training samples):
  dt <- cbind.data.frame(y=unlist(truths_training),rlist::list.rbind(res))
  
  if(length(unique(train_all[,outcome]))>2){
    #Continuous outcome
    fit <- glm(formula = y ~ ., data=dt)
  }else{
    fit <- glm(formula = y ~ ., data=dt, family=binomial())
  }
  
  #Get predictions for validation samples:
  weight_cov <- predict(fit,newdata = rlist::list.rbind(res), type="response")
  weight_residual<-unlist(truths)-weight_cov
  
  #Global weights (on validation samples):
  fit_coef <- get_weights(ps=rlist::list.cbind(preds), y=unlist(truths), convex=convex, 
                          discrete=discrete, interact=weight_cov)
  
  return(fit_coef)
  
  ###### Fit with the Global SL 
  # Results in different models, not optimal.
  
  #Initialize learners:
  #glm_fast <- make_learner(Lrnr_glm_fast)
  #stack <- make_learner(Stack, glm_fast)
  #cv_stack_weights <- Lrnr_cv$new(stack)
  
  #Grab predictions on the validation samples: (105 X 680)
  #global_SL_t_weight<-global_SL(train_all=train_all, t=t, outcome=outcome, 
  #                              stack_pool=cv_stack_weights, stack_screen=NULL,
  #                              covars=NULL, covars_wbaseline=adapt_covs,
  #                              test_size=test_size, mini_batch=mini_batch,
  #                              cv=cv, first_window=first_window, 
  #                              window_size=window_size, gap = gap_training)
  #weight_cov<-rlist::list.rbind(global_SL_t_weight$preds_globalSL_baseline)
  #weight_residual<-unlist(global_SL_t_weight$truths)-weight_cov
  
  #Global weights:
  #get_weights(ps=rlist::list.cbind(preds), y=unlist(truths), convex=convex, 
  #            discrete=discrete, interact=weight_cov)

}

########### get the ensemble:
#Fuction to make more streamlined (continuous) Super Learners
get_sl = function(predictions, weight, samples){
  if(is.list(weight)){
    pred_fin <- lapply(seq_len(length(predictions)), function(i){
      pred <- data.frame(pred = as.matrix(t(predictions[[i]])) %*% 
                           weight[[i]])
      names(pred) <- paste0("Sample",samples[i])
      pred
    })
  }else if(!is.list(weight)){
    predictions<-rlist::list.cbind(predictions)
    pred <- data.frame(pred = as.matrix(t(predictions)) %*% weight)
    
    id <- rep(as.numeric(levels(samples))[samples],each=(nrow(pred)/length(samples)))
    
    pred$id<-id
    pred_fin<-split(pred, f = pred$id)
    pred_fin<-lapply(pred_fin, function(x) x[,1])
  }else{
    #Corresponds to the case when baseline covs not specified
    pred_fin< NULL
  }
  return(pred_fin)
}

#Fuction to make more streamlined (continuous) Super Learners loss calculation
get_sl_loss = function(predictions, truths, binary=TRUE){
  if(!is.null(predictions)){
    if(binary){
      sample_loss <- unlist(lapply(seq_len(length(predictions)), function(i){
        eval_loss(t(predictions[[i]]),truths[[i]])}))
    }else{
      sample_loss <- unlist(lapply(seq_len(length(predictions)), function(i){
        eval_loss_cont(t(predictions[[i]]),truths[[i]])}))
    }
    loss <- mean(sample_loss)
  }else{
    loss=NA
    sample_loss=NA
  }
  return(list(loss=loss, sample_loss=sample_loss))
}

########### Helper prediction and truth functions
predict_fun = function(SL,task,list=TRUE,
                       cv="folds_rolling_origin",ids=NULL,test_size=1){
  #Same as SL$predict
  res=SL$predict_fold(task, "validation")
  
  #Split res into sample-specific predictions
  if(list){
    if(cv=="folds_rolling_origin" | cv=="folds_rolling_window"){
      #Switch back to sample ids, but can specify id=seq(1:n) if prefer
      #(subject id X test times) X num of folds
      id = rep(ids,each=length(SL$training_task$folds)*test_size)
    }else if(cv=="folds_vfold"){
      #More problematic, does not correspond with samples any more
      time <- nrow(task$data)/length(ids)
      id <- rep(ids,each=time)
    }
    #Arrange to match truths
    res$id<-id
    res<-split(res, f = res$id)
    res<-lapply(res, function(x) x[,-"id"])
  }
  return(res)
}

truth_fun = function(cv, folds, test_data, outcome="hypo_event", validation_truth=TRUE, 
                     ids=NULL, test_size=1){
  if(validation_truth){
    if(cv=="folds_rolling_origin" | cv=="folds_rolling_window"){
      id = rep(rep(ids,each=test_size),length(folds))
      truth<-bind_rows(lapply(folds, function(i) test_data[i$validation_set,]))
    }else if(cv=="folds_vfold"){
      #CAUTION: Not ordered correctly, however! 
      truth<-data.frame(Y=test_data[,outcome])
      id = rep(ids, each=nrow(truth)/length(unique(test_data$subject_id)))
    }
  }else{
    if(cv=="folds_rolling_origin" | cv=="folds_rolling_window"){
      #TO DO: Make general, don't use "subject_id"
      truth<-bind_rows(lapply(folds, function(i) test_data[i$training_set,]))
      id = as.numeric(levels(truth$subject_id))[truth$subject_id]
    }else if(cv=="folds_vfold"){
      #CAUTION: Not ordered correctly, however! 
      truth<-data.frame(Y=test_data[,outcome])
      id = rep(ids, each=nrow(truth)/length(unique(test_data$subject_id)))
    }
  }
  
  #Split res into sample-specific truths
  truth$id<-id
  truths<- split(truth, f = truth$id)
  truths<-lapply(truths, function(x) x[,outcome])
  return(truths)
}

#####################
### Global learner
#####################

#Good idea to have mini_batch=test_size...
global_SL = function(train_all, t, outcome, sl, stack_pool, stack_screen=NULL, 
                     covars, covars_wbaseline, test_size, mini_batch=1, V=5,
                     cv, first_window=1, window_size=1, gap=1){
  
  #### Pool across time (fitting on all data up to time t)
  
  #Here, we impose a Markov order assumption
  test_data <- train_all[train_all$time<=t,]
  ids=seq(1:length(unique(test_data$subject_id)))
  
  #Check if we have data for all samples 
  #(requirment of some time-series CVs)
  if(nrow(test_data)!=length(ids)*t){
    #Set up only V fold CV
    folds <- origami::make_folds(test_data,
                                 fold_fun = folds_vfold,
                                 V=V)
  }else{
    #Set up proper cross-validation for multiple time-series data:
    if(cv=="folds_rolling_origin"){
      folds <- origami::make_folds(test_data,
                                   fold_fun = folds_rolling_origin_pooled,
                                   t=t,
                                   first_window = first_window,
                                   validation_size = test_size, gap = gap,
                                   batch = mini_batch)
    }else if(cv=="folds_rolling_window"){
      folds <- origami::make_folds(test_data,
                                   fold_fun = folds_rolling_window_pooled,
                                   t=t,
                                   window_size = window_size,
                                   validation_size = test_size, gap = gap,
                                   batch = mini_batch)
    }else if(cv=="folds_vfold"){
      folds <- origami::make_folds(test_data,
                                   fold_fun = folds_vfold,
                                   V=V)
    }
  }
  
  # create the sl3 task with time-varying covariates:
  if(!is.null(covars)){
    task <- make_sl3_Task(
      data = test_data, covariates = covars,
      outcome = outcome, folds=folds, drop_missing_outcome = T)
    
    #Fit the global learner:
    globalSL<-stack_pool$train(task)
    globalSL_preds<-predict_fun(SL=globalSL, task=task, cv=cv, ids=ids, test_size=test_size) #per individual
  }else{
    task<-NULL
    globalSL<-NULL
    globalSL_preds<-NULL
  }
  
  # create the sl3 task with time-varying and baseline covariates:
  if(!is.null(covars_wbaseline)){
    task_baseline <- make_sl3_Task(
      data = test_data, covariates = covars_wbaseline,
      outcome = outcome, folds=folds, drop_missing_outcome = T)
    
    #Fit the global learner:
    globalSL_baseline<-stack_pool$train(task_baseline)
    globalSL_baseline_preds<-predict_fun(SL=globalSL_baseline, task=task_baseline, 
                                         cv=cv, ids=ids, test_size=test_size)
  }else{
    task_baseline<-NULL
    globalSL_baseline<-NULL
    globalSL_baseline_preds<-NULL
  }
  
  #Fit the regular Super Learner (with baseline covariates)
  #regularSL <- sl$train(task_baseline)
  regularSL <- NULL
  
  if(!is.null(stack_screen)){
    
    #Fit the global learner with screeners:
    if(!is.null(covars)){
      
      #Fit the global learner:
      globalSL_screen<-stack_screen$train(task)
      #Predict on validation samples
      globalSL_screen_preds<-predict_fun(SL=globalSL_screen, task=task, cv=cv, ids=ids, test_size=test_size)
    }else{
      globalSL_screen<-NULL
      globalSL_screen_preds<-NULL
    }
    
    # create the sl3 task with time-varying and baseline covariates
    if(!is.null(covars_wbaseline)){
     
      #Fit the global learner:
      globalSL_screen_baseline<-stack_screen$train(task_baseline)
      globalSL_screen_baseline_preds<-predict_fun(SL=globalSL_screen_baseline, task=task_baseline, 
                                                  cv=cv, ids=ids, test_size=test_size)
    }else{
      globalSL_screen_baseline<-NULL
      lobalSL_screen_baseline_preds<-NULL
    }
  
  }else{
    globalSL_screen<-NULL
    globalSL_screen_preds<-NULL
    globalSL_screen_baseline<-NULL
    globalSL_screen_baseline_preds<-NULL
  }
  
  #Get truths (validation):
  truths<-truth_fun(cv=cv, folds=folds, test_data=test_data, ids=ids, test_size=test_size, outcome=outcome)
  truths_training<-truth_fun(cv=cv, folds=folds, test_data=test_data, ids=ids, test_size=test_size, 
                             validation_truth=FALSE, outcome=outcome)
  return(list(t=t,
              #Save SL fits 
              regularSL=regularSL,
              globalSL=globalSL,
              globalSL_baseline=globalSL_baseline,
              globalSL_screen=globalSL_screen,
              globalSL_screen_baseline=globalSL_screen_baseline,
              #Save SL predictions (vfor validation samples, indexed by sample)
              preds_globalSL=globalSL_preds,
              preds_globalSL_baseline=globalSL_baseline_preds,
              preds_globalSL_screen=globalSL_screen_preds,
              preds_globalSL_screen_baseline=globalSL_screen_baseline_preds,
              #Save truths (for validation samples, indexed by sample)
              truths=truths,
              truths_training=truths_training,
              #Save tasks and folds (one for all)
              task=task, task_baseline=task_baseline, folds=folds
  ))
}

#######################
### Individual learner
#######################

individual_SL = function(train_all, t, id, cv, first_window, window_size,
                         test_size, mini_batch, covars, stack_individual, gap){
  
  #Subset to sample with subject id id
  train_one <- train_all[train_all$subject_id %in% id,]
  
  #Use data until time point until time t:
  test_data <- train_one[train_one$time<=t,]
  
  #Set up proper cross-validation for single time-series data:
  if(cv=="folds_rolling_origin"){
    folds <- origami::make_folds(test_data,
                                 fold_fun = folds_rolling_origin,
                                 first_window = first_window,
                                 validation_size = test_size, gap = gap,
                                 batch = mini_batch)
  }else if(cv=="folds_rolling_window"){
    folds <- origami::make_folds(test_data,
                                 fold_fun = folds_rolling_window,
                                 window_size = window_size,
                                 validation_size = test_size, gap = gap,
                                 batch = mini_batch)
  }
  
  truths <- truth_fun(cv=cv, folds=folds, test_data, ids=as.numeric(id), test_size=test_size)
  #truths<-lapply(folds, function(i) data.frame(test_data[i$validation, get(outcome)]))
  truths<-bind_rows(truths)[,1]
  
  # create the sl3 task
  task <- make_sl3_Task(
    data = test_data, covariates = covars, outcome = outcome, folds = folds)
  
  #Fit the regular Super Learner (with baseline covariates)
  #indregularSL <- sl$train(task) #Issue with folds for some reason
  indregularSL<-NULL
  
  #Fit the individual learner:
  individualSL <- stack_individual$train(task)
  preds_individualSL<-individualSL$predict_fold(task, "validation")
  
  return(list(t=t, 
              #Save SL fits 
              indregularSL=indregularSL,
              individualSL=individualSL, 
              #Save SL predictions (vfor validation samples, indexed by sample)
              preds_individualSL=preds_individualSL,
              #Save truths (for validation samples, indexed by sample)
              truths=truths,
              #Save tasks and folds (one for each sample)
              task=task, folds=folds
  ))
}