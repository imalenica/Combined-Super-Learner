##Evaluate the loss:
#Binary outcome: Binary Cross-Entropy / Log Loss
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

#Continuous outcome: MSE
eval_loss_cont <- function(ps, y){
  
  #Per row operation, so across all validation time-points.
  loss<-apply(ps, 1, function(p){
    #This is now sum over all h time-points
    estloss<-sum((p-y)^2, na.rm = TRUE)
    estloss})
  
  return(loss=loss)
}
#####################
### Global learner
#####################

global_SL = function(train_all, t, outcome, 
                     sl, stack_pool, stack_screen=NULL, 
                     covars, covars_wbaseline,
                     test_size=5, mini_batch=5, V=NULL,
                     cv="folds_rolling_origin", first_window=1, window_size=1){
  
  #Pool across time (fitting on all data up to time t)
  #Here, we impose a Markov order assumption
  test_data <- train_all[train_all$time<=t,]
  
  #Set up proper cross-validation for multiple time-series data:
  if(cv=="folds_rolling_origin"){
    folds <- origami::make_folds(test_data,
                                 fold_fun = folds_rolling_origin_pooled,
                                 t=t,
                                 first_window = first_window,
                                 validation_size = test_size, gap = 0,
                                 batch = mini_batch)
  }else if(cv=="folds_rolling_window"){
    folds <- origami::make_folds(test_data,
                                 fold_fun = folds_rolling_window_pooled,
                                 t=t,
                                 window_size = window_size,
                                 validation_size = test_size, gap = 0,
                                 batch = mini_batch)
  }else if(cv=="folds_vfold"){
    folds <- origami::make_folds(test_data,
                                 fold_fun = folds_vfold,
                                 V=V)
  }
  
  # create the sl3 task with time-varying covariates
  task <- make_sl3_Task(
    data = test_data, covariates = covars,
    outcome = outcome, folds=folds)
  
  # create the sl3 task:
  task_baseline <- make_sl3_Task(
    data = test_data, covariates = covars_wbaseline,
    outcome = outcome, folds=folds)
  
  #### PROBLEM: Regular sl has an issue with these folds
  #Fit the regular Super Learner (with baseline covariates)
  #regularSL <- sl$train(task_baseline)
  regularSL <- NULL
  
  #Fit the global learner:
  globalSL<-stack_pool$train(task)
  globalSL_baseline<-stack_pool$train(task_baseline)
  
  if(!is.null(stack_screen)){
    #Fit the global learner with screeners:
    globalSL_screen<-stack_screen$train(task)
    globalSL_screen_baseline<-stack_screen$train(task_baseline)
  }else{
    globalSL_screen<-NULL
    globalSL_screen_baseline<-NULL
  }
  
  return(list(t=t,
              regularSL=regularSL,
              globalSL=globalSL,
              globalSL_baseline=globalSL_baseline,
              globalSL_screen=globalSL_screen,
              globalSL_screen_baseline=globalSL_screen_baseline))
}

#######################
### Individual learner
#######################

individual_SL = function(train_all,t,id, cv="folds_rolling_origin",
                         first_window=1, window_size=5,
                         test_size,mini_batch,
                         covars_wbaseline,stack_individual){
  
  #Subset to sample with subject id id
  train_one <- train_all[train_all$subject_id %in% id,]
  
  #Use data until time point until time t:
  test_data <- train_one[train_one$time<=t,]
  
  #Set up proper cross-validation for single time-series data:
  if(cv=="folds_rolling_origin"){
    folds <- origami::make_folds(test_data,
                                 fold_fun = folds_rolling_origin,
                                 first_window = first_window,
                                 validation_size = test_size, gap = 0,
                                 batch = mini_batch)
  }else if(cv=="folds_rolling_window"){
    folds <- origami::make_folds(test_data,
                                 fold_fun = folds_rolling_window,
                                 window_size = window_size,
                                 validation_size = test_size, gap = 0,
                                 batch = mini_batch)
  }

  # create the sl3 task
  task <- make_sl3_Task(
    data = test_data, covariates = covars_wbaseline,
    outcome = outcome, folds=folds)
  
  #Fit the regular Super Learner (with baseline covariates)
  #indregularSL <- sl$train(task) #Issue with folds for some reason
  indregularSL<-NULL
  
  #Fit the global learner:
  individualSL<-stack_individual$train(task)
  
  return(list(t=t, 
              indregularSL=indregularSL,
              individualSL=individualSL))
}

#######################
### Combined learner
#######################

#The idea behind the Combined Online SL is to capitalize on both the Global and Individual Learner.
#It will train the Global SL with and without baseline covariates, and with and without screeners.
#This effectively mimics various level of smoothing over baseline covariates.
#In addition, the Combined Online SL also uses the Individual SL, which learns only from the sample we try to optimize over.
#We look at the performance over all individual samples- meaning that we have n different weights generated.
#Note that this means that our Global SL will use the same function over all n samples, 
#whereas the Individual SL will use different functions over all n samples.
#Loss generated is a sum over all the samples, evaluated only over validation time points.

### For each t, derive weights at horizon h
#train_all: the full dataset used
#outcome: name of the outcome column.
#t: time until which we train the SuperLearners
#stack_pool: sl3 stack object corresponding to the learners used for the Global Super Learner. The learners in it should not use screeners.
#stack_individual: sl3 stack object corresponding to the Individual Super Learner. The learners in it should not use screeners.
#stack_screen: sl3 stack object. The learners in it should use screeners.
#sl: sl3 SuperLearner object.
#covars: Time varying covariates
#covars_baseline: Baseline covariates
#id: Estimate Combined SuperLearner performance over a single subject. Final loss and weights will be id specific.
#cv: Time-series cross-validation used. Options are folds_rolling_origin and folds_rolling_window. 
#gap: time between the last trained time point and the first prediction time point
#h: length of horizon at which we evaluate the loss
#test_size: size of the validation set used in the training procedure (in time points).
#mini_batch: increase in the training set size from the first iteration in time points. 
#first_window: size of the training set (for example, 1 to first window). This parameter is specific to the 
#folds_rolling_origin time-series cross-validation.
#window_size: size of the training set (always of size window_size). This parameter is specific to the 
#folds_rolling_window time-series cross-validation.

combine_SL2 = function(train_all, outcome, t, 
                      stack_pool, stack_individual, stack_screen=NULL, sl,
                      covars, covars_baseline,
                      id=NULL, cv="folds_rolling_origin",
                      gap=0, h=1, test_size=5, mini_batch=5,
                      first_window=1, window_size=5){
  
  #Combine time-varying and baseline covariates
  covars_wbaseline <- c(covars, covars_baseline)
  
  #Allows for a single id combined SL
  if(is.null(id)){
    #Get all the individual samples:
    samples <- unique(train_all$subject_id)
  }else{
    #If id provided, then just the id:
    samples <- id
  }
  
  #Train the Global SLs:
  global_SL_t<-global_SL(train_all=train_all, t=t, outcome=outcome, 
                         sl=sl, stack_pool=stack_pool, stack_screen=stack_screen,
                         covars=covars, covars_wbaseline=covars_wbaseline,
                         test_size=test_size, mini_batch=mini_batch,
                         cv=cv, first_window=first_window, window_size=window_size)
  global_SL_t_reg<-global_SL(train_all=train_all, t=t, outcome=outcome, 
                         sl=sl, stack_pool=stack_pool, stack_screen=stack_screen,
                         covars=covars, covars_wbaseline=covars_wbaseline,
                         test_size=test_size, mini_batch=mini_batch,
                         cv="folds_vfold", V=10, first_window=first_window, window_size=window_size)
  
  #Create individual SLs for all samples:
  individual_SL_t <- lapply(samples, function(x){individual_SL(train_all=train_all,t=t,id=x, cv=cv,
                                                               first_window=first_window,
                                                               window_size=window_size,
                                                               test_size=test_size,
                                                               mini_batch=mini_batch,
                                                               covars_wbaseline=covars_wbaseline,
                                                               stack_individual=stack_individual)})
  
  # create the sl3 task:
  tasks <- lapply(samples, function(x){train_one <- train_all[train_all$subject_id %in% x,]
  make_sl3_Task(data = train_one[(t+1+gap):(t+gap+h),], covariates = covars, outcome = outcome)
  })
  
  # create the sl3 task with baseline covariates:
  tasks_baseline <- lapply(samples, function(x){train_one <- train_all[train_all$subject_id %in% x,]
  make_sl3_Task(data = train_one[(t+1+gap):(t+gap+h),], covariates = covars_wbaseline, outcome = outcome)
  })
  
  lrn<-length(global_SL_t$globalSL$.__enclos_env__$private$.learner_names)
  learner_names_global<-global_SL_t$globalSL$.__enclos_env__$private$.learner_names
  learner_names_individual<-individual_SL_t[[1]]$individualSL$.__enclos_env__$private$.learner_names

  ### Get all predictions:
  
  #### Regular SL: (using sl3 directly)
  #pred_regular_SL <- lapply(tasks_baseline, function(task){global_SL_t$regularSL$predict(task)})
  pred_global_SL_baseline_reg <- lapply(tasks_baseline, function(task){matrix(unlist(lapply(1:lrn, function(x){
    global_SL_t_reg$globalSL_baseline$fit_object$learner_fits[[x]]$predict(task)})),nrow=h,ncol=lrn)})
  
  #### Global SL:
  pred_global_SL <- lapply(tasks, function(task){matrix(unlist(lapply(1:lrn, function(x){
    global_SL_t$globalSL$fit_object$learner_fits[[x]]$predict(task)})),nrow=h,ncol=lrn)})
  #Global SL with baseline covariates:
  pred_global_SL_baseline <- lapply(tasks_baseline, function(task){matrix(unlist(lapply(1:lrn, function(x){
    global_SL_t$globalSL_baseline$fit_object$learner_fits[[x]]$predict(task)})),nrow=h,ncol=lrn)})
  #Global SL with screeners:
  if(!is.null(stack_screen)){
    pred_global_SL_screen <- lapply(tasks, function(task){
      p<-global_SL_t$globalSL_screen$fit_object$learner_fits[[1]]$predict(task)
      as.matrix(p)
    })
    #Global SL with screeners and baseline covariates:
    pred_global_SL_screen_baseline <- lapply(tasks_baseline, function(task){
      p<-global_SL_t$globalSL_screen_baseline$fit_object$learner_fits[[1]]$predict(task)
      as.matrix(p)
    })
  }else{
    pred_global_SL_screen<-NULL
    pred_global_SL_screen_baseline<-NULL
  }
  
  #### Individual SL:
  pred_individual_SL<-lapply(seq_along(tasks), function(i){
    matrix(unlist(lapply(1:lrn, function(x){
      individual_SL_t[[i]]$individualSL$fit_object$learner_fits[[x]]$predict(tasks_baseline[[i]])})),nrow=h,ncol=lrn)
  })
  #Regular Individual SL: (using sl3 directly)
  #pred_regular_individual_SL<-lapply(seq_along(tasks), function(i){
  #  individual_SL_t[[i]]$indregularSL$predict(tasks_baseline[[i]])})
  
  if(!is.null(stack_screen)){
    learners<-c(paste("GlobalSL", learner_names_global, sep="_"),
                paste("GlobalSL_baseline", learner_names_global, sep="_"),
                paste("GlobalSL_screen", learner_names_global, sep="_"),
                paste("GlobalSL_screenbaseline_", learner_names_global, sep="_"),
                paste("IndividualSL", learner_names_individual, sep="_")) 
    
    #For sample i, get predictions:
    preds <- lapply(seq_along(tasks), function(i){
      ps<-cbind.data.frame(pred_global_SL[[i]],
                           pred_global_SL_baseline[[i]],
                           pred_global_SL_screen[[i]],
                           pred_global_SL_screen_baseline[[i]],
                           pred_individual_SL[[i]])
      names(ps)<-learners
      t(ps)})
  }else{
    learners<-c(paste("GlobalSL", learner_names_global, sep="_"),
                paste("GlobalSL_baseline", learner_names_global, sep="_"),
                paste("IndividualSL", learner_names_individual, sep="_")) 
    
    #For sample i, get predictions:
    preds <- lapply(seq_along(tasks), function(i){
      ps<-cbind.data.frame(pred_global_SL[[i]],
                           pred_global_SL_baseline[[i]],
                           pred_individual_SL[[i]])
      names(ps)<-learners
      t(ps)})
  }
  
  #Get the truth:
  truths <- lapply(samples, function(x){
    train_one <- train_all[train_all$subject_id %in% x,]
    train_one[(t+1+gap):(t+gap+h), outcome]
  })
  
  ### Evaluate the loss for discrete learners, for each sample:
  if(length(unique(train_all[,outcome]))>2){
    #Continuous outcome
    loss <- lapply(seq_along(tasks), function(i){eval_loss_cont(preds[[i]],truths[[i]])})
  }else{
    #Binary outcome
    loss <- lapply(seq_along(tasks), function(i){eval_loss(preds[[i]],truths[[i]])})
  }
  
  #Assign weights to discrete learners:
  get_weights = function(ps,y,l){
    
    #Gives a coefficient based on ALL the predictions (NOT time specific!)
    fit_coef <- nnls::nnls(t(as.matrix(ps)), as.matrix(y))
    fit_coef <- fit_coef$x
    
    if(sum(fit_coef)==0){
      warning("All algorithms have zero weight", call. = FALSE)
      
    }
    
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
    
    return(fit_coef)
  }
  
  #### Get weights for each sample
  #Get weights with combined SuperLearner
  weights <- lapply(seq_along(tasks), function(i){
    get_weights(preds[[i]],truths[[i]],loss[[i]])
  })
  #Get weights for just the regular SuperLearner
  weights_reg_sl <- lapply(seq_along(tasks), function(i){
    get_weights(t(pred_global_SL_baseline_reg[[i]]),truths[[i]], l=0)
  })
  #Get weights for just the individual SuperLearner
  weights_ind_sl <- lapply(seq_along(tasks), function(i){
    get_weights(t(pred_individual_SL[[i]]),truths[[i]], l=0)
  })
  
  #### Get the ensemble predictions
  #Get the prediction with weighted combination of learners:
  preds_fin <- lapply(seq_along(tasks), function(i){
    #Generates weighted prediction for each validation time-point
    pred <- data.frame(pred = as.matrix(t(preds[[i]])) %*% weights[[i]])
    names(pred) <- paste0("Sample",samples[i])
    pred
  })
  #Get the regular SL prediction
  pred_regular_SL <- lapply(seq_along(tasks), function(i){
    #Generates weighted prediction for each validation time-point (regular SL)
    pred <- data.frame(pred = as.matrix((pred_global_SL_baseline_reg[[i]])) %*% weights_reg_sl[[i]])
    names(pred) <- paste0("Sample",samples[i])
    pred
  })
  #Get the individual SL prediction
  pred_individual_SL <- lapply(seq_along(tasks), function(i){
    #Generates weighted prediction for each validation time-point (regular SL)
    pred <- data.frame(pred = as.matrix((pred_individual_SL[[i]])) %*% weights_ind_sl[[i]])
    names(pred) <- paste0("Sample",samples[i])
    pred
  })
  
  #Connect weights and learners:
  fit_coef <- lapply(seq_along(tasks), function(i){
    data.frame(learners=learners,coefs=weights[[i]],convex_coefs=(weights[[i]]/sum(weights[[i]])))
  })
  
  ### Evaluate loss over final SuperLearners!
  #Note that this is the loss over validation samples. Evaluate the loss for discrete learners, for each sample:
  if(length(unique(train_all[,outcome]))>2){
    #Continuous outcome
    loss_online_SL<-unlist(lapply(seq_along(tasks), function(i){eval_loss_cont(t(preds_fin[[i]]),truths[[i]])}))
    loss_regular_SL<-unlist(lapply(seq_along(tasks), function(i){eval_loss_cont(t(pred_regular_SL[[i]]),truths[[i]])}))
    loss_individual_SL<-unlist(lapply(seq_along(tasks), function(i){eval_loss_cont(t(pred_individual_SL[[i]]),truths[[i]])}))
  }else{
    #Binary outcome
    loss_online_SL<-unlist(lapply(seq_along(tasks), function(i){eval_loss(t(preds_fin[[i]]),truths[[i]])}))
    loss_regular_SL<-unlist(lapply(seq_along(tasks), function(i){eval_loss(t(pred_regular_SL[[i]]),truths[[i]])}))
    loss_individual_SL<-unlist(lapply(seq_along(tasks), function(i){eval_loss(t(pred_individual_SL[[i]]),truths[[i]])}))
  }
  
  losses_all<-cbind.data.frame(loss_online_SL=loss_online_SL,
                               loss_regular_SL=loss_regular_SL,
                               loss_individual_SL=loss_individual_SL)

  return(list(#Final, weighted prediction.
    preds_fin=preds_fin, 
    #Final, weighted prediction for the regular SL
    pred_regular_SL=pred_regular_SL,
    #Final, weighted prediction for the one sample SL
    pred_individual_SL=pred_individual_SL,
    #t+h truth for all the samples.
    truth=truths, 
    #Predictions for each individual learner.
    preds=preds,  
    #Coefficient for each learner
    fit_coef=fit_coef, 
    #Loss for each learner.
    loss=loss, 
    #Final loss for the online, regular and individual SL
    losses_all=losses_all,
    #Learners used.
    learners=learners, 
    #Evaluation done at time t.
    t=t, 
    #Evaluation done at horizon t+h.
    h=h  
  ))
}