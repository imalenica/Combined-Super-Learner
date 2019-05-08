#Evaluate the loss:
#Binary Cross-Entropy / Log Loss
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

#####################
### Global learner
#####################
global_SL = function(t){
  
  #Pool across time (fitting on all data up to time t)
  #Here, we imposse a Markov order assumption (independence between times)
  
  test_data <- train_all[train_all$time<=t,]
  
  #TO DO: Find a smart way to pull across time and samples, and have
  #validation be samples in the next time for all samples.
  
  #Set up proper cross-validation for time-series data:
  #folds <- origami::make_folds(test_data,
  #                             fold_fun = folds_rolling_window,
  #                             window_size = training_size,
  #                             validation_size = test_size, gap = 0,
  #                             batch = mini_batch)
  
  # create the sl3 task with 
  task <- make_sl3_Task(
    data = test_data, covariates = covars,
    outcome = outcome) #,folds=folds)
  
  # create the sl3 task:
  task_baseline <- make_sl3_Task(
    data = test_data, covariates = covars_wbaseline,
    outcome = outcome) #, folds=folds)
  
  #Fit the regular Super Learner (with baseline covariates)
  regularSL <- sl$train(task_baseline)
  
  #Fit the global learner:
  globalSL<-stack$train(task)
  globalSL_baseline<-stack$train(task_baseline)
  
  #Fit the global learner with screeners:
  globalSL_screen<-stack_screen$train(task)
  globalSL_screen_baseline<-stack_screen$train(task_baseline)
  
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
individual_SL = function(t,id,first_window,test_size,mini_batch){
  
  #Subset to sample with subject id id
  train_one <- train_all[train_all$subject_id %in% id,]
  
  #Pool across all samples and all time point until time t:
  test_data <- train_one[train_one$time<=t,]
  
  #Set up proper cross-validation for time-series data:
  folds <- origami::make_folds(test_data,
                               fold_fun = folds_rolling_origin,
                               first_window = first_window,
                               validation_size = test_size, gap = 0,
                               batch = mini_batch)
  
  # create the sl3 task
  task <- make_sl3_Task(
    data = test_data, covariates = covars_wbaseline,
    outcome = outcome, folds=folds)
  
  #Fit the regular Super Learner (with baseline covariates)
  #indregularSL <- sl$train(task) #Issue with folds for some reason
  indregularSL<-NULL
  
  #Fit the global learner:
  individualSL<-stack$train(task)
  
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
#t: time until which we train the SuperLearners
#gap: time between the last trained time point and the first prediction time point
#h: horizon at which we evaluate the loss
#train_all: the full dataset used
#first_window: size of the training set (for example, 1 to first window).
#test_size: size of the validation set (in time points).
#mini_batch: by how many time points to increase the training set after the first iteration?

combine_SL = function(t,gap,h,train_all,
                      first_window, test_size,mini_batch){
  
  num_p<-(t+h)-(t)
  
  #Get all the individual samples:
  samples <- unique(train_all$subject_id)
  
  #Train the Global SLs:
  global_SL_t<-global_SL(t=t)
  
  #Create individual SLs for all samples:
  individual_SL_t <- lapply(samples, function(x){individual_SL(t=t,id=x,
                                                               first_window,test_size,mini_batch)})
  
  # create the sl3 task:
  tasks <- lapply(samples, function(x){train_one <- train_all[train_all$subject_id %in% x,]
  make_sl3_Task(data = train_one[(t+1+gap):(t+gap+h),], covariates = covars, outcome = outcome)
  })
  
  # create the sl3 task with baseline covariates:
  tasks_baseline <- lapply(samples, function(x){train_one <- train_all[train_all$subject_id %in% x,]
  make_sl3_Task(data = train_one[(t+1+gap):(t+gap+h),], covariates = covars_wbaseline, outcome = outcome)
  })
  
  lrn<-length(global_SL_t$globalSL$print())
  
  ### Get all predictions:
  
  #Regular SL:
  pred_regular_SL <- lapply(tasks_baseline, function(task){global_SL_t$regularSL$predict(task)})
  #Global SL:
  pred_global_SL <- lapply(tasks, function(task){matrix(unlist(lapply(1:lrn, function(x){
    global_SL_t$globalSL$fit_object$learner_fits[[x]]$predict(task)})),nrow=num_p,ncol=lrn)})
  #Global SL with baseline covariates:
  pred_global_SL_baseline <- lapply(tasks_baseline, function(task){matrix(unlist(lapply(1:lrn, function(x){
    global_SL_t$globalSL_baseline$fit_object$learner_fits[[x]]$predict(task)})),nrow=num_p,ncol=lrn)})
  #Global SL with screeners:
  pred_global_SL_screen <- lapply(tasks, function(task){
    p<-global_SL_t$globalSL_screen$fit_object$learner_fits[[1]]$predict(task)
    as.matrix(p)
  })
  #Global SL with screeners and baseline covariates:
  pred_global_SL_screen_baseline <- lapply(tasks_baseline, function(task){
    p<-global_SL_t$globalSL_screen_baseline$fit_object$learner_fits[[1]]$predict(task)
    as.matrix(p)
  })
  #Individual SL:
  pred_individual_SL<-lapply(seq_along(tasks), function(i){
    matrix(unlist(lapply(1:lrn, function(x){
      individual_SL_t[[i]]$individualSL$fit_object$learner_fits[[x]]$predict(tasks_baseline[[i]])})),nrow=num_p,ncol=lrn)
  })
  #Regular Individual SL:
  #pred_regular_individual_SL<-lapply(seq_along(tasks), function(i){
  #  individual_SL_t[[i]]$indregularSL$predict(tasks_baseline[[i]])})
  
  #TO DO: allow for different learners for global and individual SLs
  learners<-c(paste("GlobalSL", global_SL_t$globalSL$print(), sep="_"),
              paste("GlobalSL_baseline", global_SL_t$globalSL$print(), sep="_"),
              paste("GlobalSL_screen", global_SL_t$globalSL$print(), sep="_"),
              paste("GlobalSL_screenbaseline_", global_SL_t$globalSL$print(), sep="_"),
              paste("IndividualSL", global_SL_t$globalSL$print(), sep="_")) 
  
  #For sample i, get predictions:
  preds <- lapply(seq_along(tasks), function(i){
    ps<-cbind.data.frame(pred_global_SL[[i]],
                         pred_global_SL_baseline[[i]],
                         pred_global_SL_screen[[i]],
                         pred_global_SL_screen_baseline[[i]],
                         pred_individual_SL[[i]])
    names(ps)<-learners
    t(ps)})
  
  #Get the truth:
  truths <- lapply(samples, function(x){
    train_one <- train_all[train_all$subject_id %in% x,]
    train_one[(t+1+gap):(t+gap+h),"event"][[1]]
  })
  
  #Evaluate the loss for discrete learners, for each sample:
  loss <- lapply(seq_along(tasks), function(i){eval_loss(preds[[i]],truths[[i]])})
  
  #Assign weights to discrete learners:
  get_weights = function(ps,y,l){
    #Gives a coefficient based on ALL the predictions (NOT time specific!)
    fit_coef <- nnls::nnls(t(ps), y)
    fit_coef<-fit_coef$x
    
    #PROBLEM: If no weights given, distribute weights based the loss?
    #if(sum(fit_coef)==0){
    #This "rewards" learnes with the smallest loss, but barely.
    #fit_coef<-as.numeric(l)
    #fit_coef<-sum(fit_coef)-fit_coef
    #fit_coef<-fit_coef/sum(fit_coef)
    
    #Give all weight to the discrete learner:
    #fit_coef[which.min(l)]<-1
    #}
    
    if(sum(fit_coef)>1){
      fit_coef<-fit_coef/sum(fit_coef)
    }
    return(fit_coef)
  }
  
  #Get weights for each sample:
  weights <- lapply(seq_along(tasks), function(i){
    get_weights(preds[[i]],truths[[i]],loss[[i]])
  })
  #Get the prediction with weighted combination of learners:
  preds_fin <- lapply(seq_along(tasks), function(i){
    #Generates weighted prediction for each validation time-point
    pred <- data.frame(pred = as.matrix(t(preds[[i]])) %*% weights[[i]])
    names(pred) <- paste0("Sample",samples[i])
    pred
  })
  #Connect weights and learners:
  fit_coef <- lapply(seq_along(tasks), function(i){
    data.frame(learners=learners,coefs=weights[[i]])
  })
  
  #Evaluate loss over final SuperLearners!
  #Note that this is the loss over validation samples
  loss_online_SL<-unlist(lapply(seq_along(tasks), function(i){eval_loss(t(preds_fin[[i]]),truths[[i]])}))
  loss_regular_SL<-unlist(lapply(seq_along(tasks), function(i){eval_loss(t(pred_regular_SL[[i]]),truths[[i]])}))
  #loss_individual_SL<-unlist(lapply(seq_along(tasks), function(i){eval_loss(t(pred_regular_individual_SL[[i]]),truths[[i]])}))
  losses_all<-cbind.data.frame(loss_online_SL=loss_online_SL,
                               loss_regular_SL=loss_regular_SL)
  #loss_individual_SL=loss_individual_SL
  
  return(list(#Final, weighted prediction.
    preds_fin=preds_fin, 
    #Final, weighted prediction for the regular SL
    pred_regular_SL=pred_regular_SL,
    #Final, weighted prediction for the one sample SL
    #pred_regular_individual_SL=pred_regular_individual_SL,
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