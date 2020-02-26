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
#stack_pool: sl3 stack object corresponding to the learners used for the Global Super Learner. 
#            The learners in it should not use screeners.
#stack_individual: sl3 stack object corresponding to the Individual Super Learner. 
#                  The learners in it should not use screeners.
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
#gap_training: time between the last trained time point and the first prediction time point in the training set cross-validation
#convex: TRUE or FALSE denoting whether the NN-restricted coefficients to be further restricted to have a fixed positive sum
#discrete: TRUE or FALSE denoting whether the cross-validation selector should be only lrnr given weight 
#outcome_summary: TRUE if the outcome is a pre-specified summary of the data, over a specified time interval. 
#                 This would imply that a \code{add_summary_outcome} was used on the data prior to running the 
#                 \code{combine_SL} function, thus already specifying the test window prediction size and horizon. 

combine_SL = function(train_all, outcome, t, stack_pool, stack_individual, 
                      stack_screen = NULL, covars, covars_baseline,
                      id = NULL, cv = "folds_rolling_origin", gap = 30, 
                      h = 15, test_size = 15, mini_batch = 15, 
                      first_window = 15, window_size = 15, gap_training = 30,
                      convex = TRUE, discrete = FALSE, outcome_summary=FALSE, 
                      adapt_covs=NULL, historical_fit=NULL){
  
  #Combine time-varying and baseline covariates
  covars_wbaseline <- c(covars, covars_baseline)
  
  #Convert all elements from factor:
  if(length(which( sapply( train_all, class ) == 'factor' )[-1])>0){
    w <- which( sapply( train_all, class ) == 'factor' )[-1]
    train_all[w] <- train_all[w] %>% retype()
  }

  #Get the number of unique samples:
  n <- length(unique(train_all$subject_id))

  #Allows for a single id combined SL
  if(is.null(id)){
    #Get all the individual samples:
    samples <- unique(train_all$subject_id)
  }else{
    #If id provided, then just the id:
    samples <- unique(data.frame(train_all[,id]))[,1]
  }
  
  #Check there is enough data:
  train_one <- train_all[train_all$subject_id %in% samples[1],]
  if(all(is.na(train_one[(t+1+gap):(t+gap+h),]))){
    warning("There is not enough time-points. Try making the t parameter smaller.")
    stop()
  }

  ##############################################################
  # Train all separate learners, predictions by sample
  ############################################################## 
  
  #Data design already takes into account the test size! 
  #(as a summary over a pre-specified region)
  if(outcome_summary){
    
    ### Train the Global SLs:
    global_SL_t<-global_SL(train_all=train_all, t=t, outcome=outcome, 
                           stack_pool=stack_pool, stack_screen=stack_screen,
                           covars=covars, covars_wbaseline=covars_wbaseline,
                           test_size=1, mini_batch=mini_batch,
                           cv=cv, first_window=first_window, 
                           window_size=window_size, gap = gap_training)
    global_SL_t_reg<-global_SL(train_all=train_all, t=t, outcome=outcome, 
                               stack_pool=stack_pool, stack_screen=stack_screen,
                               covars=covars, covars_wbaseline=covars_wbaseline,
                               test_size=1, cv="folds_vfold", V=10)
    ### Create individual SLs for all samples:
    individual_SL_t <- lapply(samples, function(x){
      individual_SL(train_all=train_all,t=t,id=x,cv=cv,first_window=first_window,
                    window_size=window_size, test_size=1,
                    mini_batch=mini_batch,covars=covars,
                    stack_individual=stack_individual, gap=gap_training)})
    
  }else{
    
    ### Train the Global SLs:
    global_SL_t<-global_SL(train_all=train_all, t=t, outcome=outcome, 
                           stack_pool=stack_pool, stack_screen=stack_screen,
                           covars=covars, covars_wbaseline=covars_wbaseline,
                           test_size=test_size, mini_batch=mini_batch,
                           cv=cv, first_window=first_window, 
                           window_size=window_size, gap = gap_training)
    global_SL_t_reg<-global_SL(train_all=train_all, t=t, outcome=outcome, 
                               stack_pool=stack_pool, stack_screen=stack_screen,
                               covars=covars, covars_wbaseline=covars_wbaseline,
                               test_size=test_size, cv="folds_vfold", V=5)
    ### Create individual SLs for all samples:
    individual_SL_t <- lapply(samples, function(x){
      individual_SL(train_all=train_all,t=t,id=x,cv=cv,first_window=first_window,
                    window_size=window_size, test_size=test_size,
                    mini_batch=mini_batch,covars=covars_wbaseline,
                    stack_individual=stack_individual, gap=gap_training)})
    if(!is.null(historical_fit)){
      
      hist_preds <- list() #Eh
      for(i in 1:length(individual_SL_t)){
        
        folds<-individual_SL_t[[i]]$folds
        training_task<-individual_SL_t[[i]]$task
        
        hist_preds[[i]] <- bind_rows(lapply(folds, function(fold) {
          test_set_in_training_task <- validation_task(training_task, fold)
          hist_fits <- cbind.data.frame(lapply(historical_fit, function(fit){
            fit$predict_fold(test_set_in_training_task, "full")
          }))
        }))
      }
      
      learner_names_historical<-colnames(hist_preds[[1]])
    }
    
  }
 
  #Get all the truths (for validation samples):
  truths <- global_SL_t$truths
  truths_training <- global_SL_t$truths_training
  truths_reg <- global_SL_t_reg$truths
  
  #Get all the learner names:
  learner_names_global<-colnames(global_SL_t$preds_globalSL[[1]])
  learner_names_global_screen<-colnames(global_SL_t$preds_globalSL_screen[[1]])
  learner_names_individual<-colnames(individual_SL_t[[1]]$preds_individualSL)
  
  ##############################################################
  #Extract predictions for all learners, by sample:
  ##############################################################
  
  if(!is.null(stack_screen)){
    #No historical data, with screeners
    if(is.null(historical_fit)){
      learners<-c(paste("GlobalSL", learner_names_global, sep="_"),
                  paste("GlobalSL_baseline", learner_names_global, sep="_"),
                  paste("GlobalSL_screen", learner_names_global_screen, sep="_"),
                  paste("GlobalSL_screenbaseline", learner_names_global_screen, sep="_"),
                  paste("IndividualSL", learner_names_individual, sep="_")) 
      
      learners_reg<-c(paste("GlobalSL", learner_names_global, sep="_"),
                      paste("GlobalSL_baseline", learner_names_global, sep="_"),
                      paste("GlobalSL_screen", learner_names_global_screen, sep="_"),
                      paste("GlobalSL_screenbaseline", learner_names_global_screen, sep="_"))
      
      learners_individual<-(paste("IndividualSL", learner_names_individual, sep="_"))
      
      #For sample i, get predictions:
      preds <- lapply(seq_len(n), function(i){
        ps<-cbind.data.frame(global_SL_t$preds_globalSL[[i]],
                             global_SL_t$preds_globalSL_baseline[[i]],
                             global_SL_t$preds_globalSL_screen[[i]],
                             global_SL_t$preds_globalSL_screen_baseline[[i]],
                             individual_SL_t[[i]]$preds_individualSL)
        names(ps)<-learners
        t(ps)})
      
      preds_reg <- lapply(seq_len(n), function(i){
        ps<-cbind.data.frame(global_SL_t_reg$preds_globalSL[[i]],
                             global_SL_t_reg$preds_globalSL_baseline[[i]],
                             global_SL_t_reg$preds_globalSL_screen[[i]],
                             global_SL_t_reg$preds_globalSL_screen_baseline[[i]])
        names(ps)<-learners_reg
        t(ps)})
      
      preds_individual <- lapply(seq_len(n), function(i){
        ps<-individual_SL_t[[i]]$preds_individualSL
        names(ps)<-learners_individual
        t(ps)
      })
    #With historical data, with screeners
    }else{
      learners<-c(paste("GlobalSL", learner_names_global, sep="_"),
                  paste("GlobalSL_baseline", learner_names_global, sep="_"),
                  paste("GlobalSL_screen", learner_names_global_screen, sep="_"),
                  paste("GlobalSL_screenbaseline", learner_names_global_screen, sep="_"),
                  paste("IndividualSL", learner_names_individual, sep="_"),
                  paste("HistoricalSL", learner_names_historical, sep="_")) 
      
      learners_reg<-c(paste("GlobalSL", learner_names_global, sep="_"),
                      paste("GlobalSL_baseline", learner_names_global, sep="_"),
                      paste("GlobalSL_screen", learner_names_global_screen, sep="_"),
                      paste("GlobalSL_screenbaseline", learner_names_global_screen, sep="_"))
      
      learners_individual<-(paste("IndividualSL", learner_names_individual, sep="_"))
      
      #For sample i, get predictions:
      preds <- lapply(seq_len(n), function(i){
        ps<-cbind.data.frame(global_SL_t$preds_globalSL[[i]],
                             global_SL_t$preds_globalSL_baseline[[i]],
                             global_SL_t$preds_globalSL_screen[[i]],
                             global_SL_t$preds_globalSL_screen_baseline[[i]],
                             individual_SL_t[[i]]$preds_individualSL,
                             hist_preds[[i]])
        names(ps)<-learners
        t(ps)
        print(i)})
      
      preds_reg <- lapply(seq_len(n), function(i){
        ps<-cbind.data.frame(global_SL_t_reg$preds_globalSL[[i]],
                             global_SL_t_reg$preds_globalSL_baseline[[i]],
                             global_SL_t_reg$preds_globalSL_screen[[i]],
                             global_SL_t_reg$preds_globalSL_screen_baseline[[i]])
        names(ps)<-learners_reg
        t(ps)})
      
      preds_individual <- lapply(seq_len(n), function(i){
        ps<-individual_SL_t[[i]]$preds_individualSL
        names(ps)<-learners_individual
        t(ps)
      })
    }
  ##No screeners
  }else{
    if(is.null(historical_fit)){
      learners<-c(paste("GlobalSL", learner_names_global, sep="_"),
                  paste("GlobalSL_baseline", learner_names_global, sep="_"),
                  paste("IndividualSL", learner_names_individual, sep="_"))
      
      learners_reg<-c(paste("GlobalSL", learner_names_global, sep="_"),
                      paste("GlobalSL_baseline", learner_names_global, sep="_"))
      
      learners_individual<-(paste("IndividualSL", learner_names_individual, sep="_"))
      
      #For sample i, get predictions:
      preds <- lapply(seq_len(n), function(i){
        ps<-cbind.data.frame(global_SL_t$globalSL_preds[[i]],
                             global_SL_t$globalSL_baseline_preds[[i]],
                             individual_SL_t[[i]]$preds_individualSL)
        names(ps)<-learners
        t(ps)})
      
      preds_reg <- lapply(seq_len(n), function(i){
        ps<-cbind.data.frame(global_SL_t_reg$globalSL_preds[[i]],
                             global_SL_t_reg$globalSL_baseline_preds[[i]])
        names(ps)<-learners_reg
        t(ps)})
      
      preds_individual <- lapply(seq_len(n), function(i){
        ps<-individual_SL_t[[i]]$preds_individualSL
        names(ps)<-learners_individual
        t(ps)
      })
    }else{
      learners<-c(paste("GlobalSL", learner_names_global, sep="_"),
                  paste("GlobalSL_baseline", learner_names_global, sep="_"),
                  paste("IndividualSL", learner_names_individual, sep="_"), 
                  paste("HistoricalSL", learner_names_historical, sep="_")
                  )
      
      learners_reg<-c(paste("GlobalSL", learner_names_global, sep="_"),
                      paste("GlobalSL_baseline", learner_names_global, sep="_"))
      
      learners_individual<-(paste("IndividualSL", learner_names_individual, sep="_"))
      
      #For sample i, get predictions:
      preds <- lapply(seq_len(n), function(i){
        ps<-cbind.data.frame(global_SL_t$globalSL_preds[[i]],
                             global_SL_t$globalSL_baseline_preds[[i]],
                             individual_SL_t[[i]]$preds_individualSL,
                             hist_preds[[i]])
        names(ps)<-learners
        t(ps)})
      
      preds_reg <- lapply(seq_len(n), function(i){
        ps<-cbind.data.frame(global_SL_t_reg$globalSL_preds[[i]],
                             global_SL_t_reg$globalSL_baseline_preds[[i]])
        names(ps)<-learners_reg
        t(ps)})
      
      preds_individual <- lapply(seq_len(n), function(i){
        ps<-individual_SL_t[[i]]$preds_individualSL
        names(ps)<-learners_individual
        t(ps)
      })
    }
  }
  
  ##############################################################
  ### Evaluate the loss for discrete learners, for each sample:
  ##############################################################
  
  #TO DO: Will this work with the summarized outcome?
  if(length(unique(train_all[,outcome]))>2){
    #Continuous outcome
    loss <- lapply(seq_len(n), function(i){
      eval_loss_cont(preds[[i]],truths[[i]])})
    loss_reg <- lapply(seq_len(n), function(i){
      eval_loss_cont(preds_reg[[i]],truths_reg[[i]])})
  }else{
    #Binary outcome
    loss <- lapply(seq_len(n), function(i){
      eval_loss(preds[[i]],truths[[i]])})
    loss_reg <- lapply(seq_len(n), function(i){
      eval_loss(preds_reg[[i]],truths_reg[[i]])})
  }
  
  ##############################################################
  #### Get weights, possibly adaptively, for each sample:
  ##############################################################
  
  #Get weights with combined SuperLearner (on validation samples)
  
  ##### Stratified on the sample id:
  
  #Get weights for just the Global SuperLearner
  weights <- lapply(seq_len(n), function(i){
    get_weights(preds[[i]], truths[[i]], loss[[i]], convex, discrete)
  })
  #Get weights for just the regular SuperLearner
  weights_reg_sl <- lapply(seq_len(n), function(i){
    get_weights(preds_reg[[i]], truths_reg[[i]], loss_reg[[i]],convex, discrete)
  })
  #Get weights for just the individual SuperLearner
  weights_ind_sl <- lapply(seq_len(n), function(i){
    get_weights(t(individual_SL_t[[i]]$preds_individualSL),
                truths[[i]], l=0, convex, FALSE)
  })
  
  ##### Global weights:
  
  #Get weights for just the Global SuperLearner
  weights_all <- get_weights(ps=rlist::list.cbind(preds), y=unlist(truths), l=unlist(loss),
                             convex=convex, discrete=discrete)
  #Get weights for just the regular SuperLearner
  weights_reg_sl_all <- get_weights(ps=rlist::list.cbind(preds_reg), y=unlist(truths_reg), l=unlist(loss_reg),
                                    convex=convex, discrete=discrete)
  #Get weights for just the individual SuperLearner
  weights_ind_sl_all <- get_weights(ps=rlist::list.cbind(preds_individual), y=unlist(truths), l=0,
                                    convex=convex, discrete=FALSE)
  
  #Caution: used ONLY for the Combined SL.
  if(!is.null(adapt_covs)){
    if(any(!(adapt_covs %in% covars_baseline))){
      warning("It's only possible to adapt weights to baseline covariates.")
      stop()
    }
    
    weights_adaptive <- get_weights_adaptive(train_all=train_all, adapt_covs=adapt_covs, preds=preds,
                                             truths=truths, truths_training=truths_training, outcome=outcome)
  }else{
    weights_adaptive <- NULL
  }
  
  #"weights_all"
  weights_id <- list(weights=weights,
                     weights_reg_sl=weights_reg_sl,
                     weights_ind_sl=weights_ind_sl)
  weights_all <- list(weights=weights_all,
                      weights_reg_sl=weights_reg_sl_all,
                      weights_ind_sl=weights_ind_sl_all)
  weights_adaptive <- list(weights=weights_adaptive)
  
  weights_combined <- list(weights_id$weights, 
                           weights_all$weights,
                           weights_adaptive$weights)
  
  #Discrete predictions:
  preds_all_training <- list(preds=preds,
                             preds_reg=preds_reg,
                             preds_ind=preds_individual)
  truths_all <- list(truths=truths, truths_reg=truths_reg)  
  
  ##############################################################
  #### Pick final SL weights based on the validation samples:
  ##############################################################
  
  #Combined SLs:
  preds_fin <- get_sl(predictions=preds_all_training$preds, weight=weights_id$weights, samples=samples)
  preds_fin_all <- get_sl(predictions=preds_all_training$preds, weight=weights_all$weights, samples=samples)
  preds_fin_adapt <- get_sl(predictions=preds_all_training$preds, weight=weights_adaptive$weights, samples=samples)
  preds_combined <- list(preds_fin,preds_fin_all,preds_fin_adapt)
  
  #Regular SLs:
  pred_regular <- get_sl(predictions=preds_all_training$preds_reg, weight=weights_id$weights_reg_sl, samples=samples)
  pred_regular_all <- get_sl(predictions=preds_all_training$preds_reg, weight=weights_all$weights_reg_sl, samples=samples)
  
  #Individualized SLs:
  pred_individual <- get_sl(predictions=preds_all_training$preds_ind, weight=weights_id$weights_ind_sl, samples=samples)
  pred_individual_all <- get_sl(predictions=preds_all_training$preds_ind, weight=weights_all$weights_ind_sl, 
                                samples=samples)
  
  ### Evaluate loss over final SuperLearners!
  #Note that this is the loss over validation samples. 
  #Evaluate the loss for discrete learners, for each sample:
  if(length(unique(train_all[,outcome]))>2){
    #Continuous outcome
    loss_fin_SL <- get_sl_loss(predictions=preds_fin, truths=truths_all$truths, binary = FALSE)
    loss_fin_all_SL <- get_sl_loss(predictions=preds_fin_all, truths=truths_all$truths, binary = FALSE)
    loss_fin_adapt_SL <- get_sl_loss(predictions=preds_fin_adapt, truths=truths_all$truths, binary = FALSE)
    
    #Get only the most relevant options for Regular and Individualized:
    loss_regular_SL <- get_sl_loss(predictions=pred_regular_all, truths=truths_all$truths_reg, binary = FALSE)
    loss_individual_SL <- get_sl_loss(predictions=pred_individual, truths=truths_all$truths, binary = FALSE)
    
    winner_index <-which.min(c(loss_fin_SL$loss, loss_fin_all_SL$loss, loss_fin_adapt_SL$loss))
    preds_combined <- preds_combined[[winner_index]]
    weights_combined <- weights_combined[[winner_index]]
  }else{
    #Binary outcome
    loss_fin_SL <- get_sl_loss(predictions=preds_fin, truths=truths_all$truths)
    loss_fin_all_SL <- get_sl_loss(predictions=preds_fin_all, truths=truths_all$truths)
    loss_fin_adapt_SL <- get_sl_loss(predictions=preds_fin_adapt, truths=truths_all$truths)

    #Get only the most relevant options for Regular and Individualized:
    loss_regular_SL <- get_sl_loss(predictions=pred_regular_all, truths=truths_all$truths_reg)
    loss_individual_SL <- get_sl_loss(predictions=pred_individual, truths=truths_all$truths)
    
    winner_index <-which.min(c(loss_fin_SL$loss, loss_fin_all_SL$loss, loss_fin_adapt_SL$loss))
    preds_combined_win <- preds_combined[[winner_index]]
    weights_combined_win <- weights_combined[[winner_index]]
  }
  
  preds_cv_validation <- list(preds_fin=preds_combined_win,
                              preds_regular=pred_regular_all,
                              preds_individual=pred_individual)
  
  ##############################################################
  ### Get all predictions for validation data:
  #   1. Data-adaptively pick the weights here? 
  #      (do we need another layer of CV)
  ##############################################################
  
  # create the sl3 task:
  tasks <- lapply(samples, function(x){
    train_one <- train_all[train_all$subject_id %in% x,]
    make_sl3_Task(data = train_one[(t+1+gap):(t+gap+h),], 
                  covariates = covars, outcome = outcome, 
                  folds = 1)
  })
  
  # create the sl3 task with baseline covariates:
  tasks_baseline <- lapply(samples, function(x){
    train_one <- train_all[train_all$subject_id %in% x,]
    make_sl3_Task(data = train_one[(t+1+gap):(t+gap+h),], 
                  covariates = covars_wbaseline, outcome = outcome,
                  folds = 1)
  })
  
  #### Regular SL: (using sl3 directly)
  pred_global_SL_reg <- lapply(tasks, function(task){
    p<-global_SL_t_reg$globalSL$predict(task)
    as.matrix(p)
  })
  pred_global_SL_baseline_reg <- lapply(tasks_baseline, function(task){
    p<-global_SL_t_reg$globalSL_baseline$predict(task)
    as.matrix(p)
  })
  if(!is.null(stack_screen)){
    pred_global_SL_screen_reg <- lapply(tasks, function(task){
      p<-global_SL_t_reg$globalSL_screen$predict(task)
      as.matrix(p)
    })
    pred_global_SL_screen_baseline_reg <- lapply(tasks_baseline, function(task){
      p<-global_SL_t_reg$globalSL_screen_baseline$predict(task)
      as.matrix(p)
    })
  }else{
    pred_global_SL_screen_reg<-NULL
    pred_global_SL_screen_baseline_reg<-NULL
  }
  
  #### Global SL:
  pred_global_SL <- lapply(tasks, function(task){
    p<-global_SL_t$globalSL$predict(task)
    as.matrix(p)
  })
  pred_global_SL_baseline <- lapply(tasks_baseline, function(task){
    p<-global_SL_t$globalSL_baseline$predict(task)
    as.matrix(p)
  })
  if(!is.null(stack_screen)){
    pred_global_SL_screen <- lapply(tasks, function(task){
      p<-global_SL_t$globalSL_screen$predict(task)
      as.matrix(p)
    })
    pred_global_SL_screen_baseline <- lapply(tasks_baseline, function(task){
      p<-global_SL_t$globalSL_screen_baseline$predict(task)
      as.matrix(p)
    })
  }else{
    pred_global_SL_screen<-NULL
    pred_global_SL_screen_baseline<-NULL
  }
  
  #### Individual SL:
  pred_individual_SL<-lapply(seq_along(tasks), function(i){
    t(individual_SL_t[[i]]$individualSL$predict(tasks[[i]]))
  })
  
  #Get the truth:
  truths_validation <- lapply(samples, function(x){
    train_one <- train_all[train_all$subject_id %in% x,]
    train_one[(t+1+gap):(t+gap+h), outcome]
  })
  
  if(!is.null(stack_screen)){
    learners<-c(paste("GlobalSL", learner_names_global, sep="_"),
                paste("GlobalSL_baseline", learner_names_global, sep="_"),
                paste("GlobalSL_screen", learner_names_global_screen, sep="_"),
                paste("GlobalSL_screenbaseline",learner_names_global_screen,sep="_"),
                paste("IndividualSL", learner_names_individual, sep="_")) 
    
    learners_reg<-c(paste("GlobalSL", learner_names_global, sep="_"),
                    paste("GlobalSL_baseline", learner_names_global, sep="_"),
                    paste("GlobalSL_screen", learner_names_global_screen, sep="_"),
                    paste("GlobalSL_screenbaseline",learner_names_global_screen,sep="_"))
    
    #For sample i, get predictions:
    preds_validation <- lapply(seq_along(tasks), function(i){
      ps<-cbind.data.frame(pred_global_SL[[i]],
                           pred_global_SL_baseline[[i]],
                           pred_global_SL_screen[[i]],
                           pred_global_SL_screen_baseline[[i]],
                           pred_individual_SL[[i]])
      names(ps)<-learners
      t(ps)})
    
    preds_reg_validation <- lapply(seq_along(tasks), function(i){
      ps<-cbind.data.frame(pred_global_SL_reg[[i]],
                           pred_global_SL_baseline_reg[[i]],
                           pred_global_SL_screen_reg[[i]],
                           pred_global_SL_screen_baseline_reg[[i]])
      names(ps)<-learners_reg
      t(ps)})
    
  }else{
    learners<-c(paste("GlobalSL", learner_names_global, sep="_"),
                paste("GlobalSL_baseline", learner_names_global, sep="_"),
                paste("IndividualSL", learner_names_individual, sep="_")) 
    
    learners_reg<-c(paste("GlobalSL", learner_names_global, sep="_"),
                    paste("GlobalSL_baseline", learner_names_global, sep="_")) 
    
    #For sample i, get predictions:
    preds_validation <- lapply(seq_along(tasks), function(i){
      ps<-cbind.data.frame(pred_global_SL[[i]],
                           pred_global_SL_baseline[[i]],
                           pred_individual_SL[[i]])
      names(ps)<-learners
      t(ps)})
    
    preds_reg_validation <- lapply(seq_along(tasks), function(i){
      ps<-cbind.data.frame(pred_global_SL_reg[[i]],
                           pred_global_SL_baseline_reg[[i]])
      names(ps)<-learners_reg
      t(ps)})
  }
  
  #### Get the final ensemble predictions on an extrenal dataset
  #Get the prediction with weighted combination of learners:
  preds_fin_validation <- get_sl(predictions = preds_validation, weight=weights_combined_win, samples = samples)
  preds_regular_validation <- get_sl(predictions = preds_reg_validation, weight=weights_all$weights_reg_sl, 
                                     samples = samples)
  preds_individual_validation <- get_sl(predictions = pred_individual_SL, weight=weights_id$weights_ind_sl, 
                                     samples = samples)

  ### Evaluate loss over final SuperLearners!
  #Note that this is the loss over validation samples. Evaluate the loss for discrete learners, for each sample:
  if(length(unique(train_all[,get(outcome)]))>2){
    #Continuous outcome
    #Binary outcome
    loss_fin_SL_validation <- get_sl_loss(predictions=preds_fin_validation, truths=truths_validation, binary = FALSE)
    loss_regular_SL_validation <- get_sl_loss(predictions=preds_regular_validation, 
                                              truths=truths_validation, binary = FALSE)
    loss_individual_SL_validation <- get_sl_loss(predictions=preds_individual_validation, 
                                                 truths=truths_validation, binary = FALSE)
  }else{
    #Binary outcome
    loss_fin_SL_validation <- get_sl_loss(predictions=preds_fin_validation, truths=truths_validation)
    loss_regular_SL_validation <- get_sl_loss(predictions=preds_regular_validation, truths=truths_validation)
    loss_individual_SL_validation <- get_sl_loss(predictions=preds_individual_validation, truths=truths_validation)
  }
  
  losses_validation<-cbind.data.frame(loss_fin_SL=loss_fin_SL_validation,
                               loss_regular_SL=loss_regular_SL_validation,
                               loss_individual_SL=loss_individual_SL_validation)
  
  return(list(
    #Final, weighted prediction (outside validation samples)
    preds_fin=preds_fin_validation, 
    #Final, weighted prediction for the regular SL (outside validation samples)
    pred_regular_SL=preds_regular_validation,
    #Final, weighted prediction for the one sample SL (outside validation samples)
    pred_individual_SL=preds_individual_validation,
    
    #SuperLearner predictions from interval validation:
    preds_cv_validation=preds_cv_validation,
    #Predictions for each individual learner.
    preds_all=preds_all_training,
    
    #t+h truth for all the samples (outside validation samples)
    truth_validation = truths_validation,
    #Combined SL truths (interval validation samples):
    truth_cv_validation = global_SL_t$truths, 
    #Combined SL truths (interval training samples):
    truth_cv_training = global_SL_t$truths_training,
    #Regular SL truths (interval validation samples):
    truth_cv_regular = global_SL_t_reg$truths, 
    
    #Final loss for the online, regular and individual SL (outside validation samples)
    losses_validation=losses_validation,
    #Final loss for the online, regular and individual SL (interval validation samples)
    losses_cv_validation=cbind(SL_id=loss_fin_SL$loss, SL_all=loss_fin_all_SL$loss, SL_adapt=loss_fin_adapt_SL$loss),
    
    #Data-adaptively picked weights:
    weights_winner = weights_combined_win,
    #All considered weights:
    weights_all = weights_combined,
    
    #Learners used.
    learners=learners, 
    #Evaluation done at time t.
    t=t, 
    #Evaluation done at horizon t+h.
    h=h
  ))
}
