
# id : Subject id to simulate from
# t : Number of times we want to interrupt the individual-based-model time-series.
# fit: Object of sample_and_fit
# df:  Data
# W:   Covariates to be included in the final dataset
# var: 
# historical_true : Include the observed historical data, or simulate based on the pooled fit.
# nsim	: Number of periods for the simulated series.
# bootstrap : Do simulation using resampled errors rather than normally 
#             distributed errors or errors provided as innov.
# future : Produce sample paths that are future to and conditional on the data 
#          in object. Otherwise simulate unconditionally.
# seed : Either NULL or an integer that will be used in a call to set.seed 
#        before simulating the time series.
# noise : If true, we add extra white noise to the model.

run_simulation <- function(id=13569, t=0, fit=NULL, df=NULL, W=NULL, var="abpmean",
                           nsim = 300, seed = 4197, future = TRUE, historical_true=FALSE,
                           bootstrap = TRUE, noise=TRUE) {
  
  #Get samples used for the pooled fit:
  id_pool <- fit$fit_pooled$accuracy$subject_id
  id_pool <- id_pool[!(id_pool %in% id)]
  
  W_id <- df[df$subject_id==id,W]
  W_pool <- df[df$subject_id %in% id_pool, W]
  
  ids_pool <- df[df$subject_id %in% id_pool, "subject_id"]

  #Get the pooled model
  model_pool <- fit$fit_pooled$fit
  #Get the individual model
  model_individual <- fit$fit_individuals[[as.character(id)]]$fit

  #Simulate individual and pooled time-series:
  #PROBLEM: Oh, oh: might result in negative values...
  abpmean_truth <- fit$fit_individuals[[as.character(id)]]$plot_data$test
  abpmean_individual <- ts(simulate(model_individual, nsim = nsim, bootstrap = bootstrap, 
                                    future = future, seed = seed))
  abpmean_pooled <- ts(simulate(model_pool, nsim = nsim, bootstrap = bootstrap, 
                                future = future, seed = seed))
  
  if(historical_true){
    abpmean_historical <- data.frame(df[df$subject_id %in% id_pool, c("subject_id", "abpmean")])
    abpmean_historical$subject_id <- as.numeric(levels(abpmean_historical$subject_id))[abpmean_historical$subject_id]
    abpmean_historical<-split(abpmean_historical, f = abpmean_historical$subject_id)
    abpmean_historical<-lapply(abpmean_historical, function(x) ts(x[,2]))
  }else{
    n_hist <- length(fit$fit_pooled$plot_data)-1
    abpmean_historical <- lapply(seq(n_hist), function(i){
      ts(simulate(model_pool, nsim = nsim, bootstrap = bootstrap, 
                  future = future))
    })
  }
  
  if(noise){
    abpmean_pooled <- abpmean_pooled + 2*arima.sim(model=list(order = c(0, 0, 0)), n=nsim)
    abpmean_individual <- abpmean_individual + 2*arima.sim(model=list(order = c(0, 0, 0)), n=nsim)
    abpmean_historical <- lapply(abpmean_historical, function(abpmean_hist){
      abpmean_hist + 2*arima.sim(model=list(order = c(0, 0, 0)), n=nsim) 
    })
  }
  
  time <- seq(1, nsim, 1)
  
  #TO DO: Recode this...
  #Combine pooled and individual model into one  based on 
  #the number of interruptions
  start<-1
  abpmean_combined <- NULL
  diff<-nsim/(t+1)
  for(i in 1:(t+1)){
    if((i %% 2 != 0) & (start<nsim)){
      abpmean_combined <- c(abpmean_combined, abpmean_individual[start:(diff*i)])
      start <- start + diff
    }else if ((i %% 2 == 0) & (start<nsim)){
      abpmean_combined <- c(abpmean_combined, abpmean_pooled[start:(diff*i)])
      start <- start + diff
    }
  }
  if(length(abpmean_combined)<nsim){
    diff<-nsim-length(abpmean_combined)
    abpmean_combined<-c(abpmean_combined,abpmean_individual[(nsim-diff+1):nsim])
  }
  
  #### Create individualized dataset:
  abpmean_combined <- ts(abpmean_combined)
  res <- data.frame(subject_id=id, W_id, abpmean=abpmean_combined)
  
  #Add outcome:
  #final <- new_Y_sol1(train_all = res, cutoff = 65)  
  res$Y_15 <- lead(abpmean_combined, 15)
  res$Y_20 <- lead(abpmean_combined, 20)
  res$Y_25 <- lead(abpmean_combined, 25)
  res$Y_30 <- lead(abpmean_combined, 30)
  
  #### Create historical dataset:
  abpmean_hist <- matrix(do.call(rbind, abpmean_historical), nrow = nrow(W_pool), ncol = 1)
  res_hist <- data.frame(subject_id=ids_pool, W_pool, abpmean=abpmean_hist)
  
  #Add outcome:
  #final <- new_Y_sol1(train_all = res, cutoff = 65)  
  res_hist <- res_hist %>%
    dplyr::group_by(subject_id) %>%
    dplyr::mutate(Y_15 = lead(abpmean, 15)) %>%
    dplyr::mutate(Y_20 = lead(abpmean, 20)) %>%
    dplyr::mutate(Y_25 = lead(abpmean, 25)) %>%
    dplyr::mutate(Y_30 = lead(abpmean, 30))

  return(list(data_ind=res,
              data_hist=res_hist,
              abpmean_combined=abpmean_combined,
              abpmean_individual=abpmean_individual,
              abpmean_pooled=abpmean_pooled,
              abpmean_truth=abpmean_truth))
}

