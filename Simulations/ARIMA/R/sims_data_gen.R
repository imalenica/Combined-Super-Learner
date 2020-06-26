### Data generating processes used in simulations
suppressMessages(library(gratis))
suppressMessages(library(tsfeatures))

#ARIMA
data_gen_v0 <- function(n,t){
  sim_historical <- list() 
  ts_offset <- function(W){
    start <- 0.5*W[,"care_unit"] + 0.02*W[,"age"] + 0.5*W[,"sex"]
    return(start)
  }
  
  #Construct historical time-series
  for(i in 1:n){
    W <- cbind.data.frame(id        = rep(i, t),
                          sex       = rbinom(n = t, size = 1, prob = 0.5),
                          age       = round(runif(n = t, min = 19, max=90)),
                          care_unit = round(runif(n = t, min = 0, max=2)))
    
    sim_ts <- as.numeric(75 + 
                           arima.sim(model=list(ar=c(0.6,0.4,0.1,-0.1,-0.05)),n=t))
    sim_historical[[i]] <- cbind.data.frame(series   = sim_ts,
                                            time     = seq(1:t),
                                            lag1     = lead(sim_ts, n = 1),
                                            lag2     = lead(sim_ts, n = 2),
                                            lag3     = lead(sim_ts, n = 3),
                                            lag4     = lead(sim_ts, n = 4),
                                            lag5     = lead(sim_ts, n = 5),
                                            W,
                                            Y        = lead(sim_ts, n = 15))
    
  }
  #Remove NAs: issues with sl3
  sim_historical <- do.call(rbind,sim_historical)
  sim_hist_cc <- complete.cases(sim_historical)
  sim_historical <- sim_historical[sim_hist_cc,]
  
  #Construct individual time-series
  W <- cbind.data.frame(id       = rep(n+1, 1),
                        sex      = rbinom(n = 1, size = 1, prob = 0.5),
                        age      = round(runif(n = 1, min = 19, max=90)),
                        care_unit = round(runif(n = 1, min = 0, max=2)))
  
  sim_ts <- as.numeric(75 + 
                         arima.sim(model=list(ma=c(-0.5,-0.3,-0.1,0.2,0.5)),n=t))
  sim_individual <- cbind.data.frame(series   = sim_ts,
                                     time     = seq(1:t),
                                     lag1     = lead(sim_ts, n = 1),
                                     lag2     = lead(sim_ts, n = 2),
                                     lag3     = lead(sim_ts, n = 3),
                                     lag4     = lead(sim_ts, n = 4),
                                     lag5     = lead(sim_ts, n = 5),
                                     W,
                                     Y        = lead(sim_ts, n = 15))
  
  #Construct full data set
  sim_all <- rbind.data.frame(sim_historical, sim_individual)
  sim_all_cc <- complete.cases(sim_all)
  sim_all <- sim_all[sim_all_cc,]
  
  return(list(full=sim_all, 
              historical=sim_historical, 
              individual=sim_individual))
}

#ARIMA with W offset
data_gen_v1 <- function(n,t){
  sim_historical <- list() 
  ts_offset <- function(W){
    start <- 0.5*W[,"care_unit"] + 0.02*W[,"age"] + 0.5*W[,"sex"]
    return(start)
  }
  
  #Construct historical time-series
  for(i in 1:n){
    W <- cbind.data.frame(id       = rep(i, t),
                          sex      = rbinom(n = t, size = 1, prob = 0.5),
                          age      = round(runif(n = t, min = 19, max=90)),
                          care_unit = round(runif(n = t, min = 0, max=2)))
    
    sim_ts <- as.numeric(75 + ts_offset(W) + 
                           arima.sim(model=list(ar=c(0.6,0.4,0.1,-0.1,-0.05)),n=t))
    sim_historical[[i]] <- cbind.data.frame(series   = sim_ts,
                                            time     = seq(1:t),
                                            lag1     = lead(sim_ts, n = 1),
                                            lag2     = lead(sim_ts, n = 2),
                                            lag3     = lead(sim_ts, n = 3),
                                            lag4     = lead(sim_ts, n = 4),
                                            lag5     = lead(sim_ts, n = 5),
                                            W,
                                            Y        = lead(sim_ts, n = 15))
    
  }
  #Remove NAs: issues with sl3
  sim_historical <- do.call(rbind,sim_historical)
  sim_hist_cc <- complete.cases(sim_historical)
  sim_historical <- sim_historical[sim_hist_cc,]
  
  #Construct individual time-series
  W <- cbind.data.frame(id       = rep(n+1, 1),
                        sex      = rbinom(n = 1, size = 1, prob = 0.5),
                        age      = round(runif(n = 1, min = 19, max=90)),
                        care_unit = round(runif(n = 1, min = 0, max=2)))
  
  sim_ts <- as.numeric(75 + ts_offset(W) + 
                         arima.sim(model=list(ma=c(-0.5,-0.3,-0.1,0.2,0.5)),n=t))
  sim_individual <- cbind.data.frame(series   = sim_ts,
                                     time     = seq(1:t),
                                     lag1     = lead(sim_ts, n = 1),
                                     lag2     = lead(sim_ts, n = 2),
                                     lag3     = lead(sim_ts, n = 3),
                                     lag4     = lead(sim_ts, n = 4),
                                     lag5     = lead(sim_ts, n = 5),
                                     W,
                                     Y        = lead(sim_ts, n = 15))
  #Construct full data set
  sim_all <- rbind.data.frame(sim_historical, sim_individual)
  sim_all_cc <- complete.cases(sim_all)
  sim_all <- sim_all[sim_all_cc,]
  
  return(list(full=sim_all, 
              historical=sim_historical, 
              individual=sim_individual))
}

#ARIMA mixtures
data_gen_v2 <- function(n,t){
  sim_historical <- list() 
  ts_offset <- function(W){
    start <- 0.5*W[,"care_unit"] + 0.02*W[,"age"] + 0.5*W[,"sex"]
    return(start)
  }
  
  #Construct historical time-series
  for(i in 1:n){
    W <- cbind.data.frame(id       = rep(i, t),
                          sex      = rbinom(n = t, size = 1, prob = 0.5),
                          age      = round(runif(n = t, min = 19, max=90)),
                          care_unit = round(runif(n = t, min = 0, max=2)))
    
    sim_ts <- as.numeric(75 + ts_offset(W) + 
                           arima.sim(model=list(ar=c(0.6,0.4,0.1,-0.1,-0.05)),n=t))
    sim_historical[[i]] <- cbind.data.frame(series   = sim_ts,
                                            time     = seq(1:t),
                                            lag1     = lead(sim_ts, n = 1),
                                            lag2     = lead(sim_ts, n = 2),
                                            lag3     = lead(sim_ts, n = 3),
                                            lag4     = lead(sim_ts, n = 4),
                                            lag5     = lead(sim_ts, n = 5),
                                            W,
                                            Y        = lead(sim_ts, n = 15))
    
  }
  #Remove NAs: issues with sl3
  sim_historical <- do.call(rbind,sim_historical)
  sim_hist_cc <- complete.cases(sim_historical)
  sim_historical <- sim_historical[sim_hist_cc,]
  
  #Construct individual time-series
  W <- cbind.data.frame(id       = rep(n+1, 1),
                        sex      = rbinom(n = 1, size = 1, prob = 0.5),
                        age      = round(runif(n = 1, min = 19, max=90)),
                        care_unit = round(runif(n = 1, min = 0, max=2)))
  
  sim_ts <- as.numeric(75 +  
                         c(arima.sim(model=list(ma=c(-0.5,-0.3,-0.1,0.2,0.5)),n=(1*t/2)),
                           arima.sim(model=list(ar=c(0.6,0.4,0.1,-0.1,-0.05)),n=(1*t/2))))
  sim_individual <- cbind.data.frame(series   = sim_ts,
                                     time     = seq(1:t),
                                     lag1     = lead(sim_ts, n = 1),
                                     lag2     = lead(sim_ts, n = 2),
                                     lag3     = lead(sim_ts, n = 3),
                                     lag4     = lead(sim_ts, n = 4),
                                     lag5     = lead(sim_ts, n = 5),
                                     W,
                                     Y        = lead(sim_ts, n = 15))
  
  #Construct full data set
  sim_all <- rbind.data.frame(sim_historical, sim_individual)
  sim_all_cc <- complete.cases(sim_all)
  sim_all <- sim_all[sim_all_cc,]
  
  return(list(full=sim_all, 
              historical=sim_historical, 
              individual=sim_individual))
}

#Mixture autoregressive models (MAR)
data_gen_v3 <- function(n,t){
  sim_historical <- list() 
  ts_offset <- function(W){
    start <- 0.5*W[,"care_unit"] + 0.02*W[,"age"] + 0.5*W[,"sex"]
    return(start)
  }
  
  #Construct historical time-series
  for(i in 1:n){
    W <- cbind.data.frame(id       = rep(i, t),
                          sex      = rbinom(n = t, size = 1, prob = 0.5),
                          age      = round(runif(n = t, min = 19, max=90)),
                          care_unit = round(runif(n = t, min = 0, max=2)))
    
    sim_ts <- as.numeric(75 + ts_offset(W) + 
                           generate_ts_with_target(n = 1, ts.length = t, freq = 1, seasonal = 0,
                                                   features = c("entropy", "stability", "stl_features"),
                                                   selected.features = c("entropy", "stability", "trend"),
                                                   target = c(0.7, 1.7, 0.05)))
    
    sim_historical[[i]] <- cbind.data.frame(series   = sim_ts,
                                            time     = seq(1:t),
                                            lag1     = lead(sim_ts, n = 1),
                                            lag2     = lead(sim_ts, n = 2),
                                            lag3     = lead(sim_ts, n = 3),
                                            lag4     = lead(sim_ts, n = 4),
                                            lag5     = lead(sim_ts, n = 5),
                                            W,
                                            Y        = lead(sim_ts, n = 15))
    
  }
  #Remove NAs: issues with sl3
  sim_historical <- do.call(rbind,sim_historical)
  sim_hist_cc <- complete.cases(sim_historical)
  sim_historical <- sim_historical[sim_hist_cc,]
  
  #Construct individual time-series
  W <- cbind.data.frame(id       = rep(n+1, t),
                        sex      = rbinom(n = 1, size = 1, prob = 0.5),
                        age      = round(runif(n = 1, min = 19, max=90)),
                        care_unit = round(runif(n = 1, min = 0, max=2)))
  
  sim_ts <- as.numeric(75 + ts_offset(W) + 
                         generate_ts_with_target(n = 1, ts.length = t, freq = 1, seasonal = 0,
                                                 features = c("entropy", "stability", "stl_features"),
                                                 selected.features = c("entropy", "stability", "trend"),
                                                 target = c(0.7, 0.1, 0.05)))
  
  sim_individual <- cbind.data.frame(series   = sim_ts,
                                     time     = seq(1:t),
                                     lag1     = lead(sim_ts, n = 1),
                                     lag2     = lead(sim_ts, n = 2),
                                     lag3     = lead(sim_ts, n = 3),
                                     lag4     = lead(sim_ts, n = 4),
                                     lag5     = lead(sim_ts, n = 5),
                                     W,
                                     Y        = lead(sim_ts, n = 15))
  #Construct full data set
  sim_all <- rbind.data.frame(sim_historical, sim_individual)
  sim_all_cc <- complete.cases(sim_all)
  sim_all <- sim_all[sim_all_cc,]
  
  return(list(full=sim_all, 
              historical=sim_historical, 
              individual=sim_individual))
}