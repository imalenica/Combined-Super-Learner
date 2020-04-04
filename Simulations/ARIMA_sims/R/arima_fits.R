### Fits individual and pooled ARIMA models

sample_and_fit = function(dat=NULL, pool_size=100, ind_size=100, 
                          variable= "abpmean", seed=11){
  
  set.seed(seed)
  
  #Get samples
  samples <- unique(dat$subject_id)
  
  #Get a pooled model over a subsample:
  pool_samples <- base::sample(samples, size=pool_size)
  pool_data<-dat[dat$subject_id %in% pool_samples,]
  #Get individual samples to fit: 
  ind_samples <- base::sample(samples, size=ind_size)
  ind_data<-dat[dat$subject_id %in% ind_samples,]
  
  fit_pooled <- run_auto_arima(df=pool_data, variable = variable, pool=TRUE)
  fit_individuals<-lapply(as.list(ind_samples), function(id){
    run_auto_arima(df=ind_data, variable = variable, id=id, pool=FALSE)})
  names(fit_individuals) <- ind_samples
  
  return(list(fit_pooled=fit_pooled,
              fit_individuals=fit_individuals,
              pool_data=pool_data,
              ind_data=ind_data,
              pool_samples=pool_samples,
              ind_samples=ind_samples))
}

auto_arima <- function(df, variable, x){
  
  # 1. subset data to only contain that sample and make sure it's in order
  sub <- df[df$subject_id == x, ]
  sub_ord <- sub[order(sub$time_and_date), ]
  
  # 2. split first 80% of data into training set, and last 20% into test set
  split <- round(nrow(sub_ord) * .8)
  train <- sub_ord[1:split, ]
  train_var <- as.numeric(unlist(train[,which(colnames(train)==variable)]))
  test <- sub_ord[(split + 1):nrow(sub_ord), ]
  test_var <- as.numeric(unlist(test[,which(colnames(test) == variable)]))
  
  # 3. fit auto ARIMA
  fit <- auto.arima(y = train_var, stepwise = FALSE,
                    approximation = FALSE)
  
  # 4. obtain model coefficients                                          
  fit_coef <- fit$coef
  
  # 5. calculate residuals for subseuqent ACF and PACF plotting
  fit_resid <- residuals(fit)
  
  # 6. measure accuracy on the test set
  model_test <- Arima(test_var, model = fit)
  onestep_forecast <- fitted(model_test)
  # accuracy of the one-step ahead out of sample forecasts
  onestep_accuracy <- round(accuracy(model_test), 4)
  
  plot_info <- list(onestep_forecast = onestep_forecast,
                    test = test_var)
  
  # 7. make relevant results pretty and return them
  accuracy_onestep <- onestep_accuracy["Training set",]
  
  return(list(fit = fit,
              coefficients = fit_coef,
              accuracy = accuracy_onestep,
              plot_info = plot_info,
              residuals = fit_resid))
}

run_auto_arima <- function(df, variable = "abpmean", id=NULL, pool=TRUE){
  
  if(pool){
    
    #######################################
    ### Pooled ARIMA, across all samples
    #######################################
    
    samples <- unique(df$subject_id)
    
    train_test_samples <- lapply(samples, function(x){
      # 1. subset data to only contain that sample and make sure it's in order
      sub <- df[df$subject_id == x, ]
      sub_ord <- sub[order(sub$time_and_date), ]
      
      # 2. split first 80% of data into training set, and last 20% into test set
      split <- round(nrow(sub_ord) * .8)
      train <- sub_ord[1:split, ]
      train_var <- as.numeric(unlist(train[,which(colnames(train)==variable)]))
      test <- sub_ord[(split + 1):nrow(sub_ord), ]
      test_var <- as.numeric(unlist(test[,which(colnames(test) == variable)]))
      
      return(list(train=train, test=test,
                  train_var= train_var,
                  test_var=test_var))
    })
    
    #Pad individual time-series with a 120 min period.
    train_test_samples<-lapply(train_test_samples, function(x){
      x$train_var <- c(rep(NA, 120), x$train_var)
      return(list(train_var=x$train_var,
                  test_var=x$test_var))
    })
    
    #Combined all the train time-series into one:
    train_pooled <- ts(unlist(lapply(train_test_samples, function(x){x$train_var})))
    test_pooled <- ts(unlist(lapply(train_test_samples, function(x){x$test_var})))
    
    # 3. fit auto ARIMA
    fit <- auto.arima(y = train_pooled, stepwise = FALSE,
                      approximation = FALSE)
    
    # 4. obtain model coefficients                                          
    fit_coef <- fit$coef
    
    # 5. calculate residuals for subseuqent ACF and PACF plotting
    fit_resid <- residuals(fit)
    
    # 6. measure pooled accuracy on the test set:
    model_test <- Arima(test_pooled, model = fit)
    onestep_forecast <- fitted(model_test)
    onestep_accuracy <- round(accuracy(model_test), 4)
    plot_info_pooled <- list(onestep_forecast = onestep_forecast,
                             test = test_pooled)
    accuracy_onestep_pooled <- onestep_accuracy["Training set",]
    
    # 7. measure individual accuracy
    indiv_stats<-lapply(train_test_samples, function(x){
      model_test <- Arima(x$test, model = fit)
      onestep_forecast <- fitted(model_test)
      onestep_accuracy <- round(accuracy(model_test), 4)
      plot_info <- list(onestep_forecast = onestep_forecast,
                        test = x$test)
      accuracy_onestep <- onestep_accuracy["Training set",]
      return(list(accuracy_onestep=accuracy_onestep,
                  plot_info=plot_info))
    })
    list_acc <- lapply(indiv_stats, function(x) data.frame(as.list(x$accuracy_onestep)))
    list_acc <- list_acc[lapply(list_acc, length) > 0]
    accuracy <- rbindlist(list_acc, fill = TRUE)
    accuracy <- data.frame(subject_id = samples, accuracy)
    plot_info <- lapply(indiv_stats, function(x) x$plot_info)
    names(plot_info) <- samples
    
    return(list(fit = fit,                                    #Pooled ARIMA fit
                coefficients = fit_coef,                      #Coefs of the pooled ARIMA fit
                accuracy_pooled = accuracy_onestep_pooled,    #Stats on pooled test sets
                accuracy = accuracy,                          #Stats on individual test sets
                plot_data_pooled = plot_info_pooled,          #Pooled predictions and truth (test set)
                plot_data = plot_info,                        #Individual predictions and truth (test set)
                residuals = fit_resid                         #Residuals from the pooled training set
    ))
  }else if(!pool){
    
    #######################################
    ### Separate ARIMA for each time-series
    #######################################
    
    if(is.null(id)){
      
      #Fit on all samples in df
      samples <- unique(df$subject_id)
      
      # for each sample:
      fit_list <- lapply(samples, function(i) {auto_arima(df=df, x = i, variable=variable)})
      
      # make results pretty
      names(fit_list) <- samples
      
      list_coef <- lapply(fit_list, function(x) data.frame(as.list(x$coefficients)))
      list_coef <- list_coef[lapply(list_coef, length) > 0]
      coefficients <- rbindlist(list_coef, fill = TRUE)
      coefficients <- data.frame(subject_id = names(list_coef), coefficients)
      list_acc <- lapply(fit_list, function(x) data.frame(as.list(x$accuracy)))
      list_acc <- list_acc[lapply(list_acc, length) > 0]
      accuracy <- rbindlist(list_acc, fill = TRUE)
      accuracy <- data.frame(subject_id = names(list_acc), accuracy)
      plot_data <- lapply(fit_list, function(x) x$plot_info)
      names(plot_data) <- samples
      residuals <- lapply(fit_list, function(x) x$residuals)
      names(residuals) <- samples
      
      fit <- lapply(fit_list, function(x) x$fit)
      names(fit) <- samples
      
      return(list(fit = fit,
                  coefficients = coefficients,
                  accuracy = accuracy,
                  plot_data = plot_data,
                  residuals = residuals))
      
    }else if(!is.null(id)){
      
      # Fit only on the single specified sample
      
      fit_list <- auto_arima(df=df, x=id, variable=variable)
      
      return(list(fit = fit_list$fit,
                  coefficients = fit_list$coefficients,
                  accuracy = fit_list$accuracy,
                  plot_data = fit_list$plot_info,
                  residuals = fit_list$residuals))
    }
  }
}