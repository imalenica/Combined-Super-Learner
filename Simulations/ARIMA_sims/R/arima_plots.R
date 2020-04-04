##Plot prediction of the pooled and individualized ARIMA fit

plot_forecast <- function(fit, type=c("pooled", "individual", "time-series"), ts=NULL, id) {
  
  if(type=="pooled"){
    
    #Take the pooled prediction for sample id:
    onestep_forecast <- fit$fit_pooled$plot_data[[as.character(id)]]$onestep_forecast
    test_outcome <- fit$fit_pooled$plot_data[[as.character(id)]]$test
    
    return(ggplot() +
             geom_line(aes(x = as.numeric(time(test_outcome)),
                           y = as.numeric(test_outcome),
                           col = "blue")) +
             geom_line(aes(x = as.numeric(time(onestep_forecast)),
                           y = as.numeric(onestep_forecast),
                           col = "red")) +
             scale_color_discrete(name = "",
                                  labels = c("truth", "pooled forecast")) +
             labs(title = paste0(type, " fit for subject ", id),
                  x = "Time in test set (minutes)",
                  y = "Mean blood pressure") +
             theme(legend.position="top"))
  }else if(type=="individual"){
    
    #Take the individual prediction for sample id:
    onestep_forecast <- fit$fit_individuals[[as.character(id)]]$plot_data$onestep_forecast
    test_outcome <- fit$fit_individuals[[as.character(id)]]$plot_data$test
    
    return(ggplot() +
             geom_line(aes(x = as.numeric(time(test_outcome)),
                           y = as.numeric(test_outcome),
                           col = "blue")) +
             geom_line(aes(x = as.numeric(time(onestep_forecast)),
                           y = as.numeric(onestep_forecast),
                           col = "red")) +
             scale_color_discrete(name = "",
                                  labels = c("truth", "individual forecast")) +
             labs(title = paste0(type, " fit for subject ", id),
                  x = "Time in test set (minutes)",
                  y = "Mean blood pressure") +
             theme(legend.position="top"))
    
  }else if("time-series"){
    return(ggplot() +
             geom_line(aes(x = as.numeric(time(ts)),
                           y = as.numeric(ts),
                           col = "blue")) +
             scale_color_discrete(name = "",
                                  labels = c("time-series")) +
             labs(title = paste0("Time-Series for subject ", id),
                  x = "Time in test set (minutes)",
                  y = "Mean blood pressure") +
             theme(legend.position="top"))
  }
}
