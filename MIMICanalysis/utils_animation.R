make_plots <- function(movie_tbl, times, xmin, xmax){
 for (i in 1:length(times)) {
    current_time <- times[i]
    truth <- movie_tbl[time <= current_time, c("truth", "time")]
    
    if(any(movie_tbl$real_time <= current_time, na.rm = T)){
      forecast <- movie_tbl[real_time <= current_time, c("forecast", "time")]
      pad_time <- seq(min(truth$time), min(forecast$time)-1)
      pad_forecast <- rep(NA, length(pad_time))
      pad_forecast <- data.table(forecast=pad_forecast, time=pad_time)
      forecast <- rbind(pad_forecast, forecast)
      
      tbl <- suppressMessages(full_join(truth, forecast))
      par(mgp=c(2.5,.8,0),mar=c(5,4,4,2)+0.1)
      plot(truth$time, truth$truth, type = "l", lwd=2, cex.lab=2, cex.axis=1.5,
           xlab="Time (minutes)", ylab="MAP", xlim=c(xmin, xmax), 
           ylim = range(c(0, 120)), panel.first=grid(lty=1))
      lines(tbl$time, tbl$forecast, type = "l", lwd=2, col = "red")
      Hmisc::minor.tick()
    } else {
      par(mgp=c(2.5,.8,0),mar=c(5,4,4,2)+0.1)
      plot(truth$time, truth$truth, type = "l", lwd=2, cex.lab=2, cex.axis=1.5,
           xlab="Time (minutes)", ylab="MAP", xlim=c(xmin, xmax), 
           ylim = range(c(0, 120)), panel.first=grid(lty=1))
      Hmisc::minor.tick()
    }
  }
}

make_forecast_movie <- function(movie_tbl, name){
  xmin <- min(movie_tbl$time, na.rm = T)
  xmax <- max(movie_tbl[!is.na(truth),"time"])
  times <- seq(xmin, xmax, by=2)
  animation::saveVideo(make_plots(movie_tbl, times, xmin, xmax), 
                       nmax=length(times), video.name=name, ani.width=2000, 
                       ani.height=1200, interval=0.5)
}
