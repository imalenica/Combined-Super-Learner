######################## video one forecast alongside truth ####################
make_plots <- function(movie_tbl, times, xmin, xmax, outcome_type){
  if(outcome_type == "continuous"){
    ylims <- range(c(0, 120))
    ylabs <- "MAP"
  } else if (outcome_type == "binomial"){
    ylims <- range(c(0, 1))
    ylabs <- "AHE"
  }
  for (i in 1:length(times)) {
    current_time <- times[i]
    truth <- movie_tbl[time <= current_time, c("truth", "time")]

    if(any(movie_tbl$forecast_time <= current_time, na.rm = T)){
      forecast <- movie_tbl[forecast_time <= current_time, c("forecast", "time")]
      pad_time <- seq(min(truth$time), min(forecast$time)-1)
      pad_forecast <- rep(NA, length(pad_time))
      pad_forecast <- data.table(forecast=pad_forecast, time=pad_time)
      forecast <- rbind(pad_forecast, forecast)

      tbl <- suppressMessages(full_join(truth, forecast))
      par(mgp=c(2.5,.8,0),mar=c(5,4,4,2)+0.1)
      plot(tbl$time, tbl$truth, type = "l", lwd=2, cex.lab=2, cex.axis=1.5,
           xlab="Time (minutes)", ylab=ylabs, xlim=c(xmin, xmax),
           ylim = ylims, panel.first=grid(lty=1))
      lines(tbl$time, tbl$forecast, type = "l", lwd=2, col = "red")
      Hmisc::minor.tick()
    } else {
      par(mgp=c(2.5,.8,0),mar=c(5,4,4,2)+0.1)
      plot(truth$time, truth$truth, type = "l", lwd=2, cex.lab=2, cex.axis=1.5,
           xlab="Time (minutes)", ylab=ylabs, xlim=c(xmin, xmax),
           ylim = ylims, panel.first=grid(lty=1))
      Hmisc::minor.tick()
    }
  }
}

make_forecast_movie <- function(movie_tbl, name, outcome_type){
  xmin <- min(movie_tbl$time, na.rm = T)
  xmax <- max(movie_tbl[!is.na(movie_tbl$truth),"time"])
  times <- seq(xmin, xmax, by=2)
  animation::saveVideo(
    make_plots(movie_tbl, times, xmin, xmax, outcome_type = outcome_type),
    nmax=length(times), video.name=name, ani.width=2000, ani.height=1200,
    interval=0.3, verbose=F, autoplay=F
  )
}

##################### video multiple forecasts alongside truth #################
make_multiforecast_plots <- function(movie_tbl, cols_list,
                                     truth_line_size = 0.05,
                                     forecast_line_size = 0.1, xmin, xmax,
                                     times, outcome_type){

  all <- c(
    movie_tbl$truth,
    unlist(lapply(cols_list, function(x){
      col_idx <- which(colnames(movie_tbl) == x$forecast)
      movie_tbl[, col_idx]
    }))
  )
  if(outcome_type == "continuous"){
    ymin <- round(min(all, na.rm = T) - 10, -1)
    ymax <- round(max(all, na.rm = T) + 10, -1)
    ylabs <- "MAP"
  } else if (outcome_type == "binomial"){
    ymin <- 0
    ymax <- 1
    ylabs <- "AHE"
  }

  plots <- list()
  for (i in 1:length(times)) {
    current_time <- times[i]
    truth <- movie_tbl[movie_tbl$time <= current_time, c("truth", "time")]

    tbls <- list()
    for(j in 1:length(cols_list)){
      cols <- cols_list[[j]]
      forecast_tbl <- movie_tbl[, c("time", cols$forecast_time, cols$forecast)]
      colnames(forecast_tbl)[2] <- "forecast_time"
      if(any(forecast_tbl$forecast_time <= current_time, na.rm = T)){
        forecast <- na.omit(forecast_tbl[forecast_tbl$forecast_time <= current_time, c(3,1)])
        pad_time <- seq(min(truth$time), min(forecast$time)-1)
        pad_forecast <- rep(NA, length(pad_time))
        pad_forecast <- data.table(pad_forecast, time=pad_time)
        colnames(pad_forecast)[1] <- cols$forecast
        forecast <- rbind(pad_forecast, forecast)
        tbls[[j]] <- suppressMessages(full_join(truth, forecast))
      } else {
        tbls[[j]] <- truth
      }
    }
    tbl <- suppressMessages(data.frame(tbls %>% reduce(full_join)))
    tbl <- melt(tbl, id.vars = "time")
    tbl$mysize <- rep(forecast_line_size, nrow(tbl))
    tbl$mysize[tbl$variable=="truth"] <- truth_line_size

    if(length(unique(tbl$variable)) == 1){
      plots[[i]] <- ggplot(tbl, aes(x=time, y=value, size=mysize, color=variable)) +
        geom_line() +
        xlim(xmin, xmax) +
        ylim(ymin, ymax) +
        scale_size(range = c(forecast_line_size, truth_line_size), guide="none") +
        labs(x="Time (minutes)", y=ylabs, title="", color="Type") +
        scale_color_manual(values="#666666")
      theme(legend.title = element_blank())
    } else {
      plots[[i]] <- ggplot(tbl, aes(x=time, y=value, size=mysize, color=variable)) +
        geom_line() +
        xlim(xmin, xmax) +
        ylim(ymin, ymax) +
        scale_size(range = c(forecast_line_size, truth_line_size), guide="none") +
        labs(x="Time (minutes)", y=ylabs, title="", color="Type") +
        scale_color_manual(values=c("#666666", "#A6761D", "#E6AB02", "#66A61E",
                                    "#E7298A", "#7570B3", "#D95F02", "#1B9E77")) +
        theme(legend.title = element_blank())
    }
  }
  return(plots)
}


make_multiforecast_movie_mp4 <- function(movie_tbl, name, cols_list,
                                         truth_line_size = 0.5,
                                         forecast_line_size = 0.3,
                                         outcome_type, savio = TRUE){
  xmin <- min(movie_tbl$time, na.rm = T)
  xmax <- max(movie_tbl[!is.na(movie_tbl$truth),"time"])
  times <- seq(xmin, xmax, by = 2)
  if(savio){
    ani.options(ffmpeg = '/global/home/users/rachelvphillips/ffmpeg-git-20210501-amd64-static/ffmpeg')
  }
  animation::saveVideo(expr = {
    plots <- make_multiforecast_plots(movie_tbl, cols_list,
                                      truth_line_size = truth_line_size,
                                      forecast_line_size = forecast_line_size,
                                      xmin, xmax, times,
                                      outcome_type = outcome_type
    )
    png(ani.options("img.fmt"), width = 6, height = 4, units = "in", res = 700)
    for(i in 1:length(plots)){
      print(plots[[i]])
    }
    dev.off()
  }, use.dev = FALSE, nmax = length(times), video.name = name, ani.width = 2000,
  ani.height = 1200, interval = 0.3, autobrowse = FALSE, verbose = FALSE)
}

