library(ggplot2)
make_weights_plot <- function(weights_control, max_time = 240,
                              scale_x_breaks = c(0, 30, 60, 90, 120, 150, 180, 210, 240),
                              scale_x_labels = c(0, 30, "","","","", 180, "", 240),
                              xlab = "Lag (hours) from current time"){
  times <- seq(1:max_time)
  weights <- rep(1, length(times))
  if (!is.null(weights_control$window)) {
    window <- max(times) - weights_control$window
    weights <- weights * ifelse(times <= window, 0, 1)
  }
  if (!is.null(weights_control$rate_decay)) {
    lags <- max(times) - times
    if (!is.null(weights_control$delay_decay)) {
      lags_delayed <- lags - weights_control$delay_decay
      lags <- ifelse(lags_delayed < 0, 0, lags_delayed)
    }
    weights <- weights * (1 - weights_control$rate_decay)^lags
  }
  df <- data.frame(weights, time = abs(times-max_time))
  ggplot(data = df, aes(x=time, y=weights)) +
    geom_line() +
    theme_bw() +
    labs(x=xlab, y="Weight") +
    scale_x_continuous(breaks = scale_x_breaks, labels = scale_x_labels)
}
pdf("~/Combined-Super-Learner/MIMICIII/weights.pdf", width = 6, height = 4)
make_weights_plot(
  weights_control = list(window = NULL, delay_decay = 90, rate_decay = 0.001),
  max_time = 1440,
  scale_x_breaks = seq(from = 0, to = 1440, by = 60),
  scale_x_labels = seq(0:24)-1
)
dev.off()

make_weights_plot(
  weights_control = list(window = NULL, delay_decay = 90, rate_decay = 0.001),
  max_time = 1440,
  scale_x_breaks = seq(from = 0, to = 1440, by = 60),
  scale_x_labels = seq(0:24)-1
)
ggsave("~/Combined-Super-Learner/MIMICIII/weights.tiff", width=6, height=4, dpi=800)
