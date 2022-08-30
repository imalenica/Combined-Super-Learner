d <- read.csv("~/Downloads/performance_summary-2.csv")
d <- data.table(d)
dd <- t(data.table(t(colMeans(d[, names(d) := lapply(.SD, as.numeric)]))))[-1,]
write.csv(dd, "~/Downloads/mean_summary.csv")
dd <- df[,c(1:12),with=F]
dd <- melt(dd, id.vars = "id")

plot_osl <- ggplot(data = dd, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() +
  labs(x="", y="MSE", title="MSE by Learner") +
  coord_flip() +
  theme(legend.position="none")

dd <- df[,c(1,13:26),with=F]
dd <- melt(dd, id.vars = "id")
plot_sl <- ggplot(data = dd, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() +
  labs(x="", y="MSE", title="MSE by Learner") +
  coord_flip() +
  theme(legend.position="none")

dd <- df[,c(1,27:42),with=F]
dd <- melt(dd, id.vars = "id")
plot_historical <- ggplot(data = dd, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() +
  labs(x="", y="MSE", title="MSE by Learner") +
  coord_flip() +
  theme(legend.position="none")

dd <- df[,c(1,43:52),with=F]
dd <- melt(dd, id.vars = "id")
plot_ind <- ggplot(data = dd, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() +
  labs(x="", y="MSE", title="MSE by Learner") +
  coord_flip() +
  theme(legend.position="none")


dd <- melt(df, id.vars = "id")
plot_all <- ggplot(data = dd, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() +
  labs(x="", y="MSE", title="MSE by Learner") +
  coord_flip() +
  theme(legend.position="none")

