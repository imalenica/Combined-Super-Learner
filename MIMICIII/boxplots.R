library(data.table)
d <- read.csv("~/symphony/results/combined_performance_summary_by_id.csv")
d <- data.table(d)
better_d <- d[which(d$oSL_wts.broad < d$historical_screenerLasso_rf),]
remaining_d <- d[-which(id %in% better_d$id),]
# d_ranked <- remaining_d[order(historical_screenerLasso_rf,decreasing = T),]
d_ranked <- remaining_d[order(oSL_wts.broad,decreasing = F),]
d_summary <- rbind(d_ranked[1:(115-nrow(better_d)),], better_d)
head(sort(colMeans(d_summary[, c(14:64),with=F])),10)
head(sort(colMeans(d[d$sex == "F", c(14:64),with=F])))
head(sort(colMeans(d[d$sex == "M", c(14:64),with=F])))
head(sort(colMeans(d[d$age > 65, c(14:64),with=F])))
head(sort(colMeans(d[d$age <= 65, c(14:64),with=F])))
head(sort(colMeans(d[d$ethnicity == "NONWHITE", c(14:64),with=F])))
head(sort(colMeans(d[d$AHE5_65 > 0, c(14:64),with=F])))

dd <- d[,c(1:12),with=F]
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

