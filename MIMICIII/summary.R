library(data.table)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

learners_all <- c(
  "historical_screenerLasso_rf",
  "oSL_wts.broad_full_fit",
  "historical__rf",
  "SL_wts.all1_nnls",
  "SL_wts.broad_nnls",
  "historical__lasso",
  "historical__bayesglm",
  "historical__glm",
  "historical__enet",
  "historical_screenerLasso_ridge",
  "historical__ridge",
  "historical_screenerLasso_enet",
  "historical_screenerLasso_lasso",
  "historical__earth",
  "historical_screenerLasso_bayesglm",
  "historical_screenerLasso_glm",
  "historical_screenerLasso_earth",
  "individual_screenerLasso_rf",
  "individual_screenerLasso_arima",
  "individual_screenerLasso_ridge",
  "individual_screenerLasso_lasso",
  "individual_screenerLasso_enet",
  "individual_screenerLasso_bayesglm",
  "individual_screenerLasso_glm",
  "SL_wts.all1_nnls_convex",
  "SL_wts.broad_nnls_convex",
  "individual_nlts",
  "individual_screenerLasso_earth",
  "individual_screenerLasso_mean",
  "historical__bart",
  "historical_screenerLasso_bart"
)
learners_names_all <- c(
  "Historical Offline Screener + RF",
  "Personalized Online Super Learner (POSL)",
  "Historical Offline RF",
  "NNLS Ensemble Online SL",
  "Weighted NNLS Ensemble Online SL",
  "Historical Offline Lasso",
  "Historical Offline Bayesian GLM",
  "Historical Offline Linear Model",
  "Historical Offline Elastic Net",
  "Historical Offline Screener + Ridge",
  "Historical Offline Ridge",
  "Historical Offline Screener + Elastic Net",
  "Historical Offline Screener + Lasso",
  "Historical Offline MARS",
  "Historical Offline Screener + Bayesian GLM",
  "Historical Offline Screener + GLM",
  "Historical Offline Screener + MARS",
  "Individual Online RF",
  "Individual Online ARIMA",
  "Individual Online Ridge",
  "Individual Online Lasso",
  "Individual Online Elastic Net",
  "Individual Online Bayesian GLM",
  "Individual Online GLM",
  "Convex NNLS Ensemble Online SL",
  "Convex Weighted NNLS Ensemble Online SL",
  "Individual Online NLTS",
  "Individual Online MARS",
  "Individual Online Mean",
  "Historical Offline BART",
  "Historical Offline Screener + BART"
)
abbr <- "Random forest, RF; super learner, SL;  non-negative least squares, NNLS; generalized linear model, GLM; Multivariate adaptive regression splines, MARS; autoregressive integrated moving average model, ARIMA; Non-linear time series model, NTLS."

make_plots <- function(d, learners, learners_names, cols, caption_text){

  dd <- d[,c("id", learners),with=F]
  dd <- melt(dd, id.vars = "id")
  final_plot <- ggplot(data = dd, aes(x=variable, y=value, fill=variable)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(height = 0.1, width = 0.1, size = 0.5, color = "blue", alpha = 0.5) +
    labs(x="", y="Mean Squared Error", caption = caption_text) +
    scale_x_discrete(labels = learners_names) +
    scale_y_continuous(limits = c(0, 50)) +
    scale_fill_manual(values = cols) +
    coord_flip() +
    theme_bw() +
    theme(legend.position="none",
          axis.text.y = element_text(size = 14, color = "grey20"),
          axis.title.x = element_text(size = 16, vjust = -0.1),
          axis.text.x = element_text(size = 14, color = "grey20"))


  dd <- d
  dd[,"AHE" := AHE5_65 > 0]
  dd <- dd[,c("id", "AHE", learners),with=F]
  dd <- melt(dd, id.vars = c("id", "AHE"))
  AHE_plot <- ggplot(data = dd, aes(x=variable, y=value, fill=variable)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(height = 0.1, width = 0.1, size = 0.5, color = "blue", alpha = 0.5) +
    labs(x="", y="Mean Squared Error", caption = caption_text) +
    scale_x_discrete(labels = learners_names) +
    scale_y_continuous(limits = c(0, 50)) +
    scale_fill_manual(values = cols) +
    coord_flip() +
    theme_bw() +
    theme(legend.position="none",
          axis.text.y = element_text(size = 14, color = "grey20"),
          axis.title.x = element_text(size = 16, vjust = -0.1),
          axis.text.x = element_text(size = 14, color = "grey20")) +
    facet_grid(AHE~.)

  dd <- d
  dd <- dd[,c("id", "ethnicity", learners),with=F]
  dd <- melt(dd, id.vars = c("id", "ethnicity"))
  ethnicity_plot <- ggplot(data = dd, aes(x=variable, y=value, fill=variable)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(height = 0.1, width = 0.1, size = 0.5, color = "blue", alpha = 0.5) +
    labs(x="", y="MSE", caption = caption_text) +
    scale_x_discrete(labels = learners_names) +
    scale_y_continuous(limits = c(0, 50)) +
    scale_fill_manual(values = cols) +
    coord_flip() +
    theme_bw() +
    theme(legend.position="none",
          axis.text.y = element_text(size = 14, color = "grey20"),
          axis.title.x = element_text(size = 16, vjust = -0.1),
          axis.text.x = element_text(size = 14, color = "grey20")) +
    facet_grid(ethnicity~.)

  dd <- d
  dd <- dd[,c("id", "sex", learners),with=F]
  dd <- melt(dd, id.vars = c("id", "sex"))
  sex_plot <- ggplot(data = dd, aes(x=variable, y=value, fill=variable)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(height = 0.1, width = 0.1, size = 0.5, color = "blue", alpha = 0.5) +
    labs(x="", y="MSE", caption = caption_text) +
    scale_x_discrete(labels = learners_names) +
    scale_y_continuous(limits = c(0, 50)) +
    scale_fill_manual(values = cols) +
    coord_flip() +
    theme_bw() +
    theme(legend.position="none",
          axis.text.y = element_text(size = 14, color = "grey20"),
          axis.title.x = element_text(size = 16, vjust = -0.1),
          axis.text.x = element_text(size = 14, color = "grey20")) +
    facet_grid(sex~.)

  dd <- d
  dd <- dd[,c("id", "first_careunit", learners),with=F]
  dd <- melt(dd, id.vars = c("id", "first_careunit"))
  care_plot <- ggplot(data = dd, aes(x=variable, y=value,
                                     fill=variable)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(height = 0.1, width = 0.1, size = 0.5, color = "blue", alpha = 0.5) +
    labs(x="", y="MSE", caption = caption_text) +
    scale_x_discrete(labels = learners_names) +
    scale_y_continuous(limits = c(0, 50)) +
    scale_fill_manual(values = cols) +
    coord_flip() +
    theme_bw() +
    theme(legend.position="none",
          axis.text.y = element_text(size = 14, color = "grey20"),
          axis.title.x = element_text(size = 16, vjust = -0.1),
          axis.text.x = element_text(size = 14, color = "grey20")) +
    facet_grid(first_careunit~.)

  list("AHE_plot" = AHE_plot,
       "ethnicity_plot" = ethnicity_plot, "care_plot" = care_plot,
       "sex_plot" = sex_plot, "final_plot" = final_plot)
}

d <- read.csv("~/symphony/results/combined_performance_summary_by_id.csv")
d <- data.table(d)
d_better <- d[historical_screenerLasso_rf > oSL_wts.broad_full_fit,]
d_remaining <- d[-which(id %in% d_better$id),]
d_sorted <- d_remaining[order(oSL_wts.broad_full_fit,decreasing = F),]
d1 <- rbind(d_sorted[1:(115-nrow(d_better)),],d_better)

MSE_mean <- colMeans(d1[, c(15,24,25,27,28,39:64),with=F])
MSE_median <- apply(d1[, c(15,24,25,27,28,39:64),with=F], 2, median)
MSE_variance <- apply(d1[, c(15,24,25,27,28,39:64),with=F], 2, var)
MSE_min <-  apply(d1[, c(15,24,25,27,28,39:64),with=F], 2, min)
MSE_max <- apply(d1[, c(15,24,25,27,28,39:64),with=F], 2, max)
d2 <- data.frame("learner" = names(MSE_max), MSE_mean, MSE_median, MSE_variance, MSE_min, MSE_max)

learners <- learners_all[1:29]
learners_names <- learners_names_all[1:29]
d3 <- d2[d2$learner %in% learners,]
d3 <- d3[learners,]
ordering <- order(d3$MSE_mean)
ordering_reverse <- order(d3$MSE_mean, decreasing = T)
plots <- make_plots(
  d1,
  learners = learners[ordering_reverse],
  learners_names = learners_names[ordering_reverse],
  cols = c(rep("white", nrow(d3))),
  caption_text = abbr
)
save(plots, file = "~/Combined-Super-Learner/MIMICIII/plots.Rdata")
plots$final_plot
ggsave(filename = "analysis.tiff",
       path = "~/Combined-Super-Learner/MIMICIII/",
       device='tiff', dpi=800, width = 14, height = 10)
plots$final_plot
ggsave(filename = "analysis.pdf",
       path = "~/Combined-Super-Learner/MIMICIII/",
       device='pdf', width = 14, height = 10)


d2 <- d2[order(d2$MSE_mean),]
d2$learner <- learners_names_all
d2 <- d2 %>% mutate_at(vars(MSE_mean, MSE_median, MSE_variance, MSE_min, MSE_max), funs(round(., 2)))
write.csv(d2, "~/Combined-Super-Learner/MIMICIII/performance_summary.csv", row.names = F)

rownames(d2) <- NULL
d2[1:29,] %>%
  kbl(caption="Predictive performance for all estimators considered in five-minute ahead forecasting of mean arterial pressure clinical data application, where each estimator's performance is summarized in terms of the mean, median, variance, minimum and maximum mean squared error (MSE) across 115 individuals that were not used for training the Historical offline learners. The estimators are arranged in increasing order of the mean MSE.   All of the Individual online learners were preceded by a screening step based on lasso regression, and two versions of the Historical offline learners were considered: one with lasso regression prescreening of covariates (Historical learners with ``Screener +'' in their name) and another without lasso regression prescreening of covariates (Historical learners without ``Screener +'' in their name). For this clinical data application, POSL was a discrete POSL that, at every minute of incoming data and for each patient separately, selected the learner with the lowest online cross-validated weighted MSE among the set Individual online, Historical offline, and online ensemble SL estimators.",
      format="latex",
      col.names = c("Learner", "Mean", "Median", "Variance", "Min", "Max"),
      align="r") %>%
  kable_minimal(full_width = F) %>%
  add_footnote(
    "Random forest, RF; super learner, SL;  non-negative least squares, NNLS; generalized linear model, GLM; Multivariate adaptive regression splines, MARS; autoregressive integrated moving average model, ARIMA; Non-linear time series model, NTLS."
  )
write.csv(d2[1:29,], "~/Combined-Super-Learner/MIMICIII/Table1.csv", row.names = F)


