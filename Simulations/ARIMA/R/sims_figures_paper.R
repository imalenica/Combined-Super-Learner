suppressMessages(library(here))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(data.table))
suppressMessages(library(gridExtra))
library(RColorBrewer)
library(readr)
library(ggpubr)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

load(here("Simulations/ARIMA/Results/fit_mixture_stationary_v4_summary.Rdata"))
load(here("Simulations/ARIMA/Results/fit_stationary_v4_summary.Rdata"))
load(here("Simulations/ARIMA/Results/fit_W_offset_stationary_v4_summary.Rdata"))
load(here("Simulations/ARIMA/Results/fit_MAR_v4_summary.Rdata"))

#################################################################
### Plot majority weights (Historical vs. Individual learners)
#################################################################

sz <- 1.1
p0 <- ggplot() + 
  geom_line(aes(x = parse_number(row.names(res_nnls_convex_v0)),
                y = res_nnls_convex_v0$`Historical SL`,
                col = "blue"), size=sz) +
  geom_line(aes(x = parse_number(row.names(res_nnls_convex_v0)),
                y = res_nnls_convex_v0$`Individual SL`,
                col = "red"), size=sz) +
  xlab("") + ylab("Sum of ensemble weights") + 
  scale_color_brewer(palette = "Dark2", name="Learner", labels = c("Historical SL", "Individual SL")) +
  #scale_color_discrete(name = "Learner", labels = c("Historical SL", "Individual SL")) +
  ggtitle("(a) ARIMA process") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92"))

p1 <- ggplot() + 
  geom_line(aes(x = parse_number(row.names(res_nnls_convex_v1)),
                y = res_nnls_convex_v1$`Historical SL`,
                col = "blue"), size=sz) +
  geom_line(aes(x = parse_number(row.names(res_nnls_convex_v1)),
                y = res_nnls_convex_v1$`Individual SL`,
                col = "red"), size=sz) +
  xlab(" ") + ylab(" ") +
  scale_color_brewer(palette = "Dark2", name="Learner", labels = c("Historical SL", "Individual SL")) +
  ggtitle("(b) ARIMA process with X-dependent offset") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92"))  

p2 <- ggplot() + 
  geom_line(aes(x = parse_number(row.names(res_nnls_convex_v2)),
                y = res_nnls_convex_v2$`Historical SL`,
                col = "blue"), size=sz) +
  geom_line(aes(x = parse_number(row.names(res_nnls_convex_v2)),
                y = res_nnls_convex_v2$`Individual SL`,
                col = "red"), size=sz) +
  xlab("Time") + ylab("Sum of ensemble weights") +
  scale_color_brewer(palette = "Dark2", name="Learner", labels = c("Historical SL", "Individual SL")) +
  ggtitle("(c) Interrupted ARIMA process") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92"))

p3 <- ggplot() + 
  geom_line(aes(x = parse_number(row.names(res_nnls_convex_v3)),
                y = res_nnls_convex_v3$`Historical SL`,
                col = "blue"), size=sz) +
  geom_line(aes(x = parse_number(row.names(res_nnls_convex_v3)),
                y = res_nnls_convex_v3$`Individual SL`,
                col = "red"), size=sz) +
  xlab("Time") + ylab(" ") +
  scale_color_brewer(palette = "Dark2", name="Learner", labels = c("Historical SL", "Individual SL")) +
  ggtitle("(d) Mixture of Gaussian AR models") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92")) 

ggarrange(p0, p1, p2, p3, nrow=2, ncol=2, common.legend = TRUE, legend="right", font.label="bold")

#################################################################
###   Plot loss for different Super Learners
#################################################################

cbp <- c("#000000", "#E69F00", "#4B0082", "#009E73", "khaki")

sz2 <- 1.1
p0 <- ggplot() + 
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v0)),
                y = loss_nnls_convex_v0$`V-fold SL`,
                col = "V-fold SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v0)),
                y = loss_nnls_convex_v0$`Historical SL`,
                col = "Historical SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v0)),
                y = loss_nnls_convex_v0$`Individual SL`,
                col = "Individual SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v0)),
                y = loss_nnls_convex_v0$`Online SL`,
                col = "Online SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v0)),
                y = loss_nnls_convex_v0$`Personalized SL`,
                col = "Personalized SL"), size=sz2) +
  scale_color_manual(values = cbp, name="Super Learner") +
  ylab("MSE") + xlab(" ") + 
  labs(title = "(a) ARIMA process") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92"))

p1 <- ggplot() + 
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v1)),
                y = loss_nnls_convex_v1$`V-fold SL`,
                col = "V-fold SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v1)),
                y = loss_nnls_convex_v1$`Historical SL`,
                col = "Historical SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v1)),
                y = loss_nnls_convex_v1$`Individual SL`,
                col = "Individual SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v1)),
                y = loss_nnls_convex_v1$`Online SL`,
                col = "Online SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v1)),
                y = loss_nnls_convex_v1$`Personalized SL`,
                col = "Personalized SL"), size=sz2) +
  scale_color_manual(values = cbp, name="Super Learner") +
  ylab(" ") + xlab(" ") + 
  labs(title = "(b) ARIMA process with X-dependent offset") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92"))

p2 <- ggplot() + 
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v2)),
                y = loss_nnls_convex_v2$`V-fold SL`,
                col = "V-fold SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v2)),
                y = loss_nnls_convex_v2$`Historical SL`,
                col = "Historical SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v2)),
                y = loss_nnls_convex_v2$`Individual SL`,
                col = "Individual SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v2)),
                y = loss_nnls_convex_v2$`Online SL`,
                col = "Online SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v2)),
                y = loss_nnls_convex_v2$`Personalized SL`,
                col = "Personalized SL"), size=sz2) +
  scale_color_manual(values = cbp, name="Super Learner") +
  xlab("Time") + ylab("MSE") +
  labs(title = "(c) Interrupted ARIMA process") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92"))

p3 <- ggplot() + 
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v3)),
                y = loss_nnls_convex_v3$`V-fold SL`,
                col = "V-fold SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v3)),
                y = loss_nnls_convex_v3$`Historical SL`,
                col = "Historical SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v3)),
                y = loss_nnls_convex_v3$`Individual SL`,
                col = "Individual SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v3)),
                y = loss_nnls_convex_v3$`Online SL`,
                col = "Online SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v3)),
                y = loss_nnls_convex_v3$`Personalized SL`,
                col = "Personalized SL"), size=sz2) +
  scale_color_manual(values = cbp, name="Super Learner") +
  xlab("Time") + ylab(" ") +
  labs(title = "(d) Mixture of Gaussian AR models") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92"))

ggarrange(p0, p1, p2, p3, nrow=2, ncol=2, common.legend = TRUE, legend="right", font.label="bold")

#################################################################
###   Plot forecast predictions vs. truth for one iteration
#################################################################

load(here("Simulations/ARIMA/Results/fit_preds_extra_v4.Rdata"))

forecast_sim0 <- do.call(rbind,res_sim0$forecast)
forecast_sim1 <- do.call(rbind,res_sim1$forecast)
forecast_sim2 <- do.call(rbind,res_sim2$forecast)
forecast_sim3 <- do.call(rbind,res_sim3$forecast)

sz <- 1.1
p0 <- ggplot() + 
  geom_line(aes(x = as.numeric(row.names(forecast_sim0)),
                y = round(forecast_sim0$truth,0), col = "Truth"), size=sz) +
  geom_line(aes(x = as.numeric(row.names(forecast_sim0)),
                y = round(forecast_sim0$persSL,0), col = "POSL"), size=sz) +
  xlab("") + ylab("Outcome") + 
  scale_color_brewer(palette = "Dark2", name="Outcome", 
                     labels = c("Truth", "POSL")) +
  ggtitle("(a) ARIMA process") +
  ylim(60,90)+
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92"))

p1 <- ggplot() + 
  geom_line(aes(x = as.numeric(row.names(forecast_sim1)),
                y = round(forecast_sim1$truth,0), col = "Truth"), size=sz) +
  geom_line(aes(x = as.numeric(row.names(forecast_sim1)),
                y = round(forecast_sim1$persSL,0), col = "POSL"), size=sz) +
  xlab("") + ylab("Outcome") + 
  scale_color_brewer(palette = "Dark2", name="Outcome", 
                     labels = c("POSL", "Truth")) +
  ggtitle("(b) ARIMA process with X-dependent offset") +
  ylim(60,90)+
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92"))  

p2 <- ggplot() + 
  geom_line(aes(x = as.numeric(row.names(forecast_sim2)),
                y = round(forecast_sim2$truth,0), col = "Truth"), size=sz) +
  geom_line(aes(x = as.numeric(row.names(forecast_sim2)),
                y = round(forecast_sim2$persSL,0), col = "POSL"), size=sz) +
  xlab("") + ylab("Outcome") + 
  scale_color_brewer(palette = "Dark2", name="Outcome", 
                     labels = c("POSL", "Truth")) +
  ylim(60,90)+
  ggtitle("(c) Interrupted ARIMA process") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92"))

p3 <- ggplot() + 
  geom_line(aes(x = as.numeric(row.names(forecast_sim3)),
                y = round(forecast_sim3$truth,0), col = "Truth"), size=sz) +
  geom_line(aes(x = as.numeric(row.names(forecast_sim3)),
                y = round(forecast_sim3$persSL,0), col = "POSL"), size=sz) +
  xlab("") + ylab("Outcome") + 
  scale_color_brewer(palette = "Dark2", name="Outcome", 
                     labels = c("POSL", "Truth")) +
  ylim(60,90)+
  ggtitle("(d) Mixture of Gaussian AR models") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92")) 

ggarrange(p0, p1, p2, p3, nrow=2, ncol=2, common.legend = TRUE, legend="right", font.label="bold")

#################################################################
###   Plot CV predictions vs. truth for one iteration
#################################################################

load(here("Simulations/ARIMA/Results/fit_preds_extra_v4.Rdata"))

sz <- 1.1
p0 <- ggplot() + 
  geom_line(aes(x = as.numeric(row.names(preds_sim0)),
                y = round(preds_sim0$truth,0), col = "Truth"), size=sz) +
  geom_line(aes(x = as.numeric(row.names(preds_sim0)),
                y = round(preds_sim0$nnls_convexSL,2), col = "POSL"), size=sz) +
  xlab("") + ylab("Outcome") + 
  scale_color_brewer(palette = "Dark2", name="Outcome", 
                     labels = c("Truth", "POSL")) +
  ggtitle("(a) ARIMA process") +
  ylim(60,90)+
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92"))

p1 <- ggplot() + 
  geom_line(aes(x = as.numeric(row.names(preds_sim1)),
                y = round(preds_sim1$truth,0), col = "Truth"), size=sz) +
  geom_line(aes(x = as.numeric(row.names(preds_sim1)),
                y = round(preds_sim1$nnls_convexSL,2), col = "POSL"), size=sz) +
  xlab("") + ylab("Outcome") + 
  scale_color_brewer(palette = "Dark2", name="Outcome", 
                     labels = c("POSL", "Truth")) +
  ggtitle("(b) ARIMA process with X-dependent offset") +
  ylim(60,90)+
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92"))  

p2 <- ggplot() + 
  geom_line(aes(x = as.numeric(row.names(preds_sim2)),
                y = round(preds_sim2$truth,0), col = "Truth"), size=sz) +
  geom_line(aes(x = as.numeric(row.names(preds_sim2)),
                y = round(preds_sim2$nnls_convexSL,2), col = "POSL"), size=sz) +
  xlab("") + ylab("Outcome") + 
  scale_color_brewer(palette = "Dark2", name="Outcome", 
                     labels = c("POSL", "Truth")) +
  ylim(60,90)+
  ggtitle("(c) Interrupted ARIMA process") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92"))

p3 <- ggplot() + 
  geom_line(aes(x = as.numeric(row.names(preds_sim3)),
                y = round(preds_sim3$truth,0), col = "Truth"), size=sz) +
  geom_line(aes(x = as.numeric(row.names(preds_sim3)),
                y = round(preds_sim3$nnls_convexSL,2), col = "POSL"), size=sz) +
  xlab("") + ylab("Outcome") + 
  scale_color_brewer(palette = "Dark2", name="Outcome", 
                     labels = c("POSL", "Truth")) +
  ylim(60,90)+
  ggtitle("(d) Mixture of Gaussian AR models") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(face="bold", angle=90, hjust=1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(color="black", size=16, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(colour="black", size=16, face="bold"), 
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(colour = "grey92")) 

ggarrange(p0, p1, p2, p3, nrow=2, ncol=2, common.legend = TRUE, legend="right", font.label="bold")





