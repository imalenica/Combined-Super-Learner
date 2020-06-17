suppressMessages(library(here))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(data.table))
suppressMessages(library(gridExtra))
library(readr)
library(ggpubr)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

load(here("Simulations/ARIMA/Results/fit_mixture_stationary_v3_summary.Rdata"))
load(here("Simulations/ARIMA/Results/fit_stationary_v3_summary.Rdata"))
load(here("Simulations/ARIMA/Results/fit_W_offset_stationary_v3_summary.Rdata"))
load(here("Simulations/ARIMA/Results/fit_MAR_v3_summary.Rdata"))

#################################################################
### Plot majority weights (Historical vs. Individual learners)
#################################################################

sz <- 1.2

p0 <- ggplot() + 
  geom_line(aes(x = parse_number(row.names(res_nnls_convex_v0)),
                y = res_nnls_convex_v0$`Historical SL`,
                col = "blue"), size=sz) +
  geom_line(aes(x = parse_number(row.names(res_nnls_convex_v0)),
                y = res_nnls_convex_v0$`Individual SL`,
                col = "red"), size=sz) +
  xlab("") + ylab("Sum of ensemble weights") + 
  scale_color_discrete(name = "Learner", labels = c("Historical SL", "Individual SL")) +
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
 
sz2 <- 1.1

p0 <- ggplot() + 
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v0)),
                y = loss_nnls_convex_v0$`Pooled SL`,
                col = "Pooled SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v0)),
                y = loss_nnls_convex_v0$`Historical SL`,
                col = "Historical SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v0)),
                y = loss_nnls_convex_v0$`Individual SL`,
                col = "Individual SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v0)),
                y = loss_nnls_convex_v0$`Personalized SL`,
                col = "Personalized SL"), size=sz2) +
  scale_color_manual(values = c('Pooled SL' = 'blue', 'Historical SL' = 'red',
                                'Individual SL' = 'gray', 'Personalized SL' = 'black'),
                     name="Super Learner") +
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
                y = loss_nnls_convex_v1$`Pooled SL`,
                col = "Pooled SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v1)),
                y = loss_nnls_convex_v1$`Historical SL`,
                col = "Historical SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v1)),
                y = loss_nnls_convex_v1$`Individual SL`,
                col = "Individual SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v1)),
                y = loss_nnls_convex_v1$`Personalized SL`,
                col = "Personalized SL"), size=sz2) +
  scale_color_manual(values = c('Pooled SL' = 'blue', 'Historical SL' = 'red',
                                'Individual SL' = 'gray', 'Personalized SL' = 'black'),
                     name="Super Learner") +
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
                y = loss_nnls_convex_v2$`Pooled SL`,
                col = "Pooled SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v2)),
                y = loss_nnls_convex_v2$`Historical SL`,
                col = "Historical SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v2)),
                y = loss_nnls_convex_v2$`Individual SL`,
                col = "Individual SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v2)),
                y = loss_nnls_convex_v2$`Personalized SL`,
                col = "Personalized SL"), size=sz2) +
  scale_color_manual(values = c('Pooled SL' = 'blue', 'Historical SL' = 'red',
                                'Individual SL' = 'gray', 'Personalized SL' = 'black'),
                     name="Super Learner") +
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
                y = loss_nnls_convex_v3$`Pooled SL`,
                col = "Pooled SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v3)),
                y = loss_nnls_convex_v3$`Historical SL`,
                col = "Historical SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v3)),
                y = loss_nnls_convex_v3$`Individual SL`,
                col = "Individual SL"), size=sz2) +
  geom_line(aes(x = as.numeric(row.names(loss_nnls_convex_v3)),
                y = loss_nnls_convex_v3$`Personalized SL`,
                col = "Personalized SL"), size=sz2) +
  scale_color_manual(values = c('Pooled SL' = 'blue', 'Historical SL' = 'red',
                                'Individual SL' = 'gray', 'Personalized SL' = 'black'),
                     name="Super Learner") +
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
