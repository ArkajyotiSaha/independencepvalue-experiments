library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)

power <- function(x){mean(x < 0.05)}

p <- 100
plot_data <- data.frame()
for(ratio in c(1.1, 1.5, 2)){
  n <- p * ratio
  load(file=paste0("Simulation_results/alternative_p_",p,"_n_",n,".RData"))
  plot_data <- rbind(plot_data, result)
}
names(plot_data) <- c("Selective", "Classical", "Effect", "p", "Select.Prob")
plot_data$Method <- factor(rep(c("n = 1.1p","n = 1.5p","n = 2p"), each = 1000000), levels = c("n = 1.1p", "n = 1.5p", "n = 2p"))
plot_data <- plot_data[plot_data[,1]!= 999,]
plot_data_melt <- melt(plot_data, na.rm = FALSE, id = c("Effect", "p", "Select.Prob", "Method"), value.name = "pvalue")
names(plot_data_melt)[5] <- "Inference"
plot_data_melt$Inference <- factor(plot_data_melt$Inference, levels = c("Classical", "Selective"))
plot_data_a <- plot_data_melt[plot_data_melt$p >= 95,]
plot_data_a  <- as.data.frame(plot_data_a  %>% 
                                  group_by(Method, Inference, Bin.Select.Prob = cut(Select.Prob, breaks = seq(0, 1, by = 0.1))) %>% 
                                  filter(length(pvalue) >= 100) %>%
                                  ungroup())

p1 <- ggplot(plot_data_a, aes(Select.Prob, pvalue, color=Method, linetype = Inference)) +
  stat_summary_bin(fun = "power", geom = "line", orientation = 'x', breaks = seq(0, 1, 0.1)) + 
  stat_summary_bin(fun = "power", geom = "point", orientation = 'x',breaks = seq(0, 1, 0.1)) +
  xlab("Probability of the selection event") + ylab("Power at level \u03b1 = 0.05") + labs(color='') + labs(linetype='') +
  ggtitle(expression(paste("(a) ",r,"(",  hat(P),")" ,' \u2264 5'))) + theme(text = element_text(size = 20)) +
  theme(plot.title = element_text(hjust = 0.5))


plot_data_b <- plot_data_melt[plot_data_melt$p < 95,]
plot_data_b  <- as.data.frame(plot_data_b  %>% 
                                group_by(Method, Inference, Bin.Select.Prob = cut(Select.Prob, breaks = seq(0, 1, by = 0.1))) %>% 
                                filter(length(pvalue) >= 100) %>%
                                ungroup())


p2 <- ggplot(plot_data_b, aes(Select.Prob, pvalue, linetype = Inference, color=Method)) +
  stat_summary_bin(fun = "power", geom = "line", orientation = 'x', breaks = seq(0, 1, 0.1)) + 
  stat_summary_bin(fun = "power", geom = "point", orientation = 'x',breaks = seq(0, 1, 0.1)) +
  xlab("Probability of the selection event") + ylab("Power at level \u03b1 = 0.05") + labs(color='') + labs(linetype='') +
  ggtitle(expression(paste("(b) ",r,"(",  hat(P),")" ,' > 5')))  + theme(text = element_text(size = 20)) +
  theme(plot.title = element_text(hjust = 0.5))

png(file="Figures/Figure3.png",
    width=800, height=420)
ggpubr::ggarrange(p1, p2, ncol=2, common.legend=TRUE, legend="bottom")
dev.off()

