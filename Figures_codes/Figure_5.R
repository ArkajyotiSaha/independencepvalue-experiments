library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)

power <- function(x){mean(x < 0.05)}

p <- 100
network <- 3 # 3 for E.coli & 4 for S.cerevisiae
plot_information <- list()
eta_vector <- c(0, 0.5, 3)
for(index in 1:3){
  eta <- eta_vector[index]
  plot_data <- data.frame()
  for(ratio in c(1.2, 1.4, 1.5, 2)){
    n <- p * ratio
    load(file=paste0("Real_data_results/network_", network, "p_", p, "_n_", n, "_eta_", eta, ".RData"))
    plot_data <- rbind(plot_data, result)
  }
  plot_data$Method <- factor(rep(c("n = 1.2p", "n = 1.4p","n = 1.5p","n = 2p"), each = 10000), levels = c("n = 1.2p", "n = 1.4p","n = 1.5p","n = 2p"))
  names(plot_data)[1:3] <- c("Selective", "Classical", "Effect")
  plot_data <- plot_data[plot_data[,1]!= 999,]
  plot_data_melt <- melt(plot_data, na.rm = FALSE, id = c("Effect", "Method"), value.name = "pvalue")
  names(plot_data_melt)[3] <- "Inference"
  plot_data_melt$Inference <- factor(plot_data_melt$Inference, levels = c("Classical", "Selective"))
  plot_data_melt$Effect <- 1 - plot_data_melt$Effect
  plot_data_melt <- as.data.frame(plot_data_melt %>% 
                                   group_by(Method, Inference, Bin.Effect = cut(Effect, breaks = seq(0, 1, by = 0.1))) %>% 
                                   filter(length(pvalue) >= 100) %>%
                                   ungroup())
  plot_information[[index]] <- ggplot(plot_data_melt, aes(Effect, pvalue, color=Method, linetype = Inference)) +
    stat_summary_bin(fun = "power", geom = "line", orientation = 'x', breaks = seq(0, 1, 0.1)) + 
    stat_summary_bin(fun = "power", geom = "point", orientation = 'x',breaks = seq(0, 1, 0.1)) +
    xlab(expression(paste("Effect Size (", hat('\u0394'), ")"))) + 
    ylab("Power at level \u03b1 = 0.05") + labs(color='') + labs(linetype='') +
    theme(text = element_text(size = 30)) + theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1)
}

p1 <- plot_information[[1]] + ggtitle(paste0("(a) \u03b7 = ", eta_vector[1]))  
p2 <- plot_information[[2]] + ggtitle(paste0("(b) \u03b7 = ", eta_vector[2])) 
p3 <- plot_information[[3]] + ggtitle(paste0("(c) \u03b7 = ", eta_vector[3]))

png(file="Figures/Figure5.png",
    width=1800, height=600)
ggpubr::ggarrange(p1, p2, p3, ncol=3, common.legend=TRUE, legend="bottom")
dev.off()

