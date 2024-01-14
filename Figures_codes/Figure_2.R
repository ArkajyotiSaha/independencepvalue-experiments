library(ggplot2)
library(ggpubr)


p <- 100
plot_data <- data.frame()
for(ratio in c(1.1, 1.5, 2)){
  n <- p * ratio
  c0 <- sqrt(log(p)/n)
  load(file = paste0("Simulation_results/global_null_p_", p, "_n_", n, "_c0_", c0, ".RData"))
  plot_data <- rbind(plot_data, result)
}

names(plot_data) <- c("Selective", "Classical")
plot_data$Method <- rep(c("n = 1.1p","n = 1.5p","n = 2p"), each = 100000)
plot_data$Method <- factor(plot_data$Method, levels = c("n = 1.1p", "n = 1.5p", "n = 2p"))
plot_data <- plot_data[plot_data[,1] != 999,]

base_plot <- ggplot(data = plot_data) +
  xlab("Theoretical quantiles")  +
  ylab("Empirical quantiles") +
  geom_abline()

p1 <- base_plot + geom_qq(distribution = qunif,
              aes(sample = Classical, color = Method)) + 
  ggtitle(paste0("(a) Classical Inference")) + 
  theme(text = element_text(size = 20)) + 
  theme(plot.title = element_text(hjust = 0.5)) + labs(color = '')

p2 <- base_plot + geom_qq(distribution = qunif,
                          aes(sample = Selective, color = Method)) +
  ggtitle(paste0("(b) Selective Inference")) + 
  theme(text = element_text(size = 20)) + 
  theme(plot.title = element_text(hjust = 0.5)) + labs(color = '')

png(file = "Figures/Figure2.png",
    width = 800, height = 420)
ggpubr::ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "bottom")
dev.off()

