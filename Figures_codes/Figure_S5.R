library(ggplot2)
library(reshape2)

p <- 10
ratio <- 1.5
n <- p * ratio
c0 <- sqrt(log(p)/n)

load(file = paste0("Simulation_results/SimulationS5_global_null_p_", p, "_n_", n, "_c0_", c0, ".RData"))
result <- as.data.frame(result)
names(result) <- c("p_0", "p_LRT", "p", "p'")
plot_data <- melt(result, variable.name = "Method", value.name = "pvalue")
plot_data$Method <- factor(plot_data$Method, levels = c("p_LRT", "p_0", "p'", "p"))
plot_data <- plot_data[plot_data$pvalue != 999,]

png(file = "Figures/FigureS5.png",
    width = 400, height = 440)
ggplot(data = plot_data) +
  geom_abline() + geom_qq(distribution = qunif, aes(sample = pvalue, color = Method)) + 
  xlab("Theoretical quantiles") + ylab("Empirical quantiles")  + 
  theme(text = element_text(size = 20)) + theme(plot.title = element_text(hjust = 0.5)) + labs(color = '') + 
  scale_color_discrete(labels = c(expression(paste(p[LRT], ' (6)')), expression(paste(p[est], ' (Section S11.1.2)')), expression(paste(p, "' (22)")), expression(paste(p, " (11)")))) + 
  theme(legend.position = "bottom") + theme(legend.text = element_text(size = 20), legend.text.align = 0) + guides(color = guide_legend(nrow = 2, byrow = TRUE))
dev.off()


