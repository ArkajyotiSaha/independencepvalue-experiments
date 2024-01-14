library(ggplot2)
library(reshape2)

p <- 10
ratio <- 1.5
n <- p * ratio
c0 <- sqrt(log(p)/n)

load(file = paste0("Simulation_results/histogram_S4_global_null_p_", p, "_n_", n, "_c0_", c0, ".RData"))
result <- as.data.frame(result[result[,1] != 999,])
names(result) <- c("Setup 2", "Setup 1")
plot_data <- melt(result, variable.name = "Method", value.name = "pvalue")
plot_data$Method <- factor(plot_data$Method, levels = c("Setup 1", "Setup 2"))

png(file = "Figures/FigureS4.png",
    width = 1200, height = 800)
ggplot2::ggplot(plot_data, aes(x = pvalue, fill = Method)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), color = "#e9ecef", alpha = 0.6, position = 'identity') + scale_fill_manual(values = c("#69b3a2", "#404080")) + labs(fill = "") + xlab("p-values") + ylab("Fraction") + theme(text = element_text(size = 15)) +
  theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size = 40)) + theme(legend.position = "bottom") + theme(
    plot.title = element_text(face = "bold"), axis.title.x = element_text(size = 50), axis.title.y  = element_text(size = 50))
dev.off()