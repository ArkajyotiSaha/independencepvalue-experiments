library(ggplot2)
library(reshape2)

p <- 6
c0 <- 0.5
ratio <- 1.5
n <- p * ratio
load(file = paste0("Simulation_results/Simulation1c_p_", p, "_n_", n, "_c0_", c0, ".RData"))
result <- as.data.frame(result[result[,1] != 999,])
names(result) <- c("Selective", "Classical")
plot_data <- melt(result, variable.name = "Inference", value.name = "pvalue")
plot_data$Inference <- factor(plot_data$Inference, levels = c("Classical", "Selective"))

png(file = "Figures/Figure1c.png",
    width = 1200, height = 800)
ggplot2::ggplot(plot_data, aes(x = pvalue, fill = Inference)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), color = "#e9ecef", alpha = 0.6, position = 'identity') + scale_fill_manual(values = c("#69b3a2", "#404080")) + labs(fill="") + xlab("p-values") + ylab("Fraction") + theme(text = element_text(size = 15)) + ggtitle("(c) Histogram of p-values") +
  theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size = 40)) + theme(legend.position = "top" ) + theme(
    plot.title = element_text(face = "bold"), axis.title.x = element_text(size = 50), axis.title.y  = element_text(size = 50))
dev.off()