library(ggplot2)
library(reshape2)
library(dplyr)

power <- function(x){mean(x < 0.05)}

p <- 3
ratio <- 1.5
n <- floor(p * ratio)

load(file = paste0("Simulation_results/SimulationS6_alternative_p_", p, "_n_", n, ".RData"))
result <- as.data.frame(result)
names(result) <- c("p'", "p", "Effect", "zeta")
plot_data <- melt(result, variable.name = "Method", value.name = "pvalue", id = c("Effect", "zeta"))
plot_data$Method <- factor(plot_data$Method, levels = c("p'", "p"))
plot_data <- plot_data[plot_data$pvalue != 999,]
plot_data$Effect <- 1 - plot_data$Effect

plot_data_a <- plot_data[plot_data$zeta < 0.05,]
plot_data_a  <- as.data.frame(plot_data_a  %>% 
                                group_by(Method, Bin.Effect = cut(Effect, breaks = seq(0, 1, by = 0.1))) %>% 
                                filter(length(pvalue) >= 100) %>%
                                ungroup())

p1 <- ggplot(plot_data_a, aes(Effect, pvalue, color = Method)) +
  stat_summary_bin(fun = "power", geom = "line", orientation = 'x', breaks = seq(0, 1, 0.1)) + 
  stat_summary_bin(fun = "power", geom = "point", orientation = 'x', breaks = seq(0, 1, 0.1)) +
  xlab("Effect Size (\u0394)") + ylab("Power at level \u03b1 = 0.05") + labs(color = '') + scale_colour_manual(values = c("#00BFC4", "#C77CFF"), labels = c( expression(paste(p, "' (S23)")), expression(paste(p, " (11)")))) +
  ggtitle(expression(paste("(a) ",zeta,"(",  x,")",'  < 0.05'))) + theme(text = element_text(size = 20)) +
  theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "bottom")


plot_data_b <- plot_data[plot_data$zeta > 0.45,]
plot_data_b  <- as.data.frame(plot_data_b  %>% 
                                group_by(Method, Bin.Effect = cut(Effect, breaks = seq(0, 1, by = 0.1))) %>% 
                                filter(length(pvalue) >= 100) %>%
                                ungroup())

p2 <- ggplot(plot_data_b, aes(Effect, pvalue, color = Method)) +
  stat_summary_bin(fun = "power", geom = "line", orientation = 'x', breaks = seq(0, 1, 0.1)) + 
  stat_summary_bin(fun = "power", geom = "point", orientation = 'x', breaks = seq(0, 1, 0.1)) +
  xlab("Effect Size (\u0394)") + ylab("Power at level \u03b1 = 0.05") + labs(color = '') + scale_colour_manual(values = c("#00BFC4", "#C77CFF"), labels = c( expression(paste(p, "' (S23)")), expression(paste(p, " (11)")))) +
  ggtitle(expression(paste("(b) ",zeta,"(",  x,")",'  > 0.45'))) + theme(text = element_text(size = 20)) +
  theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "bottom")

png(file = "Figures/FigureS6.png",
    width = 800, height = 420)
ggpubr::ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "bottom")
dev.off()
