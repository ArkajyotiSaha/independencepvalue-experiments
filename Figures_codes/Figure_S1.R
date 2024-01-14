library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)

p <- 500

plot_list <- list()
filtering_methods =  c("Variance", "Mean")
for(index in 1:2){
  plot_data <- data.frame()
  for(ratio in c(1.1, 1.5, 2)){
    n <- p * ratio/5
    c0 <- sqrt(log(p/5)/n)
    load(file = paste0("Simulation_results/Variable_selection_",filtering_methods[index],"_filtering_global_null_p_", p, "_n_", n, "_c0_", c0, ".RData"))
    plot_data <- rbind(plot_data, result)
  }
  
  names(plot_data) <- c("Selective","Classical")
  plot_data$Method <- rep(c("n = 1.1p/5","n = 1.5p/5","n = 2p/5"), each = 10000)
  plot_data$Method <- factor(plot_data$Method, levels = c("n = 1.1p/5", "n = 1.5p/5", "n = 2p/5"))
  plot_data <- plot_data[plot_data[,1]!= 999,]
  
  base_plot <- ggplot(data = plot_data) +
    xlab("Theoretical quantiles")  +
    ylab("Empirical quantiles") +
    geom_abline()
  plot_list[[index]] <- list()
  plot_list[[index]][[1]] <- base_plot + geom_qq(distribution = qunif,
                                                 aes(sample = Classical, color = Method)) + 
    ggtitle(paste0("Classical Inference")) + 
    theme(text = element_text(size = 30)) + 
    theme(plot.title = element_text(hjust = 0.5)) + labs(color = '')
  
  plot_list[[index]][[2]] <- base_plot + geom_qq(distribution = qunif,
                                                 aes(sample = Selective, color = Method)) +
    ggtitle(paste0("Selective Inference")) + 
    theme(text = element_text(size = 30)) + 
    theme(plot.title = element_text(hjust = 0.5)) + labs(color = '')
}

variance_filter_plot <- ggpubr::ggarrange(plot_list[[1]][[1]], plot_list[[1]][[2]], ncol = 2, legend = "none")
variance_filter_plot_heading <- annotate_figure(variance_filter_plot, top = text_grob("(a) Variance filtering", 
                                                                                      size = 40))

mean_filter_plot <- ggpubr::ggarrange(plot_list[[2]][[1]], plot_list[[2]][[2]], ncol = 2, legend = "none")
mean_filter_plot_heading <- annotate_figure(mean_filter_plot, top = text_grob("(b) Mean filtering", 
                                                                              size = 40))
common_legend <- get_legend(plot_list[[1]][[1]] + theme(text = element_text(size = 40), legend.position = "bottom"))

png(file = "Figures/FigureS1.png",
    width = 1050, height = 1150)
grid.arrange(arrangeGrob(variance_filter_plot_heading, 
                         mean_filter_plot_heading, nrow = 2), 
             common_legend, heights = c(15, 1))
dev.off()

