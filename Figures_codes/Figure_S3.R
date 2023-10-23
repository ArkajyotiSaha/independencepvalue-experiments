library(ggplot2)
library(ggpubr)


p <- 100

plot_list <- list()
mu =  c(2, 10, 30)
for(index in 1:3){
  plot_data <- data.frame()
  for(ratio in c(1.1, 1.5, 2)){
    n <- p * ratio
    c0 <- sqrt(log(p)/n)
    load(file=paste0("Simulation_results/Non_Gaussian_Poisson_distribution_",mu[index],"_global_null_p_", p, "_n_", n, "_c0_", c0, ".RData"))
    plot_data <- rbind(plot_data, result)
  }
  names(plot_data) <- c("Selective","Classical")
  plot_data$Method <- rep(c("n = 1.1p","n = 1.5p","n = 2p"), each = 10000)
  plot_data$Method <- factor(plot_data$Method, levels = c("n = 1.1p", "n = 1.5p", "n = 2p"))
  plot_data <- plot_data[plot_data[,1]!= 999,]
  
  base_plot <- ggplot(data = plot_data) +
    xlab("Theoretical quantiles")  +
    ylab("Empirical quantiles") +
    geom_abline()
  plot_list[[index]] <- list()
  plot_list[[index]][[1]] <- base_plot + geom_qq(distribution = qunif,
                                                 aes(sample = Classical, color = Method)) + 
    ggtitle(bquote(mu~"="~.(mu[index])))+ theme(text = element_text(size = 30)) + 
    theme(plot.title = element_text(hjust = 0.5)) + labs(color='')
  
  plot_list[[index]][[2]] <- base_plot + geom_qq(distribution = qunif,
                                                 aes(sample = Selective, color = Method)) +
    ggtitle(bquote(mu~"="~.(mu[index])))+ theme(text = element_text(size = 30)) + 
    theme(plot.title = element_text(hjust = 0.5)) + labs(color='')
}

classical_inference_plot <- ggpubr::ggarrange(plot_list[[1]][[1]], plot_list[[2]][[1]], plot_list[[3]][[1]], ncol=3, legend="none")
classical_inference_plot_heading <- annotate_figure(classical_inference_plot, top = text_grob("(a) Classical Inference", 
                                                                                              size = 40))

selective_inference_plot <- ggpubr::ggarrange(plot_list[[1]][[2]], plot_list[[2]][[2]], plot_list[[3]][[2]], ncol=3, legend="none")
selective_inference_plot_heading <- annotate_figure(selective_inference_plot, top = text_grob("(b) Selective Inference", 
                                                                                              size = 40))
common_legend <- get_legend(plot_list[[1]][[1]] + theme(text = element_text(size = 40), legend.position = "bottom"))

png(file="Figures/FigureS3.png",
    width=1500, height=1150)
grid.arrange(arrangeGrob(classical_inference_plot_heading, 
                         selective_inference_plot_heading, nrow=2), 
             common_legend, heights=c(15, 1))
dev.off()

