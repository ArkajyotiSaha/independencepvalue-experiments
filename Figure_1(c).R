library(ggplot2)


p <- 6
c0 <- 0.5
ratio <- 1.5
p_val_list <- list()
for(indica in 1){
  n <- p * ratio
  zip <- list.files(pattern = paste0("p_",p,"_ratio_",ratio,"_c0_",c0,"_starting_"), path = "Simulation_results/histogram", full.names = T)
  mat_res <- rep(0,0)
  for(j in 1:length(zip)){
    load(zip[j])
    mat_res <- c(mat_res, dip_res[2,][dip_res[1,]!=999])
  }
  p_val_list[[indica]] <- mat_res
}
total_mat <- cbind(sort(p_val_list[[1]]), 1:length(sort(p_val_list[[1]]))/length(sort(p_val_list[[1]])))
Method_names <- c(rep("Classical", each = length(sort(p_val_list[[1]]))))
total_mat <- as.data.frame(total_mat)
total_mat$Method <-  Method_names
for(indica in 1){
  n <- p * ratio
  zip <- list.files(pattern = paste0("p_",p,"_ratio_",ratio,"_c0_",c0,"_starting_"), path = "Simulation_results/histogram", full.names = T)
  mat_res <- rep(0,0)
  for(j in 1:length(zip)){
    load(zip[j])
    mat_res <- c(mat_res, dip_res[1,][dip_res[1,]!=999])
  }
  p_val_list[[indica]] <- mat_res
}
total_mat_new <- cbind(sort(p_val_list[[1]]), 1:length(sort(p_val_list[[1]]))/length(sort(p_val_list[[1]])))
Method_names <- c(rep("Selective", each = length(sort(p_val_list[[1]]))))
total_mat_new  <- as.data.frame(total_mat_new )
total_mat_new$Method <-  Method_names
total_mat_hist <- rbind(total_mat, total_mat_new)
p1 <- ggplot2::ggplot(total_mat_hist, aes(x = V1, fill=Method)) +
  geom_histogram(aes(y = stat(count / sum(count))), color="#e9ecef", alpha=0.6, position = 'identity') + scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + xlab("p-values") + ylab("Fraction") + theme(text = element_text(size = 15)) + ggtitle("(c) Histogram of p-values") +
  theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size = 40)) + theme(legend.position="top") + theme(
    plot.title = element_text(face="bold"), axis.title.x = element_text( size=50), axis.title.y  = element_text( size=50)
  )



png(file="Figures/Figure1c.png",
    width=1200, height=800)
p1
dev.off()