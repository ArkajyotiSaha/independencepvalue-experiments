library(ggplot2)
library(ggpubr)

unpack_zip_null <- function(files){
  mat_Res <- matrix(0, 0, 3)
  for(j in 1:length(files)){
    load(files[j])
    mat_Res <- rbind(mat_Res, t(dip_res[1:3,]))
  }
  return(mat_Res)
}

plot_list <- list()


Coeffmatrix<-expand.grid(ratio = c(1.1, 1.5,  2))
p <- 100
c0 <- 0.2
p_val_list <- list()
for(indica in 1:3){
  ratio = Coeffmatrix$ratio[indica]
  n <- p * ratio
  zip <- list.files(pattern = paste0("p_",p,"_ratio_",ratio,"_c0_",c0,"_starting_"), path = "Simulation_results/global_null", full.names = T)
  sim <- unpack_zip_null(zip)
  sim  <- sim[sim[,1] != 999, ]
  p_val_list[[indica]] <- sim[,2]
}
total_mat <- cbind(sort(p_val_list[[1]]), 1:length(sort(p_val_list[[1]]))/length(sort(p_val_list[[1]])))

for(i in 2:3){
  total_mat <- rbind(total_mat, cbind(sort(p_val_list[[i]]), 1:length(sort(p_val_list[[i]]))/length(sort(p_val_list[[i]]))))
}

Method_names <- c(rep("n = 1.1p", each = length(sort(p_val_list[[1]]))), rep("n = 1.5p", each = length(sort(p_val_list[[2]]))), rep("n = 2p", each = length(sort(p_val_list[[3]]))))

Method_names <- factor(Method_names, levels = c("n = 1.1p", "n = 1.5p", "n = 2p"))
total_mat <- as.data.frame(total_mat)
total_mat$Method <-  Method_names
p1 <- ggplot(total_mat, aes(x=V2, y=V1, color=Method)) + geom_point()+ geom_line() + geom_abline(intercept = 0, slope = 1, color="black") + xlab("Theoretical quantiles") + ylab("Empirical quantiles") + ggtitle(paste0("(a) Classical Inference")) + theme(text = element_text(size = 20)) + theme(plot.title = element_text(hjust = 0.5)) + labs(color='')




plot_list <- list()
Coeffmatrix<-expand.grid(ratio = c(1.1, 1.5,  2))
p <- 100
c0 <- 0.2
p_val_list <- list()
for(indica in 1:3){
  ratio = Coeffmatrix$ratio[indica]
  n <- p * ratio
  zip <- list.files(pattern = paste0("p_",p,"_ratio_",ratio,"_c0_",c0,"_starting_"), path = "Simulation_results/global_null", full.names = T)
  sim <- unpack_zip_null(zip)
  sim  <- sim[sim[,1] != 999, ]
  p_val_list[[indica]] <- sim[,1]
}


total_mat <- cbind(sort(p_val_list[[1]]), 1:length(sort(p_val_list[[1]]))/length(sort(p_val_list[[1]])))

for(i in 2:3){
  total_mat <- rbind(total_mat, cbind(sort(p_val_list[[i]]), 1:length(sort(p_val_list[[i]]))/length(sort(p_val_list[[i]]))))
}

Method_names <- c(rep("n = 1.1p", each = length(sort(p_val_list[[1]]))), rep("n = 1.5p", each = length(sort(p_val_list[[2]]))), rep("n = 2p", each = length(sort(p_val_list[[3]]))))

Method_names <- factor(Method_names, levels = c("n = 1.1p", "n = 1.5p", "n = 2p"))
total_mat <- as.data.frame(total_mat)

total_mat$Method <-  Method_names

p2 <- ggplot2::ggplot(total_mat, aes(x=V2, y=V1, color=Method)) + geom_point()+ geom_line() + geom_abline(intercept = 0, slope = 1, color="black") + xlab("Theoretical quantiles") + ylab("Empirical quantiles") + ggtitle(paste0("(b) Selective Inference")) + theme(text = element_text(size = 20)) + theme(plot.title = element_text(hjust = 0.5)) + labs(color='')


png(file="Figures/Figure2.png",
    width=1200, height=800)
ggpubr::ggarrange(p1, p2, ncol=2, common.legend=TRUE, legend="bottom")
dev.off()

