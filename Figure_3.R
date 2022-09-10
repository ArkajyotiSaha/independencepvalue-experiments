unpack_zip_res_new <- function(files){
  mat_Res <- matrix(0, 0, 4)
  for(j in 1:length(files)){
    load(files[j])
    mat_Res <- rbind(mat_Res, t(dip_res[1:4,]))
  }
  return(mat_Res)
}
library(ggplot2)
library(ggpubr)
Coeffmatrix_part<-expand.grid(part <- seq(0.1, 1, 0.1), ratio = c(1.1, 1.5, 2))
Coeffmatrix_part<- cbind(Coeffmatrix_part, rep(0, nrow(Coeffmatrix_part)))
Coeffmatrix_part<- cbind(Coeffmatrix_part, rep(0, nrow(Coeffmatrix_part)))
Coeffmatrix<-expand.grid(ratio = c(1.1, 1.5, 2))
p <- 100
for(i in 1:nrow(Coeffmatrix)){
  ratio <- Coeffmatrix$ratio[i]
  zip <- list.files(pattern = paste0("p_",p,"_ratio_",ratio,"_starting_"), path = "Simulation_results/alternative", full.names = T)
  sim_res <- unpack_zip_res_new(zip)
  tot_result_comb_1 <- sim_res[sim_res[,1] != 999, ]
  tot_result_comb <- tot_result_comb_1[tot_result_comb_1[,4] >= 95,]
  tot_result_comb[,3] <- 1 - tot_result_comb[,3]
  for(jip in 1:10){
    gh_selec <- tot_result_comb[tot_result_comb[,3] <= (jip * 0.1) & tot_result_comb[,3] > ((jip-1) * 0.1),]
    if(is.vector(gh_selec)){
      gh_selec <- t(as.matrix(gh_selec))
    }
    if(nrow(gh_selec)>100){
      Coeffmatrix_part[((i-1)*10 + jip), 3:4] <- c(mean(gh_selec[,1]< 0.05, na.rm = T), mean(gh_selec[,2]< 0.05, na.rm = T))
    }
    if(nrow(gh_selec) <= 100){
      Coeffmatrix_part[((i-1)*10 + jip), 3:4] <- NA
    }
  }
}
Coeffmatrix_part <- Coeffmatrix_part[complete.cases(Coeffmatrix_part),]
Coeffmatrix_part_melt <- as.data.frame(rbind(Coeffmatrix_part[,1:3], Coeffmatrix_part[,c(1:2, 4)]))
names(Coeffmatrix_part_melt ) <- c("part", "n", "Power")
Coeffmatrix_part_melt$Inference <- rep(c("Selective", "Classical"), each = nrow(Coeffmatrix_part))
Coeffmatrix_part_melt$n <- paste0("n = ", Coeffmatrix_part_melt$n, "p")
Coeffmatrix_part_melt$n <- as.factor(Coeffmatrix_part_melt$n)


p1 <- ggplot(Coeffmatrix_part_melt, aes(x=part, y=Power, color=n, linetype = Inference)) + geom_point()+ geom_line() + 
  xlab("Effect Size (\u0394)") + ylab("Power at level \u03b1 = 0.05") + labs(color='') + 
  ggtitle(expression(paste("(a) ",r,"(",  hat(P),")" ,' \u2264 5'))) + theme(text = element_text(size = 20)) + 
  theme(plot.title = element_text(hjust = 0.5))  + ylim(0, 0.9)




Coeffmatrix_part<-expand.grid(part <- seq(0.1, 1, 0.1), ratio = c(1.1, 1.5, 2))
Coeffmatrix_part<- cbind(Coeffmatrix_part, rep(0, nrow(Coeffmatrix_part)))
Coeffmatrix_part<- cbind(Coeffmatrix_part, rep(0, nrow(Coeffmatrix_part)))
Coeffmatrix<-expand.grid(ratio = c(1.1, 1.5, 2))
p <- 100
for(i in 1:nrow(Coeffmatrix)){
  ratio <- Coeffmatrix$ratio[i]
  zip <- list.files(pattern = paste0("p_",p,"_ratio_",ratio,"_starting_"), path = "Simulation_results/alternative", full.names = T)
  sim_res <- unpack_zip_res_new(zip)
  tot_result_comb_1 <- sim_res[sim_res[,1] != 999, ]
  tot_result_comb <- tot_result_comb_1[tot_result_comb_1[,4] < 95,]
  tot_result_comb[,3] <- 1 - tot_result_comb[,3]
  for(jip in 1:10){
    gh_selec <- tot_result_comb[tot_result_comb[,3] <= (jip * 0.1) & tot_result_comb[,3] > ((jip-1) * 0.1),]
    if(is.vector(gh_selec)){
      gh_selec <- t(as.matrix(gh_selec))
    }
    if(nrow(gh_selec)>100){
      Coeffmatrix_part[((i-1)*10 + jip), 3:4] <- c(mean(gh_selec[,1]< 0.05, na.rm = T), mean(gh_selec[,2]< 0.05, na.rm = T))
    }
    if(nrow(gh_selec) <= 100){
      Coeffmatrix_part[((i-1)*10 + jip), 3:4] <- NA
    }
  }
}
Coeffmatrix_part <- Coeffmatrix_part[complete.cases(Coeffmatrix_part),]
Coeffmatrix_part_melt <- as.data.frame(rbind(Coeffmatrix_part[,1:3], Coeffmatrix_part[,c(1:2, 4)]))
names(Coeffmatrix_part_melt ) <- c("part", "n", "Power")
Coeffmatrix_part_melt$Inference <- rep(c("Selective", "Classical"), each = nrow(Coeffmatrix_part))
Coeffmatrix_part_melt$n <- paste0("n = ", Coeffmatrix_part_melt$n, "p")
Coeffmatrix_part_melt$n <- as.factor(Coeffmatrix_part_melt$n)

p2 <- ggplot(Coeffmatrix_part_melt, aes(x=part, y=Power, color=n, linetype = Inference)) + geom_point()+ geom_line() + 
  xlab("Effect Size (\u0394)") + ylab("Power at level \u03b1 = 0.05") + labs(color='') + 
  ggtitle(expression(paste("(b) ",r,"(",  hat(P),")" ,' > 5'))) + theme(text = element_text(size = 20)) + 
  theme(plot.title = element_text(hjust = 0.5)) + ylim(0, 0.9)


png(file="Figures/Figure3.png",
    width=1200, height=800)
ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend="bottom")
dev.off()
