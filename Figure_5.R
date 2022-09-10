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
network <- 4 # 3 for E.coli & 4 for S.cerevisiae
Coeffmatrix_part<-expand.grid(part <- seq(0.1, 1, 0.1), ratio = c(1.2, 1.4, 1.5, 2))
Coeffmatrix_part<- cbind(Coeffmatrix_part, rep(0, nrow(Coeffmatrix_part)))
Coeffmatrix_part<- cbind(Coeffmatrix_part, rep(0, nrow(Coeffmatrix_part)))
Coeffmatrix<-expand.grid( ratio = c(1.2, 1.4, 1.5, 2))
p <- 100
noise <- 0
for(i in 1:nrow(Coeffmatrix)){
  ratio <- Coeffmatrix$ratio[i]
  zip <- list.files(path = "Real_data_results", pattern = paste0("network_",network,"p_",p,"_ratio_",ratio,"_noise_",noise,"_starting_"), full.names = T)
  sim_res <- unpack_zip_res_new(zip)
  tot_result_comb_1 <- sim_res[sim_res[,1] != 999, ]
  tot_result_comb <- tot_result_comb_1[tot_result_comb_1[,4] >= 50,]
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

p1 <- ggplot(Coeffmatrix_part_melt, aes(x=part, y=Power, color=n, linetype = Inference)) + geom_point()+ 
  geom_line() + xlab(expression(paste("Effect Size (", hat('\u0394'), ")"))) + 
  ylab("Power at level \u03b1 = 0.05") + labs(color='') + ggtitle(paste0("(a) \u03b7 = ",noise)) + 
  theme(text = element_text(size = 20)) + theme(plot.title = element_text(hjust = 0.5))+ ylim(0,1.01)

Coeffmatrix_part<-expand.grid(part <- seq(0.1, 1, 0.1), ratio = c(1.2, 1.4, 1.5, 2))
Coeffmatrix_part<- cbind(Coeffmatrix_part, rep(0, nrow(Coeffmatrix_part)))
Coeffmatrix_part<- cbind(Coeffmatrix_part, rep(0, nrow(Coeffmatrix_part)))
Coeffmatrix<-expand.grid( ratio = c(1.2, 1.4, 1.5, 2))
noise <- 0.5
for(i in 1:nrow(Coeffmatrix)){
  ratio <- Coeffmatrix$ratio[i]
  zip <- list.files(path = "Real_data_results", pattern = paste0("network_",network,"p_",p,"_ratio_",ratio,"_noise_",noise,"_starting_"), full.names = T)
  sim_res <- unpack_zip_res_new(zip)
  tot_result_comb_1 <- sim_res[sim_res[,1] != 999, ]
  tot_result_comb <- tot_result_comb_1[tot_result_comb_1[,4] >= 50,]
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

p2 <- ggplot(Coeffmatrix_part_melt, aes(x=part, y=Power, color=n, linetype = Inference)) + geom_point()+ 
  geom_line() + xlab(expression(paste("Effect Size (", hat('\u0394'), ")"))) + 
  ylab("Power at level \u03b1 = 0.05") + labs(color='') + ggtitle(paste0("(b) \u03b7 = ",noise)) + 
  theme(text = element_text(size = 20)) + theme(plot.title = element_text(hjust = 0.5))+ ylim(0,1.01) 

Coeffmatrix_part<-expand.grid(part <- seq(0.1, 1, 0.1), ratio = c(1.2, 1.4, 1.5, 2))
Coeffmatrix_part<- cbind(Coeffmatrix_part, rep(0, nrow(Coeffmatrix_part)))
Coeffmatrix_part<- cbind(Coeffmatrix_part, rep(0, nrow(Coeffmatrix_part)))
Coeffmatrix<-expand.grid( ratio = c(1.2, 1.4, 1.5, 2))
noise <- 3
for(i in 1:nrow(Coeffmatrix)){
  ratio <- Coeffmatrix$ratio[i]
  zip <- list.files(path = "Real_data_results", pattern = paste0("network_",network,"p_",p,"_ratio_",ratio,"_noise_",noise,"_starting_"), full.names = T)
  sim_res <- unpack_zip_res_new(zip)
  tot_result_comb_1 <- sim_res[sim_res[,1] != 999, ]
  tot_result_comb <- tot_result_comb_1[tot_result_comb_1[,4] >=50,]
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

p3 <- ggplot(Coeffmatrix_part_melt, aes(x=part, y=Power, color=n, linetype = Inference)) + geom_point()+ 
  geom_line() + xlab(expression(paste("Effect Size (", hat('\u0394'), ")"))) + 
  ylab("Power at level \u03b1 = 0.05") + labs(color='') + ggtitle(paste0("(c) \u03b7 = ",noise)) + 
  theme(text = element_text(size = 20)) + theme(plot.title = element_text(hjust = 0.5))+ ylim(0,1.01)

png(file="Figures/Figure5.png",
    width=1800, height=600)
ggarrange(p1, p2, p3, ncol=3, common.legend = TRUE, legend="bottom")
dev.off()