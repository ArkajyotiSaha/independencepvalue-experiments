library(independencepvalue)
library(clusterGeneration)
library(MASS)
Coeffmatrix<-expand.grid(p = c(100), ratio = c(1.1, 1.5, 2), starting = seq(1, 1000000, by = 10000 ))
indica <- as.numeric(Sys.getenv('SGE_TASK_ID'))

p = Coeffmatrix$p[indica]
ratio = Coeffmatrix$ratio[indica]
n = p * ratio
starting = Coeffmatrix$starting[indica]

effect_size <- function(Sigma, CP, k){
  S <- Sigma
  p <- nrow(S)
  ptemp <- sum(CP==k)
  if(2*ptemp >= p){
    p1 <- ptemp
    S11_x <- as.matrix(S[which(CP==k), which(CP==k)])#S11
    S22_x <- as.matrix(S[which(CP!=k), which(CP!=k)])#S22
    S12_x <- as.matrix(S[which(CP==k), which(CP!=k)])#S_12
  }
  if(2*ptemp < p){
    p1 = p - ptemp
    S11_x <- as.matrix(S[which(CP!=k), which(CP!=k)])#S11
    S22_x <- as.matrix(S[which(CP==k), which(CP==k)])#S22
    S12_x <- as.matrix(S[which(CP!=k), which(CP==k)])#S_12
  }
  #ensures that p2 <= p1.
  p2 <- p - p1
  if(p1 > p2){rp <- p2}#Define r(P)
  if(p1 <= p2){rp <- p1}
  inv_S11_x <- solve(S11_x)#inv(S11)
  inv_S11_x_half <- amen::mhalf(inv_S11_x)#S_11^{-1/2}
  inv_S22_x <- solve(S22_x)#inv(S22)
  inv_S22_x_half <- amen::mhalf(inv_S22_x)#S_22^{-1/2}
  tilde_S12_x <- inv_S11_x_half %*% S12_x %*% inv_S22_x_half#S_12^W = covariance matrix of whitened X_1 and X_2
  svdecom <- svd(tilde_S12_x)#compact SVD
  singular_values <- svdecom$d#Lambda
  test_stat <- prod(1-singular_values^2)
  return(test_stat)
}

eval_total <- function(i0, p, n) {
  set.seed(i0)
  Sigma_cov <- clusterGeneration::genPositiveDefMat(p)$Sigma + matrix(runif(1), p, p)
  Sigma_cor <- Sigma_cov
  for(i in 1:p){
    for(j in 1:p){
      Sigma_cor[i,j] <- Sigma_cov[i,j]/sqrt(Sigma_cov[i,i] * Sigma_cov[j,j])
    }
  }
  Sigma <- Sigma_cor
  
  Var_S <- Sigma
  diag_S <- diag(1/sqrt(diag(Var_S)))
  Cor_S <- diag_S %*% Var_S %*% diag_S
  dis_S <- 1 - abs(Cor_S)
  test <- as.dist(dis_S, diag = TRUE)
  clust_result <- hclust(test, method = "single")
  c0 <- 1 - clust_result$height[p - 1] + sqrt(log(p)/n)
  c0 <- min(c0, 0.99999999)
  
  set.seed(i0)
  X <- MASS::mvrnorm(n = n, rep(0, p), Sigma)
  block_diag_structure <- independencepvalue::block_diag(cor(X), c= c0)
  if(length(unique(block_diag_structure))> 1){
    set.seed(i0)
    k1 <- sample(unique(block_diag_structure), 1)
    set.seed(i0)
    t1 <- independencepvalue::selective_p_val(S=cov(X), n=n, CP=block_diag_structure, c=c0, k=k1, d0=5, mc_iter=1000)
    set.seed(i0)
    t2 <- independencepvalue::classical_p_val(S=cov(X), n=n, CP=block_diag_structure, k=k1, mc_iter=1000)
    t3 <- effect_size(Sigma, CP=block_diag_structure, k=k1)
    p1 <- sum(block_diag_structure==k1)
    if(2*p1 < p){p1 <- p - sum(block_diag_structure==k1)}
    return(c(t1, t2, t3, p1))
  }
  else{
    return(c(999, 999, 999, 999))
  }
}

t1 <- proc.time()
dip_res <- sapply(starting:(starting + 9999), eval_total, p, n)
t2 <- proc.time()

time_tot <- t2 - t1

save(dip_res, time_tot, file = paste0("Simulation_results/alternative/p_",p,"_ratio_",ratio,"_starting_",starting,".RData"))

