library(MASS)
library(independencepvalue)
library(clusterGeneration)

source("Simulation_codes/Effect_size.R")
source("Simulation_codes/Selection_probability.R")
p <- 100
#Code for simulation under alternative. 
#input: i0 = Seed, p = # of variables, n = # of samples. 
#output: A vector of seletive and classical p-values, popoulation effect size, number of nonzero canonical correlations, and probability of the selection event.
eval_total <- function(i0, p, n) {
  set.seed(i0)
  Sigma_cov <- clusterGeneration::genPositiveDefMat(p)$Sigma + matrix(runif(1), p, p)#generate a dense p.d. matrix. 
  Sigma_cor <- Sigma_cov
  for(i in 1:p){
    for(j in 1:p){
      Sigma_cor[i, j] <- Sigma_cov[i, j]/sqrt(Sigma_cov[i, i] * Sigma_cov[j, j])#Standardize the variance and set it as population coariance matrix as in S6.1 of the Supplementaary Materials.
    }
  }
  Sigma <- Sigma_cor
  
  Var_S <- Sigma
  diag_S <- diag(1/sqrt(diag(Var_S)))
  Cor_S <- diag_S %*% Var_S %*% diag_S
  dis_S <- 1 - abs(Cor_S)
  test <- as.dist(dis_S, diag = TRUE)
  clust_result <- hclust(test, method = "single")
  c0 <- 1 - clust_result$height[p - 1] + sqrt(log(p)/n)#choose a threshold based on the population covariance matrix as in S6.2 of the Supplementary Materials. 
  c0 <- min(c0, 0.99999999)
  
  set.seed(i0)
  X <- MASS::mvrnorm(n = n, rep(0, p), Sigma)
  block_diag_structure <- independencepvalue::block_diag(cor(X), c = c0)
  if(length(unique(block_diag_structure))> 1){
    set.seed(i0)
    k1 <- sample(unique(block_diag_structure), 1)
    set.seed(i0)
    t1 <- independencepvalue::selective_p_val(S = cov(X), n = n, CP = block_diag_structure, c = c0, k = k1, d0 = 5, mc_iter = 1000, maxeval = 10000)
    selectprob <- selection(S = cov(X), n = n, CP = block_diag_structure, c = c0, k = k1, d0 = 5, mc_iter = 1000, maxeval = 10000)
    set.seed(i0)
    t2 <- independencepvalue::classical_p_val(S = cov(X), n = n, CP = block_diag_structure, k = k1, mc_iter = 1000)
    t3 <- effect_size(Sigma, CP = block_diag_structure, k = k1)
    p1 <- sum(block_diag_structure == k1)
    if(2*p1 < p){p1 <- p - sum(block_diag_structure == k1)}
    return(c(t1, t2, t3, p1, selectprob))
  }
  else{
    return(c(999, 999, 999, 999, 999))
  }
}


for(ratio in c(1.1, 1.5, 2)){
  n <- p * ratio
  result <- t(sapply(1:1000000, eval_total, p, n))
  save(result, file = paste0("Simulation_results/alternative_p_", p, "_n_", n, ".RData"))
}

