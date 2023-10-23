library(MASS)
library(independencepvalue)

source("Simulation_codes/Effect_size.R")
source("Simulation_codes/Utils.R")

eval_total <- function(i0, p, n) {
  set.seed(i0)
  Sigma_cov <- clusterGeneration::genPositiveDefMat(p)$Sigma + matrix(runif(1), p, p)
  Sigma_cor <- Sigma_cov
  for(i in 1:p){
    for(j in 1:p){
      Sigma_cor[i, j] <- Sigma_cov[i, j]/sqrt(Sigma_cov[i, i] * Sigma_cov[j, j])
    }
  }
  Sigma <- Sigma_cor
  
  Var_S <- Sigma
  diag_S <- diag(1/sqrt(diag(Var_S)))
  Cor_S <- diag_S %*% Var_S %*% diag_S
  dis_S <- 1 - abs(Cor_S)
  test <- as.dist(dis_S, diag=TRUE)
  clust_result <- hclust(test, method="single")
  c0 <- 1 - clust_result$height[p - 1] + sqrt(log(p)/n)
  c0 <- min(c0, 0.99999999)
  
  set.seed(i0)
  X <- MASS::mvrnorm(n=n, rep(0, p), Sigma)
  block_diag_structure <- independencepvalue::block_diag(cor(X), c=c0)
  if(length(unique(block_diag_structure))> 1){
    set.seed(i0)
    k1 <- sample(unique(block_diag_structure), 1)
    test_hyp <- independencepvalue:::test_stat_CCA(cov(X), block_diag_structure, k1)
    sip <- sapply(1:1000, sample_simul, i0, p, nrow(test_hyp$S11), n, test_hyp$S11, test_hyp$S22_x, c0)
    zp <- sum(as.numeric(sip[2, ]))
    print(zp)
    if(zp < 100){
      sip <- sapply(1:(min((1000 * 100/zp), 1e+05)), sample_simul, i0, p, nrow(test_hyp$S11), n, test_hyp$S11, test_hyp$S22_x, c0)
    }
    t0 <- mean(test_hyp$statistic >= sip[1, ][sip[2, ] == TRUE])
    set.seed(i0)
    t1 <- independencepvalue::selective_p_val(S=cov(X), n=n, CP=block_diag_structure, c=c0, k=k1, d0=1, mc_iter=1000)
    t2 <- effect_size(Sigma, CP=block_diag_structure, k=k1)
    t3 <- min(svd(test_stat_CCA_new(cov(X), block_diag_structure, k1)$S12w)$u^2)
    return(c(t0, t1, t2, t3))
  }
  else{
    return(c(999, 999, 999, 999))
  }
}


p <- 3
ratio <- 1.5
n <- floor(p * ratio)

result <- t(sapply(1:10000, eval_total, p, n))

save(result, file=paste0("Simulation_results/SimulationS6_alternative_p_", p, "_n_", n, ".RData"))


