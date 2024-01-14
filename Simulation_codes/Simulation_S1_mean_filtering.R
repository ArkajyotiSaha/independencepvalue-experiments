library(MASS)
library(independencepvalue)
library(matrixStats)

#Code for simulation under global null after mean filtering. 
#input: i0 = Seed, p = # of variables, n = # of samples, c0 = threshold. 
#output: A vector of seletive and classical p-values.
eval_total <- function(i0, p, n, c0) {
  Sigma <- diag(p)
  set.seed(i0)
  X <- mvrnorm(n = n, rep(0, p), Sigma)#simulate data from identity matrix
  colVar_X <- colMeans(X)
  orderedcolVar_X <- order(colVar_X, decreasing = T)
  X <- X[, orderedcolVar_X[1:(p/5)]]#mean filtering, keep the variables with highest means
  block_diag_structure <- block_diag(cor(X), c = c0)#theshold to obtain the groups
  if(length(unique(block_diag_structure))> 1){
    set.seed(i0)
    k1 <- sample(unique(block_diag_structure), 1)#randomly choose a group for independence testing
    set.seed(i0)
    t1 <- selective_p_val(S = cov(X), n = n, CP = block_diag_structure, c = c0, k = k1, d0 = 5, mc_iter = 1000) 
    set.seed(i0)
    t2 <- classical_p_val(S = cov(X), n = n, CP = block_diag_structure, k = k1, mc_iter = 1000)
    return(c(t1, t2))
  }
  else{
    return(c(999, 999))
  }
}



p <- 500
for(ratio in c(1.1, 1.5, 2)){
  n <- p * ratio/5
  c0 <- sqrt(log(p/5)/n)
  result <- t(sapply(1:10000, eval_total, p, n, c0))
  save(result, file = paste0("Simulation_results/Variable_selection_Mean_filtering_global_null_p_", p, "_n_", n, "_c0_", c0, ".RData"))
}

