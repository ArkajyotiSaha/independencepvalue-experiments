library(MASS)
library(independencepvalue)

#Code for simulation under global null with a given threshold. 
#input: i0 = Seed, p = #of variables, n = #of samples, c0 = threshold. 
#output: A vector of seletive and classical p-values.
eval_total <- function(i0, p, n, c0) {
  Sigma <- diag(p)
  set.seed(i0)
  X <- MASS::mvrnorm(n = n, rep(0, p), Sigma)#simulate data from identity matrix
  block_diag_structure <- independencepvalue::block_diag(cor(X), c = c0)#threshold to obtain the groups
  if(length(unique(block_diag_structure))> 1){
    set.seed(i0)
    k1 <- sample(unique(block_diag_structure), 1)#randomly choose a group for independence testing
    set.seed(i0)
    t1 <- independencepvalue::selective_p_val(S = cov(X), n = n, CP = block_diag_structure, c = c0, k = k1, d0 = 5, mc_iter = 1000)
    set.seed(i0)
    t2 <- independencepvalue::classical_p_val(S = cov(X), n = n, CP = block_diag_structure, k = k1, mc_iter = 1000)
    return(c(t1, t2))
  }
  else{
    return(c(999, 999))
  }
}

p <- 100

for(ratio in c(1.1, 1.5, 2)){
  n <- p * ratio
  c0 <- sqrt(log(p)/n)
  result <- t(sapply(1:100000, eval_total, p, n, c0))
  save(result, file = paste0("Simulation_results/global_null_p_", p, "_n_", n, "_c0_", c0, ".RData"))
}


