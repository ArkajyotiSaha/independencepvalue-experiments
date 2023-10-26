library(MASS)
library(independencepvalue)


eval_total <- function(i0, p, n, c0, nu) {
  set.seed(i0)
  X <- matrix(rt(n*p, df=nu), n, p)
  block_diag_structure <- independencepvalue::block_diag(cor(X), c=c0)
  if(length(unique(block_diag_structure))> 1){
    set.seed(i0)
    k1 <- sample(unique(block_diag_structure), 1)
    set.seed(i0)
    t1 <- independencepvalue::selective_p_val(S=cov(X), n=n, CP=block_diag_structure, c=c0, k=k1, d0=5, mc_iter=1000)
    set.seed(i0)
    t2 <- independencepvalue::classical_p_val(S=cov(X), n=n, CP=block_diag_structure, k=k1, mc_iter=1000)
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
  for(nu in c(2, 10, 30)){
    result <- t(sapply(1:10000, eval_total, p, n, c0, nu))
    save(result, file=paste0("Simulation_results/Non_Gaussian_t_distribution_",nu,"_global_null_p_", p, "_n_", n, "_c0_", c0, ".RData"))
  }
}

