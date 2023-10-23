library(Matrix)
library(MASS)
library(independencepvalue)

p <- 6
ratio <- 1.5
n <- p * ratio
c0 <- 0.5

Sigma <- create_example(p = p, a = 0.6, b = 0.3)


eval_total <- function(i0, p, n, c0, Sigma) {
  set.seed(i0)
  X <- mvrnorm(n=n, rep(0, p), Sigma)
  block_diag_structure <- independencepvalue::block_diag(cor(X), c=c0)
  if(length(unique(block_diag_structure))> 1){
    set.seed(i0)
    k1 <- sample(unique(block_diag_structure), 1)
    set.seed(i0)
    t1 <- independencepvalue::selective_p_val(S=cov(X), n=n, CP=block_diag_structure, c=c0, k=k1, d0=5, mc_iter=1000, maxeval = 10000)
    set.seed(i0)
    t2 <- independencepvalue::classical_p_val(S=cov(X), n=n, CP=block_diag_structure, k=k1, mc_iter=1000)
    return(c(t1, t2))
  }
  else{
    return(c(999, 999))
  }
}

result <- t(sapply(1:10000, eval_total, p, n, c0, Sigma))

save(result, file=paste0("Simulation_results/Simulation1c_p_", p, "_n_", n, "_c0_", c0, ".RData"))


