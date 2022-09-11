library(MASS)
library(independencepvalue)
Coeffmatrix<-expand.grid(p = c(100), ratio = c(1.1, 1.5, 2), c0 = 0.2, starting = seq(1, 100000, by = 1000))

indica <- as.numeric(Sys.getenv('SGE_TASK_ID'))

p = Coeffmatrix$p[indica]
ratio = Coeffmatrix$ratio[indica]
n = p * ratio
c0 = Coeffmatrix$c0[indica]
starting = Coeffmatrix$starting[indica]

eval_total <- function(i0, p, n, c0) {
  Sigma <- diag(p)
  set.seed(i0)
  X <- mvrnorm(n = n, rep(0, p), Sigma)
  block_diag_structure <- independencepvalue::block_diag(cor(X), c= c0)
  if(length(unique(block_diag_structure))> 1){
    set.seed(i0)
    k1 <- sample(unique(block_diag_structure), 1)
    set.seed(i0)
    t1 <- independencepvalue::selective_p_val(S=cov(X), n=n, CP=block_diag_structure, c=c0, k=k1, d0=5, mc_iter=1000)
    set.seed(i0)
    t2 <- independencepvalue::classical_p_val(S=cov(X), n=n, CP=block_diag_structure, k=k1, mc_iter=1000)
    p1 <- sum(block_diag_structure==k1)
    if(2*p1 < p){p1 <- p - sum(block_diag_structure==k1)}
    return(c(t1, t2, p1))
  }
  else{
    return(c(999, 999, 999))
  }
}


t1 <- proc.time()
dip_res <- sapply(starting:(starting + 999), eval_total, p, n, c0)
t2 <- proc.time()

time_tot <- t2 - t1

save(dip_res, time_tot, file = paste0("Simulation_results/global_null/p_",p,"_ratio_",ratio,"_c0_",c0,"_starting_",starting,".RData"))



