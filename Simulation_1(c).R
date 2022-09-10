Coeffmatrix<-expand.grid(p = c(6), ratio = c(1.5), c0 = 0.5, starting = seq(1, 9901, by = 100 ))
library(Matrix)
library(MASS)
library(QRM)
library(independencepvalue)
indica <- as.numeric(Sys.getenv('SGE_TASK_ID'))

p = Coeffmatrix$p[indica]
ratio = Coeffmatrix$ratio[indica]
n = p * ratio
a = 0.6
b = 0.3
c0 <- 0.5
starting = Coeffmatrix$starting[indica]


Sigma_11 <- QRM::equicorr(p/2, a)
Sigma_22 <- QRM::equicorr(p/2, a)

Sigma <- as.matrix(Matrix::bdiag(Sigma_11, Sigma_22))
Sigma[((p/2+1):p),(1:p/2)] <- b
Sigma[(1:p/2), ((p/2+1):p)] <- b


eval_total <- function(i0, p, n, c0, Sigma) {
  set.seed(i0)
  X <- mvrnorm(n = n, rep(0, p), Sigma)
  block_diag_structure <- block_diag(cor(X), c= c0)
  if(length(unique(block_diag_structure))> 1){
    set.seed(i0)
    k1 <- sample(unique(block_diag_structure), 1)
    set.seed(i0)
    t1 <- selective_p_val(S=cov(X), n=n, CP=block_diag_structure, c=c0, k=k1, d0=5, mc_iter=1000)
    set.seed(i0)
    t2 <- classical_p_val(S=cov(X), n=n, CP=block_diag_structure, k=k1, mc_iter=1000)
    p1 <- sum(block_diag_structure==k1)
    if(2*p1 < p){p1 <- p - sum(block_diag_structure==k1)}
    return(c(t1, t2, p1))
  }
  else{
    return(c(999, 999, 999))
  }
}

t1 <- proc.time()
dip_res <- sapply(starting:(starting + 99), eval_total, p, n, c0)
t2 <- proc.time()

time_tot <- t2 - t1

save(dip_res, time_tot, file = paste0("Simulation_results/histogram/p_",p,"_ratio_",ratio,"_c0_",c0,"_starting_",starting,".RData"))


