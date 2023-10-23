library(MASS)
library(independencepvalue)

source("Simulation_codes/Utils.R")

eval_total <- function(i0, p, n, c0) {
  Sigma <- diag(p)
  set.seed(i0)
  X <- MASS::mvrnorm(n=n, rep(0, p), Sigma)
  block_diag_structure <- independencepvalue::block_diag(cor(X), c=c0)
  if(length(unique(block_diag_structure))> 1){
    set.seed(i0)
    k1 <- sample(unique(block_diag_structure), 1)
    Sigma_hat <- matrix(a, p, p)
    diag(Sigma_hat) <- 1
    Sigma_hat[which(block_diag_structure == k1), which(block_diag_structure != k1)] <- 0
    Sigma_hat[which(block_diag_structure != k1), which(block_diag_structure == k1)] <- 0
    sip <- sapply(1:1000, MC_unconditional, Sigma_hat, block_diag_structure, k1, n, c0)
    test_hyp <- independencepvalue:::test_stat_CCA(cov(X), block_diag_structure, k1)
    zp <- sum(as.numeric(sip[2, ]))
    if(zp < 100){
      sip <- sapply(1:(min((1000 * 100/zp), 1e+05)), MC_unconditional, Sigma_hat, block_diag_structure, k1, n, c0)
    }
    t1 <- mean(test_hyp$statistic >= sip[1, ][sip[2, ] == TRUE])
    Sigma_hat <- diag(p)
    sip <- sapply(1:1000, MC_unconditional, Sigma_hat, block_diag_structure, k1, n, c0)
    test_hyp <- independencepvalue:::test_stat_CCA(cov(X), block_diag_structure, k1)
    zp <- sum(as.numeric(sip[2, ]))
    if(zp < 100){
      sip <- sapply(1:(min((1000 * 100/zp), 1e+05)), MC_unconditional, Sigma_hat, block_diag_structure, k1, n, c0)
    }
    t2 <- mean(test_hyp$statistic >= sip[1, ][sip[2, ] == TRUE])
    return(c(t1, t2))
  }
  else{
    return(c(999, 999))
  }
}


p <- 10
ratio <- 1.5
n <- p * ratio
c0 <- sqrt(log(p)/n)
a <- 0.99

result <- t(sapply(1:10000, eval_total, p, n, c0))

save(result, file=paste0("Simulation_results/histogram_S4_global_null_p_", p, "_n_", n, "_c0_", c0, ".RData"))

