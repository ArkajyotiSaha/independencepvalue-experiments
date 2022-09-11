library(Matrix)
library(MASS)
library(independencepvalue)
library(plot.matrix)
library(QRM)


p <- 6
n <- 9
a <- 0.6
b <- 0.3
c0 <- 0.5

Sigma_11 <- QRM::equicorr(p/2, a)
Sigma_22 <- QRM::equicorr(p/2, a)

Sigma <- as.matrix(Matrix::bdiag(Sigma_11, Sigma_22))
Sigma[((p/2+1):p), (1:p/2)] <- b
Sigma[(1:p/2), ((p/2+1):p)] <- b


eval_total <- function(i0, p, n, c0, Sigma) {
  set.seed(i0)
  X <- MASS::mvrnorm(n=n, rep(0, p), Sigma)
  block_diag_structure <- independencepvalue::block_diag(cor(X), c=c0)
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
i0 <- 9768
eval_total(i0, p, n, c0, Sigma)

set.seed(i0)
X <- mvrnorm(n=n, rep(0, p), Sigma)
png(file="Figures/Figure1ab.png",
    width=1600, height=800)
par(fig=c(0.01, 0.5, 0, 0.8))
plot(Sigma, breaks=c(0, 1), main="(a) Absolute population correlation", xlab=NA, ylab=NA, col=rev(heat.colors(10)), cex.lab=1.5, cex.axis=2.5, cex.main=3.75, cex.sub=2.5, key=NULL)
par(fig=c(0.51, 1, 0, 0.8),new=TRUE)
plot(abs(cor(X)), breaks=c(0, 1), main="(b) Absolute sample correlation", xlab=NA, ylab=NA, col=rev(heat.colors(10)), cex.lab=1.5, cex.axis=2.5, cex.main=3.75, cex.sub=2.5, key=NULL)
par(fig=c(0, 1, 0, 1),new=TRUE)
rect(
  head(seq(0.85, 5.85, 5/10), -1),
  6.25,
  tail(seq(0.85, 5.85, 5/10), -1),
  6.5,
  col=rev(heat.colors(10))
)
mtext((1:10)/10, side=3, at=tail(seq(0.85, 5.85, 5/10), -1)-0.25, cex=2.5)
dev.off()



