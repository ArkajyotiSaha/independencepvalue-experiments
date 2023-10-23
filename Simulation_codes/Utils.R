library(amen)
library(MASS)
library(Matrix)
library(independencepvalue)

test_stat_CCA_new <- function (S, CP, k)
{
  p <- nrow(S)
  ptemp <- sum(CP == k)
  if (2 * ptemp >= p) {
    S11_x <- as.matrix(S[which(CP == k), which(CP == k)])
    S22_x <- as.matrix(S[which(CP != k), which(CP != k)])
    S12_x <- as.matrix(S[which(CP == k), which(CP != k)])
  }
  if (2 * ptemp < p) {
    S11_x <- as.matrix(S[which(CP != k), which(CP != k)])
    S22_x <- as.matrix(S[which(CP == k), which(CP == k)])
    S12_x <- as.matrix(S[which(CP != k), which(CP == k)])
  }
  S11_x_half <- amen::mhalf(S11_x)
  inv_S11_x_half <- solve(S11_x_half)
  S22_x_half <- amen::mhalf(S22_x)
  inv_S22_x_half <- solve(S22_x_half)
  tilde_S12_x <- inv_S11_x_half %*% S12_x %*% inv_S22_x_half
  svdecom <- svd(tilde_S12_x)
  singular_values <- svdecom$d
  test_stat <- prod(1 - singular_values^2)
  return(list(statistic = test_stat, S12w = tilde_S12_x))
}


sample_simul <- function(i, i0, p, p1, n, S11, S22, c0){
  set.seed(i+i0)
  X <- MASS::mvrnorm(n=n, rep(0, p), diag(p))
  test_inf <- test_stat_CCA_new(cov(X), c(rep(1, p1), rep(2, p-p1)), 1)
  test_info <- list()
  test_info[[1]] <- test_inf$statistic
  S_X_new_cross <- amen::mhalf(as.matrix(S11)) %*% test_inf$S12w %*% amen::mhalf(as.matrix(S22))
  S_X_new <- as.matrix(Matrix::bdiag(S11, S22))
  S_X_new[(1:p1),(p1+1):p] <- S_X_new_cross
  S_X_new[(p1+1):p,(1:p1)] <- t(S_X_new_cross)
  cor_X_new <- solve(diag(sqrt(diag(S_X_new)))) %*% S_X_new %*% solve(diag(sqrt(diag(S_X_new))))
  test_info[[2]] <- as.numeric(all(abs(cor_X_new[(1:p1),(p1+1):p]) < c0))
  return(test_info)
}

MC_unconditional <- function(i, Sigma_hat, CP, k, n, c0){
  p <- nrow(Sigma_hat)
  set.seed(i)
  X_hat <- MASS::mvrnorm(n=n, rep(0, p), Sigma_hat)
  test_stat <- independencepvalue:::test_stat_CCA(cov(X_hat), CP, k)$statistic
  check <- all(abs(cor(X_hat)[which(CP == k), which(CP != k)]) <= c0)
  return(list(statistic = test_stat, status = check))
}