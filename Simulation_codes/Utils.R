library(amen)
library(MASS)
library(Matrix)
library(independencepvalue)

# Given a sample covariance matrix S and a group of variables, returns the cross-covariance matrix between the whitened variables and the test statistics as a function of canonical correlations.
# Modify independencepvalue:::test_stat_CCA() to return to desired components. 
# Reuse the code from independencepvalue:::test_stat_CCA().
# Please refer to the documentation of independencepvalue:::test_stat_CCA() for details.
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

# Function for Monte Carlo simulation for selective inference, when we only condition on the diagonal blocks. 
# input: i, i0 = Integers, used to set the seed and keep track of the iterations, p = # of variables, p1 = # of nonzero canonical correlations,
#        n = # of samples, S11 and S22 = two fixed diagonal blocks, c0 = threshold. 
# output: A list containing, a) test statistic corresponding to the simulation, and b) a logical scalar indicating if the simulation resulted in satisfying the selection event.
sample_simul <- function(i, i0, p, p1, n, S11, S22, c0){
  set.seed(i+i0)
  X <- MASS::mvrnorm(n = n, rep(0, p), diag(p))
  test_inf <- test_stat_CCA_new(cov(X), c(rep(1, p1), rep(2, p-p1)), 1)#simulate the whitened cross-covariance from identity, since it is independent of the diagonal blocks under the null (Lemma 3 of the Supplementary Materials S11.2)
  test_info <- list()
  test_info[[1]] <- test_inf$statistic
  S_X_new_cross <- amen::mhalf(as.matrix(S11)) %*% test_inf$S12w %*% amen::mhalf(as.matrix(S22))#new S_12
  S_X_new <- as.matrix(Matrix::bdiag(S11, S22))
  S_X_new[(1:p1),(p1+1):p] <- S_X_new_cross
  S_X_new[(p1+1):p,(1:p1)] <- t(S_X_new_cross)#new S
  cor_X_new <- solve(diag(sqrt(diag(S_X_new)))) %*% S_X_new %*% solve(diag(sqrt(diag(S_X_new))))
  test_info[[2]] <- as.numeric(all(abs(cor_X_new[(1:p1),(p1+1):p]) < c0))#check if this new S satisfies the selection event
  return(test_info)
}

# Function for Monte Carlo simulation for selective inference, given specific values of the nuisance parameters in S11.1. 
# input: i, i0 = Integers, Sigma_hat = a block diagonal matrix with the nuisance parameters as the diagonal blocks, 
#        CP = A vector indicating the group membership corresponding to the block diagnal structure of Sigma_hat, k = the group to be tested for independence with the remaining variables, 
#        n = # of samples, c0 = threshold. 
# output: A list containing, a) test statistic corresponding to the simulation, and b) a logical scalar indicating if the simulation resulted in satisfying the selection event.
MC_unconditional <- function(i, Sigma_hat, CP, k, n, c0){
  p <- nrow(Sigma_hat)
  set.seed(i)
  X_hat <- MASS::mvrnorm(n = n, rep(0, p), Sigma_hat)#simulate data using Sigma_hat as the population covariance matrix
  test_stat <- independencepvalue:::test_stat_CCA(cov(X_hat), CP, k)$statistic#compute the test statistic for the simulated data
  check <- all(abs(cor(X_hat)[which(CP == k), which(CP != k)]) <= c0)#check if the simulated data satisfies the selection event
  return(list(statistic = test_stat, status = check))
}