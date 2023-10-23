library(amen)

effect_size <- function(Sigma, CP, k){
  S <- Sigma
  p <- nrow(S)
  ptemp <- sum(CP==k)
  if(2*ptemp>=p){
    p1 <- ptemp
    S11_x <- as.matrix(S[which(CP==k), which(CP==k)])#S11
    S22_x <- as.matrix(S[which(CP!=k), which(CP!=k)])#S22
    S12_x <- as.matrix(S[which(CP==k), which(CP!=k)])#S_12
  }
  if(2*ptemp < p){
    p1 <- p - ptemp
    S11_x <- as.matrix(S[which(CP!=k), which(CP!=k)])#S11
    S22_x <- as.matrix(S[which(CP==k), which(CP==k)])#S22
    S12_x <- as.matrix(S[which(CP!=k), which(CP==k)])#S_12
  }
  #ensures that p2<=p1.
  p2 <- p - p1
  if(p1 > p2){rp <- p2}#Define r(P)
  if(p1<=p2){rp <- p1}
  inv_S11_x <- solve(S11_x)#inv(S11)
  inv_S11_x_half <- amen::mhalf(inv_S11_x)#S_11^{-1/2}
  inv_S22_x <- solve(S22_x)#inv(S22)
  inv_S22_x_half <- amen::mhalf(inv_S22_x)#S_22^{-1/2}
  tilde_S12_x <- inv_S11_x_half %*% S12_x %*% inv_S22_x_half#S_12^W=covariance matrix of whitened X_1 and X_2
  svdecom <- svd(tilde_S12_x)#compact SVD
  singular_values <- svdecom$d#Lambda
  test_stat <- prod(1-singular_values^2)
  return(test_stat)
}