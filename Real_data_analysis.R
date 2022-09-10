Coeffmatrix<-expand.grid(network = c(3, 4), p = c(100), ratio = c(1.2, 1.4, 1.5, 2), sigma = c(0, 0.5, 3), starting = seq(1, 10000, by = 1000))
library(independencepvalue)
indica <- as.numeric(Sys.getenv('SGE_TASK_ID'))
network <- Coeffmatrix$network[indica]
p <- Coeffmatrix$p[indica]
ratio <- Coeffmatrix$ratio[indica]
sigma <- Coeffmatrix$sigma[indica]
starting <- Coeffmatrix$starting[indica]
n <- p * ratio


library(data.table)
fip <- fread(paste0("DREAM5_data/gold_standard_edges_only/DREAM5_NetworkInference_Edges_Network",network,".tsv"))
dip <- fread(paste0("DREAM5_data/Training/net",network,"_expression_data.tsv"))
add_mat <- matrix(0, ncol(dip), ncol(dip))
fip1 <- as.numeric(gsub(".*?([0-9]+).*", "\\1", fip[,1]$V1))
fip2 <- as.numeric(gsub(".*?([0-9]+).*", "\\1", fip[,2]$V2))

for(i in 1:nrow(fip)){
  add_mat[fip1[i], fip2[i]] <- fip[,3]$V3[i]
  add_mat[fip2[i], fip1[i]] <- fip[,3]$V3[i]
}

library(Corbi)
library(igraph)
dip <- as.matrix(dip)

effect_size <- function(Sigma, CP, k){
  S <- Sigma
  p <- nrow(S)
  ptemp <- sum(CP==k)
  if(2*ptemp >= p){
    p1 <- ptemp
    S11_x <- as.matrix(S[which(CP==k), which(CP==k)])#S11
    S22_x <- as.matrix(S[which(CP!=k), which(CP!=k)])#S22
    S12_x <- as.matrix(S[which(CP==k), which(CP!=k)])#S_12
  }
  if(2*ptemp < p){
    p1 = p - ptemp
    S11_x <- as.matrix(S[which(CP!=k), which(CP!=k)])#S11
    S22_x <- as.matrix(S[which(CP==k), which(CP==k)])#S22
    S12_x <- as.matrix(S[which(CP!=k), which(CP==k)])#S_12
  }
  #ensures that p2 <= p1.
  p2 <- p - p1
  if(p1 > p2){rp <- p2}#Define r(P)
  if(p1 <= p2){rp <- p1}
  inv_S11_x <- solve(S11_x)#inv(S11)
  inv_S11_x_half <- amen::mhalf(inv_S11_x)#S_11^{-1/2}
  inv_S22_x <- solve(S22_x)#inv(S22)
  inv_S22_x_half <- amen::mhalf(inv_S22_x)#S_22^{-1/2}
  tilde_S12_x <- inv_S11_x_half %*% S12_x %*% inv_S22_x_half#S_12^W = covariance matrix of whitened X_1 and X_2
  svdecom <- svd(tilde_S12_x)#compact SVD
  singular_values <- svdecom$d#Lambda
  test_stat <- prod(1-singular_values^2)
  return(test_stat)
}


eval_total <- function(i0, dip, add_mat, p, n) {
  length_list <- sort(rowSums(add_mat), decreasing = T)
  length_list_index <- order(rowSums(add_mat), decreasing = T)
  select_index_row <- length_list_index[1]
  set.seed(i0+1)
  select_index <- c(sample(which(add_mat[select_index_row,] == 1), (p - 1)), select_index_row)
  
  
  #choose n
  set.seed(i0+2)
  row_index <- sample(1:nrow(dip), n)
  
  row_index_complmnt <- setdiff(1:nrow(dip), row_index)
  
  
  #Reference data
  
  dip_short <- dip[row_index_complmnt,select_index]
  X_ref <- dip_short
  tx_ref <- X_ref
  set.seed(i0+3)
  X_ref <- tx_ref + sigma * sqrt(diag(cov(tx_ref)))* matrix(rnorm(length(row_index_complmnt) * length(select_index)), length(row_index_complmnt), length(select_index))
  Sigma <- cov(X_ref)
  
  Var_S <- Sigma
  diag_S <- diag(1/sqrt(diag(Var_S)))
  Cor_S <- diag_S %*% Var_S %*% diag_S
  dis_S <- 1 - abs(Cor_S)
  test <- as.dist(dis_S, diag = TRUE)
  clust_result <- hclust(test, method = "single")
  c0 <- 1 - clust_result$height[p - 1] + sqrt(log(p)/n)
  c0 <- min(c0, 0.99999999)
  
  ##t data
  
  dip_short <- dip[row_index,select_index]
  X <- dip_short
  tx <- X
  
  set.seed(i0+4)
  X <- tx + sigma * sqrt(diag(cov(tx)))* matrix(rnorm(length(row_index) * length(select_index)), length(row_index), length(select_index))
  block_diag_structure <- block_diag(cor(X), c= c0)
  if(length(unique(block_diag_structure))> 1){
    set.seed(i0)
    k1 <- sample(unique(block_diag_structure), 1)
    set.seed(i0)
    t1 <- selective_p_val(S=cov(X), n=n, CP=block_diag_structure, c=c0, k=k1, d0=5, mc_iter=1000)
    set.seed(i0)
    t2 <- classical_p_val(S=cov(X), n=n, CP=block_diag_structure, k=k1, mc_iter=1000)
    t3 <- effect_size(Sigma, CP=block_diag_structure, k=k1)
    p1 <- sum(block_diag_structure==k1)
    if(2*p1 < p){p1 <- p - sum(block_diag_structure==k1)}
    return(c(t1, t2, t3, p1))
  }
  else{
    return(c(999, 999, 999, 999))
  }
}


t1 <- proc.time()
dip_res <- sapply(starting:(starting + 999), eval_total, dip, add_mat, p, n)
t2 <- proc.time()

time_tot <- t2 - t1

save(dip_res, time_tot, file = paste0("Real_data_results/network_",network,"p_",p,"_ratio_",ratio,"_noise_",sigma,"_starting_",starting,".RData"))


