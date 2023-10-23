library(independencepvalue)

source("Simulation_codes/Effect_size.R")

eval_total <- function(i0, dip, add_mat, p, n, eta) {
  length_list <- sort(rowSums(add_mat), decreasing=T)
  length_list_index <- order(rowSums(add_mat), decreasing=T)
  select_index_row <- length_list_index[1]
  set.seed(i0+1)
  select_index <- c(sample(which(add_mat[select_index_row,]==1), (p - 1)), select_index_row)
  
  
  #choose n
  set.seed(i0+2)
  row_index <- sample(1:nrow(dip), n)
  row_index_complmnt <- setdiff(1:nrow(dip), row_index)
  
  
  #Reference data
  dip_short <- dip[row_index_complmnt, select_index]
  X_ref <- dip_short
  tx_ref <- X_ref
  set.seed(i0+3)
  X_ref <- tx_ref + eta * sqrt(diag(cov(tx_ref)))* matrix(rnorm(length(row_index_complmnt) * length(select_index)), length(row_index_complmnt), length(select_index))
  Sigma <- cov(X_ref)
  
  Var_S <- Sigma
  diag_S <- diag(1/sqrt(diag(Var_S)))
  Cor_S <- diag_S %*% Var_S %*% diag_S
  dis_S <- 1 - abs(Cor_S)
  test <- as.dist(dis_S, diag=TRUE)
  clust_result <- hclust(test, method="single")
  c0 <- 1 - clust_result$height[p - 1] + sqrt(log(p)/n)
  c0 <- min(c0, 0.99999999)
  
  dip_short <- dip[row_index, select_index]
  X <- dip_short
  tx <- X
  
  set.seed(i0+4)
  X <- tx + eta * sqrt(diag(cov(tx)))* matrix(rnorm(length(row_index) * length(select_index)), length(row_index), length(select_index))
  block_diag_structure <- independencepvalue::block_diag(cor(X), c=c0)
  if(length(unique(block_diag_structure))> 1){
    set.seed(i0)
    k1 <- sample(unique(block_diag_structure), 1)
    set.seed(i0)
    t1 <- independencepvalue::selective_p_val(S=cov(X), n=n, CP=block_diag_structure, c=c0, k=k1, d0=5, mc_iter=1000, maxeval = 10000)
    set.seed(i0)
    t2 <- independencepvalue::classical_p_val(S=cov(X), n=n, CP=block_diag_structure, k=k1, mc_iter=1000)
    t3 <- effect_size(Sigma, CP=block_diag_structure, k=k1)
    return(c(t1, t2, t3))
  }
  else{
    return(c(999, 999, 999))
  }
}
p <- 100
for(network in c(3:4)){
  fip <- data.table::fread(paste0("Real_data_results/DREAM5_data/gold_standard_edges_only/DREAM5_NetworkInference_Edges_Network", network, ".tsv"))
  dip <- data.table::fread(paste0("Real_data_results/DREAM5_data/Training/net", network, "_expression_data.tsv"))
  add_mat <- matrix(0, ncol(dip), ncol(dip))
  fip1 <- as.numeric(gsub(".*?([0-9]+).*", "\\1", fip[, 1]$V1))
  fip2 <- as.numeric(gsub(".*?([0-9]+).*", "\\1", fip[, 2]$V2))
  
  for(i in 1:nrow(fip)){
    add_mat[fip1[i], fip2[i]] <- fip[, 3]$V3[i]
    add_mat[fip2[i], fip1[i]] <- fip[, 3]$V3[i]
  }
  dip <- as.matrix(dip)
  for(ratio in c(1.2, 1.4, 1.5, 2)){
    n <- p * ratio
    for(eta in c(0, 0.5, 3)){
      result <- t(sapply(1:10000, eval_total, dip, add_mat, p, n, eta))
      save(result, file = paste0("Real_data_results/network_", network, "p_", p, "_n_", n, "_eta_", eta, ".RData"))
    }
  }
}


