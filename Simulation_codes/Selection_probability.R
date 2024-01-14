library(rcdd)
library(volesti)
library(geometry)
library(CholWishart)
library(independencepvalue)
library(SimplicialCubature)
library(future.apply)

# Compute the probability of the selection event using numerical integration. 
# This is the denominator of the ratio computed in independencepvalue:::selective_p_val_integrate().
# Reuse the code from independencepvalue:::selective_p_val_integrate().
# Please refer to the documentation of independencepvalue:::selective_p_val_integrate() for details.
selection_integrate <- function(n, L, g, test_hyp, tol = 1e-05, maxeval = 1e4) {
  p1 <- nrow(test_hyp$S11)
  p2 <- nrow(test_hyp$S22)
  p <- p1 + p2
  
  du <- 0 # initialize du
  P <- rcdd::makeH(L, g, x = NULL) # create the Half space representation of the polytope
  PV_d <- rcdd::scdd(P) # compute the convex hull representation of the polytope.
  V_d <- as.matrix(PV_d$output[,-c(1, 2)])
  Pi <- volesti::Vpolytope(V = V_d) # create a convenient Vertex representation from the convex hull. 
  triang = try(geometry::delaunayn(Pi@V), silent = TRUE) # try the Delaunay triangulation
  if (!inherits(triang, 'try-error')) { # if the triangulation is successful
    prod_res <- CholWishart::lmvgamma(n/2, p2) - ((CholWishart::lmvgamma((n-p+p2)/2, p2) + CholWishart::lmvgamma(p2/2, p) + CholWishart::lmvgamma((p - p2)/2, p)))
    alpha <- log(pi ^ (p2 / 2) * 2 ^ p) + prod_res # compute the constant
    f_tot <- function(x) {
      dt <- independencepvalue:::dCCA(n, p, p2, alpha, x)
      return(dt) # compute the function in the numerator
    }
    part_int <- function(i, Pi, triang, f_tot) { #code to perform integration
      if (stats::var(round(Pi@V[triang[i, ],][, 1], 5)) == 0) {
        Pi@V[triang[i, ],][, 1][length(Pi@V[triang[i, ],][, 1])] <-
          Pi@V[triang[i, ], ][, 1][length(Pi@V[triang[i, ],][, 1])] * (1 - 10 ^ (-5)) # Account for numerical instabilities on the triangulation by adding negligble deviation
      }
      return(SimplicialCubature::adaptIntegrateSimplex(f_tot, t(Pi@V[triang[i, ],]), fDim = 2, tol = tol, maxEvals = maxeval)$integral) # perform numerical integration on the simplices
    }
    par_I_tot_list <- future.apply::future_sapply(1:nrow(triang), part_int, Pi, triang, f_tot, future.seed = TRUE)
    du <- sum(par_I_tot_list[1, ])   
    if (sum(par_I_tot_list[1, ]) == 0) { # if the integral in the denominator is effectively zero, set du to 0.
      du <- 0
    }
  }
  return(du)
}

# Compute the probability of the selection event directly from Beta distribution. 
# This is the denominator of the ratio computed in independencepvalue:::selective_p_val_beta().
# Reuse the code from independencepvalue:::selective_p_val_beta().
# Please refer to the documentation of independencepvalue:::selective_p_val_beta() for details.
selection_beta <- function(S, CP, k, n, c, test_hyp) {
  p <- nrow(S)
  diag_S <- diag(1 / sqrt(diag(S)))
  R <- diag_S %*% S %*% diag_S
  g_u <- min(1, c ^ 2 * (1 - test_hyp$statistic) / max(abs(R[CP == k, CP != k])) ^ 2)
  I_denom_tot <- stats::pbeta(c(0, g_u), (p - 1) / 2, (n - p) / 2)
  return((I_denom_tot[2] - I_denom_tot[1]))
}

# Compute the probability of the selection event using Monte Carlo. 
# This is the denominator of the ratio computed in independencepvalue:::selective_p_val_MC().
# Reuse the code from independencepvalue:::selective_p_val_MC().
# Please refer to the documentation of independencepvalue:::selective_p_val_MC() for details.
selection_MC <- function(n, L, g, test_hyp, mc_iter) {
  p1 <- nrow(test_hyp$S11)
  p2 <- nrow(test_hyp$S22)
  p <- p1 + p2
  sip <- future.apply::future_sapply(1:mc_iter, function(i)
      independencepvalue:::MC_function_selective(p, p2, n, L, g), future.seed = TRUE)
  zp <- sum(as.numeric(sip[2, ]))
  if (zp < 100) {
    sip <- future.apply::future_sapply(1:(min((mc_iter * 100 / zp), 100000)), function(i)
        independencepvalue:::MC_function_selective(p, p2, n, L, g), future.seed = TRUE)
  }
  return( mean(sip[2, ] == TRUE))
}

# Compute the probability of the selection event in the selective p value
selection <- function(S, CP, k, n, c, d0 = 5, tol = 1e-05, maxeval = 1e5, mc_iter = 1000){
  test_hyp <- independencepvalue:::test_stat_CCA(S, CP, k)
  p1 <- nrow(test_hyp$S11)
  p2 <- nrow(test_hyp$S22)
  if (p2 == 1 & p2 <= d0) {
    du <- selection_beta(S, CP, k, n, c, test_hyp)
  }
  else {
    du <- 0
    L <- independencepvalue:::form_L(test_hyp)
    g <- c(rep(c, 2 * p1 * p2), rep(0, p2), rep(1, p2))
    if (p2 <= d0) {
      du <- selection_integrate(n, L, g, test_hyp, tol, maxeval)
    }
    if (du[1] <= 0 || du[1] >= 1 || p2 > d0) {
      # use Monte Carlo approach
      du <- selection_MC(n, L, g, test_hyp, mc_iter)
    }
  }
  return(du)
}







