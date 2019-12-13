# Code for semiparametric IV project when Sigma is
# homoskedastic and does not depend on Z.

library(numDeriv)
library(rootSolve)
library(mvtnorm)

# Compute the efficient score
efficient_score_2 <- function(beta, g_hat_i, Z_i, T_i, Y_i){
  comp_1 = g_hat_i*(Y_i - beta[1]*T_i - beta[2])
  comp_2 = Y_i - beta[1]*T_i - beta[2]
  res = c(comp_1, comp_2)
  return(res)
}

# Evalute the estimating equation at beta = c(beta, beta_0)
est_equation_2 <- function(beta, g_hat, Z, T, Y){
  n = length(Z)
  est_eqn = 0
  for (i in 1:n)
    est_eqn = est_eqn + efficient_score_2(beta, g_hat[i],
                                          Z[i], T[i], Y[i])
  return(est_eqn)

}

# compute the robust sanwich estimator of variance
sandwich_variance_2 <- function(param_0, g_hat, Z, T, Y){

  A_m = 0
  B_m = 0
  m = length(Z)
  for (i in 1:m){

    phi_i = efficient_score_2(param_0, g_hat[i], Z[i], T[i], Y[i])


    B_i = phi_i %*% t(phi_i)
    A_i = -jacobian(efficient_score_2, param_0, g_hat_i = g_hat[i],
                    Z_i = Z[i], T_i = T[i], Y_i = Y[i])
    A_m = A_m + A_i
    B_m = B_m + B_i
  }
  A_m = A_m/m
  B_m = B_m/m
  var_cov_mat = (solve(A_m)%*%B_m%*%(t(solve(A_m))))/m
  return(var_cov_mat)
}


##############################################################################
# No baseline covariates.
# Estimate beta and beta_0.
# Estimate gamma(Z) using NW-regression.

estimate_beta <- function(dataset, start_0){
  colnames(dataset) <- c('Z', 'T', 'Y')
  dataset = dataset[order(dataset$Z),]

  n = dim(dataset)[1]
  T = dataset$T
  Y = dataset$Y
  Z = dataset$Z

  # Estimate E[T | Z], E[Y | Z], E[T^2 | Z], E[Y^2 | Z], E[TY | Z]
  T_given_Z = ksmooth(Z, T, n.points = n, x.points = Z)$y

  # Feed into the estimating equation and solve for (beta, beta_0)
  point_est = multiroot(est_equation_2, start = start_0, g_hat = T_given_Z,
                        Z = Z, T = T, Y = Y)$root

  var_cov_matrix = sandwich_variance_2(point_est, T_given_Z,
                                       Z, T, Y)
  sd_beta = sqrt(var_cov_matrix[1,1])
  est_beta = point_est[1]
  return(c(est_beta, sd_beta))

}
