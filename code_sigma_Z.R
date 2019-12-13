# Code for semiparametric IV project when Sigma is
# heteroskedastic and depends on Z.

library(numDeriv)
library(rootSolve)
library(mvtnorm)

# Compute the efficient score
efficient_score <- function(beta, g_hat_i, Z_i, T_i, Y_i,
                            sigma_sq_given_z_i, cov_given_z_i,
                            tau_sq_given_z_i){
  # Var(Y - beta T | Z)
  var_given_z_i = sigma_sq_given_z_i - 2*beta[1]*cov_given_z_i + beta[1]^2*tau_sq_given_z_i

  # If Var(Y - beta T | Z) evaluated at the data point is 0,
  # return (0, 0)
  if (var_given_z_i == 0) return(c(0,0))
  else {
    comp_1 = g_hat_i*(Y_i - beta[1]*T_i - beta[2])/var_given_z_i
    comp_2 = (Y_i - beta[1]*T_i - beta[2])/var_given_z_i
    res = c(comp_1, comp_2) # beta and intercept
    return(res)
  }
}

# Evaluate the estimating equation at beta = c(beta, beta_0)
est_equation <- function(beta, g_hat, Z, T, Y,
                         sigma_sq_given_z, cov_given_z,
                         tau_sq_given_z){

  n = length(Z)
  est_eqn = 0
  for (i in 1:n)
    est_eqn = est_eqn + efficient_score(beta, g_hat[i],
                                        Z[i], T[i], Y[i],
                                        sigma_sq_given_z[i],
                                        cov_given_z[i],
                                        tau_sq_given_z[i])

  return(est_eqn)

}

# Compute the robust variance estimator
sandwich_variance <- function(param_0, g_hat, Z, T, Y,
                              sigma_sq_given_z, cov_given_z,
                              tau_sq_given_z){

  A_m = 0
  B_m = 0
  m = dim(dataset)[1]
  for (i in 1:m){

    phi_i = efficient_score(param_0, g_hat[i], Z[i], T[i], Y[i],
             sigma_sq_given_z[i], cov_given_z[i], tau_sq_given_z[i])


    B_i = phi_i %*% t(phi_i)
    A_i = -jacobian(efficient_score, param_0, g_hat_i = g_hat[i],
                    Z_i = Z[i], T_i = T[i], Y_i = Y[i],
                    sigma_sq_given_z_i = sigma_sq_given_z[i],
                    cov_given_z_i = cov_given_z[i],
                    tau_sq_given_z_i = tau_sq_given_z[i])
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
# Estimate gamma(Z) and Sigma(Z) using NW-regression.

# Return the estimted treatment effect and its SE

estimate_beta_estimate_var <- function(dataset, start_0){
  colnames(dataset) <- c('Z', 'T', 'Y')
  dataset = dataset[order(dataset$Z),]

  n = dim(dataset)[1]
  T = dataset$T
  Y = dataset$Y
  Z = dataset$Z

  # Estimate E[T | Z], E[Y | Z], E[T^2 | Z], E[Y^2 | Z], E[TY | Z]
  T_given_Z = ksmooth(Z, T, n.points = n, x.points = Z)$y
  Y_given_Z = ksmooth(Z, Y, n.points = n, x.points = Z)$y
  T_sq_given_Z = ksmooth(Z, T^2, n.points = n, x.points = Z)$y
  Y_sq_given_Z = ksmooth(Z, Y^2, n.points = n, x.points = Z)$y
  TY_given_Z = ksmooth(Z, T*Y, n.points = n, x.points = Z)$y

  # Compute Var(T | Z), Var(Y | Z), Cov(TY | Z)
  tau_sq_given_z  = T_sq_given_Z - T_given_Z^2
  sigma_sq_given_z = Y_sq_given_Z - Y_given_Z^2
  cov_given_z = TY_given_Z - T_given_Z*Y_given_Z

  # Feed into the estimating equation and solve for (beta, beta_0)
  point_est = multiroot(est_equation, start = start_0, g_hat = T_given_Z,
                        Z = Z, T = T, Y = Y,
                        sigma_sq_given_z = sigma_sq_given_z,
                        cov_given_z = cov_given_z,
                        tau_sq_given_z = tau_sq_given_z)$root

  var_cov_matrix = sandwich_variance(point_est, T_given_Z,
                                     Z, T, Y,
                                     sigma_sq_given_z, cov_given_z,
                                     tau_sq_given_z)
  sd_beta = sqrt(var_cov_matrix[1,1])
  est_beta = point_est[1]
  return(c(est_beta, sd_beta))

}




