
intercept_inits <- function(){
  to_return <- list(
    psi = rbeta(1,1,1),
    rho = rbeta(1,1,1),
    gamma = rbeta(1,1,1),
    phi = rbeta(1,1,1),
    delta = rbeta(1,1,1),
    z = matrix(
      1,
      ncol = constant_list$nseason,
      nrow = constant_list$nsite
    )
  )
  return(to_return)
}

covariate_inits <- function(){
  to_return <- list(
    psi_beta = rnorm(constant_list$ncovar_psi),
    rho_beta = rnorm(constant_list$ncovar_rho),
    gamma_beta = rnorm(constant_list$ncovar_gamma),
    phi_beta = rnorm(constant_list$ncovar_phi),
    delta_beta = rnorm(constant_list$ncovar_delta),
    z = matrix(
      1,
      ncol = constant_list$nseason,
      nrow = constant_list$nsite
    )
  )
  return(to_return)
}
