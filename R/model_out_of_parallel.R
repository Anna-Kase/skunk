
# packages
library(nimble)

# source data prep and model
source("./R/parallel_data_prep.R") 
source("./nimble/three_cov_model.R")

my_inits <- function(){
  to_return <- list(
    delta_beta = rnorm(constant_list$ncovar_delta),
    gamma_beta = rnorm(constant_list$ncovar_gamma),
    phi_beta = rnorm(constant_list$ncovar_phi),
    psi_beta = rnorm(constant_list$ncovar_psi),
    rho_beta = rnorm(constant_list$ncovar_rho),
    z = matrix(
      1,
      ncol = constant_list$nseason,
      nrow = constant_list$nsite
    )
  )
  return(to_return)
}


# fit model
longest_shot <- nimble::nimbleMCMC(
  code=three_cov_model,
  constants = constant_list,
  data = data_list,
  monitors = c("psi_beta", "rho_beta", "delta_beta", "gamma_beta", "phi_beta"),
  inits = my_inits,
  nburnin = 5000,
  nchains = 1,
  niter = 20000
)



plot(longest_shot[,"psi_beta[1]"], type="l")
plot(longest_shot[,"rho_beta[1]"], type="l")
plot(longest_shot[,"phi_beta[1]"], type="l")
plot(longest_shot[,"gamma_beta[1]"], type="l")
plot(longest_shot[,"delta_beta[1]"], type="l")
