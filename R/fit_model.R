

library(nimble)
library(sf)
library(dplyr)
library(MCMCvis)

sf::sf_use_s2(FALSE)

source("./R/create_m_n.R")
source("./R/init_functions.R")

# intercept only model
model_type <- "intercept"
source("./nimble/intercept_only.R")
source("./R/prep_data.R")

long_shot <- nimble::nimbleMCMC(
  code=my_model,
  constants = constant_list,
  data = data_list,
  monitors = c("psi", "rho", "delta", "gamma", "phi"),
  inits = intercept_inits()
)


MCMCvis::MCMCsummary(
  long_shot,
  digits=2
)

# covariate model
model_type <- "covariate"
source("./R/create_m_n.R") 
source("./nimble/cov_model_1.R")
source("./R/prep_data.R")

nimbleOptions(verbose=FALSE)
longer_shot <- nimble::nimbleMCMC(
  code=cov_model,
  constants = constant_list,
  data = data_list,
  monitors = c("psi_beta", "rho_beta", "delta_beta", "gamma_beta", "phi_beta"),
  nburnin = 5000,
  nchains = 2,
  niter = 15000
)


MCMCvis::MCMCsummary(
  longer_shot,
  digits=2
)