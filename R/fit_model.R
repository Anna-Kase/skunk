

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


#saveRDS(long_shot, "../intercept_only.rds")




# covariate model
model_type <- "covariate"
source("./R/create_m_n.R") 
source("./nimble/cov_model_1.R")
source("./R/prep_data.R")

longer_shot <- nimble::nimbleMCMC(
  code=cov_model,
  constants = constant_list,
  data = data_list,
  monitors = c("psi_beta", "rho_beta", "delta_beta", "gamma_beta", "phi_beta"),
  nburnin = 5000,
  nchains = 2,
  niter = 30000
)

# model started running at 8:20 am 29 Sept 2022
# 8:50 check-in -- model is ~1/8 through first chain
# 9:24 (I missed 9:20) ~1 hour and the model is just over
#     a quarter of the way through the first chain
# 10:22 model is one dash away from being halfway through the first chain
# 11:35 model is ~5/8 through first chain
# 12:23 first chain is 3 dashes away from complete
# 12:35 we made it to the second chain!
# 2:26 a little less than 3/8 the way through the second chain
# 3:20 we are halfway through the second chain!
# 4:34 - closed the laptop at ~3:50 to grab an Uber to get to the 
# 4:45 train. The model appears to have kept running and we are
#     at ~5/8 through the second chain
# 5:20 - model is just over 7/8 through the second chain
# 5:30 on the dot the model finished!



MCMCvis::MCMCsummary(
  longer_shot,
  digits=2
)


#                mean    sd    2.5%    50%  97.5% Rhat n.eff
# delta_beta[1] -1.500 0.330 -2.2000 -1.500 -0.880    1  4879
# delta_beta[2] -0.740 0.770 -2.4000 -0.710  0.660    1  3093
# gamma_beta[1] -2.400 0.210 -2.8000 -2.400 -2.000    1  5819
# gamma_beta[2] -0.230 0.140 -0.5200 -0.230  0.048    1  5990
# phi_beta[1]    0.094 0.250 -0.4000  0.095  0.580    1  4173
# phi_beta[2]   -0.450 0.230 -0.9000 -0.440 -0.011    1  3441
# psi_beta[1]   -1.600 0.320 -2.2000 -1.600 -0.990    1  5244
# psi_beta[2]   -0.570 0.240 -1.1000 -0.550 -0.120    1  4831
# rho_beta[1]   -2.000 0.061 -2.1000 -2.000 -1.900    1  7127
# rho_beta[2]    0.096 0.052 -0.0058  0.097  0.200    1  6902



# Model with 3 Covariates (urb, water_dist, open_dev) with NO interactions
# covariate model
model_type <- "covariate3"
source("./R/create_m_n.R") 
source("./nimble/three_cov_model.R")
source("./R/prep_data.R")

longest_shot <- nimble::nimbleMCMC(
  code=three_cov_model,
  constants = constant_list,
  data = data_list,
  monitors = c("psi_beta", "rho_beta", "delta_beta", "gamma_beta", "phi_beta"),
  nburnin = 5000,
  nchains = 2,
  niter = 30000
)




