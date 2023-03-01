

library(dplyr)


output <- readRDS("./skunk_rds/spatial_covariates_fall2.RDS")

output <- do.call("rbind", output)

head(output)

# still subsampling the output because the object size is too large for
# R to handle otherwise

set.seed(1995)

my_samples2 <- sample(
  1:nrow(output),
  2500
)

cov_mod_sub <- output[my_samples2,]



# "Extract" MCMC posterior
source("./R/mcmc_functions.R")
mc <- split_mcmc(cov_mod_sub)

head(mc)

source("./R/spatial_covariates_fall_data_prep.R")

z <- z_prob <- array(
  NA, 
  dim=c(
    length(my_samples2),
    constant_list$nsite, 
    constant_list$nseason
  )
)

for(i in 1:nsite){
  if(data_list$y[i,1] > 0 & !is.na(data_list$y[i,1])){
    z_prob[,i,1] <- 1
  } else{
    psi <- plogis(mc$psi_beta %*% constant_list$X_psi[i,])
    rho <- plogis(mc$rho_beta %*% constant_list$X_rho[i,])
    num <- psi * (1-rho)^constant_list$J[i,1]
    den <- num + (1-psi)
    z_prob[,i,1] <- num/den
  }
}

z[,,1] <- rbinom(
  prod(dim(z)[1:2]),
  size = 1,
  prob = z_prob[,,1]
)


delta_bar <- zeta <- array(dim=dim(z_prob))

# Create d_vec (aka neighborhood colonization covariate) matrix

d_vec <- array(NA, dim=c(length(my_samples2), nsite, nsite, nseason))
for(i in 1:nsite){
  for(ti in 1:nseason){
    d_vec[,i,,ti] <- (
      mc$delta_beta %*%
        t(constant_list$delta_array[i,,])) + 
      as.numeric(mc$delta_fall) * constant_list$season_vec[ti]
  }
}
d_vec <- plogis(d_vec)

# Now fill in zm, zeta, and delta_bar
pb <- txtProgressBar(max = nsite)   # Setting up progress bar
for(i in 1:nsite){
  setTxtProgressBar(pb,i)
  # multiplying our z array (every MCMC step, for every individual site, at t=1),
  # by our m matrix (how many "nearby" neighbors each site has) to fill in
  # our zm array
  zm <- sweep(             
    z[,1:nsite,1],
    2,
    constant_list$m[i, 1:nsite],
    FUN = "*"
  )
  # These lines creates the tmp_zeta object that holds the total number of
  # neighbors for each site. Then it creates the zeta array which holds
  # the number of neighbors but as a numeric variable rather than a factor,
  # then we work through multiplying our zm matrix (aka the nearby neighbors)
  # by the log(1-d_vec) (our nearby colonization covariate), then summing
  # those probabilities, convert the log probabilities back to 
  # numeric probabilities, and then fill those into the delta_bar array
  tmp_zeta <- rowSums(zm)
  zeta[,i,1] <- as.numeric(tmp_zeta>0)
  zm <- zm * log(1 - d_vec[,i,,1])
  zm <- rowSums(zm)
  zm <- 1 - exp(zm)
  delta_bar[,i,1] <- zm
}

# Now keep filling in the z_prob matrix for the rest of the seasons
# where if the species was not observed then the probability of z is 
# hard coded to be 1, but if the species was observed then we calculate
# that value based on the linear predictors from the dynamic occupancy
# model formula: 
# (z[i,t] * theta[i,t]) + ((1 - z[i,t]) * I) + ((1 - z[i,t]) * (1 - I) * gamma)
pb <- txtProgressBar(max = nseason, min = 2)
for(t in 2:nseason){
  setTxtProgressBar(pb, t)
  for(i in 1:nsite){
    
    if(data_list$y[i,t] > 0 & !is.na(data_list$y[i,t])){
      z_prob[,i,t] <- 1
    } else{
      phi <- plogis(mc$phi_beta %*% constant_list$X_phi[i,])
      gamma <- plogis(mc$gamma_beta %*% constant_list$X_gamma[i,] +
                        mc$gamma_fall * constant_list$season_vec[t])
      rho <- plogis(mc$rho_beta %*% constant_list$X_rho[i,])
      num <- (z[,i,t-1] * phi) + 
        ((1-z[,i,t-1]) * zeta[,i,t-1] * delta_bar[,i,t-1]) +
        ((1-z[,i,t-1]) * (1-zeta[,i,t-1])*gamma)
      dnm1 <- 1 - num
      # add rho
      num <- num * (1 - rho)^constant_list$J[i,t]
      z_prob[,i,t] <- num / (num + dnm1)
    }
    # Then following the calculation structure from the first season above, 
    # we fill in the z matrix using the z probabilities we just calculated
  }
  z[,,t] <- rbinom(
    prod(dim(z)[1:2]),
    size = 1,
    prob = z_prob[,,t]
  )
  # Then create the zm matrix with the neighborhood colonization to 
  # use to calculate delta_bar for all of the seasons
  for(i in 1:nsite){
    zm <- sweep(
      z[,1:nsite,t],
      2,
      constant_list$m[i, 1:nsite],
      FUN = "*"
    )
    # Now calculate and fill in delta_bar for all of the seasons
    tmp_zeta <- rowSums(zm)
    zeta[,i,t] <- as.numeric(tmp_zeta>0)
    zm <- zm * log(1 - d_vec[,i,,t])
    zm <- rowSums(zm)
    zm <- 1 - exp(zm)
    delta_bar[,i,t] <- zm
  }
}





# parameter estimations

# psi = z*phi + (1-z)*I*d_vec +
# (1-z)*(1-I)*gamma

psi_est <- (z[,i,t-1] * phi) + 
  ((1-z[,i,t-1]) * zeta[,i,t-1] * delta_bar[,i,t-1]) +
  ((1-z[,i,t-1]) * (1-zeta[,i,t-1])*gamma)

mean(psi_est)

# logit(gamma) = B0 + B1*X1 + B2*X2 + etc.

gamma_est <- plogis(mc$gamma_beta %*% constant_list$X_gamma[i,] +
                  mc$gamma_fall * constant_list$season_vec[t])

mean(gamma_est)

# logit(theta) = B0 + B1*X1 + B2*X2 + etc.

phi_est <- plogis(mc$phi_beta %*% constant_list$X_phi[i,])

mean(phi_est)

# rho

rho_est <- plogis(mc$rho_beta %*% constant_list$X_rho[i,])

mean(rho_est)