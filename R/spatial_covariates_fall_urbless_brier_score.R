
library(dplyr)


# Spatial covariates with fall and urbless term brier scores

# read in saved RDS file if not already loaded
chain_output <- readRDS("../skunk/skunk_rds/spatial_covariates_fall_urbless.rds")

chain_output <- do.call("rbind", chain_output)

# Sub-sample the posterior so these arrays are not as huge
set.seed(1995)

my_samples2 <- sample(
  1:nrow(chain_output),
  50
)

cov_mod_sub <- chain_output[my_samples2,]

# "Extract" MCMC posterior
source("./R/mcmc_functions.R")
mc <- split_mcmc(cov_mod_sub)



# Create nsite and nseason objects from original data prep for
# ease of access
source("./R/spatial_covariates_fall_urbless_data_prep.R")
nsite <- constant_list$nsite
nseason <- constant_list$nseason


# Create empty arrays to house calculated z values and the probabilities needed
# to calculate those z values (aka indicator function for species presence)
z <- z_prob <- array(
  NA, 
  dim=c(
    length(my_samples2), 
    constant_list$nsite,
    constant_list$nseason))


# Calculating z probabilities for the first season of data (t=1)
# and store them in z_prob
# In the original data we gave the model, when the species was not
# detected, the value in the array will be 1. If the species was 
# detected at a site, then the z probability will be:
# ((psi * (1-rho))^J) / (((psi * (1-rho))^J) + (1-psi) )
# CHECK WHAT IS GOING ON WITH NA SITES
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



# Fill in the z array for the first season based on z probabilities 
# This is array holds the indicator term for species presence
z[,,1] <- rbinom(
  prod(dim(z)[1:2]),
  size = 1,
  prob = z_prob[,,1]
)

rm(psi, rho, chain_output)
gc()

# Fill in delta_bar (aka neighborhood colonization probabilities) at t=1


# Create empty delta_bar and zeta arrays
delta_bar <- zeta <- array(dim=dim(z_prob))

# Create d_vec (aka neighborhood colonization covariate) matrix

d_vec <- array(NA, dim=c(length(my_samples2), nsite, nsite, nseason))
for(i in 1:nsite){
  for(ti in 1:nseason){
  d_vec[,i,,ti] <- (mc$delta_beta %*% t(constant_list$delta_array[i,,,ti])) +
    (as.numeric(mc$delta_fall) * constant_list$season_vec[ti]) +
    (as.numeric(mc$delta_year) * constant_list$year_vec[ti])
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
      gamma <- plogis(mc$gamma_beta %*% constant_list$X_gamma[i,,t] +
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





# Make psi_prob aka forecast probabilities (f[i,t]) for Brier Scores

# Create empty matrix
fpsi_prob <- array(
  NA, 
  dim=c(length(
    my_samples2), 
    constant_list$nsite, 
    5))

# # Fill in psi values at t = 1 (just the straight psi values from the posterior)
# for(i in 1:nsite){
#   psi <- plogis(mc$psi_beta %*% constant_list$X_psi[i,])
#   fpsi_prob[,i,1] <- psi
# }

# Now fill in the the psi_prob array for the rest of the seasons
# Here we do not hard code anything, rather if the species was not present
# we leave the psi probability as the psi value from the MCMC posterior.
# If the species was present, then the psi probability is calculated based
# on the linear predictors from the dynamic occupancy model formula: 
# (z[i,t] * theta[i,t]) + ((1 - z[i,t]) * I) + ((1 - z[i,t]) * (1 - I) * gamma)
fz <- array(NA, dim=c(length(my_samples2), nsite, 5))
fz[,,1] <- z[,,nseason]

fzeta <- array(NA, dim=c(length(my_samples2), nsite, 5))
fzeta[,,1] <- zeta[,,nseason]

fdelta_bar <- array(NA, dim=c(length(my_samples2), nsite, 5))
fdelta_bar[,,1] <- delta_bar[,,nseason]

rm(delta_bar)
rm(d_vec)
gc()

fd_vec <- array(NA, dim=c(length(my_samples2), nsite, nsite, 5))
for(i in 1:nsite){
  for(ti in 1:5){
    fd_vec[,i,,ti] <- (mc$delta_beta %*% t(constant_list$delta_array[i,,,ti])) +
      (as.numeric(mc$delta_fall) * c(0,0,0,0,1)[ti]) +
      (as.numeric(mc$delta_year) * c(0,0,0,0,1)[ti])
  }
}
fd_vec <- plogis(fd_vec)

# find number of detection days nsite, 5 matrix - add to line 193
for(t in 2:5){
  for(i in 1:nsite){
    phi <- plogis(mc$phi_beta %*% constant_list$X_phi[i,])
    gamma <- plogis(mc$gamma_beta %*% constant_list$X_gamma[i,,t] +
                      (as.numeric(mc$gamma_fall)) * c(0,0,0,0,1)[t])
    fpsi_prob[,i,t] <- (fz[,i,t-1] * phi) +
      ((1-fz[,i,t-1]) * fzeta[,i,t-1] * fdelta_bar[,i,t-1]) +
      ((1-fz[,i,t-1]) * (1-fzeta[,i,t-1])*gamma)
    fz[,,t] <- rbinom(
      prod(
        dim(fpsi_prob)[1:2]
      ),
      1,
      prob = fpsi_prob[,,t]
    )
    for(i in 1:nsite){
      zm <- sweep(
        z[,1:nsite,t],
        2,
        constant_list$m[i, 1:nsite],
        FUN = "*"
      )
      # Now calculate and fill in delta_bar for all of the seasons
      tmp_zeta <- rowSums(zm)
      fzeta[,i,t] <- as.numeric(tmp_zeta>0)
      zm <- zm * log(1 - fd_vec[,i,,t])
      zm <- rowSums(zm)
      zm <- 1 - exp(zm)
      fdelta_bar[,i,t] <- zm
    }
  }
}



# bring in hold out data

fyd <- read.csv("./data/complete_data.csv")
fyd <- fyd[grep("21", fyd$Season),]
fyd$Season <- factor(fyd$Season, levels = unique(fyd$Season))
fyd <- fyd[order(fyd$Season, fyd$Site), ]
fynsite <- n_distinct(fyd$Site)
fynseason <- n_distinct(fyd$Season)

fy <- matrix(
  fyd$Y,
  ncol=fynseason,
  nrow=fynsite
)

# fy should be a matrix that is nsite by 4.
# make binary
fy[fy>0] <- 1

# get number of mcmc samples
nsim <- dim(fz)[1]
brier_post <- array(
  NA,
  dim = c(
    length(my_samples2),
    fynsite,
    4
  )
)



fj <- matrix(
  fyd$J,
  ncol=fynseason,
  nrow=fynsite
)

fj[fj==0] <- NA

# drop the 'known' season at this point so we are just
#  left with the forecasts.
fz <- fz[,,-1]
fpsi_prob <- fpsi_prob[,,-1]

# get number of mcmc samples
nsim <- dim(fz)[1]
brier_post <- array(
  NA,
  dim = c(
    length(my_samples2),
    fynsite,
    ncol(fj)
  )
)


for(i in 1:nsite){
  rho <- plogis(mc$rho_beta %*% constant_list$X_rho[i,])
  for(t in 1:4){
    # skip if no data
    if(is.na(fj[i,t])){
      next
    }
    # species observed
    sp_det <- fpsi_prob[,i,t] * (1 - (1 - rho)^fj[i,t])
    # species not observed
    sp_not_det <- (1 - fpsi_prob[,i,t]) +
      fpsi_prob[,i,t] * (1 - rho)^fj[i,t]
    eval_prob <- (sp_det * fy[i,t]) +
      (sp_not_det * (1 - fy[i,t]))
    
    
    brier_post[,i,t] <- (fy[i,t] - eval_prob)^2
  }
}
# total accuracy
mean(brier_post, na.rm = TRUE)


