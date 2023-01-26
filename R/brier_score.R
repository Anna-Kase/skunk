

# Brier Scores

# Source original data and data list, MCMC functions, 
# and read in the saved model
model_type <- "intercept"
source("./R/create_m_n.R")
source("./R/prep_data.R")
source("./R/mcmc_functions.R")
long_shot <- readRDS("../intercept_only.rds")


# Sub-sample the posterior so these arrays are not as huge
set.seed(1995)
my_samples <- sample(
  1:nrow(long_shot),
  5000
)

long_shot <- long_shot[my_samples,]

# "Extract" MCMC posterior
mc <- split_mcmc(long_shot)

# Create nsite and nseason objects from original data prep for
# ease of access
nsite <- constant_list$nsite
nseason <- constant_list$nseason


# Create empty arrays to house calculated z values and the probabilities needed
# to calculate those z values (aka indicator function for species presence)
z <- z_prob <- array(NA, dim=c(length(my_samples), 
                               constant_list$nsite, 
                               constant_list$nseason))


# Calculating z probabilities for the first season of data (t=1)
# and store them in z_prob
# In the original data we gave the model, when the species was not
# detected, the value in the array will be 1. If the species was 
# detected at a site, then the z probability will be:
# ((psi * (1-rho))^J) / (((psi * (1-rho))^J) + (1-psi) )
for(i in 1:nsite){
  if(data_list$y[i,1] > 0 & !is.na(data_list$y[i,1])){
    z_prob[,i,1] <- 1
  } else{
    num <- mc$psi * (1-mc$rho)^constant_list$J[i,1]
    den <- num + (1-mc$psi)
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


# Fill in delta_bar (aka neighborhood colonization probabilities) at t=1


# Create empty delta_bar and zeta arrays
delta_bar <- zeta <- array(dim=dim(z_prob))

# Create d_vec (aka neighborhood colonization covariate) matrix
d_vec <- matrix(mc$d, nrow=dim(z_prob)[1], ncol=nsite)

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
  zm <- zm * log(1 - d_vec)
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
      
      num <- (z[,i,t-1] * mc$phi) + 
        ((1-z[,i,t-1]) * zeta[,i,t-1] * delta_bar[,i,t-1]) +
        ((1-z[,i,t-1]) * (1-zeta[,i,t-1])*mc$gamma)
      dnm1 <- 1 - num
      # add rho
      num <- num * (1 - mc$rho)^constant_list$J[i,t]
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
    zm <- zm * log(1 - d_vec)
    zm <- rowSums(zm)
    zm <- 1 - exp(zm)
    delta_bar[,i,t] <- zm
  }
}


# Make psi_prob aka forecast probabilities (f[i,t]) for Brier Scores

# Create empty matrix
psi_prob <- array(NA, dim=c(length(my_samples), constant_list$nsite, constant_list$nseason+4))

# Fill in psi values at t = 1 (just the straight psi values from the posterior)
for(i in 1:nsite){
  psi_prob[,i,1] <- mc$psi
}

# Now fill in the the psi_prob array for the rest of the seasons
# Here we do not hard code anything, rather if the species was not present
# we leave the psi probability as the psi value from the MCMC posterior.
# If the species was present, then the psi probability is calculated based
# on the linear predictors from the dynamic occupancy model formula: 
# (z[i,t] * theta[i,t]) + ((1 - z[i,t]) * I) + ((1 - z[i,t]) * (1 - I) * gamma)
for(t in 2:nseason){
  for(i in 1:nsite){
      num <- (z[,i,t-1] * mc$phi) + 
        ((1-z[,i,t-1]) * zeta[,i,t-1] * delta_bar[,i,t-1]) +
        ((1-z[,i,t-1]) * (1-zeta[,i,t-1])*mc$gamma)
      psi_prob[,i,t] <- num
    }
  }


# Calculate Brier Scores based on the difference of mean squares of the
# psi_prob array (the occurrence probability
# based on observed species presence, presence in the neighborhood, and
# colonization from other) and the z array (species occurrence)
BS <- (psi_prob - z)^2
mean(BS)

# hard code future data collection as median number of days
#  sampled n(i.e., 28).
pb <- txtProgressBar(min = nseason+1,max = nseason+4)
for(t in (nseason+1):(nseason+4)){
  setTxtProgressBar(pb, t)
  for(i in 1:nsite){
      num <- (z[,i,t-1] * mc$phi) + 
        ((1-z[,i,t-1]) * zeta[,i,t-1] * delta_bar[,i,t-1]) +
        ((1-z[,i,t-1]) * (1-zeta[,i,t-1])*mc$gamma)
      dnm1 <- 1 - num
      # add rho
      num <- num * (1 - mc$rho)^median(constant_list$J)
      z_prob[,i,t] <- num / (num + dnm1)
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
    zm <- zm * log(1 - d_vec)
    zm <- rowSums(zm)
    zm <- 1 - exp(zm)
    delta_bar[,i,t] <- zm
  }
}

for(t in (nseason+1):(nseason+4)){
  for(i in 1:nsite){
    num <- (z[,i,t-1] * mc$phi) + 
      ((1-z[,i,t-1]) * zeta[,i,t-1] * delta_bar[,i,t-1]) +
      ((1-z[,i,t-1]) * (1-zeta[,i,t-1])*mc$gamma)
    psi_prob[,i,t] <- num
  }
}


