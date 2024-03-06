
library(nimble)
library(dplyr)


# Intercept only model brier score


# read in saved RDS file if not already loaded
chain_output <- readRDS("./skunk_rds/intercept_only.rds")


# Sub-sample the posterior so these arrays are not as huge
set.seed(1995)

my_samples2 <- sample(
  1:nrow(chain_output),
  2500
)

mod_sub <- chain_output[my_samples2,]

# "Extract" MCMC posterior
source("./R/mcmc_functions.R")
mc <- split_mcmc(mod_sub)

# Create nsite and nseason objects from original data prep for
# ease of access
source("./R/intercept_data_prep.R")
nsite <- constant_list$nsite
nseason <- constant_list$nseason


# Create empty arrays to house calculated z values and the probabilities needed
# to calculate those z values (aka indicator function for species presence)
z <- z_prob <- array(NA, dim=c(length(my_samples2), constant_list$nsite, 
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
fpsi_prob <- array(NA, dim=c(length(my_samples2), 
                            constant_list$nsite, constant_list$nseason))


# forecasted z matrix
fz <- array(NA, dim=c(length(my_samples2), nsite, 5))
fz[,,1] <- z[,,nseason]

# forecasted z_prob 
fz_prob <- array(NA, dim=c(length(my_samples2), nsite, 5))
fz_prob[,,1] <- z_prob[,,nseason]

# forecasted site neighbors
fzeta <- array(NA, dim=c(length(my_samples2), nsite, 5))
fzeta[,,1] <- zeta[,,nseason]

# forecasted zm*(log(1-d_vec))
fdelta_bar <- array(NA, dim=c(length(my_samples2), nsite, 5))
fdelta_bar[,,1] <- delta_bar[,,nseason]


for(t in 2:5){
  for(i in 1:nsite){
    fpsi_prob[,i,t] <- (fz[,i,1] * mc$phi) +
      ((1-fz[,i,1]) * fzeta[,i,1] * fdelta_bar[,i,1]) +
      ((1-fz[,i,1]) * (1-fzeta[,i,1])*mc$gamma)
  }
  # now simulating the forecasted z matrix
  fz[,,t] <- rbinom(
    prod(dim(fz)[1:2]),
    size = 1,
    prob = fpsi_prob[,,t]
  )
}




# This generates fz. From there, we need to get
# 1. The detection / non-detection data we held out.
# 2. The number of days each camera trap was operational
#      on those deployments.
# 3. The detection probability for each site

# for 1, let's call it fy. for 2, I'll call it fj.


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


fj <- matrix(
  fyd$J,
  ncol=fynseason,
  nrow=fynsite
)

fj[fj==0] <- NA


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



# drop the 'known' season at this point so we are just
#  left with the forecasts.
fz <- fz[,,-1]
fpsi_prob <- fpsi_prob[,,-1]

for(i in 1:nsite){
  for(t in 1:4){
    # skip if no data
    if(is.na(fj[i,t])){
      next
    }
    # species observed
    sp_det <- fz_prob[,i,t] * (1 - (1 - mc$rho)^fj[i,t])
    # species not observed
    sp_not_det <- (1 - fz_prob[,i,t]) +
      fz_prob[,i,t] * (1 - mc$rho)^fj[i,t]
    eval_prob <- (sp_det * fy[i,t]) +
      (sp_not_det * (1 - fy[i,t]))
    
    
    brier_post[,i,t] <- (fy[i,t] - eval_prob)^2
  }
}
# total accuracy
mean(brier_post, na.rm = TRUE)


