
library(dplyr)


locs <- readRDS("./data/simulation_neighbors.RDS")
locs[[1]]


covs <- read.csv("./data/scaled_simulation_covariates.csv")

output <- readRDS("./skunk_rds/spatial_covariates_fall2.RDS")
source("./R/mcmc_functions.R")
output <- do.call("rbind", output)
set.seed(89)
nsamp <- 2500
output <- output[sample(1:nrow(output), nsamp),]
mc <- split_mcmc(output)


dm_covs <- cbind(1, covs$urb, covs$water_dist, covs$open_dev)
dm_delta <- vector("list", length = nrow(covs))
for(i in 1:nrow(covs)){
  nn <- length(locs[[i]])
  tmp <- array(1, dim = c(nn, 4))
  tmp[,2] <-dm_covs[i,2] - dm_covs[locs[[i]],2]  
  tmp[,3] <- dm_covs[i,3] - dm_covs[locs[[i]],3]
  tmp[,4] <- dm_covs[i,4] - dm_covs[locs[[i]],4]
  dm_delta[[i]] <- tmp
}





# number of simulations from mcmc
nsim <- 100

nsite <- nrow(dm_covs)
nseason <- 12
fall_vec <- c(0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0,
              0, 1, 0, 0, 0, 1)

oz <- array(0, dim = c(nsite, nseason))

pb <- txtProgressBar(max = nsamp)
for(s in 1:nsamp){
  setTxtProgressBar(pb, i)
# simulation model
  tmp_psi <- dm_covs %*% mc$psi_beta[s,]
  tmp_psi <- plogis(tmp_psi)
  z <- matrix(NA, nrow=nsite, ncol=nseason)

for(i in 1:nsim){  
  # initial occupancy
  
  # sample across the landscape
  z[,1] <- rbinom(
    nrow(tmp_psi),
    1,
    tmp_psi
  )
  
  for(j in 2:12){
    my_prob <- rep(NA, nsite)
    for(site in 1:nsite){
      nn <- locs[[site]]
      has_n <- as.numeric(sum(z[nn,j-1]) > 0)
      if(z[site,j-1] == 1){
        tmp_phi <- dm_covs[site,] %*% mc$phi_beta[s,]
        my_prob[site] <- plogis(tmp_phi)
      }
      if(z[site,j-1] == 0 & !has_n){
        tmp_gamma <- dm_covs[site,] %*% mc$gamma_beta[s,] + 
          mc$gamma_fall[s,] * fall_vec[j]
        my_prob[site] <- plogis(tmp_gamma)
      }
      if(z[site,j-1] == 0 & has_n){
        tmp_delta <- dm_delta[[site]] %*% mc$delta_beta[s,] +
          mc$delta_fall[s,] * fall_vec[j]
        tmp_delta <- plogis(tmp_delta)
        my_prob[site] <- 1 - exp(
          t(log(1 - tmp_delta)) %*%
          (z[locs[[site]],j-1]) 
        )
      }
    }
      
    z[,j] <- rbinom(
      length(my_prob),
      1,
      my_prob
    )
    }  
}
  
  oz <- oz + z
}

write.csv(oz, "./data/original_z_sim_values.csv")

mu_z <- oz/(nsim*nsamp)

tmp_sd <- rowSums(oz)
sd_mat <- matrix(0, nrow = nsite, ncol = nsim*nsamp*nseason)
for(i in 1:nrow(sd_mat)){
  sd_mat[i,1:tmp_sd[i]] <- 1
}
site_sd <- apply(sd_mat, 1, sd)
