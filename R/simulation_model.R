library(dplyr)

locs <- readRDS("./data/simulation_neighbors.RDS")

covs <- read.csv("./data/scaled_simulation_covariates.csv")

output <- readRDS("./skunk_rds/spatial_covariates_fall2.RDS")
source("./R/mcmc_functions.R")
output <- do.call("rbind", output)
set.seed(89)
nsamp <- 1000
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


nsite <- nrow(dm_covs)
nseason <- 12
# this is from our model object, we are not
#  going to use all of it.
fall_vec <- c(0, 0, 1, 0, 0, 0, 1, 0, 0,
              0, 1, 0, 0, 0, 1, 0, 0, 0,
              1, 0, 0,
              0, 1, 0, 0, 0, 1)

oz <- array(0, dim = c(nsite, nseason))
plot_array <- array(0, dim=c(nsamp, 3, nseason-1))

pb <- txtProgressBar(max = nsamp)
for(s in 1:nsamp){
  setTxtProgressBar(pb, s)
  
# simulation model
  tmp_psi <- dm_covs %*% mc$psi_beta[s,]
  tmp_psi <- plogis(tmp_psi)
  z <- matrix(NA, nrow=nsite, ncol=nseason)
  # initial occupancy
  # sample across the landscape
  z[,1] <- rbinom(
    nrow(tmp_psi),
    1,
    tmp_psi
  )
  
  for(j in 2:12){
    my_prob <- rep(NA, nsite)
    has_n <- rep(NA, nsite)
    
    for(site in 1:nsite){
      has_n[site] <- as.numeric(sum(z[locs[[site]],j-1]) > 0)
    }
    # phi sites
    phi_sites <- which(z[,j-1] == 1)
    if(length(phi_sites)>0){
      tmp_phi <- dm_covs[phi_sites,] %*% mc$phi_beta[s,]
      my_prob[phi_sites] <- plogis(tmp_phi)
      plot_array[s,1,j-1] <- length(phi_sites)
    }
    gamma_sites <- which(z[,j-1] == 0 & !has_n)
      if(length(gamma_sites)>0){
        tmp_gamma <- dm_covs[gamma_sites,] %*% mc$gamma_beta[s,] + 
          mc$gamma_fall[s,] * fall_vec[j]
        my_prob[gamma_sites] <- plogis(tmp_gamma)
        plot_array[s,2,j-1] <- length(gamma_sites)
      }
    
    delta_sites <- which(z[,j-1] == 0 & has_n)
      if(length(delta_sites)>0){
        X <- dm_delta[delta_sites]
        Y <- locs[delta_sites]
        tmp_delta <- sapply(
              1:length(delta_sites),
              function(idx){
                td <- ( X[[idx]] %*% mc$delta_beta[s,] + 
                         mc$delta_fall[s,] * fall_vec[j])
                td <- log(1 - plogis(td))
                td <- 1 - exp(
                  t(td) %*% z[Y[[idx]],j-1]
                )
                td
              }
        )
        my_prob[delta_sites] <- tmp_delta
        plot_array[s,3,j-1] <- length(delta_sites)
      }
      
    z[,j] <- rbinom(
      length(my_prob),
      1,
      my_prob
    )


    
    }  
  oz <- oz + z

}

plot_array <- apply(
  plot_array / nrow(covs), 
  c(2,3),
  quantile,
  probs = c(0.025,0.5,0.975)
)


#write.csv(oz, "./data/original_z_sim_values.csv")
saveRDS(plot_array, "./data/phi_gamma_delta_by_season.RDS")
