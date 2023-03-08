
library(dplyr)
library(parallel)

locs <- readRDS("./data/simulation_neighbors.RDS")

covs <- read.csv("./data/scaled_simulation_covariates.csv")

output <- readRDS("./skunk_rds/spatial_covariates_fall2.RDS")
source("./R/mcmc_functions.R")
output <- do.call("rbind", output)
set.seed(89)
nsamp <- 10
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
plot_array <- array(0, dim=c(nsamp, 3, nseason))

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
      
    }
    gamma_sites <- which(z[,j-1] == 0 & !has_n)
      if(length(gamma_sites)>0){
        tmp_gamma <- dm_covs[gamma_sites,] %*% mc$gamma_beta[s,] + 
          mc$gamma_fall[s,] * fall_vec[j]
        my_prob[gamma_sites] <- plogis(tmp_gamma)
      
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
      
      }
      
    z[,j] <- rbinom(
      length(my_prob),
      1,
      my_prob
    )
    
    plot_array[s,1,j] <- length(my_prob[phi_sites])
    plot_array[s,2,j] <- length(my_prob[gamma_sites])
    plot_array[s,3,j] <- length(my_prob[delta_sites])
    
    }  
  oz <- oz + z

}

#write.csv(oz, "./data/original_z_sim_values.csv")

# average occupancy each site and season
mu_z <- oz/(nsamp)

# overall occupancy throughout study
gmu_z <- rowMeans(mu_z)

# variation in occupancy across seasons
occ_sd <- apply(mu_z, 1, sd)

gr <- read.csv(
  "./data/raw_sim_covariates/point_locs.csv"
)

library(sf)
sf::sf_use_s2(FALSE)
gr <- sf::st_as_sf(
  gr,
  coords = c("x","y"),
  crs = 32616
)
gr$occ_prob <- gmu_z
gr$occ_sd <- occ_sd

plot(gr["occ_prob"], pch = 15)
plot(gr["occ_sd"], pch = 15)



# Changes in phi, gamma, delta over simulated time

library(bbplot)

windows(7, 3)

{
  par(mar = c(5, 2, 0.5, 0.5), oma = c(0, 5, 0, 0), lend = 1)
  bbplot::blank(xlim = c(1,12), ylim = c(0, 1), bty = "l")
  
  bbplot::axis_blank(1, at = seq(1, 12, by = 1))
  bbplot::axis_blank(2)
  bbplot::axis_text(text = seq(1, 12, 1), side = 1, line = 0.9,
                    at = seq(1, 12, 1))
  bbplot::axis_text(side = 2, las = 1, line = 0.4)
  bbplot::axis_text("Simulated Season", side = 1, line = 2.5)
  bbplot::axis_text("Proportion Sites", side = 2, outer = TRUE, at = 0.6, 
                    line = 2)
  
  points(
    x = c(1:12)-0.125,
    y = ((sum(plot_array[1:nsamp,1,1:12]))/nsite)/nsamp,
    pch = 18,
    col = "black",
    cex = 2
  )
}




  