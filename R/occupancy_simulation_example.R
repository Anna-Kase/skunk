
psi <- c(-2.27, -0.17)

gamma <- c(2.58, 0.59)

phi <- c(0.68, 0.09)

# make a grid of locations, # site by site by year

locs <- array(
  NA,
  dim = c(100,100, 12)
)
loc_count <- array(
  0,
  dim = c(100,100, 12)
)


one_loc <- locs[,,1]

my_neighs <- vector("list", length = prod(dim(one_loc)))

for(i in 1:nrow(one_loc)){
  for(j in 1:ncol(one_loc)){
    tmp_neighs <- matrix(
      c(
        i + 1, j,
        i, j-1,
        i - 1, j,
        i, j+1
      ),
      ncol = 2,
      byrow = TRUE
    )
    to_go <- which(
      tmp_neighs < 1 | tmp_neighs > 100,
      arr.ind = TRUE
    )
    if(nrow(to_go)>0){
      tmp_neighs <- tmp_neighs[-to_go[,1],]
    }
    tmp_neighs <- c(i - 1, i +1, j - 1, j +1 )
    tmp_neigh_locs
  }
}



for(i in 1:nsim){  
  # initial occupancy
  tmp_psi <- dm %*% psi
  tmp_psi <- plogis(tmp_psi)
  # sample across the landscape
  locs[,,1] <- rbinom(
    nrow(tmp_psi),
    1,
    tmp_psi
  )
  
  for(j in 2:12){
    zeta[i,j-1] <- 
      locs[,,j-1] %*% tmp_neighs
    
    tmp_gamma <- plogis(dm %*% gamma)
    tmp_phi <- plogis(dm %*% phi)
    
    tmp_delta <- 1 - exp(
      (locs[,,j-1] * tmp_neighs) %*%
        log(1 - (plogis(dm %*% delta))) #***
    )
  }  
}

# ***This is the part where I keep going back and forth if we need 
  # more specificity on # of neighbors. Originally, the indicator term
  # took on a zero or a one that would turn this entire delta term
  # on and off. So it would make sense that we push this through the
  # same equation from the original model where it wasn't explicit how
  # many neighbors a site had, and we used the delta output in 
  # conjunction with a design matrix to obtain the delta values used
  # for plotting. Therefore, it seems like that process should be repeated
  # here for the simulations.
  # One of my original ideas was that we would need to have a bigger
  # if loop that would use a specific delta probability based on the
  # number of neighbors a site had. But that no longer seems correct
  # the more I have thought about this. 








# covariates -- this case we are using the 
# intercept and open space
dm <- cbind(1, rnorm(100*100))
dm[,2] <- sort(dm[,2])

# number of simulations from mcmc
nsim <- 100



pb <- txtProgressBar(max = nsim)
for(i in 1:nsim){
  setTxtProgressBar(pb, i)
  # initial occupancy
  tmp_psi <- dm %*% psi
  tmp_psi <- plogis(tmp_psi)
  # sample across the landscape
  locs[,,1] <- rbinom(
    nrow(tmp_psi),
    1,
    tmp_psi
  )
  for(j in 2:12){
    zeta[i,j-1] <- 
      locs[,,j-1] %*% tmp_neighs
    
    tmp_gamma <- plogis(dm %*% gamma)
    tmp_phi <- plogis(dm %*% phi)
    
    tmp_delta <- 1 - exp(
      (locs[,,j-1] * tmp_neighs) %*%
        log(1 - (plogis(dm %*% delta)))
    )
    
    my_pred <- (tmp_phi[,1] * locs[,,j-1]) +
      ((1 - locs[,,j-1]) * zeta[, j-1] * tmp_delta[,j-1]) +
      (tmp_gamma[,1] * (1 - locs[,,j-1]) * zeta[, j-1])
    
    locs[,,j] <- rbinom(
      nrow(tmp_gamma),
      1,
      my_pred
    )
  }
  loc_count <- loc_count + locs
}

overall_occupancy <- loc_count / nsim

my_mean <- apply(
  overall_occupancy,
  c(1,2),
  mean
)
image(my_mean)

my_sd <- apply(
  overall_occupancy,
  c(1,2),
  sd
)
image(my_sd)

one_loc <- locs[,,1]

my_neighs <- vector("list", length = prod(dim(one_loc)))

for(i in 1:nrow(one_loc)){
  for(j in 1:ncol(one_loc)){
    tmp_neighs <- matrix(
      c(
        i + 1, j,
        i, j-1,
        i - 1, j,
        i, j+1
      ),
      ncol = 2,
      byrow = TRUE
    )
    to_go <- which(
      tmp_neighs < 1 | tmp_neighs > 100,
      arr.ind = TRUE
    )
    if(nrow(to_go)>0){
      tmp_neighs <- tmp_neighs[-to_go[,1],]
    }
    tmp_neighs <- c(i - 1, i +1, j - 1, j +1 )
    tmp_neigh_locs
  }
}





# get delta neighbor 1
dprob <- rbeta(1e5, 4, 6)
dprob <- matrix(
  dprob,
  ncol = 1
)

# phi matrix
phi <- matrix(
  rbeta(1e5 * 200, 9,5),
  ncol = 200,
  nrow = 1e5
)

dprob <- dprob[,rep(1, ncol(phi))]

# my occupancy
occ1 <- dprob / (dprob + (1 - phi))
occ_summary <- apply(
  occ1,
  2,
  quantile,
  probs = c(0.025,0.5,0.975)
)

