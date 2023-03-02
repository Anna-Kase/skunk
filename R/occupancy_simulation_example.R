
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

# covariates
dm <- cbind(1, rnorm(100*100))
dm[,2] <- sort(dm[,2])

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
    tmp_gamma <- plogis(dm %*% gamma)
    tmp_phi <- plogis(dm %*% phi)
    my_pred <- (tmp_gamma[,1] * (1 - locs[,,j-1])) +
      (tmp_phi[,1] * locs[,,j-1])
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

