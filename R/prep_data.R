



# intercept only model
if(model_type == "intercept"){
  #data list
  data_list <- list(
    y=y
  )
  
  #constant list
  constant_list <- list(
    J=J,
    nsite=nsite,
    nseason=nseason,
    m=site_m
  ) 
  
  dat <- dat[grep("18|19", dat$Season),]
  
  dat$Season <- factor(dat$Season, levels=unique(dat$Season))
  
  dat <- dat[order(dat$Season, dat$Site), ]
  
  nsite <- n_distinct(dat$Site)
  nseason <- n_distinct(dat$Season)
  
  # y matrix
  y <- matrix(
    dat$Y,
    ncol=nseason,
    nrow=nsite
  )
  
  # J matrix
  J <- matrix(
    dat$J,
    ncol=nseason,
    nrow=nsite
  )
  
  y[J==0] <- NA
  
}

# covariate model 1 (urb)
if(model_type == "covariate"){
  
  dat <- dat[grep("18|19", dat$Season),]
  
  dat$Season <- factor(dat$Season, levels=unique(dat$Season))
  
  dat <- dat[order(dat$Season, dat$Site), ]
  
  nsite <- n_distinct(dat$Site)
  nseason <- n_distinct(dat$Season)
  
  # y matrix
  y <- matrix(
    dat$Y,
    ncol=nseason,
    nrow=nsite
  )
  
  # J matrix
  J <- matrix(
    dat$J,
    ncol=nseason,
    nrow=nsite
  )
  
  y[J==0] <- NA
  
  
  covs <- read.csv("./data/site_covariates.csv")
  
  dm <- cbind(1, covs$urb)
  site_km <- site_dist
  units(site_km) <- "km"
  diag(site_km) <- 1
  delta1 <- matrix(
    NA,
    nrow = nsite,
    ncol = nsite
  )
  for(i in 1:nsite){
    delta1[i,] <- covs$urb[i] - covs$urb
  }
  delta1 <- delta1/units::drop_units(site_km)
  
  #data list
  data_list <- list(
    y=y
  )
  
  #constant list
  constant_list <- list(
    J=J,
    nsite=nsite,
    nseason=nseason,
    m=site_m,
    X_psi = dm,
    X_rho = dm,
    X_gamma = dm,
    X_phi = dm,
    X_delta1 = delta1,
    ncovar_psi = ncol(dm),
    ncovar_rho = ncol(dm),
    ncovar_gamma = ncol(dm),
    ncovar_phi = ncol(dm),
    ncovar_delta = 2
  ) 
  
}


<<<<<<< Updated upstream
=======
# covariate model 3 (urb, dist_water, and open_dev)

if(model_type == "covariate3"){
  
  dat <- dat[grep("18|19", dat$Season),]
  
  dat$Season <- factor(dat$Season, levels=unique(dat$Season))
  
  dat <- dat[order(dat$Season, dat$Site), ]
  
  nsite <- n_distinct(dat$Site)
  nseason <- n_distinct(dat$Season)
  
  # y matrix
  y <- matrix(
    dat$Y,
    ncol=nseason,
    nrow=nsite
  )
  
  # J matrix
  J <- matrix(
    dat$J,
    ncol=nseason,
    nrow=nsite
  )
  
  y[J==0] <- NA
  
  
  covs <- read.csv("./data/site_covariates.csv")
  
  dm <- cbind(1, covs$urb, covs$water_dist, covs$open_dev)
  site_km <- site_dist
  units(site_km) <- "km"
  diag(site_km) <- 1
  
  site_km_array <- array(NA,dim=c(nsite,nsite,3))
  site_km_array[,,1:3] <- c(site_km, site_km, site_km) 
  
  delta_array <- array(NA, dim=c(nsite,nsite,4))
  for(i in 1:nsite){
    delta_array[,,1] <- 1
    delta_array[i,,2] <- covs$urb[i] - covs$urb
    delta_array[i,,3] <- covs$water_dist[i] - covs$water_dist
    delta_array[i,,4] <- covs$open_dev[i] - covs$open_dev
  }
  delta_array[,,2:4] <- delta_array[,,2:4]/site_km_array
  
  
  #data list
  data_list <- list(
    y=y
  )
  
  #constant list
  constant_list <- list(
    J=J,
    nsite=nsite,
    nseason=nseason,
    m=site_m,
    X_psi = dm,
    X_rho = dm,
    X_gamma = dm,
    X_phi = dm,
    ncovar_psi = ncol(dm),
    ncovar_rho = ncol(dm),
    ncovar_gamma = ncol(dm),
    ncovar_phi = ncol(dm),
    ncovar_delta = dim(delta_array)
  ) 
  
}


>>>>>>> Stashed changes
