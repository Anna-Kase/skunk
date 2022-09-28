



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
  
  dat <- dat[grep("18", dat$Season),]
  
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


