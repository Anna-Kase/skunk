
library(sf)
library(dplyr)
sf::sf_use_s2(
  FALSE
)


# our distance
l <- 2500

# read in data
dat <- read.csv(
  "./data/full_capture_history.csv"
)

# make a spatial data frame
dat <- sf::st_as_sf(
  dat,
  coords = c(
    "Long",
    "Lat"
  ),
  crs = 4326
)

# convert to UTMs
dat <- sf::st_transform(
  dat,
  crs = 32616
)

# get all the unique sites, then compute distances

unq <- dplyr::distinct(
  dat[,"Site"]
)

# double check to make sure all sites are unique,
# This should be 0
if(
  sum(
    duplicated(unq$Site)
  ) > 0
){
  stop("You have duplicate site names, fix these!")
}


site_dist <- sf::st_distance(
  unq
)
diag(site_dist) <- -1

# make the M matrix, but we don't know max(n) yet.
m_list <- vector(
  "list",
  length = nrow(
    site_dist
  )
)
names(m_list) <- unq$Site

for(i in 1:nrow(site_dist)){
  nearby <- which(
    as.numeric(
      site_dist[i,]
    ) < l
  )
  m_list[[i]] <- nearby[-which(nearby == i)]
}


# check max
n <- sapply(
  m_list,
  length
)

max_n <- max(n)

# create m
m <- matrix(
  0,
  ncol = max_n,
  nrow = nrow(site_dist)
)

# fill it in
for(i in 1:nrow(site_dist)){
  if(
    length(m_list[[i]]) == 0
  ){
    next
  }
  m[i,1:n[i]] <- m_list[[i]]
}

# change m so that can accommodate sites with no neighbors

site_m <- site_dist
site_m <- units::drop_units(site_m)
site_m[site_m <= l] <- 1
site_m[site_m > l] <- 0
diag(site_m) <- 0




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

dm <- cbind(1, covs$urb, scale(covs$water_dist), scale(covs$open_dev))
site_km <- site_dist
units(site_km) <- "km"
diag(site_km) <- 1

site_km_array <- array(NA,dim=c(nsite,nsite,3))
site_km_array[,,1:3] <- c(site_km, site_km, site_km) 

delta_array <- array(NA, dim=c(nsite,nsite,4))
for(i in 1:nsite){
  delta_array[,,1] <- 1
  delta_array[i,,2] <-dm[i,2] - dm[,2]  
  delta_array[i,,3] <- dm[i,3] - dm[,3]
  delta_array[i,,4] <- dm[i,4] - dm[,4]
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
  ncovar_delta = dim(delta_array)[3],
  delta_array = delta_array
) 

