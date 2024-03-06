


# Spatial covariates with fall term and 
# urbanization tolerance interaction data prep

# read in data
complete <- read.csv("./data/complete_data.csv")


# Make coordinates spatial
complete <- sf::st_as_sf(
  complete,
  coords = c(
    "Long",
    "Lat"
  ),
  crs = 4326
)

# convert to UTMs
complete <- sf::st_transform(
  complete,
  crs = 32616
)

unq <- dplyr::distinct(
  complete[,"Site"]
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

# take out the 2021 data to use as hold out/validation data
complete <- complete %>% 
  dplyr::filter(!Season %in% c("JA21", "AP21", "JU21", "OC21"))


# NEW STUFF

# get all the unique years
unq_years <- unique(
  substr(complete$Season, 3, 4)
)


# and the unique seasons, just hard coding these
unq_seasons <- c("JA", "AP", "JU", "OC")


# combine
my_seasons <- expand.grid(
  unq_seasons, unq_years
)
# paste them together
my_seasons <- apply(
  my_seasons,
  1,
  paste0,
  collapse = ""
)
# drop JA14
my_seasons <- my_seasons[-1]

# get all the unique sites
unq_sites <- unique(complete$Site)


# combine seasons and sites

full_history <- expand.grid(
  Season = my_seasons,
  Site = unq_sites
)
full_history$Season <- factor(
  full_history$Season,
  levels = my_seasons
)
full_history <- full_history[
  order(full_history$Season, full_history$Site),
]

# get just the unique sites
site_locs <- dplyr::distinct(
  complete[, c("Site", "geometry")]
)

# join on site locations
tmp_complete <- dplyr::inner_join(
  full_history,
  site_locs,
  by = "Site"
)

# join to complete
complete <- dplyr::left_join(
  full_history,
  complete[,c("Season", "Site", "Crs", "Y", "J")],
  by = c("Season", "Site")
)

# turn NA to 0's
complete$J[is.na(complete$J)] <- 0
complete$Y[is.na(complete$Y)] <- 0
complete$Crs <- 4326

# BACK TO OLD STUFF


# Create m matrix

# set distance
l <- 2500

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



# Data and Constant lists and covariates 

# SOME SMALL CHANGES HERE - MOSTLY FORMATTING

complete$Season <- factor(
  complete$Season, 
  levels = my_seasons)
complete <- complete[order(complete$Season, complete$Site), ]
nsite <- dplyr::n_distinct(complete$Site)
nseason <- dplyr::n_distinct(complete$Season)

# y matrix
y <- matrix(
  complete$Y,
  ncol=nseason,
  nrow=nsite
)

# J matrix
J <- matrix(
  complete$J,
  ncol=nseason,
  nrow=nsite
)

y[J==0] <- NA


# creating binary seasonal vector that will "turn on" when it is fall
# to account for a shift in colonization due to offspring dispersal, and 
# off in every other season -- HOWEVER, the way in which the model is coded,
# delta_bar and gamma (our colonization parameters) are time lagged because
# their value depends upon skunk presence in the previous time step,
# therefore, to have this seasonal effect "turn on" at the correct time,
# we have to code it so that it actually turns on in the season preceding
# the season of interest to account for that lag. In our case, the season 
# of interest is the fall, so we turn on the effect during the summer.

# GOING BACK TO USING OCTOBER AS THE FALL TURN ON POINT

my_seasons <- unique(complete$Season)
season_vec <- rep(0, length(my_seasons))
season_vec[grep("OC", my_seasons)] <- 1



# Year just needs to be a numeric vector, not a factor. We also
#  want it to switch during the fall. Easiest way to spin this
#  up is to take the cumulative sum of season_vec and then add 1.
year_vec <- cumsum(season_vec) + 1


# Covariate scripts - can be sourced individually,
# but one file is written that contains all three spatial
# covariates -- to avoid bugs becuase they seem to randomly
# appear in these scripts and these take a hot minute to run,
# just load in the final .csv file
# source("./R/creating_spatial_points.R")
# source("./R/scaled_covariates.R")
# source("./R/dist_water_cov.R")
# source("./R/managed_lawn_cov.R")

# read in the covariate data created by sourced code above
covs <- read.csv("./data/site_covariates.csv")

# intial design matrix, going to modify for gamma and delta
dm <- cbind(1, covs$urb, scale(covs$water_dist), scale(covs$open_dev))

dmt <- array(
  NA,
  dim = c(nrow(dm), ncol(dm)+2, nseason)
)
# scale year_vec
year_vec <- scale(year_vec)
for(i in 1:nseason){
  dmt[,1:4,i] <- dm
  dmt[,5,i] <- as.numeric(year_vec[i])
  dmt[,6,i] <- covs$urb * as.numeric(year_vec[i])
}


site_km <- site_dist
units(site_km) <- "km"
diag(site_km) <- 1

# ADDING A BIT HERE
site_km <- apply(
  site_km,
  2,
  function(x) as.numeric(x)
)

site_km_array <- array(NA,dim=c(nsite,nsite,4))

for(i in 1:dim(site_km_array)[3]){
  site_km_array[,,i] <- site_km
}


# also needs to vary through time now. What is funny about this
#  is that because these are differences the year term just
#  drops out (we'll need to add it in a different way).
delta_array <- array(NA, dim=c(nsite,nsite,5, nseason))
for(i in 1:nsite){
  for(j in 1:nseason){
    delta_array[,,1,j] <- 1
    delta_array[i,,2,j] <-dm[i,2] - dm[,2]  
    delta_array[i,,3,j] <- dm[i,3] - dm[,3]
    delta_array[i,,4,j] <- dm[i,4] - dm[,4]
    delta_array[i,,5,j] <- dmt[i,6,j] - dmt[,6,j] 
  }
}
for(j in 1:nseason){
  delta_array[,,2:5,j] <- delta_array[,,2:5,j]/site_km_array
}



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
  X_gamma = dmt, 
  X_phi = dm,
  ncovar_psi = ncol(dm),
  ncovar_rho = ncol(dm),
  ncovar_gamma = dim(dmt)[2],
  ncovar_phi = ncol(dm),
  ncovar_delta = dim(delta_array)[3],
  delta_array = delta_array,
  season_vec = season_vec,
  year_vec = year_vec
) 
