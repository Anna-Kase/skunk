


# Skunk Data Prep

# read in old and new datasets
og <- read.csv("./data/full_capture_history.csv")

updated <- read.csv("../../../Downloads/full_capture_history.csv")

# remove two columns that the new dataset had that the old one did not
updated <- updated %>% dplyr::select(-c(Start, End))

# full join both datasets
both <- dplyr::left_join(updated, og)

# remove sites that were too close to other sites
both <- both %>% dplyr::filter(!Site %in% c("C05-EGG1", "C05-BMP2", "D02-MOP1",
                                            "D05-BRP1", "D07-NOC1", "R06-SLT2",
                                            "R06-SLT3", "S07-SAG2", "S07-STJ1",
                                            "S07-WGL0"))

# # copy the J column so one can be the original and the new can
# # have zeros converted to NAs
# both <- both %>% 
#   mutate(naJ = J)
# 
# # convert zeros in new J column to be NAs
# both["naJ"][both["naJ"] == 0] <- NA
# 
# # calculate the mean number of NAs in the new J column - 
# # these will later become NAs in the detection data (Y) as 
# # well - but they will be estimated by the model
# perc_na <- both %>% 
#   group_by(Season) %>% 
#   summarize(
#     perc_na <- mean(is.na(naJ))
#   )

# remove a few early seasons of data to create a complete seasonal
# detection history
# JU16 also was removed because it was not behaving

complete <- both %>% 
  dplyr::filter(!Season %in% c("JA11", "JA13", "JU13", "JA14", "JU16"))

# save this dataset as a .csv because it needs to be used in the 
# scripts that create the covariate document so the sites are correct
# and everything matches
write.csv(complete, "./data/complete_data.csv")


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

complete$Season <- factor(complete$Season, levels=unique(complete$Season))
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

my_seasons <- unique(complete$Season)
season_vec <- rep(0, length(my_seasons))
season_vec[grep("JU", my_seasons)] <- 1



# source all of the covariate scripts
source("./R/creating_spatial_points.R")
source("./R/scaled_covariates.R")
source("./R/dist_water_cov.R")
source("./R/managed_lawn_cov.R")

# read in the covariate data created by sourced code above
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



