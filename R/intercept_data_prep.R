
# Intercept only model data prep

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
