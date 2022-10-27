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


