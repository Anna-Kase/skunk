
# Creating spatial points and projecting them into UTM to be 
# sourced for use in covariate creation


site_coords <- read.csv(
  "./data/full_capture_history.csv"
)

#unique sites
site_coords <- distinct(site_coords[,c("Site","Long", "Lat")])


# Create spatial points
# You must put the correct CRS for respective city
sites <- sf::st_as_sf(
  site_coords,
  coords = c("Long", "Lat"),
  crs = 4326
)

#reproject to UTM
sites <- sf::st_transform(
  sites,
  crs=32616
)
