

##############################################
#
# Example of extracting from landuse/landcover
#
# Written by T. Gallo and M. Fidino
#
##############################################

library(uwinspatialtools)
library(dplyr)


# Need a table of site locations
# Here I use a random sample of 10 Chicago sites as an example
# Data avaliable in same github repo
# Columns: LocationName, UTM_E, UTM_N, UTMZone, City
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

# We will use the High-res Landcover for NE Illinois for this example
#  Data can be downloaded from:
#  browseURL("https://datahub.cmap.illinois.gov/dataset/high-resolution-land-cover-ne-illinois-and-nw-indiana-2010")
#  Load iLULC map
# REPLACE FILE PATH WITH LOCAL FILE PATH

my_raster_path <-
  "../../GIS/cmap/landcover_2010_chicagoregion.img"



# read it in
my_map <- raster::raster(my_raster_path)

#  For this example we will extract the proportion canopy cover (lulc class 1)
#    and create our own 'impervious cover' value, which is the sum of multiple
#    lulc classes.
lulc_prop <- extract_raster_prop(
  my_points = sites,
  location_column = "Site",
  my_buffer = 1000,
  my_raster_data = my_map,
  lulc_cats = list(
    "tree" = 1,
    "imperv" = 5:7
  )
)

rm(my_map)
gc()

# Load 2010 statewide census data for housing and population
# Data can be downloaded from:
# browseURL("http://silvis.forest.wisc.edu/data/housing-block-change/")
# REPLACE FILE PATH AND LAYER NAME WITH LOCAL FILE PATH
pop_data <- sf::st_read(
  "../../GIS/housing_density",
  layer = "il_blk10_Census_change_1990_2010_PLA2"
)

# fix any potential issules with the vectors before trying
#  to summarise.
pop_data <- sf::st_make_valid(pop_data)

# Run function to calculate housing units, housing density, population
#  and population density.  For this example we extract population data
#  within a 1km radius buffer.

population_data <- extract_polygon(
  my_points = sites,
  location_column = "Site",
  my_buffer = 1000,
  my_shape = pop_data,
  layers = c("HU10")
)


my_covariates <- inner_join(
  x= population_data,
  y= lulc_prop,
  by = "Site"
)

my_covariates <- select(my_covariates, -LocationName.geometry)
my_covariates$HU10 <- units::drop_units(my_covariates$HU10)


#making scaled covariates
scaled_covariates <- my_covariates

scaled_covariates <- scaled_covariates %>% 
  summarize_if(
    is.numeric, .funs = scale
  ) %>% 
  data.frame()

my_pca <- prcomp(scaled_covariates, center=FALSE)

my_covariates$urb <- my_pca$x[,1]

loadings <- my_pca$rotation[,1]

my_example <- as.matrix(scaled_covariates) %*% loadings

write.csv(my_covariates, "./data/site_covariates.csv", row.names = FALSE)

