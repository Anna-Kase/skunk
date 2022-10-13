



library(uwinspatialtools)
library(dplyr)



# Code from: 

##############################################
#
# Example of extracting from landuse/landcover
#
# Written by T. Gallo and M. Fidino
#
##############################################

# with proper modifications

# source spatial points file
source("./R/creating_spatial_points.R")


#  High-res Landcover for NE Illinois 
#  Data can be downloaded from:
#  browseURL("https://datahub.cmap.illinois.gov/dataset/high-resolution-land-cover-ne-illinois-and-nw-indiana-2010")
#  Load LULC map

my_raster_path <-
  "../../GIS/cmap/landcover_2010_chicagoregion.img"


# read it in
my_map <- raster::raster(my_raster_path)

#  Extract the proportion canopy cover (lulc class 1)
#  and create 'impervious cover' value, which is the sum of multiple
#  lulc classes

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

# clean up
rm(my_map)
gc()



# Load 2010 statewide census data for housing and population
# Data can be downloaded from:
# browseURL("http://silvis.forest.wisc.edu/data/housing-block-change/")

pop_data <- sf::st_read(
  "../../GIS/housing_density",
  layer = "il_blk10_Census_change_1990_2010_PLA2"
)

# fix any potential issues with the vectors before trying
#  to summarise
pop_data <- sf::st_make_valid(pop_data)

# Run function to calculate housing units, housing density, population
#  and population density.  Extract population data
#  within a 1km radius buffer.

population_data <- extract_polygon(
  my_points = sites,
  location_column = "Site",
  my_buffer = 1000,
  my_shape = pop_data,
  layers = c("HU10")
)


# Now join the landcover and the population data into one object

my_covariates <- inner_join(
  x= population_data,
  y= lulc_prop,
  by = "Site"
)

# remove unecessary column(s) and drop units
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

