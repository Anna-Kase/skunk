
# Creating a distance to stream covariate
library(dplyr)
library(sf)

# source spatial points file
source("./R/creating_spatial_points.R")


# Dowload data from:
 browseURL("https://clearinghouse.isgs.illinois.edu/data/hydrology/streams-and-shorelines")
# load water data

water_path <- "../../GIS/water"
water <- sf::st_read(
  water_path,
  layer = "IL_Streams_From_100K_DLG_Ln"
)

# fix any potential issues
water_data <- sf::st_make_valid(
  water
)

# reproject into the UTM
water_data <- sf::st_transform(
  water_data,
  crs=32616
)

# cut down water data to the extent of the study area
# this is the same extent as the sites shapefile that was 
# sourced in
small_water <- sf::st_crop(
  water_data,
  xmin=395087,
  xmax=456196.7,
  ymin=4604022,
  ymax=4678480
)
# make valid once more after cropping
small_water <- sf::st_make_valid(
  small_water
)

# Plot it to make sure it looks correct
plot(small_water$geometry)
plot(sites$geometry, add=TRUE, pch=19, col="forestgreen")

# Selecting one column of data to work with that provides 
# information on if a linestring is a stream and selecting
# only objects that have the designation for being a stream.
# Interpretation of values on the metadata tab at same URL the data
# was downloaded from under "Entity_and_Attribute_Information"
smaller_water <- dplyr::select(
  small_water,
  ENR
)

smallest_water <- smaller_water %>% 
  dplyr::filter(ENR >= 400000) 

# Get the index of the closest stream to each site

closest_streams_index <- sf::st_nearest_feature(
  sites, 
  smallest_water
)

# Do element wise calculations, indexing
#  smallest_water by closest_streams_index
site_dist <- sf::st_distance(
  sites,
  smallest_water[closest_streams_index,],
  by_element = TRUE
)

# Remove the units (it's meters) from the distances
site_dist <- units::drop_units(
  site_dist
)

# Make it a dataframe
dist_df <- data.frame(
  water_dist = site_dist
)

# Add site designations back in so this is not just a list
# of distance values and it can be joined with the other covariates
dist_df$Site <- sites$Site

# Add the distance to water covariate to the other
# scaled covariates

if(
  file.exists(
    "./data/site_covariates.csv"
  )
){
  my_covariates <- read.csv(
    "./data/site_covariates.csv"
  )
} else {
  source("./R/scaled_covariates.R")
}

# check if dist column on here, if not add it
if(!"water_dist" %in% colnames(my_covariates)){
  more_covariates <- inner_join(
    x= my_covariates,
    y= dist_df,
    by = "Site"
  )
  write.csv(
    more_covariates,
    "./data/site_covariates.csv",
    row.names = FALSE
  )
}

# add water dist onto here
sites$dist <- dist_df$water_dist

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(
  c('blue','red')
)

#This adds a column of color values
# based on the y values
sites$col <- rbPal(100)[as.numeric(cut(sites$dist,breaks = 100))]

plot(small_water$geometry)
plot(sites$geometry, pch = 20, col = sites$col, cex=2, add=TRUE)



