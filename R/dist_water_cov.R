
# Creating a distance to stream covariate

# source spatial points file
source("./R/creating_spatial_points.R")


# Dowload data from:
# browseURL("https://clearinghouse.isgs.illinois.edu/data/hydrology/streams-and-shorelines")
# load water data

water <- sf::st_read(
  "../../GIS/water",
  layer = "IL_Streams_From_100K_DLG_Ln"
)

# fix any potential issues
water_data <- sf::st_make_valid(water)

# reproject into the UTM
water_data <- sf::st_transform(
  water_data,
  crs=32616
)

# cut down water data to the extent of the study area
# this is the same extent as the sites shapefile that was 
# sourced in
small_water <- sf::st_crop(water_data, xmin=395087, xmax=456196.7,
                           ymin=4604022, ymax=4678480)


# Plot it to make sure it looks correct
plot(small_water$geometry)
plot(sites$geometry, add=TRUE, pch=19, col="forestgreen")

# Selecting one column of data to work with that provides 
# information on if a linestring is a stream and selecting
# only objects that have the designation for being a stream.
# Interpretation of values on the metadata tab at same URL the data
# was downloaded from under "Entity_and_Attribute_Information"
smaller_water <- select(small_water, ENR)

smallest_water <- smaller_water %>% 
  filter(ENR >= 400000) 

# Calculate the distance of each site to the streams
# This creates a matrix where each row is a site and each
# column is the distance of that site to the 1026 
# different streams
dist <- sf::st_distance(sites, smallest_water)

# Remove the units (it's meters) from the distances
dist <- units::drop_units(dist)

# Make it a dataframe
dist_df <- as.data.frame(dist)

# Create a column that is the minimum value of each row
# allowing us to isolate the distance of each site to the 
# nearest stream
dist_df$water_dist <- apply(dist,1,FUN=min)

# Select out just the column of minimum distance to stream
# values for each site
select_dist <- select(dist_df, water_dist)

# Add site designations back in so this is not just a list
# of distance values and it can be joined with the other covariates
select_dist$Site <- sites$Site

# Add the distance to water covariate to the other
# scaled covariates

source("./R/scaled_covariates.R")

more_covariates <- inner_join(
  x= my_covariates,
  y= select_dist,
  by = "Site"
)

View(more_covariates)

# Save as .csv file -- this line will override 
# the covariate file made in the scaled_covariates.R script
# write.csv(more_covariates, "./data/site_covariates.csv", row.names = FALSE)


