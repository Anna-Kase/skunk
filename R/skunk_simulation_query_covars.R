library(sf)
sf::sf_use_s2(FALSE)
library(dplyr)
library(uwinspatialtools)

# read in data and convert to UTMs
dat <- read.csv(
  "./data/full_capture_history.csv"
)
dat <- dplyr::distinct(
  dat[,c("Site", "Long", "Lat")]
)

dat <- sf::st_as_sf(
  dat,
  coords = c("Long", "Lat"),
  crs = 4326
)

dat <- sf::st_transform(
  dat,
  crs = 32616
)

# get bounding box
my_box <- sf::st_as_sfc(
  sf::st_bbox(
    dat
  )
)
# add a buffer of 1000m buffer
my_box <- sf::st_buffer(
  my_box,
  1000
)

# and then make a grid
my_grid <- sf::st_make_grid(
  my_box,
  1000,
  what = "centers"
)

# make into spatial data.frame
my_grid <- sf::st_as_sf(
  data.frame(
    sf::st_coordinates(
      my_grid
    )
  ),
  coords = c("X", "Y"),
  crs = 32616
)

# read in county layers
county <- sf::st_read(
  "D:/GIS/illinois_county",
  layer = "IL_BNDY_County_Py"
)

# Query down to the counties we need. We first filter down to 
#  county name and then select the single COUNTY_NAM (county name)
#  column from the data.
county <- county %>%
  filter(
    .,
    COUNTY_NAM %in% c(
      "MCHENRY",
      "LAKE",
      "KANE",
      "DUPAGE",
      "COOK",
      "KENDALL",
      "WILL"
    )
  ) %>% 
  select(
    COUNTY_NAM
  )

county <- sf::st_transform(
  county,
  crs = 32616
)

# remove points outside of our study area
my_grid <- sf::st_intersection(
  my_grid,
  county
)

# and reduce the county lines to the study area

county_crop <-  sf::st_crop(
  county,
  my_grid
) %>% 
  sf::st_buffer(., 0)



# get stuff for urbanization pca

my_raster_path <-
  "D:/GIS/cmap/landcover_2010_chicagoregion.img"


# read it in
my_map <- raster::raster(my_raster_path)

#  Extract the proportion canopy cover (lulc class 1)
#  and create 'impervious cover' value, which is the sum of multiple
#  lulc classes

my_grid$Site <- paste0(
  "Site_", stringr::str_pad(
    1:nrow(my_grid), 
    width = 4,
    pad = "0"
  )
)

# save the grid coordinates
tmp_coords <- sf::st_coordinates(my_grid)
tmp_df <- data.frame(
  Site = my_grid$Site,
  x = tmp_coords[,1],
  y = tmp_coords[,2],
  crs = 32616
)
write.csv(
  tmp_df,
  "./data/raw_sim_covariates/point_locs.csv",
  row.names = FALSE
)
# Split my_grid into pieces,
#  RAM allocations last resort.

my_grid <- split(
  my_grid,
  factor(floor(1:nrow(my_grid)/50))
)
lulc_prop <- vector(
  "list",
  length = length(my_grid)
)

for(i in 1:length(lulc_prop)){
  cat(
    paste(i,"of", length(lulc_prop),"\n\n")
  )
  lulc_prop[[i]] <- extract_raster_prop(
    my_points = my_grid[[i]],
    location_column = "Site",
    my_buffer = 1000,
    my_raster_data = my_map,
    lulc_cats = list(
      "tree" = 1,
      "imperv" = 5:7
    )
  )
  gc()
}

lulc_prop <- dplyr::bind_rows(
  lulc_prop
)


# clean up
rm(my_map)
gc()




# Load 2010 statewide census data for housing and population
# Data can be downloaded from:
# browseURL("http://silvis.forest.wisc.edu/data/housing-block-change/")

pop_data <- sf::st_read(
  "D:/GIS/housing_density",
  layer = "il_blk10_Census_change_1990_2010_PLA2"
)

# fix any potential issues with the vectors before trying
#  to summarise
pop_data <- sf::st_make_valid(
  pop_data
)

# Run function to calculate housing units, housing density, population
#  and population density.  Extract population data
#  within a 1km radius buffer.

population_data <- vector(
  "list",
  length = length(my_grid)
)

for(i in 1:length(my_grid)){
  cat(
    paste(i,"of", length(population_data),"\n\n")
  )
  population_data[[i]] <- extract_polygon(
    my_points = my_grid[[i]],
    location_column = "Site",
    my_buffer = 1000,
    my_shape = pop_data,
    layers = c("HU10")
  )
}
population_data <- dplyr::bind_rows(
  population_data
)
# Now join the landcover and the population data into one object
my_covariates <- dplyr::inner_join(
  population_data,
  lulc_prop,
  by = "Site"
)




# remove unecessary column(s) and drop units
my_covariates <- select(my_covariates, -LocationName.geometry)
my_covariates$HU10 <- units::drop_units(my_covariates$HU10)

write.csv(
  my_covariates,
  "./data/raw_sim_covariates/urb_covars.csv",
  row.names = FALSE
)

# get distance to nearest watersource

if(class(my_grid) == "list"){
  my_grid <- dplyr::bind_rows(
    my_grid
  )
}
water_path <- "D:/GIS/IL_streams"
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
  my_grid
)
# make valid once more after cropping
small_water <- sf::st_make_valid(
  small_water
)

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
  my_grid, 
  smallest_water
)

# Do element wise calculations, indexing
#  smallest_water by closest_streams_index
site_dist <- sf::st_distance(
  my_grid,
  smallest_water[closest_streams_index,],
  by_element = TRUE
)

# Remove the units (it's meters) from the distances
site_dist <- units::drop_units(
  site_dist
)

# Make it a dataframe
dist_df <- data.frame(
  Site = my_grid$Site,
  water_dist = site_dist
)

write.csv(
  dist_df,
  "./data/raw_sim_covariates/dist2water.csv",
  row.names = FALSE
)

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





# lawn covariate
# library(dplyr)
# library(sf)
# library(uwinspatialtools)

# source spatial points file
# source("./R/creating_spatial_points.R")


# Dowload data from:
#browseURL("https://www.mrlc.gov/data/nlcd-2019-land-cover-conus")
# load water data

my_lawn_path <-
  "D:/GIS/IL_nlcd/illinois_nlcd.tiff"


# read it in
my_lawn <- raster::raster(my_lawn_path)

#  Extract the Developed, Open Space 
#  ("Includes areas with a mixture of some constructed materials, 
#  but mostly vegetation in the form of lawn grasses. Impervious 
#  surfaces account for less than 20 percent of total cover. These 
#  areas most commonly include large-lot single-family housing units, 
#  parks, golf courses, and vegetation planted in developed settings 
#  for recreation, erosion control, or aesthetic purposes.") 
#  (lulc class 21)

# split grid back up again



my_grid <- split(
  my_grid,
  factor(floor(1:nrow(my_grid)/50))
)
lawn_prop <- vector(
  "list",
  length = length(my_grid)
)

for(i in 1:length(lawn_prop)){
  
  cat(
    paste(i,"of", length(lawn_prop),"\n\n")
  )
  
  lawn_prop[[i]] <- extract_raster_prop(
    my_points = my_grid[[i]],
    location_column = "Site",
    my_buffer = 1000,
    my_raster_data = my_lawn,
    lulc_cats = list(
      "open_dev" = 21
    ) 
  )
  
}
lawn_prop <- dplyr::bind_rows(
  lawn_prop
)

# remove unecessary column(s)
lawn_prop <- select(lawn_prop, -LocationName.geometry)

write.csv(
  lawn_prop,
  "./data/raw_sim_covariates/urban_openspace.csv",
  row.names = FALSE
)
