

# lawn covariate
library(dplyr)
library(sf)
library(uwinspatialtools)

# source spatial points file
source("./R/creating_spatial_points.R")



# Dowload data from:
# browseURL("https://www.mrlc.gov/data/nlcd-2019-land-cover-conus")
# load water data

my_lawn_path <-
  "../../GIS/mrlc_lc//nlcd_2019_land_cover_l48_20210604.img"


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


lawn_prop <- extract_raster_prop(
  my_points = sites,
  location_column = "Site",
  my_buffer = 1000,
  my_raster_data = my_lawn,
  lulc_cats = list(
    "open_dev" = 21
  ) 
)

# clean up
rm(my_lawn)
gc()

# remove unecessary column(s)
lawn_prop <- select(lawn_prop, -LocationName.geometry)


# Add the managed lawn covariate to the other
# covariates
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

# check if lawn column on here, if not add it
if(!"open_dev" %in% colnames(my_covariates)){
  more_covariates <- inner_join(
    x= my_covariates,
    y= lawn_prop,
    by = "Site"
  )
  write.csv(
    more_covariates,
    "./data/site_covariates.csv",
    row.names = FALSE
  )
}


