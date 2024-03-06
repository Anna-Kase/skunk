library(sf)
sf::sf_use_s2(FALSE)
library(dplyr)

# read in points, and then figure out neighbors
my_grid <- read.csv(
  "./data/raw_sim_covariates/point_locs.csv"
)
my_grid <- sf::st_as_sf(
  my_grid,
  coords = c("x","y"),
  crs = 32616
)

# buffer out 1001m
buff_grid <- sf::st_buffer(
  my_grid,
  1001
)

# get nearest neighbors
my_nn <- st_intersects(
  buff_grid,
  my_grid,
  sparse = TRUE
)
# remove the point of interest from the neighbor pool
for(i in 1:length(my_nn)){
  tmp <- my_nn[[i]]
  if(i %in% tmp){
    tmp <- tmp[-which(tmp == i)]
  }
  my_nn[[i]] <- tmp
}

# read in the OG covariates
og_covs <- read.csv(
  "./data/site_covariates.csv"
)
# remove urb
og_covs <- og_covs[,-which(colnames(og_covs) == "urb")]

cov_mu <- og_covs %>% 
  dplyr::summarise(
    dplyr::across(
      tidyselect:::where(
        is.numeric
      ),mean
    )
  )
cov_sd <- og_covs %>% 
  dplyr::summarise(
    dplyr::across(
      tidyselect:::where(
        is.numeric
      ),sd
    )
  )

# now read in the covariates and join them
c1 <- read.csv(
  "./data/raw_sim_covariates/urb_covars.csv"
)
c2 <- read.csv(
  "./data/raw_sim_covariates/dist2water.csv"
)
c3 <- read.csv(
  "./data/raw_sim_covariates/urban_openspace.csv"
)
grid_covs <- dplyr::inner_join(
  c1,
  c2,
  by = "Site"
)
grid_covs <- dplyr::inner_join(
  grid_covs,
  c3,
  by = "Site"
)
# scale as needed
for(i in 1:5){
  grid_covs[,i+1] <- (grid_covs[,i+1] - cov_mu[,i]) / cov_sd[,i]
}

# do the pca on the OG covs


scale_covs <- og_covs %>% 
  summarize_if(
    is.numeric, .funs = scale
  ) %>% 
  data.frame()

my_pca <- prcomp(scale_covs[,c("HU10", "tree", "imperv")], center=FALSE)

my_covariates$urb <- my_pca$x[,1]

loadings <- my_pca$rotation[,1]

study_area_pca <- as.matrix(
  grid_covs[,c("HU10","tree","imperv")]
) %*% loadings

grid_df <- read.csv(
  "./data/raw_sim_covariates/point_locs.csv"
)
grid_df$urb <- study_area_pca
grid_df$water_dist <- grid_covs$water_dist 
grid_df$open_dev <- grid_covs$open_dev

write.csv(
  grid_df,
  "./data/scaled_simulation_covariates.csv",
  row.names = FALSE
)

saveRDS(
  my_nn,
  "./data/simulation_neighbors.RDS"
)
