


move_map <- read.csv(
  "./data/full_capture_history.csv"
)


# Create spatial points
# You must put the correct CRS for respective city
move_map <- sf::st_as_sf(
  move_map,
  coords = c("Long", "Lat"),
  crs = 4326
)

#reproject to UTM
move_map <- sf::st_transform(
  move_map,
  crs=32616
)

plot(move_map)

move_map <- move_map %>% 
  mutate(color = case_when(Y > 0 ~ "green",
                           TRUE ~ "black"))

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


move_map_16 <- move_map[move_map$Season %in% c("JU16","OC16"),]
move_map_17 <- move_map[move_map$Season %in% c("JA17", "AP17", "JU17","OC17"),]
move_map_18 <- move_map[move_map$Season %in% c("JA18", "AP18", "JU18","OC18"),]
move_map_19 <- move_map[move_map$Season %in% c("JA19", "AP19", "JU19","OC19"),]
move_map_20 <- move_map[move_map$Season %in% c("JA20", "AP20", "JU20","OC20"),]
move_map_21 <- move_map[move_map$Season %in% c("JA21", "AP21", "JU21","OC21"),]

move_map_16_obs <- move_map_16[move_map_16$Y > 0,]
move_map_16_not <- move_map_16[move_map_16$Y <= 0,]

move_map_17_obs <- move_map_17[move_map_17$Y > 0,]
move_map_17_not <- move_map_17[move_map_17$Y <= 0,]

move_map_18_obs <- move_map_18[move_map_18$Y > 0,]
move_map_18_not <- move_map_18[move_map_18$Y <= 0,]

move_map_19_obs <- move_map_19[move_map_19$Y > 0,]
move_map_19_not <- move_map_19[move_map_19$Y <= 0,]

move_map_20_obs <- move_map_20[move_map_20$Y > 0,]
move_map_20_not <- move_map_20[move_map_20$Y <= 0,]

move_map_21_obs <- move_map_21[move_map_21$Y > 0,]
move_map_21_not <- move_map_21[move_map_21$Y <= 0,]


line = -2
cex = 2
adj  = 0.000025
par(mfrow = c(2, 3))
par(cex=0.6)
par(bg="white")
plot(small_water$geometry)
plot(move_map_16_not$geometry, pch=20, col = move_map_16_not$color, cex=2, add=T)
plot(move_map_16_obs$geometry, pch=20, col = move_map_16_obs$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'16",cex.main=cex,col.main="red",font=2,line=line)

plot(small_water$geometry)
plot(move_map_17_not$geometry, pch=20, col = move_map_17_not$color, cex=2, add=T)
plot(move_map_17_obs$geometry, pch=20, col = move_map_17_obs$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'17",cex.main=cex,col.main="red",font=2,line=line)

plot(small_water$geometry)
plot(move_map_18_not$geometry, pch=20, col = move_map_18_not$color, cex=2, add=T)
plot(move_map_18_obs$geometry, pch=20, col = move_map_18_obs$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'18",cex.main=cex,col.main="red",font=2,line=line)

plot(small_water$geometry)
plot(move_map_19_not$geometry, pch=20, col = move_map_19_not$color, cex=2, add=T)
plot(move_map_19_obs$geometry, pch=20, col = move_map_19_obs$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'19",cex.main=cex,col.main="red",font=2,line=line)

plot(small_water$geometry)
plot(move_map_20_not$geometry, pch=20, col = move_map_20_not$color, cex=2, add=T)
plot(move_map_20_obs$geometry, pch=20, col = move_map_20_obs$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'20",cex.main=cex,col.main="red",font=2,line=line)

plot(small_water$geometry)
plot(move_map_21_not$geometry, pch=20, col = move_map_21_not$color, cex=2, add=T)
plot(move_map_21_obs$geometry, pch=20, col = move_map_21_obs$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'21",cex.main=cex,col.main="red",font=2,line=line)





