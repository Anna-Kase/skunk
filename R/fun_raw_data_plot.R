


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

move_map_16 <- move_map[move_map$Season %in% c("JU16","OC16"),]
move_map_17 <- move_map[move_map$Season %in% c("JA17", "AP17", "JU17","OC17"),]
move_map_18 <- move_map[move_map$Season %in% c("JA18", "AP18", "JU18","OC18"),]
move_map_19 <- move_map[move_map$Season %in% c("JA19", "AP19", "JU19","OC19"),]
move_map_20 <- move_map[move_map$Season %in% c("JA20", "AP20", "JU20","OC20"),]
move_map_21 <- move_map[move_map$Season %in% c("JA21", "AP21", "JU21","OC21"),]


line = -2
cex = 2
adj  = 0.000025
par(mfrow = c(2, 3))
par(cex=0.6)
par(bg="white")
plot(small_water$geometry)
plot(move_map_16$geometry, pch=20, col = move_map_16$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'16",cex.main=cex,col.main="red",font=2,line=line)

plot(small_water$geometry)
plot(move_map_17$geometry, pch=20, col = move_map_17$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'17",cex.main=cex,col.main="red",font=2,line=line)

plot(small_water$geometry)
plot(move_map_18$geometry, pch=20, col = move_map_18$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'18",cex.main=cex,col.main="red",font=2,line=line)

plot(small_water$geometry)
plot(move_map_19$geometry, pch=20, col = move_map_19$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'19",cex.main=cex,col.main="red",font=2,line=line)

plot(small_water$geometry)
plot(move_map_20$geometry, pch=20, col = move_map_20$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'20",cex.main=cex,col.main="red",font=2,line=line)

plot(small_water$geometry)
plot(move_map_21$geometry, pch=20, col = move_map_21$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'21",cex.main=cex,col.main="red",font=2,line=line)





