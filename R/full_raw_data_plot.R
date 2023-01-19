
new_skunks <- readxl::read_excel("./data/working_full_skunk_data.xlsx")


skunks <- new_skunks[complete.cases(new_skunks),]


# Create spatial points
sites <- sf::st_as_sf(
  skunks,
  coords = c("Long", "Lat"),
  crs = 4326
)

#reproject to UTM
sites <- sf::st_transform(
  sites,
  crs=32616
)


sites <- sites %>% 
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





move_map_10 <- sites[sites$Season %in% c("AP10", "JU10","OC10"),]
move_map_11 <- sites[sites$Season %in% c("JA11", "AP11", "JU11","OC11"),]
move_map_12 <- sites[sites$Season %in% c("JA12", "AP12", "JU12","OC12"),]
move_map_13 <- sites[sites$Season %in% c("JA13", "AP13", "JU13","OC13"),]
move_map_14 <- sites[sites$Season %in% c("JA14", "AP14", "JU14","OC14"),]
move_map_15 <- sites[sites$Season %in% c("JA15", "AP15", "JU15","OC15"),]
move_map_16 <- sites[sites$Season %in% c("JU16","OC16"),]
move_map_17 <- sites[sites$Season %in% c("JA17", "AP17", "JU17","OC17"),]
move_map_18 <- sites[sites$Season %in% c("JA18", "AP18", "JU18","OC18"),]
move_map_19 <- sites[sites$Season %in% c("JA19", "AP19", "JU19","OC19"),]
move_map_20 <- sites[sites$Season %in% c("JA20", "AP20", "JU20","OC20"),]
move_map_21 <- sites[sites$Season %in% c("JA21", "AP21", "JU21","OC21"),]


move_map_10_obs <- move_map_10[move_map_10$Y > 0,]
move_map_10_not <- move_map_10[move_map_10$Y <= 0,]

move_map_11_obs <- move_map_11[move_map_11$Y > 0,]
move_map_11_not <- move_map_11[move_map_11$Y <= 0,]

move_map_12_obs <- move_map_12[move_map_12$Y > 0,]
move_map_12_not <- move_map_12[move_map_12$Y <= 0,]

move_map_13_obs <- move_map_13[move_map_13$Y > 0,]
move_map_13_not <- move_map_13[move_map_13$Y <= 0,]

move_map_14_obs <- move_map_14[move_map_14$Y > 0,]
move_map_14_not <- move_map_14[move_map_14$Y <= 0,]

move_map_15_obs <- move_map_15[move_map_15$Y > 0,]
move_map_15_not <- move_map_15[move_map_15$Y <= 0,]

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



# Make the plot
line = -2
cex = 2
adj  = 0.25
par(mfrow = c(4, 3))
par(cex=0.6)
par(bg="white")

plot(small_water$geometry)
plot(move_map_10_not$geometry, pch=20, col = move_map_10_not$color, cex=2, add=T)
plot(move_map_10_obs$geometry, pch=20, col = move_map_10_obs$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'10",cex.main=cex,col.main="red",font=2,line=line)

plot(small_water$geometry)
plot(move_map_11_not$geometry, pch=20, col = move_map_11_not$color, cex=2, add=T)
plot(move_map_11_obs$geometry, pch=20, col = move_map_11_obs$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'11",cex.main=cex,col.main="red",font=2,line=line)

plot(small_water$geometry)
plot(move_map_12_not$geometry, pch=20, col = move_map_12_not$color, cex=2, add=T)
plot(move_map_12_obs$geometry, pch=20, col = move_map_12_obs$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'12",cex.main=cex,col.main="red",font=2,line=line)

plot(small_water$geometry)
plot(move_map_13_not$geometry, pch=20, col = move_map_13_not$color, cex=2, add=T)
plot(move_map_13_obs$geometry, pch=20, col = move_map_13_obs$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'13",cex.main=cex,col.main="red",font=2,line=line)

plot(small_water$geometry)
plot(move_map_14_not$geometry, pch=20, col = move_map_14_not$color, cex=2, add=T)
plot(move_map_14_obs$geometry, pch=20, col = move_map_14_obs$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'14",cex.main=cex,col.main="red",font=2,line=line)

plot(small_water$geometry)
plot(move_map_15_not$geometry, pch=20, col = move_map_15_not$color, cex=2, add=T)
plot(move_map_15_obs$geometry, pch=20, col = move_map_15_obs$color, cex=2, add=T)
title(outer=outer,adj=adj,main="'15",cex.main=cex,col.main="red",font=2,line=line)

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
