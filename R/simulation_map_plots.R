

library(sf)
library(prettymapr)
sf::sf_use_s2(FALSE)
library(bbplot)



utm_crs <- 32616

oz <- read.csv(
  "./data/original_z_sim_values.csv"
)
nsamp <- 10000

# average occupancy each site and season
mu_z <- oz/(nsamp)

# overall occupancy throughout study
gmu_z <- rowMeans(mu_z)

# variation in occupancy across seasons
occ_sd <- apply(mu_z, 1, sd)

gr <- read.csv(
  "./data/raw_sim_covariates/point_locs.csv"
)

gr <- sf::st_as_sf(
  gr,
  coords = c("x","y"),
  crs = utm_crs
)
gr$occ_prob <- gmu_z
gr$occ_sd <- occ_sd


# read in county maps

county_path <- "../../GIS/county_maps/illinois_county/"
county <- sf::st_read(
  county_path,
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
  crs = utm_crs
)


# and reduce the county lines to the study area

county_crop <-  sf::st_crop(
  county,
  gr
)


# read in stream data

# Dowload data from:
# browseURL("https://clearinghouse.isgs.illinois.edu/data/hydrology/streams-and-shorelines")
# load water data



water <- sf::read_sf("../../GIS/water/IL_Streams_From_100K_DLG_Ln.shp")



# reproject into the UTM
water <- sf::st_transform(
  water,
  crs=utm_crs
)



# Selecting one column of data to work with that provides 
# information on if a linestring is a stream and selecting
# only objects that have the designation for being a stream.
# Interpretation of values on the metadata tab at same URL the data
# was downloaded from under "Entity_and_Attribute_Information"
water <- dplyr::select(
  water,
  ENR
)

water <- water %>% 
  dplyr::filter(ENR >= 400000) 


water_crop <-  sf::st_intersection(
  water,
  county_crop
) 


# make a buffer to plot on maps to make a nice smooth edge
# around the study area
all_counties <- sf::st_union(county_crop)
county_buff <- sf::st_buffer(all_counties, 2000)
county_donut <- sf::st_sym_difference(county_buff, all_counties)



#### occupancy map ####


windows(3.5, 3.5)

svg(
  "./plots/mean_occ.svg",
  height = 3.5,
  width = 3.5
)

par(mar = c(1,1,1,1))


my_pal <- pals::ocean.speed(11)

plot(
  gr["occ_prob"],
  pch = 19,
  cex = 0.8,
  pal = my_pal,
  nbreaks = 10,
  reset = FALSE,
  main = ""
)

plot(st_geometry(water_crop),
     add=TRUE,
     col = scales::alpha("white", 0.5),
     lwd = 1.4
)
plot(st_geometry(county_donut),
     col = "white",
     border = NA,
     add = TRUE
)

plot(st_geometry(county_crop),
     col = "NA",
     add = TRUE,
     lwd = 2
)

par(xpd = NA)
addscalebar(plotepsg = utm_crs, style = "ticks",
            lwd = 2, padin = c(1.07,-0.175), label.cex = 1)
addnortharrow(pos = "topright", padin = c(0.1,0.19), scale = 0.7)



dev.off()



#### standard deviation map ####


windows(3.5, 3.5)

svg(
  "./plots/sd_occ.svg",
  height = 3.5,
  width = 3.5
)

par(mar = c(1,1,1,1))


my_pal <- pals::brewer.blues(9)

plot(
  gr["occ_sd"],
  pch = 19,
  cex = 0.8,
  pal = my_pal,
  nbreaks = 8,
  reset = FALSE,
  main = ""
)

plot(st_geometry(water_crop),
     add=TRUE,
     col = scales::alpha("white", 0.5),
     lwd = 1.4
)
plot(st_geometry(county_donut),
     col = "white",
     border = NA,
     add = TRUE
)

plot(st_geometry(county_crop),
     col = "NA",
     add = TRUE,
     lwd = 2
)

par(xpd = NA)
addscalebar(plotepsg = utm_crs, style = "ticks",
            lwd = 2, padin = c(1.07,-0.175), label.cex = 1)
addnortharrow(pos = "topright", padin = c(0.1,0.19), scale = 0.7)



dev.off()






#### mixture map ####


bb <- gr %>% 
  dplyr::mutate(
    op_cat = cut(occ_prob, breaks = 3),
    os_cat = cut(occ_sd, breaks = 3)
  )

# number of categories
ncat <- length(
  levels(
    bb$op_cat
  )
)

# names of the different scales (to match
#  colors)
scale_names <- paste0(
  rep(1:ncat, each = ncat), " - ", rep(1:ncat, ncat)
)

# the colors
scale_vals <- pals::brewer.seqseq2()

# make group names to match with scale_names
bb$group <- paste0(
  as.numeric(bb$op_cat)," - ", as.numeric(bb$os_cat)
)
# give them the scale_names as levels
bb$group <- factor(
  bb$group,
  levels = scale_names
)

# convert to numeric to match to colors
bb$color <- scale_vals[
  as.numeric(bb$group)
]



windows(3.5, 3.5)

svg(
  "./plots/bivariate.svg",
  height = 3.5,
  width = 3.5
)

par(mar = c(1,1,1,1))

plot(st_geometry(county_crop))

plot(
  st_geometry(bb),
  pch = 19,
  add = TRUE,
  cex = 0.8,
  col = bb$color
)

  plot(st_geometry(water_crop), 
     add=TRUE,
     col = scales::alpha("white", 0.5),
     lwd = 1.4
)
plot(st_geometry(county_donut),
     col = "white",
     border = NA,
     add = TRUE
)

plot(st_geometry(county_crop), 
     col = "NA",
     add = TRUE,
     lwd = 2
)

par(xpd = NA)
addscalebar(plotepsg = utm_crs, style = "ticks",
            lwd = 2, padin = c(1.07,-0.175), label.cex = 1)
addnortharrow(pos = "topright", padin = c(0.1,0.19), scale = 0.7)



dev.off()



#### county reference lines ####


windows(3.5, 3.5)

svg(
  "./plots/county_reference.svg",
  height = 3.5,
  width = 3.5
)

par(mar = c(1,1,1,1))



plot(st_geometry(water_crop),
     col = scales::alpha("white", 0.5),
     lwd = 1.4
)
plot(st_geometry(county_donut),
     col = "white",
     border = NA,
     add = TRUE
)

plot(st_geometry(county_crop),
     col = "NA",
     add = TRUE,
     lwd = 2
)

par(xpd = NA)
addscalebar(plotepsg = utm_crs, style = "ticks",
            lwd = 2, padin = c(1.07,-0.175), label.cex = 1)
addnortharrow(pos = "topright", padin = c(0.1,0.19), scale = 0.7)



dev.off()




