library(sf)
library(prettymapr)
sf::sf_use_s2(FALSE)
library(bbplot)
library(dplyr)
library(FedData)
library(terra)
library(pals)

sites <- read.csv(
  "./data/complete_data.csv"
)

sites <- sites %>% 
  dplyr::select(Site, Long, Lat, Crs) %>% 
  dplyr::distinct(.)

# 4326 = wgs84
sites <- sf::st_as_sf(
  sites,
  coords = c("Long", "Lat"),
  crs = 4326
)

# 32616 = utm for Chicago area
sites <- sf::st_transform(
  sites,
  crs = 32616
)

county_path <- "D:/GIS/counties/cb_2020_us_county_500k/cb_2020_us_county_500k.shp"
county <- sf::st_read(
  county_path
)

county <- county[county$STATE_NAME == "Illinois",]

# Query down to the counties we need.
county <- county %>%
  filter(
    .,
    NAME %in% c(
      "McHenry",
      "Lake",
      "Kane",
      "Dupage",
      "Cook",
      "Kendall",
      "Will"
    )
  ) %>% 
  select(
    NAME
  )

# reproject to utms as well
county <- sf::st_transform(
  county,
  crs = sf::st_crs(sites)
)

# add a little bit of space around the sites and then
#  get a bounding box
site_buff <- sf::st_buffer(
  sites,
  1500
)

imperv <- FedData::get_nlcd(
  site_buff,
  year = 2019,
  dataset = "impervious",
  label = "Chicago_imperv",
  force.redo = TRUE
)

# reproject everything to imperv crs
county <- sf::st_transform(
  county,
  crs = sf::st_crs(imperv)
)

sites <- sf::st_transform(
  sites,
  crs = sf::st_crs(imperv)
)

# crop counties down to imperv layer
county_crop <- sf::st_crop(
  county,
  imperv
)

# get a white to black palette for plotting purposes.
imperv_pal <- data.frame(
  value = seq(0,99,1),
  color = pals::brewer.greys(100)
)




windows(4,4)

svg(
  "./plots/study_area_map_raw.svg",
  width = 4,
  height = 4
)
terra::plot(
  imperv,
  col = imperv_pal,
  legend = FALSE,
  box = FALSE,
  axes = FALSE,
  reset = FALSE
)


plot(
  sf::st_geometry(county_crop),
  add = TRUE,
  lwd = 3
)

points(
  sf::st_coordinates(sites),
  pch = 21,
  bg = "goldenrod",
  cex = 1
)
par(xpd = NA)
addscalebar(plotepsg = 5070, style = "ticks",
            lwd = 2, padin = c(0.34,-0.175), label.cex = 1)
addnortharrow(pos = "topright", padin = c(0.1,0.0), scale = 0.7)
#dev.off()

legend.col(imperv_pal$color, c(0,100))

windows(4,4)
pdf("./plots/study_area_legend_raw.pdf", height = 4, width = 4)
par(mar = c(1,5.5,1,5.5))
par(lend = 2)
legend_image <- as.raster(matrix(rev(imperv_pal$color), ncol=1))
plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
for(i in 1:5){
  lines(
    x = c(0.57, 0.6),
    y = rep(seq(0,1,length=5)[i],2)
  )
}
rasterImage(legend_image, 0.43, 0, 0.57,1)
rect(0.43, 0, 0.57,1, border = "black", lwd = 1)
text(x=0.59, y = seq(0,1,l=5), labels = seq(0,1,l=5), pos = 4)
dev.off()

