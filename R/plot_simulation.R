library(sf)
library(prettymapr)
sf::sf_use_s2(FALSE)
library(bbplot)

plot_array <- readRDS("./data/phi_gamma_delta_by_season.RDS")

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
  crs = 32616
)
gr$occ_prob <- gmu_z
gr$occ_sd <- occ_sd






# Changes in phi, gamma, delta over simulated time


windows(7, 3)

svg(
  "./plots/colonization_changes.svg",
  height = 3,
  width = 7
)

{
  par(mar = c(3, 4, 0.5, 10), oma = c(0,0,0,0), lend = 1)
  bbplot::blank(xlim = c(1,11), ylim = c(0, 0.6), bty = "l")
  
  bbplot::axis_blank(1, at = seq(1, 11, by = 1))
  bbplot::axis_blank(2)
  bbplot::axis_text(text = seq(1, 11, 1), side = 1, line = 0.3,
                    at = seq(1, 11, 1))
  bbplot::axis_text(side = 2, las = 1, line = 0.4)
  bbplot::axis_text("Season", side = 1, line = 1.5, cex = 1.2)
  bbplot::axis_text("Proportion of sites", side = 2, at = 0.3, 
                    line = 2.25, cex = 1.2)
  my_col <- c("goldenrod", "black", "lightblue4")
  my_jiggle <- c(-0.125,0, 0.125)
  my_pch <- c(21,22,23)
  for(i in 1:3){
    for(t in 1:11){
      lines(
        x = rep(t + my_jiggle[i],2),
        y = plot_array[-2,i,t],
        lwd = 2,
        col = my_col[i]
      )
    }
    points(
      x = c(1:11) + my_jiggle[i],
      y = plot_array[2,i,],
      bg = my_col[i],
      cex = 1.2,
      pch = my_pch[i]
    )
    
  }
  par(xpd = NA)
  legend(
    x = 11.5,
    y = 0.5,
    legend = c(
      "Persistence",
      "Random colonization",
      "Local colonization"
    ),
    pch = my_pch,
    pt.bg = my_col,
    pt.cex = 1.2,
    bty = "n",
    y.intersp = 1.5
  )
  
}

dev.off()



