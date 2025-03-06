

library(dplyr)
library(bbplot)
library(MCMCvis)

output <- readRDS("./skunk_rds/spatial_covariates_fall2.RDS")

source("./R/spatial_covariates_fall_data_prep.R")

MCMCsummary(output, round=2)

source("./R/mcmc_functions.R")
output <- do.call("rbind", output)
mc <- split_mcmc(output)


sc <- read.csv("./data/site_covariates.csv")

##### urbanization ####

ubounds <- range(sc$urb)
ubounds
pretty_urb <- seq(-3, 4, length.out = 200)
urb_dm <- rbind(
  1,
  seq(-3, 4, length.out = 200),
  0,
  0
)

urb_phi <- mc$phi_beta %*% urb_dm
urb_phi <- plogis(urb_phi)

urb_phi <- apply(
  urb_phi,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
)

urb_phi <- t(urb_phi)




#### wateranization ####

ubounds <- range(sc$water_dist)
ubounds
pretty_water <- seq(0, 14000, length.out = 1000)
water_dm <- rbind(
  1,
  0,
  (pretty_water - mean(sc$water_dist))/sd(sc$water_dist),
  0
)

water_phi <- mc$phi_beta %*% water_dm
water_phi <- plogis(water_phi)

water_phi <- apply(
  water_phi,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
)

water_phi <- t(water_phi)




#### openanization ####

ubounds <- range(sc$open_dev)
ubounds
pretty_open <- seq(0, 0.4, length.out = 100)
open_dm <- rbind(
  1,
  0,
  0,
  (pretty_open - mean(sc$open_dev))/sd(sc$open_dev)
)

open_phi <- mc$phi_beta %*% open_dm
open_phi <- plogis(open_phi)

open_phi <- apply(
  open_phi,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
)

open_phi <- t(open_phi)



#### plotting ####

windows(7, 3)

tiff(
  "./plots/phi_figure2.tiff",
  height = 3,
  width = 7,
  units = "in",
  res = 600,
  compression = "lzw"
)
m <- matrix(1:3, ncol = 3)
layout(m)

{
  par(mar = c(5, 2, 0.5, 0.5), oma = c(0, 5, 0, 0), lend = 1)
  bbplot::blank(xlim = range(pretty_urb), ylim = c(0, 1), bty = "l")
  u <- par("usr")
  
  text(
    x = u[1] +0.1,
    y = u[4] - 0.1,
    labels = "A)",
    cex = 1.2,
    pos = 4
  )
  bbplot::axis_blank(1, at = seq(-3, 4, by = 1))
  bbplot::axis_blank(2)
  bbplot::axis_text(text = seq(-3, 4, 1), side = 1, line = 0.9,
                    at = seq(-3, 4, 1))
  bbplot::axis_text(side = 2, las = 1, line = 0.4)
  bbplot::axis_text("Urbanization", side = 1, line = 2.5)
  bbplot::axis_text("Pr(Persistence)", side = 2, outer = TRUE, at = 0.6, 
                    line = 2)
  bbplot::ribbon(x = pretty_urb,
                 y = urb_phi[,-2], col = "goldenrod", alpha = 0.5)
  lines(x = pretty_urb, y = urb_phi[,2], col = "goldenrod", lwd = 3)


  bbplot::blank(xlim = range(pretty_water), ylim = c(0, 1), bty = "l")
  u <- par("usr")
  
  text(
    x = u[1] + (u[1] * (1-0.93)),
    y = u[4] - 0.1,
    labels = "B)",
    cex = 1.2,
    pos = 4
  )
  bbplot::axis_blank(1, at = c(0, 4000, 9000, 14000))
  bbplot::axis_blank(2)
  bbplot::axis_text(c(0, 4, 9, 14), side = 1, line = 0.9, at = c(0, 4000, 9000, 14000)) 
  bbplot::axis_text("Distance to stream\nor river (km)", side = 1, line = 4)
  bbplot::axis_text("Pr(Persistence)", side = 2, outer = TRUE, at = 0.6, line = 2)
  bbplot::ribbon(x = pretty_water,
                 y = water_phi[,-2], col = "goldenrod", alpha = 0.5)
  lines(x = pretty_water, y = water_phi[,2], col = "goldenrod", lwd = 3)

  
  bbplot::blank(xlim = range(pretty_open), ylim = c(0, 1), bty = "l")
  u <- par("usr")
  
  text(
    x = u[1] + (u[1] * (1-0.93)),
    y = u[4] - 0.1,
    labels = "C)",
    cex = 1.2,
    pos = 4
  )
  bbplot::axis_blank(1, at = seq(0, 0.4, by = 0.1))
  bbplot::axis_blank(2)
  bbplot::axis_text(text = seq(0, 0.4, 0.1), side = 1, line = 0.9,
                    at = seq(0, 0.4, 0.1))
  bbplot::axis_text("Urban open space", side = 1, line = 2.5)
  bbplot::axis_text("Pr(Persistence)", side = 2, outer = TRUE, at = 0.6, 
                    line = 2)
  bbplot::ribbon(x = pretty_open,
                 y = open_phi[,-2], col = "goldenrod", alpha = 0.5)
  lines(x = pretty_open, y = open_phi[,2], col = "goldenrod", lwd = 3)
}
dev.off()


