
library(dplyr)
library(bbplot)

output <- readRDS("./skunk_rds/spatial_covariates_fall2.RDS")

source("./R/spatial_covariates_fall_data_prep.R")

constant_list$X_gamma

library(MCMCvis)

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

urb_gamma <- mc$gamma_beta %*% urb_dm
urb_gamma <- plogis(urb_gamma)

urb_gamma <- apply(
  urb_gamma,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
)

urb_gamma <- t(urb_gamma)


urb_gamma_fall <- (mc$gamma_beta %*% urb_dm) + mc$gamma_fall[,1]
urb_gamma_fall <- plogis(urb_gamma_fall)

urb_gamma_fall <- apply(
  urb_gamma_fall,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
)

urb_gamma_fall <- t(urb_gamma_fall)


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

water_gamma <- mc$gamma_beta %*% water_dm
water_gamma <- plogis(water_gamma)

water_gamma <- apply(
  water_gamma,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
)

water_gamma <- t(water_gamma)


water_gamma_fall <- (mc$gamma_beta %*% water_dm) + mc$gamma_fall[,1]
water_gamma_fall <- plogis(water_gamma_fall)

water_gamma_fall <- apply(
  water_gamma_fall,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
)

water_gamma_fall <- t(water_gamma_fall)



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

open_gamma <- mc$gamma_beta %*% open_dm
open_gamma <- plogis(open_gamma)

open_gamma <- apply(
  open_gamma,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
)

open_gamma <- t(open_gamma)


open_gamma_fall <- (mc$gamma_beta %*% open_dm) + mc$gamma_fall[,1]
open_gamma_fall <- plogis(open_gamma_fall)

open_gamma_fall <- apply(
  open_gamma_fall,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
)

open_gamma_fall <- t(open_gamma_fall)


#### plotting ####

windows(7, 3)

tiff(
  "./plots/gamma_figure.tiff",
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
bbplot::axis_text("Pr(Colonization)", side = 2, outer = TRUE, at = 0.6, 
                  line = 2)
bbplot::ribbon(x = pretty_urb,
               y = urb_gamma[,-2], col = "lightblue4", alpha = 0.5)
bbplot::ribbon(x = pretty_urb, 
               y = urb_gamma_fall[,-2], col = "goldenrod", alpha = 0.5)
lines(x = pretty_urb, y = urb_gamma[,2], col = "lightblue4", lwd = 3)
lines(x = pretty_urb, y = urb_gamma_fall[,2],
      col = "goldenrod", lwd = 3, lty = 2)
legend(
   x= u[1],
   y = 0.87,
  legend = c("Fall", "Spring, Summer, Winter"),
  lty = c(2,1),
  lwd = 3,
  col = c("goldenrod", "lightblue4"),
  bty = "n",
  seg.len = 3.25
)

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
bbplot::axis_text("Distance to water (km)", side = 1, line = 2.5)
bbplot::axis_text("Pr(Colonization)", side = 2, outer = TRUE, at = 0.6, line = 2)
bbplot::ribbon(x = pretty_water,
               y = water_gamma[,-2], col = "lightblue4", alpha = 0.5)
bbplot::ribbon(x = pretty_water, 
               y = water_gamma_fall[,-2], col = "goldenrod", alpha = 0.5)
lines(x = pretty_water, y = water_gamma[,2], col = "lightblue4", lwd = 3)
lines(x = pretty_water, y = water_gamma_fall[,2],
      col = "goldenrod", lwd = 3, lty = 2)



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
bbplot::axis_text("Pr(Colonization)", side = 2, outer = TRUE, at = 0.6, 
                  line = 2)
bbplot::ribbon(x = pretty_open,
               y = open_gamma[,-2], col = "lightblue4", alpha = 0.5)
bbplot::ribbon(x = pretty_open, 
               y = open_gamma_fall[,-2], col = "goldenrod", alpha = 0.5)
lines(x = pretty_open, y = open_gamma[,2], col = "lightblue4", lwd = 3)
lines(x = pretty_open, y = open_gamma_fall[,2],
      col = "goldenrod", lwd = 3, lty = 2)
}
dev.off()
