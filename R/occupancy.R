

library(dplyr)
library(bbplot)
library(MCMCvis)
library(nimble)


output <- readRDS("./skunk_rds/spatial_covariates_fall2.RDS")

source("./R/spatial_covariates_fall_data_prep.R")


MCMCsummary(output, round=2)

source("./R/mcmc_functions.R")
output <- do.call("rbind", output)
mc <- split_mcmc(output)



#### get matrix of Pr(colonization) for each number of neighbors ####


our_m <- matrix(
  c(1,1,1,0,1,1,0,0,1),
  ncol = 3,
  nrow = 3,
  byrow=TRUE
)

delta_vec_fall <-log(1- plogis(mc$delta_beta[,c(1,1,1)] + mc$delta_fall[,1]))

delta_vec_fall <- delta_vec_fall %*% our_m
delta_prob_fall <- 1 - exp(delta_vec_fall)



delta_vec <-log(1- plogis(mc$delta_beta[,c(1,1,1)]))

delta_vec <- delta_vec %*% our_m
delta_prob <- 1 - exp(delta_vec)



phi <- log(1 - plogis(mc$phi_beta))


dprob_neigh1 <- delta_prob[,rep(1, ncol(phi))]
dprob_neigh2 <- delta_prob[,rep(2, ncol(phi))]
dprob_neigh3 <- delta_prob[,rep(3, ncol(phi))]

dprob_neigh1_fall <- delta_prob_fall[,rep(1, ncol(phi))]
dprob_neigh2_fall <- delta_prob_fall[,rep(2, ncol(phi))]
dprob_neigh3_fall <- delta_prob_fall[,rep(3, ncol(phi))]

occ1 <- dprob_neigh1 / (dprob_neigh1 + (1 - phi))
occ1_summary <- apply(
  occ1,
  2,
  quantile,
  probs = c(0.025,0.5,0.975)
)


occ2 <- dprob_neigh2 / (dprob_neigh2 + (1 - phi))
occ2_summary <- apply(
  occ2,
  2,
  quantile,
  probs = c(0.025,0.5,0.975)
)


occ3 <- dprob_neigh3 / (dprob_neigh3 + (1 - phi))
occ3_summary <- apply(
  occ3,
  2,
  quantile,
  probs = c(0.025,0.5,0.975)
)


occ1_fall <- dprob_neigh1_fall / (dprob_neigh1_fall + (1 - phi))
occ1_summary_fall <- apply(
  occ1_fall,
  2,
  quantile,
  probs = c(0.025,0.5,0.975)
)


occ2_fall <- dprob_neigh2_fall / (dprob_neigh2_fall + (1 - phi))
occ2_summary_fall <- apply(
  occ2_fall,
  2,
  quantile,
  probs = c(0.025,0.5,0.975)
)


occ3_fall <- dprob_neigh3_fall / (dprob_neigh3_fall + (1 - phi))
occ3_summary_fall <- apply(
  occ3_fall,
  2,
  quantile,
  probs = c(0.025,0.5,0.975)
)



# plot

windows(4, 3) 

# tiff(
#   "./plots/delta_figure.tiff",
#   height = 3,
#   width = 4,
#   units = "in",
#   res = 600,
#   compression = "lzw"
# )


{ 
  # blank plot
  par(mar = c(5, 2, 0.5, 0.5), oma = c(0, 4, 0, 0), lend = 1)
  bbplot::blank(xlim = c(0.5, 3.5), ylim = c(0, 0.8), bty = "l")
  
  # axis hash marks
  bbplot::axis_blank(1, at = seq(1,3, by = 1))
  bbplot::axis_blank(2)
  
  # axis hash mark labels
  bbplot::axis_text(text = seq(1,3,by=1), line = 0.9,
                    side = 1, at = seq(1,3,1))
  bbplot::axis_text(side = 2, las = 1, line = 0.4)
  
  # axis labels
  bbplot::axis_text("Number of neighboring sites \n with striped skunk", 
                    side = 1, line = 3)
  bbplot::axis_text("Pr(Occupancy)", side = 2, outer = TRUE, at = 0.65)
  
  # layer on data
  
  for(i in 1:3){
    lines(
      x = rep(i-0.125,2),
      y = occ1_summary[-2,i],
      col ="purple",
      lwd = 3
    )
  }
  points(
    x = c(1:3)-0.125,
    y = occ1_summary[2,],
    pch = 18,
    col = "black",
    cex = 2
  )
  points(
    x = c(1:3)-0.125,
    y = delta_prob_fall[2,],
    pch = 18,
    col = "purple",
    cex = 1.5
  )
  
  for(i in 1:3){
    lines(
      x = rep(i+0.125,2),
      y = delta_prob[-2,i],
      col ="darkgreen",
      lwd = 3
    )
  }
  points(
    x = c(1:3)+0.125,
    y = delta_prob[2,],
    pch = 19,
    col = "black",
    cex = 1.5
  )
  points(
    x = c(1:3)+0.125,
    y = delta_prob[2,],
    pch = 19,
    col = "darkgreen",
    cex = 1
  )
  
  # legend
  legend(
    x= 0.4,
    y = 0.85,
    legend = c("Fall", "Spring, Summer, Winter"),
    pch = c(18,19),
    col = c("purple", "darkgreen"),
    bty = "n",
    cex = 0.65,
    pt.cex = c(1.5,1)
  )
  
}

dev.off()













# isolate columns of Pr(colonization) by # of neighbors

neigh1 <- delta_prob[,1]
neigh1_fall <- delta_prob_fall[,1]

neigh2 <- delta_prob[,2]
neigh2_fall <- delta_prob_fall[,2]

neigh3 <- delta_prob[,3]
neigh3_fall <- delta_prob_fall[,3]



#### get Pr(persistence) matrix for each spatial covariate ####

##### phi urbanization ####

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






#### phi wateranization ####

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





#### phi openanization ####

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




#### Turn # neighbor vectors into phi sized matrices ####

dprob_neigh1 <- delta_prob[,rep(1, ncol(urb_phi))]
dprob_neigh2 <- delta_prob[,rep(2, ncol(urb_phi))]
dprob_neigh3 <- delta_prob[,rep(3, ncol(urb_phi))]

dprob_neigh


dprob <- dprob[,rep(1, ncol(phi))]

# my occupancy
occ1 <- dprob / (dprob + (1 - phi))
occ_summary <- apply(
  occ1,
  2,
  quantile,
  probs = c(0.025,0.5,0.975)
)

