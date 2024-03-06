

library(dplyr)
library(bbplot)
library(MCMCvis)


output <- readRDS("./skunk_rds/spatial_covariates_fall2.RDS")

source("./R/spatial_covariates_fall_data_prep.R")


MCMCsummary(output, round=2)

source("./R/mcmc_functions.R")
output <- do.call("rbind", output)
mc <- split_mcmc(output)



our_m <- matrix(
  c(1,1,1,0,1,1,0,0,1),
  ncol = 3,
  nrow = 3,
  byrow=TRUE
)

 delta_vec_fall <-log(1- plogis(mc$delta_beta[,c(1,1,1)] + mc$delta_fall[,1]))

 delta_vec_fall <- delta_vec_fall %*% our_m
 delta_prob_fall <- 1 - exp(delta_vec_fall)

 delta_prob_fall <- apply(delta_prob_fall, 2,
                     quantile,
                     probs = c(0.025,0.5,0.975))
 
 
 delta_vec <-log(1- plogis(mc$delta_beta[,c(1,1,1)]))
 
 delta_vec <- delta_vec %*% our_m
 delta_prob <- 1 - exp(delta_vec)
 
 delta_prob <- apply(delta_prob, 2,
                     quantile,
                     probs = c(0.025,0.5,0.975))

 
 
 
 #### plotting ####
 
windows(4, 3) 

tiff(
  "./plots/delta_figure.tiff",
  height = 3,
  width = 4,
  units = "in",
  res = 600,
  compression = "lzw"
)
 

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
bbplot::axis_text("Pr(Colonization)", side = 2, outer = TRUE, at = 0.65)

# layer on data

for(i in 1:3){
  lines(
    x = rep(i-0.125,2),
    y = delta_prob_fall[-2,i],
    col ="goldenrod",
    lwd = 3
  )
}
points(
  x = c(1:3)-0.125,
  y = delta_prob_fall[2,],
  pch = 20,
  col = "black",
  cex = 2
)
points(
  x = c(1:3)-0.125,
  y = delta_prob_fall[2,],
  pch = 20,
  col = "goldenrod",
  cex = 1.5
)

for(i in 1:3){
  lines(
    x = rep(i+0.125,2),
    y = delta_prob[-2,i],
    col ="lightblue4",
    lwd = 3
  )
}
points(
  x = c(1:3)+0.125,
  y = delta_prob[2,],
  pch = 18,
  col = "black",
  cex = 1.5
)
points(
  x = c(1:3)+0.125,
  y = delta_prob[2,],
  pch = 18,
  col = "lightblue4",
  cex = 1
)

# legend
legend(
  x= 0.4,
  y = 0.85,
  legend = c("Fall", "Spring, Summer, Winter"),
  pch = c(20,18),
  col = c("goldenrod", "lightblue4"),
  bty = "n",
  cex = 0.65,
  pt.cex = c(1.5,1)
)

}

dev.off()
