
# how to calculate the probability of local colonizaiton using the output
# of the nimble dynamic occupancy model and depending on how many of the 
# "nearby" sites were previously occupied 

z <- matrix(
  c(1, 0, 0,
    1, 1, 0,
    1, 1, 1),
  byrow = TRUE,
  ncol = 3,
  nrow = 3
)

d_vec <- rep(0.18, 3)

1 - exp(z %*% log(1 - d_vec))
