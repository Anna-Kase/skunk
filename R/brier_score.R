


# Brier Score

# Using the formula from the Dietze Ecological Forecasting book;
# I can produce a matrix that contains the Brier Score (maybe)
# for each site and season (y matrix) for the intercept only model.

BS <- (((mean(long_shot[,4]) - y[])^2) * 2)/nsite

# This works but it doesn't feel right...
# These scores are supposed to be between 0 and 1 and I have
# a handfull in the matrix that are greater than 1.
# I also tried to automate this process using for loops but that would
# produce a null object. I don't know if I am supposed to be using the 
# mean of psi as the "forecast probability", but I think there should
# be a forecast probability for each site? Like, somehow the original
# dataset should be able to be incorporated into the posterior, but now
# that I am typing that out it also doesn't feel right... 

# Definitely needs more work...









