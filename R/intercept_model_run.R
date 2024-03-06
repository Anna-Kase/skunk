

# Intercept only model

library(scales)
library(parallel)
library(dplyr)
library(nimble)
library(MCMCvis)

source("./R/intercept_data_prep.R") 
source("./nimble/intercept_model.R")



set.seed(7)


run_MCMC_allcode <- function(seed, data, cons) {
  
  # packages
  library(nimble)
  
  source("./nimble/intercept_model.R",local=TRUE)
  
  
  set.seed(seed)
  
  
  # fitting model mcmc
  long_shot <- nimble::nimbleMCMC(
    code=my_model,
    constants = cons,
    data = data,
    monitors = c("psi", "rho", "d", "gamma", "phi")
  )
  
  CmyModel <- compileNimble(long_shot)
  
  myMCMC <- buildMCMC(CmyModel)
  CmyMCMC <- compileNimble(myMCMC)
  
  results <- runMCMC(CmyMCMC, niter = 50000, nburnin = 10000, nchains = 1,
                     setSeed = seed)
  
  return(results)
  
}



my_cluster <- makeCluster(4)

chain_output <- parLapply(cl = my_cluster, X = 1:4, 
                          fun = run_MCMC_allcode, 
                          data = data_list, cons = constant_list)



# It's good practice to close the cluster when you're done with it.
stopCluster(my_cluster)

# save the model output 
saveRDS(chain_output, file = "../skunk_rds/intercept.rds")


int <- readRDS("../skunk_rds/intercept.rds")

head(int[[1]])

plot(int[[1]][,2], type = "l")

# plot all this out with all chains

npar <- ncol(int[[1]])
pdf("./fuzzy_plots/intercept.pdf")
for(i in 1:npar){
  print(i)
  my_parm <- colnames(int[[1]])[i]
  my_range <- range(
    sapply(
      int,
      function(x) x[,i]
    )
  )
  
  plot(
    int[[1]][,i],
    ylim = my_range,
    type = "l",
    main = my_parm
  )
  my_col <- c("red", "green", "blue")
  for(j in 2:4){
    lines(
      int[[j]][,i],
      col = my_col[j-1]
    )
  }
}
dev.off()


# model summary

MCMCvis::MCMCsummary(
  int,
  digits=2
)

