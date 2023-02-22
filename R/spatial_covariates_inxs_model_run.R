library(scales)
library(parallel)
library(dplyr)
library(nimble)
library(MCMCvis)

source("./R/spatial_covariates_inxs_data_prep.R") 
source("./nimble/spatial_covariates_inxs.R")



set.seed(47)


run_MCMC_allcode <- function(seed, data, cons) {
  
  # packages
  library(nimble)
  
  source("./nimble/spatial_covariates_inxs.R",local=TRUE)
  
  my_inits <- function(){
    to_return <- list(
      delta_beta = rnorm(cons$ncovar_delta),
      gamma_beta = rnorm(cons$ncovar_gamma),
      phi_beta = rnorm(cons$ncovar_phi),
      psi_beta = rnorm(cons$ncovar_psi),
      rho_beta = rnorm(cons$ncovar_rho),
      delta_urbless = rnorm(1),
      gamma_urbless = rnorm(1),
      z = matrix(
        1,
        ncol = cons$nseason,
        nrow = cons$nsite
      )
    )
    return(to_return)
  }
  
  set.seed(seed)
  
  # get my inits list object
  core_inits <- my_inits()
  
  # fitting model mcmc
  longest_shot <- nimble::nimbleModel(
    code=spatial_covariates_inxs,
    constants = cons,
    data = data,
    inits = core_inits
  )
  
  CmyModel <- compileNimble(longest_shot)
  
  myMCMC <- buildMCMC(CmyModel)
  CmyMCMC <- compileNimble(myMCMC)
  
  results <- runMCMC(CmyMCMC, niter = 50000, nburnin = 10000, nchains = 1,
                     setSeed = seed, inits = core_inits)
  
  return(results)
  
}



my_cluster <- makeCluster(4)

chain_output <- parLapply(cl = my_cluster, X = 1:4, 
                          fun = run_MCMC_allcode, 
                          data = data_list, cons = constant_list)



# It's good practice to close the cluster when you're done with it.
stopCluster(my_cluster)

# save the model output 
saveRDS(chain_output, file = "../skunk_rds/spatial_covariates_inxs.rds")


sci <- readRDS("../skunk/spatial_covariates_inxs.rds")

head(sci[[1]])

plot(sci[[1]][,2], type = "l")

# plot all this out with all chains

npar <- ncol(sci[[1]])
pdf("./fuzzy_plots/spatial_covariates_inxs.pdf")
for(i in 1:npar){
  print(i)
  my_parm <- colnames(sci[[1]])[i]
  my_range <- range(
    sapply(
      sci,
      function(x) x[,i]
    )
  )
  
  
  plot(
    sci[[1]][,i],
    ylim = my_range,
    type = "l",
    main = my_parm
  )
  my_col <- c("red", "green", "blue")
  for(j in 2:4){
    lines(
      sci[[j]][,i],
      col = my_col[j-1]
    )
  }
}
dev.off()




MCMCvis::MCMCsummary(
  sci,
  digits=2
)






