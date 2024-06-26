
# Spatial covariates model run
library(scales)

# The One Script to Run Them All
source("./R/parallel_data_prep.R") 
source("./nimble/three_cov_model.R")



library(scales)
library(parallel)
library(dplyr)
library(nimble)
library(MCMCvis)

source("./R/spatial_covariates_data_prep.R") 
source("./nimble/spatial_covariates_model.R")

run_MCMC_allcode <- function(seed, data, cons) {
  
  # packages
  library(nimble)

  source("./nimble/spatial_covariates_model.R",local=TRUE)
  
  my_inits <- function(){
    to_return <- list(
      delta_beta = rnorm(cons$ncovar_delta),
      gamma_beta = rnorm(cons$ncovar_gamma),
      phi_beta = rnorm(cons$ncovar_phi),
      psi_beta = rnorm(cons$ncovar_psi),
      rho_beta = rnorm(cons$ncovar_rho),
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
    code=spatial_covariates,
    constants = cons,
    data = data,
    inits = core_inits
  )
  conf <- configureMCMC(
    longest_shot, 
    monitors = c(
      "delta_beta",
      "gamma_beta",
      "phi_beta",
      "psi_beta",
      "rho_beta"), 
    thin = 1
  )
  Rmcmc <- buildMCMC(conf)
  compiledList <- compileNimble(longest_shot, Rmcmc)
  Cmcmc <- compiledList$Rmcmc
  
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

# save the model output because this took a long time to run
saveRDS(chain_output, file = "../skunk_rds/spatial_covariates.rds")


sc <- readRDS("../skunk_rds/spatial_covariates.rds")

head(sc[[1]])

plot(sc[[1]][,2], type = "l")

# plot all this out with all chains

npar <- ncol(rds[[1]])
pdf("./fuzzy_plots/spatial_covariates.pdf")
for(i in 1:npar){
  print(i)
  my_parm <- colnames(sc[[1]])[i]
  my_range <- range(
    sapply(
      sc,
      function(x) x[,i]
    )
  )
  
  plot(
    sc[[1]][,i],
    ylim = my_range,
    type = "l",
    main = my_parm
    )
  my_col <- c("red", "green", "blue")
  for(j in 2:4){
    lines(
      sc[[j]][,i],
      col = my_col[j-1]
    )
  }
}
dev.off()


# model summary
MCMCvis::MCMCsummary(
  sc,
  digits=2
)





