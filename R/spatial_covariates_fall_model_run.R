library(scales)
library(parallel)
library(dplyr)
library(nimble)
library(MCMCvis)

source("./R/spatial_covariates_fall_data_prep.R") 
source("./nimble/spatial_covariates_fall_model.R")



set.seed(21)


run_MCMC_allcode <- function(seed, data, cons) {
  
  # packages
  library(nimble)
  
  source("./nimble/spatial_covariates_fall_model.R",local=TRUE)
  
  my_inits <- function(){
    to_return <- list(
      delta_beta = rnorm(cons$ncovar_delta),
      gamma_beta = rnorm(cons$ncovar_gamma),
      phi_beta = rnorm(cons$ncovar_phi),
      psi_beta = rnorm(cons$ncovar_psi),
      rho_beta = rnorm(cons$ncovar_rho),
      delta_fall = rnorm(1),
      gamma_fall = rnorm(1),
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
    code=spatial_covariates_fall,
    constants = cons,
    data = data,
    inits = core_inits
  )
  
  CmyModel <- nimble::compileNimble(longest_shot)
  
  myMCMC <- nimble::buildMCMC(CmyModel)
  CmyMCMC <- nimble::compileNimble(myMCMC)
  

  results <- nimble::runMCMC(CmyMCMC, niter = 50000, nburnin = 10000, nchains = 1,
                     setSeed = seed, inits = core_inits,
                     progressBar = getNimbleOption("MCMCprogressBar"))

  
  return(results)
  
}



my_cluster <- makeCluster(4)

chain_output <- parLapply(cl = my_cluster, X = 1:4, 
                          fun = run_MCMC_allcode, 
                          data = data_list, cons = constant_list)



# It's good practice to close the cluster when you're done with it.
stopCluster(my_cluster)

# save the model output because this took a long time to run
saveRDS(chain_output, file = "../../GitHub/spatial_covariates_fall.RDS")


scf <- readRDS("../../GitHub/skunk_rds/spatial_covariates_fall2.RDS")

head(scf[[1]])

plot(scf[[1]][,2], type = "l")

# plot all this out with all chains

npar <- ncol(scf[[1]])
pdf("./fuzzy_plots/spatial_covariates_fall.pdf")
for(i in 1:npar){
  print(i)
  my_parm <- colnames(scf[[1]])[i]
  my_range <- range(
    sapply(
      scf,
      function(x) x[,i]
    )
  )
  
  plot(
    scf[[1]][,i],
    ylim = my_range,
    type = "l",
    main = my_parm
  )
  my_col <- c("red", "green", "blue")
  for(j in 2:4){
    lines(
      scf[[j]][,i],
      col = my_col[j-1]
    )
  }
}
dev.off()



MCMCvis::MCMCsummary(
  scf,
  digits=2
)





