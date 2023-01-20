



library(scales)
library(parallel)
library(dplyr)

source("./R/full_data_prep.R") 
source("./nimble/spatial_covariates_model.R")



set.seed(7)


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
  
  CmyModel <- compileNimble(longest_shot)
  
  myMCMC <- buildMCMC(CmyModel)
  CmyMCMC <- compileNimble(myMCMC)
  
  results <- runMCMC(CmyMCMC, niter = 20000, setSeed = seed, inits = core_inits)
  
  return(results)

}



my_cluster <- makeCluster(4)

chain_output <- parLapply(cl = my_cluster, X = 1:4, 
                          fun = run_MCMC_allcode, 
                          data = data_list, cons = constant_list)



# It's good practice to close the cluster when you're done with it.
stopCluster(my_cluster)

# save the model output because this took a long time to run
saveRDS(chain_output, file = "../skunk/spatial_covariates.RDS")


rds <- readRDS("../skunk/full_model_output.RDS")

head(rds[[1]])

plot(rds[[1]][,2], type = "l")

# plot all this out with all chains

npar <- ncol(rds[[1]])
pdf("./fuzzy_plots/first_run.pdf")
for(i in 1:npar){
  print(i)
  my_parm <- colnames(rds[[1]])[i]
  my_range <- range(
    sapply(
      rds,
      function(x) x[,i]
    )
  )
  
  plot(
    rds[[1]][,i],
    ylim = my_range,
    type = "l",
    main = my_parm
    )
  my_col <- c("red", "green", "blue")
  for(j in 2:4){
    lines(
      rds[[j]][,i],
      col = my_col[j-1]
    )
  }
}
dev.off()



MCMCvis::MCMCsummary(
  rds,
  digits=2
)




