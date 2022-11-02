library(scales)

# The One Script to Run Them All
source("./R/parallel_data_prep.R") 
source("./nimble/three_cov_model.R")


library(parallel)




run_MCMC_allcode <- function(seed, data, cons) {
  
  # packages
  library(nimble)

  source("./nimble/three_cov_model.R",local=TRUE)
  
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
    code=three_cov_model,
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
  
  results <- runMCMC(Cmcmc, niter = 10000, setSeed = seed, inits = core_inits)
  
  return(results)

}


source("./R/parallel_data_prep.R") 

my_cluster <- makeCluster(4)

chain_output <- parLapply(cl = my_cluster, X = 1:4, 
                          fun = run_MCMC_allcode, 
                          data = data_list, cons = constant_list)

# It's good practice to close the cluster when you're done with it.
stopCluster(my_cluster)


# chain plots (caterpillars)
par(mfrow = c(2,2))

pdf("my_output.pdf")
for(i in 1:ncol(chain_output[[1]])){
  my_mcmc <- sapply(
    chain_output,
    function(x){
      x[,i]
    }
  )
  
  my_range <- range(my_mcmc)
  
  plot(my_mcmc[,1], ylim = my_range, main = colnames(chain_output[[1]])[i],
       col = alpha('black', 0.5), type = 'l', ylab = 'value', xlab = "step")
  lines(my_mcmc[,2], col = alpha('red', 0.5))
  lines(my_mcmc[,3], col = alpha('green', 0.5))
  lines(my_mcmc[,4], col = alpha('purple', 0.5))
}
dev.off()

