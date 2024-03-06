

# Fall term only model run


library(scales)
library(parallel)
library(dplyr)
library(nimble)
library(MCMCvis)

source("./R/fall_term_data_prep.R") 
source("./nimble/fall_term_model.R")



set.seed(28)


run_MCMC_allcode <- function(seed, data, cons) {
  
  # packages
  library(nimble)
  
  source("./nimble/fall_term_model.R",local=TRUE)
  
  my_inits <- function(){
    to_return <- list(
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
    code=fall_only,
    constants = cons,
    data = data,
    inits = core_inits
  )
  
  CmyModel <- compileNimble(longest_shot)
  
  myMCMC <- buildMCMC(CmyModel)
  CmyMCMC <- compileNimble(myMCMC)
  
  results <- runMCMC(CmyMCMC, niter = 20000, nchains = 1,
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
saveRDS(chain_output, file = "../skunk_rds/fall_term.rds")


ft <- readRDS("../skunk_rds/fall_term.rds")

head(ft[[1]])

plot(ft[[1]][,2], type = "l")

# plot all this out with all chains

npar <- ncol(ft[[1]])
pdf("./fuzzy_plots/fall_term.pdf")
for(i in 1:npar){
  print(i)
  my_parm <- colnames(ft[[1]])[i]
  my_range <- range(
    sapply(
      ft,
      function(x) x[,i]
    )
  )
  
  plot(
    ft[[1]][,i],
    ylim = my_range,
    type = "l",
    main = my_parm
  )
  my_col <- c("red", "green", "blue")
  for(j in 2:4){
    lines(
      ft[[j]][,i],
      col = my_col[j-1]
    )
  }
}
dev.off()



MCMCvis::MCMCsummary(
  ft,
  digits=2
)


