

my_inits <- function(){
  to_return <- list(
    delta_beta = rnorm(constant_list$ncovar_delta),
    gamma_beta = rnorm(constant_list$ncovar_gamma),
    phi_beta = rnorm(constant_list$ncovar_phi),
    psi_beta = rnorm(constant_list$ncovar_psi),
    rho_beta = rnorm(constant_list$ncovar_rho),
    delta_fall = rnorm(1),
    gamma_fall = rnorm(1),
    delta_year = rnorm(1),
    z = matrix(
      1,
      ncol = constant_list$nseason,
      nrow = constant_list$nsite
    )
  )
  return(to_return)
}



# fitting model mcmc
longest_shot <- nimble::nimbleMCMC(
  code=spatial_covariates_fall_urbless,
  constants = constant_list,
  data = data_list,
  inits = my_inits,
  niter = 200,
  nchains = 2
)


MCMCvis::MCMCsummary(
  longest_shot,
  digits=2
)


plot(longest_shot[[1]][,2], type = "l")

npar <- ncol(longest_shot[[1]])
pdf("./fuzzy_plots/spatial_covariates_fall_urbless.pdf")
for(i in 1:npar){
  print(i)
  my_parm <- colnames(longest_shot[[1]])[i]
  my_range <- range(
    sapply(
      longest_shot,
      function(x) x[,i]
    )
  )
  
  plot(
    longest_shot[[1]][,i],
    ylim = my_range,
    type = "l",
    main = my_parm
  )
  my_col <- c("red", "green", "blue")
  for(j in 2:4){
    lines(
      longest_shot[[j]][,i],
      col = my_col[j-1]
    )
  }
}
dev.off()

