



spatial_covariates_fall <- nimble::nimbleCode({
  for(i in 1:nsite){
    logit(psi[i,1]) <- inprod(
      psi_beta[1:ncovar_psi],
      X_psi[i, 1:ncovar_psi]
    )
    z[i,1] ~ dbern(psi[i,1])
    delta_bar[i,1] <- (1 - exp(
      inprod(
        (z[1:nsite, 1]*m[i,1:nsite]),
        log(1-d_vec[i,1:nsite,1]) 
      )
    ))
    for(t in 2:nseason){
      # indicator function for neighborhood (I)
      zeta[i,t-1] <- step(
        inprod(
          z[1:nsite, t-1],
          m[i, 1:nsite])-1
      )
      delta_bar[i,t] <- (1 - exp(
        inprod(
          (z[1:nsite,t]*m[i,1:nsite]),
          log(1-d_vec[i,1:nsite,t])
        )
      ))
      logit(phi[i, t-1]) <- inprod(
        phi_beta[1:ncovar_phi],
        X_phi[i, 1:ncovar_phi]
      )
      logit(gamma[i, t-1]) <- inprod(
        gamma_beta[1:ncovar_gamma],
        X_gamma[i, 1:ncovar_gamma]
      ) + gamma_fall*season_vec[t]
      lin_pred[i, t-1] <- ( 
        (z[i, t-1] * phi[i, t-1]) + 
          ((1-z[i, t-1])*zeta[i, t-1]*delta_bar[i, t-1]) + 
          ((1-z[i, t-1])*(1-zeta[i, t-1])*gamma[i, t-1]))
      z[i,t] ~ dbern(lin_pred[i, t-1])
    }
  }
  #observation model
  for(i in 1:nsite){
    for(tt in 1:nseason){
      logit(rho[i,tt]) <- inprod(
        rho_beta[1:ncovar_rho],
        X_rho[i,1:ncovar_rho]
      )
      y[i,tt] ~ dbin(rho[i,tt]*z[i,tt], J[i,tt])
    }
  }
  #priors
  for(psii in 1:ncovar_psi){
    psi_beta[psii] ~ dlogis(0,1)
  }
  for(gammai in 1:ncovar_gamma){
    gamma_beta[gammai] ~ dlogis(0,1)
  }
  for(phii in 1:ncovar_phi){
    phi_beta[phii] ~ dlogis(0,1)
  }
  for(rhoi in 1:ncovar_rho){
    rho_beta[rhoi] ~ dlogis(0,1)
  }
  for(deltai in 1:ncovar_delta){
    delta_beta[deltai] ~ dlogis(0,1)
  }
  for(i in 1:nsite){ 
    for(ii in 1:nsite){
      for(t in 1:nseason){
      logit(d_vec[i, ii, t]) <-inprod(
        delta_beta[1:ncovar_delta],
        delta_array[i,ii,1:ncovar_delta]
      ) + delta_fall*season_vec[t]
    }
    }
  }
  delta_fall ~ dlogis(0,1) 
  gamma_fall ~ dlogis(0,1)
})
