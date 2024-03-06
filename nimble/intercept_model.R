

my_model <- nimble::nimbleCode({
  for(i in 1:nsite){
    z[i,1] ~ dbern(psi)
    delta_bar[i,1] <- (1 - exp(
      inprod(
        (z[1:nsite,1]*m[i,1:nsite]),
        log(1-d_vec[1:nsite])
      )
    ))
    for(t in 2:nseason){
      # indicator function for neighborhood (I)
      zeta[i, t-1] <- step(
        inprod(
          z[1:nsite, t-1], 
          m[i, 1:nsite]
        ) -1
      )
      delta_bar[i,t] <- (1 - exp(
        inprod(
          (z[1:nsite,t]*m[i,1:nsite]),
          log(1-d_vec[1:nsite])
        )
      ))
      lin_pred[i, t-1] <- (z[i, t-1] * phi) + #persistence
       ((1-z[i, t-1])*zeta[i, t-1]*delta_bar[i, t-1]) + #neighborhood colonization
        ((1-z[i, t-1])*(1-zeta[i, t-1])*gamma) #out of neighborhood colonization
      z[i,t] ~ dbern(lin_pred[i,t-1])
    }
  }
  #observation model
  for(i in 1:nsite){
    for(tt in 1:nseason){
      y[i,tt] ~ dbin(rho*z[i,tt], J[i,tt])
    }
  }
  #priors
  psi ~ dbeta(1,1)
  rho ~ dbeta(1,1)
  phi ~ dbeta(1,1)
  gamma ~ dbeta(1,1)
  d ~ dbeta(1,1)
  for(i in 1:nsite){
    d_vec[i] <- d
  }
})

