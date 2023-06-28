#############
# Binary

ABWish<-function(){
  for (i in 1:Narm){
    y[i] ~ dbinom(mean[i],n[i])
    logit(mean[i]) <- mu[drug[i]] + v[study[i],drug[i]]
  }
  for (j in 1:Nstudy) { 
    v[j,1:Ndrug] ~ dmnorm(zero.AB[1:Ndrug], invR[1:Ndrug, 1:Ndrug]) }
  invR[1:Ndrug, 1:Ndrug] ~ dwish(Omega[1:Ndrug,1:Ndrug], Ndrug)
  R[1:Ndrug, 1:Ndrug] <- inverse(invR[ , ])
  for (k in 1:Ndrug){
    tau[k] <- sqrt(R[k,k])
  }
  for (j in 1:Ndrug){
    for (k in (j+1):Ndrug){
      rho[j,k] <- R[j,k]/(tau[j]*tau[k])
    }
  }
  for (k in 1:Ndrug) { mu[k] ~ dnorm(0, 0.001) }  
  for (k in 1:Ndrug) { lor[k] <- mu[k] - mu[1] }
  # ranking
  for (k in 1:Ndrug) { G[k] <- exp(mu[k])/(1+exp(mu[k])) }
  T.rank <- rank(G)
  for (k in 1:Ndrug) {
    rk[k] <- T.rank[k]
    best1[k] <- equals(rk[k],1)
    best2[k] <- equals(rk[k],2)
    best3[k]<-equals(rk[k],3)
    best12[k] <- best1[k] + best2[k]
  }
}


LARE<-function(){
  for (i in 1:Nstudy){
    for (j in 1:Narm[i]){
      y[i,j] ~ dbinom(mean[i,j],n[i,j])
      logit(mean[i,j]) <- mu[i] + lor[i,j]
    }
  }
  for (i in 1:Nstudy){
    w[i,1] <- 0
    lor[i,1] <- 0
    for (j in 2:Narm[i]){
      lor[i,j] ~ dnorm(md[i,j], inv[i,j])
      md[i,j] <- d[drug[i,j]] - d[drug[i,1]] + sw[i,j]
      #Incorporating information from overall treatment effect back 
      #into study specific treatment effect
      w[i,j] <- lor[i,j] - d[drug[i,j]] + d[drug[i,1]]
      sw[i,j] <- sum(w[i,1:(j-1)])/(j-1)
      #contrast-specific variance
      inv[i,j] <- inv.d*2*(j-1)/j
    }
  }
  for (j in 1:Nstudy) { 
    mu[j] ~ dnorm(0, 0.01) }
  d[1] <- 0
  for (k in 2:Ndrug) { 
    d[k] ~ dnorm(0, 0.01)
  }
  tau ~ dunif(0.01, 10)
  inv.d <- 1/pow(tau, 2)
  # ranking
  mp <- mean(mu[])
  for (k in 1:Ndrug) { G[k] <- exp(mp + d[k])/(1+exp(mp + d[k])) }
  T.rank <- rank(G)
  for (k in 1:Ndrug) {
    rk[k] <- T.rank[k]
    best1[k] <- equals(rk[k],1)
    best2[k] <- equals(rk[k],2)
    best3[k] <- equals(rk[k],3)
    best12[k] <- best1[k] + best2[k]
  }
}


CBWish<-function(){
  for (i in 1:Narm){
    y[i] ~ dbinom(mean[i],n[i])
    logit(mean[i]) <- mu[study[i]] + delta[study[i],drug[i]]*(1-equals(drug[i],1))
  }
  for (j in 1:Nstudy){
    delta[j,1:Ndrug] ~ dmnorm(d[1:Ndrug], invR[1:Ndrug, 1:Ndrug])
  }
  invR[1:Ndrug, 1:Ndrug] ~ dwish(Omega[1:Ndrug,1:Ndrug], Ndrug)
  R[1:Ndrug, 1:Ndrug] <- inverse(invR[ , ])
  for (k in 1:Ndrug){
    tau[k] <- sqrt(R[k,k])
  }
  for (j in 1:Ndrug){
    for (k in (j+1):Ndrug){
      rho[j,k] <- R[j,k]/(tau[j]*tau[k])
    }
  }
  for (j in 1:Nstudy) { mu[j] ~ dnorm(0, 0.01) }
  d[1] <- 0
  for (k in 2:Ndrug) { d[k] ~ dnorm(0, 0.01) }
  # ranking
  mp <- mean(mu[])
  for (k in 1:Ndrug) { T[k] <- exp(mp + d[k])/(1+exp(mp + d[k])) }
  T.rank <- rank(T)
  for (k in 1:Ndrug) {
    rk[k] <- T.rank[k]
    best1[k] <- equals(rk[k],1)
    best2[k] <- equals(rk[k],2)
    best3[k]<-equals(rk[k],3)
    best12[k] <- best1[k] + best2[k]
  }
}


#############
# Continuous case


ABWish_C<-function(){
  for (i in 1:Narm){
    # y[i] ~ dbinom(mean[i],n[i])
    # logit(mean[i]) <- mu[drug[i]] + v[study[i],drug[i]]
    y[i] ~ dnorm(mean[i], 1/pow(sigma[i]/sqrt(n[i]), 2))
    mean[i] <- mu[drug[i]] + v[study[i],drug[i]]
  }
  for (j in 1:Nstudy) { 
    v[j,1:Ndrug] ~ dmnorm(zero.AB[1:Ndrug], invR[1:Ndrug, 1:Ndrug]) }
  invR[1:Ndrug, 1:Ndrug] ~ dwish(Omega[1:Ndrug,1:Ndrug], Ndrug)
  R[1:Ndrug, 1:Ndrug] <- inverse(invR[ , ])
  for (k in 1:Ndrug){
    tau[k] <- sqrt(R[k,k])
  }
  for (j in 1:Ndrug){
    for (k in (j+1):Ndrug){
      rho[j,k] <- R[j,k]/(tau[j]*tau[k])
    }
  }
  for (k in 1:Ndrug) { mu[k] ~ dnorm(0, 0.001) }  
  for (k in 1:Ndrug) { lor[k] <- mu[k] - mu[1] }
  # # ranking
  # for (k in 1:Ndrug) { G[k] <- exp(mu[k])/(1+exp(mu[k])) }
  # T.rank <- rank(G)
  # for (k in 1:Ndrug) {
  #   rk[k] <- T.rank[k]
  #   best1[k] <- equals(rk[k],1)
  #   best2[k] <- equals(rk[k],2)
  #   best3[k]<-equals(rk[k],3)
  #   best12[k] <- best1[k] + best2[k]
  # }
}


LARE_C <-function(){
  for (i in 1:Nstudy){
    for (j in 1:Narm[i]){
      # y[i,j] ~ dbinom(mean[i,j],n[i,j])
      # logit(mean[i,j]) <- mu[i] + lor[i,j]
      y[i,j] ~ dnorm(mean[i,j], 1/pow(sigma[i,j]/sqrt(n[i,j]), 2))
      mean[i,j] <- mu[i] + lor[i,j]
    }
  }
  for (i in 1:Nstudy){
    w[i,1] <- 0
    lor[i,1] <- 0
    for (j in 2:Narm[i]){
      lor[i,j] ~ dnorm(md[i,j], inv[i,j])
      md[i,j] <- d[drug[i,j]] - d[drug[i,1]] + sw[i,j]
      #Incorporating information from overall treatment effect back 
      #into study specific treatment effect
      w[i,j] <- lor[i,j] - d[drug[i,j]] + d[drug[i,1]]
      sw[i,j] <- sum(w[i,1:(j-1)])/(j-1)
      #contrast-specific variance
      inv[i,j] <- inv.d*2*(j-1)/j
    }
  }
  for (j in 1:Nstudy) { 
    mu[j] ~ dnorm(0, 0.01) }
  d[1] <- 0
  for (k in 2:Ndrug) { 
    d[k] ~ dnorm(0, 0.01)
  }
  tau ~ dunif(0.01, 10)
  inv.d <- 1/pow(tau, 2)
  
  # # ranking
  # mp <- mean(mu[])
  # for (k in 1:Ndrug) { G[k] <- exp(mp + d[k])/(1+exp(mp + d[k])) }
  # T.rank <- rank(G)
  # for (k in 1:Ndrug) {
  #   rk[k] <- T.rank[k]
  #   best1[k] <- equals(rk[k],1)
  #   best2[k] <- equals(rk[k],2)
  #   best3[k] <- equals(rk[k],3)
  #   best12[k] <- best1[k] + best2[k]
  #}
}


CBWish_C <-function(){
  for (i in 1:Narm){
    # y[i] ~ dbinom(mean[i],n[i])
    # logit(mean[i]) <- mu[study[i]] + delta[study[i],drug[i]]*(1-equals(drug[i],1))
    y[i] ~ dnorm(mean[i], 1/pow(sigma[i]/sqrt(n[i]), 2))
    mean[i] <- mu[study[i]] + delta[study[i],drug[i]]*(1-equals(drug[i],1))
  }
  for (j in 1:Nstudy){
    delta[j,1:Ndrug] ~ dmnorm(d[1:Ndrug], invR[1:Ndrug, 1:Ndrug])
  }
  invR[1:Ndrug, 1:Ndrug] ~ dwish(Omega[1:Ndrug,1:Ndrug], Ndrug)
  R[1:Ndrug, 1:Ndrug] <- inverse(invR[ , ])
  for (k in 1:Ndrug){
    tau[k] <- sqrt(R[k,k])
  }
  for (j in 1:Ndrug){
    for (k in (j+1):Ndrug){
      rho[j,k] <- R[j,k]/(tau[j]*tau[k])
    }
  }
  for (j in 1:Nstudy) { mu[j] ~ dnorm(0, 0.01) }
  d[1] <- 0
  for (k in 2:Ndrug) { d[k] ~ dnorm(0, 0.01) }
  # ranking
  # mp <- mean(mu[])
  # for (k in 1:Ndrug) { T[k] <- exp(mp + d[k])/(1+exp(mp + d[k])) }
  # T.rank <- rank(T)
  # for (k in 1:Ndrug) {
  #   rk[k] <- T.rank[k]
  #   best1[k] <- equals(rk[k],1)
  #   best2[k] <- equals(rk[k],2)
  #   best3[k]<-equals(rk[k],3)
  #   best12[k] <- best1[k] + best2[k]
  # }
}