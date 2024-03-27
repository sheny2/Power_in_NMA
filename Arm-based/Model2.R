ABWish.het.cor<-function(){
  for (i in 1:Narm){
    y[i] ~ dbinom(mean[i],n[i])
    logit(mean[i]) <- mu[drug[i]] + v[study[i],drug[i]]
  }
  
  for (j in 1:Nstudy) { 
    v[j,1:Ndrug] ~ dmnorm(zero.AB[1:Ndrug], invR[1:Ndrug, 1:Ndrug]) }
  invR[1:Ndrug, 1:Ndrug] ~ dwish(Omega[1:Ndrug,1:Ndrug], Ndrug + 1)
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
  for (k in 1:Ndrug) { lor[k] <- mu[k] - mu[1] } # record LOR
  for (k in 1:Ndrug) { or[k] <- exp(lor[k])} # record OR

}



ABWish.het.eqcor <-function(){
  for (i in 1:Narm){
    y[i] ~ dbinom(mean[i],n[i])
    logit(mean[i]) <- mu[drug[i]] + v[study[i],drug[i]]
  }
  for (j in 1:Nstudy) { 
    v[j,1:Ndrug] ~ dmnorm(zero.AB[1:Ndrug], COV_mat[1:Ndrug, 1:Ndrug]) }
  
  for(j in 1:Ndrug){
    sigma[j] <- 1/sqrt(inv.sig.sq[j])
    inv.sig.sq[j] ~ dgamma(2, 1)
  }
  
  for(j in 1:Ndrug){
    for(k in 1:Ndrug){ 
      COV_mat[j,k] <- 1/sigma[j]*1/sigma[k]*ifelse(j == k, diag, offdiag)
    }
  }
  
  diag <- (1 + (Ndrug - 2)*rho)/(1 + (Ndrug - 2)*rho - (Ndrug - 1)*rho^2)
  offdiag <- (-rho/(1 + (Ndrug - 2)*rho - (Ndrug - 1)*rho^2))
  rho ~ dunif(-1/(Ndrug - 1), 0.9999)
  
  for (k in 1:Ndrug) { mu[k] ~ dnorm(0, 0.001) }  
  for (k in 1:Ndrug) { lor[k] <- mu[k] - mu[1] }
  for (k in 1:Ndrug) { or[k] <- exp(lor[k])}
}




NMA.het.eqcor <-function(){
  for (i in 1:Narm){
    y[i] ~ dbinom(mean[i],n[i])
    logit(mean[i]) <- mu[drug[i]] + v[study[i],drug[i]]
  }
  for (j in 1:Nstudy) { 
    v[j,1:Ndrug] ~ dmnorm(zero.AB[1:Ndrug], COV_mat[1:Ndrug, 1:Ndrug]) }
  
  
  for(j in 1:Ndrug){
    sigma[j] <- 1/sqrt(inv.sig.sq[j])
    inv.sig.sq[j] ~ dgamma(0.001, 0.001)
  }
  
  for(j in 1:Ndrug){
    for(k in 1:Ndrug){ 
      COV_mat[j,k] <- 1/sigma[j]*1/sigma[k]*ifelse(j == k, diag, offdiag)
    }
  }
  
  diag <- (1 + (Ndrug - 2)*rho)/(1 + (Ndrug - 2)*rho - (Ndrug - 1)*rho^2)
  offdiag <- (-rho/(1 + (Ndrug - 2)*rho - (Ndrug - 1)*rho^2))
  rho ~ dunif(-1/(Ndrug - 1), 0.9999)
  
  for (k in 1:Ndrug) { mu[k] ~ dnorm(0, 0.001) }  
  for (k in 1:Ndrug) { lor[k] <- mu[k] - mu[1] } # record OR
  for (k in 1:Ndrug) { or[k] <- exp(lor[k])} # record LOR
}




NMA.homo.eqcor <-function(){
  for (i in 1:Narm){
    y[i] ~ dbinom(mean[i],n[i])
    logit(mean[i]) <- mu[drug[i]] + v[study[i],drug[i]]
  }
  for (j in 1:Nstudy) { 
    v[j,1:Ndrug] ~ dmnorm(zero.AB[1:Ndrug], COV_mat[1:Ndrug, 1:Ndrug]) }
  
  sigma <- 1/sqrt(inv.sig.sq)
  inv.sig.sq ~ dgamma(0.001, 0.001)
  
  for(j in 1:Ndrug){
    for(k in 1:Ndrug){ 
      COV_mat[j,k] <- 1/sigma^2*ifelse(j == k, diag, offdiag)
    }
  }
  
  diag <- (1 + (Ndrug - 2)*rho)/(1 + (Ndrug - 2)*rho - (Ndrug - 1)*rho^2)
  offdiag <- (-rho/(1 + (Ndrug - 2)*rho - (Ndrug - 1)*rho^2))
  rho ~ dunif(-1/(Ndrug - 1), 0.9999)
  
  for (k in 1:Ndrug) { mu[k] ~ dnorm(0, 0.001) }  
  for (k in 1:Ndrug) { lor[k] <- mu[k] - mu[1] }
  for (k in 1:Ndrug) { or[k] <- exp(lor[k])}
}


##############################################################################





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
  
  for (k in 1:Ndrug) { or[k] <- exp(d[k])}
  
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





pcnetmeta_ABWish = function(){
  for(i in 1:len){
    p[i] <- phi(mu[t[i]] + vi[s[i], t[i]])
    r[i] ~ dbin(p[i], totaln[i])
    rhat[i] <- p[i]*totaln[i]
    dev[i] <- 2*(r[i]*(log(r[i]) - log(rhat[i])) +
                   (totaln[i] - r[i])*(log(totaln[i] - r[i]) - log(totaln[i] - rhat[i])))
  }
  totresdev <- sum(dev[])
  for(j in 1:nstudy){
    vi[j, 1:ntrt] ~ dmnorm(zeros[1:ntrt], T[1:ntrt, 1:ntrt])
  }
  for(j in 1:ntrt){
    AR[j] <- phi(mu[j]/sqrt(1 + invT[j,j]))
    mu[j] ~ dnorm(0, 0.001)
    sigma[j] <- sqrt(invT[j,j])
  }
  invT[1:ntrt, 1:ntrt] <- inverse(T[,])
  T[1:ntrt, 1:ntrt] ~ dwish(I[1:ntrt, 1:ntrt], ntrt + 1)
  for(j in 1:ntrt){        
    for(k in 1:ntrt){
      LRR[j,k] <- log(RR[j,k])
      LOR[j,k] <- log(OR[j,k])
      RR[j,k] <- AR[j]/AR[k]
      RD[j,k] <- AR[j] - AR[k]
      OR[j,k] <- AR[j]/(1 - AR[j])/AR[k]*(1 - AR[k])
    }
  }
}



pcnetmeta_ABhet = function(){
  for(i in 1:len){
    p[i] <- phi(mu[t[i]] + vi[s[i], t[i]])
    r[i] ~ dbin(p[i], totaln[i])
    rhat[i] <- p[i]*totaln[i]
    dev[i] <- 2*(r[i]*(log(r[i]) - log(rhat[i])) +
                   (totaln[i] - r[i])*(log(totaln[i] - r[i]) - log(totaln[i] - rhat[i])))
  }
  totresdev <- sum(dev[])
  for(j in 1:nstudy){
    vi[j, 1:ntrt] ~ dmnorm(zeros[1:ntrt], T[1:ntrt, 1:ntrt])
  }
  for(j in 1:ntrt){
    AR[j] <- phi(mu[j]/sqrt(1 + pow(sigma[j], 2)))
    mu[j] ~ dnorm(0, 0.001)
    sigma[j] <- 1/sqrt(inv.sig.sq[j])
    inv.sig.sq[j] ~ dgamma(a, b)
  }
  for(j in 1:ntrt){
    for(k in 1:ntrt){ 
      T[j,k] <- 1/sigma[j]*1/sigma[k]*ifelse(j == k, diag, offdiag)
    }
  }
  diag <- (1 + (ntrt - 2)*rho)/(1 + (ntrt - 2)*rho - (ntrt - 1)*rho^2)
  offdiag <- (-rho/(1 + (ntrt - 2)*rho - (ntrt - 1)*rho^2))
  rho ~ dunif(-1/(ntrt - 1), 0.9999)
  for(j in 1:ntrt){
    for(k in 1:ntrt){
      LRR[j,k] <- log(RR[j,k])
      LOR[j,k] <- log(OR[j,k])
      RR[j,k] <- AR[j]/AR[k]
      RD[j,k] <- AR[j] - AR[k]
      OR[j,k] <- AR[j]/(1 - AR[j])/AR[k]*(1 - AR[k])
    }
  }
}



