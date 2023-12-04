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


ABWish2<-function(){
  for(i in 1:Narm){
    p[i] <- phi(mu[t[i]] + vi[s[i], t[i]])
    r[i] ~ dbin(p[i], n[i])
    rhat[i] <- p[i]*n[i]
    dev[i] <- 2*(r[i]*(log(r[i]) - log(rhat[i])) +
                   (n[i] - r[i])*(log(n[i] - r[i]) - log(n[i] - rhat[i])))
  }
  totresdev <- sum(dev[])
  for(j in 1:Nstudy){
    vi[j, 1:Ndrug] ~ dmnorm(zeros[1:Ndrug], T[1:Ndrug, 1:Ndrug])
  }
  for(j in 1:Ndrug){
    AR[j] <- phi(mu[j]/sqrt(1 + invT[j,j]))
    mu[j] ~ dnorm(0, 0.001)
    sigma[j] <- sqrt(invT[j,j])
  }
  invT[1:Ndrug, 1:Ndrug] <- inverse(T[,])
  T[1:Ndrug, 1:Ndrug] ~ dwish(I[1:Ndrug, 1:Ndrug], Ndrug + 1)
  for(j in 1:Ndrug){        
    for(k in 1:Ndrug){
      LRR[j,k] <- log(RR[j,k])
      LOR[j,k] <- log(OR[j,k])
      RR[j,k] <- AR[j]/AR[k]
      RD[j,k] <- AR[j] - AR[k]
      OR[j,k] <- AR[j]/(1 - AR[j])/AR[k]*(1 - AR[k])
    }
  }
}




pcnetmeta_ABWish = model{
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


