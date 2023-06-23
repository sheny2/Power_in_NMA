library(tidyverse)
library(gemtc)

source("ifplot.fun.R")
source("BayesDiagnos.fun.R")
source("ranko_sucra.R")

set.seed(2023)

######## Import data
bleed <- read.table("Data_Bleeding.csv", sep=",", header=TRUE)
mace <- read.table("Data_MACE.csv", sep=",", header=TRUE)


######## Row 7 and 8 need to be combined
bleed[7,3:7] = bleed[7,3:7] + bleed[8,3:7]
bleed <- bleed[-8,]
mace[7,3:10] = mace[7,3:10] + mace[8,3:10]
mace <- mace[-8,]


Bleed_data <- bleed[,c("study", "treatment", "n", "TIMImajor")]
Bleed_data = Bleed_data[1:(nrow(Bleed_data)-2), ]

MACE_data <- mace[,c("study", "treatment", "n", "MACE")]        
MACE_data = MACE_data[1:(nrow(MACE_data)-2), ]


## Note: A is reference

trt <- c("A","B", "C","D")
# trt <- c("VKA + DAPT", "VKA + P2Y12", "NOAC + DAPT", "NOAC + P2Y12")

trts <- read.table(textConnection('id description
                                  A "VKA + DAPT"
                                  B "VKA + P2Y12"
                                  C "NOAC + DAPT"
                                  D "NOAC + P2Y12"'), header=TRUE)



Bleed_data$treatment <- trt[Bleed_data$treatment]
colnames(Bleed_data) = c("study","treatment", "sampleSize","responders")
MACE_data$treatment <- trt[MACE_data$treatment]
colnames(MACE_data) = c("study","treatment", "sampleSize","responders")


Bleed_data


######## Fit AB model

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

########
N = nrow(Bleed_data)
NS = n_distinct(Bleed_data$study)
NT = n_distinct(Bleed_data$treatment)
s = Bleed_data$study
t = Bleed_data$treatment
t = c(1,2,1,3,4,1,4,1,2,3,4)
y = Bleed_data$responders
n = Bleed_data$sampleSize

data_AB <- list('Narm'=N, 'Nstudy'=NS, 
                'Ndrug'=NT, 'study'= s, 'drug'=t, 
                'y'=y, 'n'=n ,'Omega'=diag(rep(0.2,times=4)),
                'zero.AB' = (rep(0, times=4)))
inits_AB<- list(list(mu=rep(0,4)),
                list(mu=rep(0,4)))
para_AB<-c( "lor", "tau", "best1", "best2", "best3")
fit_AB<-jags(data=data_AB, inits=inits_AB, para_AB,
             n.iter=20000, n.burnin = 1500, n.chains = 2, n.thin = 1,
             DIC=TRUE, model.file=ABWish)
#output data 
fit_AB$BUGSoutput$summary[,c(1, 3, 7)]
#saving treatment effect output
AB_trt_results<-data.frame(fit_AB$BUGSoutput$summary[,c(1, 3, 7)])
AB_trt_results <- tibble::rownames_to_column(AB_trt_results, "drug_list")
AB_trt_results<-AB_trt_results%>%
  filter(drug_list %in% c("lor[1]", "lor[2]", "lor[3]", "lor[4]"))


ABresults<-AB_trt_results%>%
  mutate(LL = as.numeric(X2.5.), 
         UL = as.numeric(X97.5.), 
         mean = as.numeric(mean))%>%
  filter(!(drug_list==1))
ggplot(ABresults, aes(y = drug_list, x =mean )) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25)+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  labs(y = "", x="", title = "Arm-Based NMA Treatments LOR")

########



# # simulation
# library(foreach)
# library(doParallel)
# # Initialize parallel backend
# cl <- makeCluster(12)
# 
# # Register the parallel backend
# registerDoParallel(cl)
# 
# S = 50

