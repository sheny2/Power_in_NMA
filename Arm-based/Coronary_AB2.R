library(tidyverse)
library(R2jags)
library(gemtc)
library(tibble)
library(gemtc)

source('Models.R')
source('Helpers.R')

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


network <- mtc.network(data.ab=Bleed_data, treatments=trts)
cons.model <- gemtc::mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                               hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                               re.prior.sd = 10)
cons.out <- mtc.run(cons.model, n.adapt=20000, n.iter=10000, thin=1)
summary(cons.out)
gemtc::forest(cons.out)


########
## simulation
library(foreach)
library(doParallel)
# Initialize parallel backend
cl <- makeCluster(detectCores())

# Register the parallel backend
registerDoParallel(cl)

S = 1000


### Major Bleeding
# pi_a = sum(bleed[bleed$treatment==1,]$TIMImajor) / sum(bleed[bleed$treatment==1,]$n)

pi_a = sum(Bleed_data[Bleed_data$treatment=="A",]$responders) / sum(Bleed_data[Bleed_data$treatment=="A",]$sampleSize)
n_study = n_distinct(Bleed_data$study)

tau = 0.24

# True OR (from paper)
OR_ab = 0.58
OR_ac = 0.70
OR_ad = 0.49
Log_OR_dat = log(c(1, OR_ab,OR_ac,OR_ad))

# pi_b =  calculate_pi(OR_ab, pi_a)
# pi_c =  calculate_pi(OR_ac, pi_a)
# pi_d =  calculate_pi(OR_ad, pi_a)


result_bleeding_all = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  Log_OR_dat = log(c(1, 0.5654593,0.6905949,0.5051204))
  tau = 0.3773
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = Bleed_data[Bleed_data$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    for (k in 1:nrow(study_data)){
      if (k == 1) {
        logit_p_ik = alpha_iB 
        p_ik = inverse_logit(logit_p_ik)
        y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
      } else {
        d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
        delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
        logit_p_ik = alpha_iB + delta_iBk
        p_ik = inverse_logit(logit_p_ik)
        y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
      }
    }
  }
  
  sim_dat = Bleed_data
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat, treatments=trts)
  cons.model <- gemtc::mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                                 hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                                 re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.iter=5000, n.adapt = 2000, thin=1)
  
  res = summary(gemtc::relative.effect(cons.out,"A", c("B","C","D")))
  
  AB_lower_95 = res$summaries$quantiles[1,1]
  AB_upper_95 = res$summaries$quantiles[1,5]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)
  
  AC_lower_95 = res$summaries$quantiles[2,1]
  AC_upper_95 = res$summaries$quantiles[2,5]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)
  
  AD_lower_95 = res$summaries$quantiles[3,1]
  AD_upper_95 = res$summaries$quantiles[3,5]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)
  
  res = summary(gemtc::relative.effect(cons.out,c("B","B","C"), c("C","D","D")))
  
  BC_lower_95 = res$summaries$quantiles[1,1]
  BC_upper_95 = res$summaries$quantiles[1,5]
  BC_reject = 1-as.numeric(0 >= BC_lower_95 & 0 <= BC_upper_95)
  
  BD_lower_95 = res$summaries$quantiles[2,1]
  BD_upper_95 = res$summaries$quantiles[2,5]
  BD_reject = 1-as.numeric(0 >= BD_lower_95 & 0 <= BD_upper_95)
  
  CD_lower_95 = res$summaries$quantiles[3,1]
  CD_upper_95 = res$summaries$quantiles[3,5]
  CD_reject = 1-as.numeric(0 >= CD_lower_95 & 0 <= CD_upper_95)
  
  reject_correct_gemtc = c(AB_reject, AC_reject, AD_reject)
  
  
  
  
  ### Data pre for jags, AB model first
  
  Log_OR_dat = log(c(1, 0.5666299,0.5335303,0.4539045))
  tau_AB= c(0.415619241, 0.495172675, 0.468347076, 0.426510282)
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = Bleed_data[Bleed_data$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    for (k in 1:nrow(study_data)){
      if (k == 1) {
        logit_p_ik = alpha_iB 
        p_ik = inverse_logit(logit_p_ik)
        y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
      } else {
        d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
        tau = tau_AB[which(trt == study_data$treatment[k])]
        delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
        logit_p_ik = alpha_iB + delta_iBk
        p_ik = inverse_logit(logit_p_ik)
        y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
      }
    }
  }
  
  sim_dat = Bleed_data
  sim_dat$responders = y_ik
  
  
  NS = 4
  NT = 4
  N = nrow(sim_dat)
  s = sim_dat$study
  t = as.integer(factor(sim_dat$treatment, levels = c("A", "B", "C", "D"), labels = c(1, 2, 3, 4)))
  y = sim_dat$responders
  n = sim_dat$sampleSize
  drug_list <- unique(sim_dat$treatment)
  Narm <- as.numeric(table(sim_dat$study))
  n.obs <- matrix(NA, nrow = NS, ncol = max(Narm))
  n.eve <- matrix(NA, nrow = NS, ncol = max(Narm))
  dr <- matrix(NA, nrow = NS, ncol = max(Narm))
  study<-unique(sim_dat$study)
  
  
  # AB models
  data_AB <- list('Narm'=N, 'Nstudy'=NS, 
                  'Ndrug'=NT, 'study'= s, 'drug'=t, 
                  'y'=y, 'n'=n ,'Omega'=diag(rep(0.2,times=4)),
                  'zero.AB' = (rep(0, times=4)))
  inits_AB<- list(list(mu=rep(0,4)),
                  list(mu=rep(0,4)))
  para_AB<-c( "lor", "tau", "best1", "best2", "best3")
  fit_AB<- jags(data=data_AB, inits=inits_AB, para_AB,
                n.iter=10000, n.burnin = 2000, n.chains = 2, n.thin = 1,
               DIC=TRUE, model.file=ABWish)
   
  # #saving treatment effect output
  AB_trt_results <- data.frame(fit_AB$BUGSoutput$summary[,c(1, 3, 7)])
  AB_trt_results <- tibble::rownames_to_column(AB_trt_results, "drug_list")
  AB_trt_results <- AB_trt_results %>%
     filter(drug_list %in% c("lor[1]", "lor[2]", "lor[3]", "lor[4]"))

  ABresults <- AB_trt_results%>%
    mutate(LL = as.numeric(X2.5.),
           UL = as.numeric(X97.5.),
           mean = as.numeric(mean))

  AB_lower_95 = ABresults[2,3]
  AB_upper_95 = ABresults[2,4]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)

  AC_lower_95 = ABresults[3,3]
  AC_upper_95 = ABresults[3,4]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)

  AD_lower_95 = ABresults[4,3]
  AD_upper_95 = ABresults[4,4]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)

  reject_correct_AB = c(AB_reject, AC_reject, AD_reject)
  
  
  

  # LA model
  Log_OR_dat = log(c(1, 0.5614587,0.6850401,0.5023263))
  tau = 0.307250329
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = Bleed_data[Bleed_data$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    for (k in 1:nrow(study_data)){
      if (k == 1) {
        logit_p_ik = alpha_iB 
        p_ik = inverse_logit(logit_p_ik)
        y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
      } else {
        d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
        delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
        logit_p_ik = alpha_iB + delta_iBk
        p_ik = inverse_logit(logit_p_ik)
        y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
      }
    }
  }
  
  sim_dat = Bleed_data
  sim_dat$responders = y_ik
  
  NS = 4
  NT = 4
  N = nrow(sim_dat)
  s = sim_dat$study
  t = as.integer(factor(sim_dat$treatment, levels = c("A", "B", "C", "D"), labels = c(1, 2, 3, 4)))
  y = sim_dat$responders
  n = sim_dat$sampleSize
  drug_list <- unique(sim_dat$treatment)
  Narm <- as.numeric(table(sim_dat$study))
  n.obs <- matrix(NA, nrow = NS, ncol = max(Narm))
  n.eve <- matrix(NA, nrow = NS, ncol = max(Narm))
  dr <- matrix(NA, nrow = NS, ncol = max(Narm))
  study<-unique(sim_dat$study)
  
  study<-unique(sim_dat$study)
  for (i in 1:NS){
    n.obs[i,1:Narm[i]] <- sim_dat$sampleSize[sim_dat$study==study[i]]
    n.eve[i,1:Narm[i]] <- sim_dat$responders[sim_dat$study==study[i]]
    dr[i,1:Narm[i]] <- match(sim_dat$treatment[sim_dat$study==study[i]],drug_list)
  }
  
  ##putting data into list form
  data_LA <- list('Narm'=Narm, 'Nstudy'=NS,'Ndrug'=NT, 'drug'=dr,'y'=n.eve,'n'=n.obs) 
  init_LA <- list(list(mu=rep(0,max(NS)), d=c(NA,rep(0,max(t)-1))),
                  list(mu=rep(0,max(NS)), d=c(NA,rep(0,max(t)-1))))
  para_LA <- c('d','tau','best1', 'best2', 'best3')
  fit_LA <- jags(data=data_LA, inits=init_LA, para_LA,
                 n.iter=10000, n.burnin = 2000, n.chains = 2, n.thin = 1,
                 DIC=TRUE, model.file=LARE)
  
  #saving treatment effect output
  LA_trt_results<-data.frame(fit_LA$BUGSoutput$summary[,c(1, 3, 7)])
  LA_trt_results <- tibble::rownames_to_column(LA_trt_results, "drug_list")
  LA_trt_results<-LA_trt_results%>%
    filter(drug_list %in% c("d[1]", "d[2]", "d[3]", "d[4]"))
  
  LAresults <- LA_trt_results%>%
    mutate(LL = as.numeric(X2.5.),
           UL = as.numeric(X97.5.),
           mean = as.numeric(mean))
  
  AB_lower_95 = LAresults[2,3]
  AB_upper_95 = LAresults[2,4]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)
  
  AC_lower_95 = LAresults[3,3]
  AC_upper_95 = LAresults[3,4]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)
  
  AD_lower_95 = LAresults[4,3]
  AD_upper_95 = LAresults[4,4]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)
  
  reject_correct_LA = c(AB_reject, AC_reject, AD_reject) 
  

  
  
  # CB model
  Log_OR_dat = log(c(1, 0.5721020,0.7047791,0.5032106))
  tau_CB = c(0.580856953,0.485949485,0.508666132,0.438989182)
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = Bleed_data[Bleed_data$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    for (k in 1:nrow(study_data)){
      if (k == 1) {
        logit_p_ik = alpha_iB 
        p_ik = inverse_logit(logit_p_ik)
        y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
      } else {
        d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
        tau = tau_CB[which(trt == study_data$treatment[k])]
        delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
        logit_p_ik = alpha_iB + delta_iBk
        p_ik = inverse_logit(logit_p_ik)
        y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
      }
    }
  }
  
  sim_dat = Bleed_data
  sim_dat$responders = y_ik
  
  NS = 4
  NT = 4
  N = nrow(sim_dat)
  s = sim_dat$study
  t = as.integer(factor(sim_dat$treatment, levels = c("A", "B", "C", "D"), labels = c(1, 2, 3, 4)))
  y = sim_dat$responders
  n = sim_dat$sampleSize
  drug_list <- unique(sim_dat$treatment)
  Narm <- as.numeric(table(sim_dat$study))
  n.obs <- matrix(NA, nrow = NS, ncol = max(Narm))
  n.eve <- matrix(NA, nrow = NS, ncol = max(Narm))
  dr <- matrix(NA, nrow = NS, ncol = max(Narm))
  study<-unique(sim_dat$study)
  
  data_CB <- list('Narm'=N, 'Nstudy'=NS, 
                  'Ndrug'=NT, 'study'= s, 'drug'=t, 
                  'y'=y, 'n'=n ,'Omega'=diag(rep(0.2,times=4)))
  inits_CB <- list(list(mu=rep(0,max(NS)), d=c(NA,rep(1,max(t)-1))),
                   list(mu=rep(0,max(NS)), d=c(NA,rep(0,max(t)-1))))
  
  para_CB<-c( "d", "tau", "best1", 'best2', 'best3')
  fit_CB<-jags(data=data_CB, inits=inits_CB, para_CB,
               n.iter=10000, n.burnin = 2000, n.chains = 2, n.thin = 1,
               DIC=TRUE, model.file=CBWish)
  
  #saving treatment effect output
  CB_trt_results<-data.frame(fit_CB$BUGSoutput$summary[,c(1, 3, 7)])
  CB_trt_results <- tibble::rownames_to_column(CB_trt_results, "drug_list")
  CB_trt_results<-CB_trt_results%>%
    filter(drug_list %in% c("d[1]", "d[2]", "d[3]", "d[4]"))
  
  CBresults <- CB_trt_results%>%
    mutate(LL = as.numeric(X2.5.),
           UL = as.numeric(X97.5.),
           mean = as.numeric(mean))
  
  AB_lower_95 = CBresults[2,3]
  AB_upper_95 = CBresults[2,4]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)
  
  AC_lower_95 = CBresults[3,3]
  AC_upper_95 = CBresults[3,4]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)
  
  AD_lower_95 = CBresults[4,3]
  AD_upper_95 = CBresults[4,4]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)
  
  reject_correct_CB = c(AB_reject, AC_reject, AD_reject) 
  
  return(rbind(reject_correct_gemtc, reject_correct_LA, reject_correct_CB, reject_correct_AB))
}


result_bleeding_all = result_bleeding_all/S

result_bleeding_all = matrix(result_bleeding_all, nrow = 4)

colnames(result_bleeding_all) <- c("Power AB","Power AC","Power AD")
result_bleeding_all = result_bleeding_all %>% as.data.frame() 
result_bleeding_all$Model = c("GEMTC", "LA", "CB", "AB")
result_bleeding_all

save(result_bleeding_all, file = "result_bleeding_all2.RData")

stopCluster(cl)

