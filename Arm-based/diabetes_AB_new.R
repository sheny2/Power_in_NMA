library(tidyverse)
library(R2jags)
library(gemtc)
library(tibble)
library(gemtc)
library(pcnetmeta)
data(diabetes)

source('Models.R')
source('Helpers.R')

source("ifplot.fun.R")
source("BayesDiagnos.fun.R")
source("ranko_sucra.R")

set.seed(2023)

######## Import data

diabetes_ab = data.frame(study = factor(diabetes$s.id), treatment = factor(diabetes$t.id), 
                         sampleSize = diabetes$n, responders = diabetes$r)

diabetes_ab$study = as.numeric(diabetes_ab$study)

network <- mtc.network(diabetes_ab)

cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                        hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                        re.prior.sd = 10)
cons.out <- mtc.run(cons.model, n.adapt=20000, n.iter=5000, thin=1)
gemtc::forest(cons.out)

estimates <- summary(cons.out)
estimates


trt = unique(diabetes_ab$treatment)

########
## simulation
library(foreach)
library(doParallel)

# Register the parallel backend
cl <- makeCluster(detectCores())

# Register the parallel backend
registerDoParallel(cl)

S = 1000

Log_OR_dat = c(log(1), estimates$summaries$statistics[,1][1:5])
OR_dat = exp(Log_OR_dat)


OR_12 = exp(estimates$summaries$statistics[,1][1])
OR_13 = exp(estimates$summaries$statistics[,1][2])
OR_14 = exp(estimates$summaries$statistics[,1][3])
OR_15 = exp(estimates$summaries$statistics[,1][4])
OR_16 = exp(estimates$summaries$statistics[,1][5])

pi_a = sum(diabetes_ab[diabetes_ab$treatment=="1",]$responders) / 
  sum(diabetes_ab[diabetes_ab$treatment=="1",]$sampleSize)

n_study = n_distinct(diabetes_ab$study)

tau_est = estimates$summaries$statistics[,1][6]

# pi_b =  calculate_pi(OR_ab, pi_a)
# pi_c =  calculate_pi(OR_ac, pi_a)
# pi_d =  calculate_pi(OR_ad, pi_a)


result_diabetes_all = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  ### Data pre for jags, AB model first
  Log_OR_dat = c(0, -0.23834311, 0.02570665, -0.13496070, -0.30540779, -0.49401317)
  tau_AB = c(0.568964683, 0.767673995, 0.630069285, 0.660131261, 0.685801466, 0.825769537)
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = diabetes_ab[diabetes_ab$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    if (study_data$treatment[1] == "1") {
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
    } else {for (k in 1:nrow(study_data)){
      d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
      tau = tau_AB[which(trt == study_data$treatment[k])]
      delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
      logit_p_ik = alpha_iB + delta_iBk
      p_ik = inverse_logit(logit_p_ik)
      y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
    }
    }
  }
  
  sim_dat = diabetes_ab
  sim_dat$responders = y_ik
  
  NS = 22
  NT = 6
  N = nrow(sim_dat)
  s = as.numeric(sim_dat$study)
  t = as.numeric(sim_dat$treatment)
  # t = as.integer(factor(sim_dat$treatment, levels = c("A", "B", "C", "D"), labels = c(1, 2, 3, 4)))
  y = sim_dat$responders
  n = sim_dat$sampleSize
  # drug_list<-unique(sim_dat$treatment)
  drug_list<- c(1,2,3,4,5,6)
  Narm <- as.numeric(table(sim_dat$study))
  n.obs <- matrix(NA,nrow=NS, ncol=max(Narm))
  n.eve <- matrix(NA,nrow=NS, ncol=max(Narm))
  dr <- matrix(NA,nrow=NS, ncol=max(Narm))
  
  # AB models
  data_AB <- list('Narm'=N, 'Nstudy'=NS, 
                  'Ndrug'=NT, 'study'= s, 'drug'=t, 
                  'y'=y, 'n'=n ,'Omega'=diag(rep(0.2,times=6)),
                  'zero.AB' = (rep(0, times=6)))
  inits_AB<- list(list(mu=rep(0,6)),
                  list(mu=rep(0,6)))
  para_AB<-c( "lor", "tau", "best1", "best2", "best3")
  fit_AB<-jags(data=data_AB, inits=inits_AB, para_AB,
               n.iter=20000, n.burnin = 5000, n.chains = 2, n.thin = 1,
               DIC=TRUE, model.file=ABWish)
  #output data 
  fit_AB$BUGSoutput$summary[,c(1, 3, 7)]
  #saving treatment effect output
  AB_trt_results<-data.frame(fit_AB$BUGSoutput$summary[,c(1, 3, 7)])
  AB_trt_results <- tibble::rownames_to_column(AB_trt_results, "drug_list")
  AB_trt_results<-AB_trt_results%>%
    filter(drug_list %in% c("lor[1]", "lor[2]", "lor[3]", "lor[4]", "lor[5]", "lor[6]"))
  
  
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
  
  AE_lower_95 = ABresults[5,3]
  AE_upper_95 = ABresults[5,4]
  AE_reject = 1-as.numeric(0 >= AE_lower_95 & 0 <= AE_upper_95)
  
  AF_lower_95 = ABresults[6,3]
  AF_upper_95 = ABresults[6,4]
  AF_reject = 1-as.numeric(0 >= AF_lower_95 & 0 <= AF_upper_95)
  
  reject_correct_AB = c(AB_reject, AC_reject, AD_reject, AE_reject, AF_reject)
  
  
  
  
  
  
  # LA model
  Log_OR_dat = c(0, -0.2934780,-0.0712749,-0.2427150, -0.4125047, -0.4914017)
  tau = 0.1345282
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = diabetes_ab[diabetes_ab$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    if (study_data$treatment[1] == "1") {
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
    } else {for (k in 1:nrow(study_data)){
      d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
      delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
      logit_p_ik = alpha_iB + delta_iBk
      p_ik = inverse_logit(logit_p_ik)
      y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
    }
    }
  }
  
  sim_dat = diabetes_ab
  sim_dat$responders = y_ik
  
  NS = 22
  NT = 6
  N = nrow(sim_dat)
  s = as.numeric(sim_dat$study)
  t = as.numeric(sim_dat$treatment)
  # t = as.integer(factor(sim_dat$treatment, levels = c("A", "B", "C", "D"), labels = c(1, 2, 3, 4)))
  y = sim_dat$responders
  n = sim_dat$sampleSize
  # drug_list<-unique(sim_dat$treatment)
  drug_list<- c(1,2,3,4,5,6)
  Narm <- as.numeric(table(sim_dat$study))
  n.obs <- matrix(NA,nrow=NS, ncol=max(Narm))
  n.eve <- matrix(NA,nrow=NS, ncol=max(Narm))
  dr <- matrix(NA,nrow=NS, ncol=max(Narm))

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
                 n.iter=20000, n.burnin = 5000, n.chains = 2, n.thin = 1,
                 DIC=TRUE, model.file=LARE)
  #output data 
  fit_LA$BUGSoutput$summary[,c(1, 3, 7)]
  
  #saving treatment effect output
  LA_trt_results<-data.frame(fit_LA$BUGSoutput$summary[,c(1, 3, 7)])
  LA_trt_results <- tibble::rownames_to_column(LA_trt_results, "drug_list")
  LA_trt_results<-LA_trt_results%>%
    filter(drug_list %in% c("d[1]", "d[2]", "d[3]", "d[4]", "d[5]", "d[6]"))
  
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
  
  AE_lower_95 = LAresults[5,3]
  AE_upper_95 = LAresults[5,4]
  AE_reject = 1-as.numeric(0 >= AE_lower_95 & 0 <= AE_upper_95)
  
  AF_lower_95 = LAresults[6,3]
  AF_upper_95 = LAresults[6,4]
  AF_reject = 1-as.numeric(0 >= AF_lower_95 & 0 <= AF_upper_95)
  
  reject_correct_LA = c(AB_reject, AC_reject, AD_reject, AE_reject, AF_reject)
  
  

  # return power results
  return(rbind(reject_correct_LA, reject_correct_AB))
}


result_diabetes_all = result_diabetes_all/S

result_diabetes_all = matrix(result_diabetes_all, nrow = 2)

colnames(result_diabetes_all) <- c("Power AB","Power AC","Power AD", "Power AE", "Power AF")
result_diabetes_all = result_diabetes_all %>% as.data.frame() 
result_diabetes_all$Model = c("LA", "AB")
result_diabetes_all

save(result_diabetes_all, file = "result_diabetes_all.RData")

stopCluster(cl)



