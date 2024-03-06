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
load("smokingcessation_ab.RData")

smokingcessation_ab$study = as.numeric(smokingcessation_ab$study)

network <- mtc.network(smokingcessation_ab)
cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                        hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                        re.prior.sd = 10)
cons.out <- mtc.run(cons.model, n.adapt=20000, n.iter=5000, thin=1)
gemtc::forest(cons.out)

estimates <- summary(cons.out)
estimates


########
## simulation
library(foreach)
library(doParallel)

# Register the parallel backend
cl <- makeCluster(detectCores())

# Register the parallel backend
registerDoParallel(cl)

S = 1000


Log_OR_dat = c(log(1), estimates$summaries$statistics[,1][1:3])
OR_dat = exp(Log_OR_dat)


OR_ab = exp(estimates$summaries$statistics[,1][1])
OR_ac = exp(estimates$summaries$statistics[,1][2])
OR_ad = exp(estimates$summaries$statistics[,1][3])

pi_a = sum(smokingcessation_ab[smokingcessation_ab$treatment=="A",]$responders) / 
  sum(smokingcessation_ab[smokingcessation_ab$treatment=="A",]$sampleSize)

n_study = n_distinct(smokingcessation_ab$study)

tau = estimates$summaries$statistics[,1][4]

# pi_b =  calculate_pi(OR_ab, pi_a)
# pi_c =  calculate_pi(OR_ac, pi_a)
# pi_d =  calculate_pi(OR_ad, pi_a)



result_smoke_all = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  Log_OR_dat = c(log(1), estimates$summaries$statistics[,1][1:3])
  tau = 0.8432
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = smokingcessation_ab[smokingcessation_ab$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    if (study_data$treatment[1] == "A") {
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
  
  sim_dat = smokingcessation_ab
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- gemtc::mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                                 hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                                 re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.iter=5000, n.adapt = 2000, thin=1)
  
  prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
  prob <- round(prob, digits=3)
  sucra <- gemtc::sucra(prob)
  rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
  A_last = as.numeric(rank_order[4] == "A")
  Correct_order = as.numeric(identical(rank_order,c("D","C","B", "A")))
  
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
  
  reject_correct_gemtc = c(AB_reject, AC_reject, AD_reject)
  
  
  ### Data pre for jags, AB model first
  Log_OR_dat = log(c(1, 1.796733, 2.132154, 2.687028))
  tau_AB = c(0.4811976, 0.8334492, 0.7683358, 0.9263793)
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = smokingcessation_ab[smokingcessation_ab$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    if (study_data$treatment[1] == "A") {
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
  
  sim_dat = smokingcessation_ab
  sim_dat$responders = y_ik
  
  NS = 24
  NT = 4
  N = nrow(sim_dat)
  s = sim_dat$study
  # t = sim_dat$treatment
  t = as.integer(factor(sim_dat$treatment, levels = c("A", "B", "C", "D"), labels = c(1, 2, 3, 4)))
  y = sim_dat$responders
  n = sim_dat$sampleSize
  drug_list<-unique(sim_dat$treatment)
  drug_list<- c("A","B","C","D")
  Narm <- as.numeric(table(sim_dat$study))
  n.obs <- matrix(NA,nrow=NS, ncol=max(Narm))
  n.eve <- matrix(NA,nrow=NS, ncol=max(Narm))
  dr <- matrix(NA,nrow=NS, ncol=max(Narm))
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
               n.iter=5000, n.burnin = 2000, n.chains = 2, n.thin = 1,
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
  Log_OR_dat = log(c(1, 1.635225, 2.321211, 2.987180))
  tau = 0.8438631
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = smokingcessation_ab[smokingcessation_ab$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    if (study_data$treatment[1] == "A") {
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
  sim_dat = smokingcessation_ab
  sim_dat$responders = y_ik
  
  NS = 24
  NT = 4
  N = nrow(sim_dat)
  s = sim_dat$study
  # t = sim_dat$treatment
  t = as.integer(factor(sim_dat$treatment, levels = c("A", "B", "C", "D"), labels = c(1, 2, 3, 4)))
  y = sim_dat$responders
  n = sim_dat$sampleSize
  drug_list<-unique(sim_dat$treatment)
  drug_list<- c("A","B","C","D")
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
                 n.iter=5000, n.burnin = 2000, n.chains = 2, n.thin = 1,
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
  Log_OR_dat = log(c(1, 1.423690, 2.163300, 2.573106))
  tau_CB = c(0.4811976, 0.8334492, 0.7683358, 0.9263793)
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = smokingcessation_ab[smokingcessation_ab$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    if (study_data$treatment[1] == "A") {
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
    } else {for (k in 1:nrow(study_data)){
      d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
      tau = tau_CB[which(trt == study_data$treatment[k])]
      delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
      logit_p_ik = alpha_iB + delta_iBk
      p_ik = inverse_logit(logit_p_ik)
      y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
    }
    }
  }
  
  sim_dat = smokingcessation_ab
  sim_dat$responders = y_ik
  
  NS = 24
  NT = 4
  N = nrow(sim_dat)
  s = sim_dat$study
  # t = sim_dat$treatment
  t = as.integer(factor(sim_dat$treatment, levels = c("A", "B", "C", "D"), labels = c(1, 2, 3, 4)))
  y = sim_dat$responders
  n = sim_dat$sampleSize
  drug_list<-unique(sim_dat$treatment)
  drug_list<- c("A","B","C","D")
  Narm <- as.numeric(table(sim_dat$study))
  n.obs <- matrix(NA,nrow=NS, ncol=max(Narm))
  n.eve <- matrix(NA,nrow=NS, ncol=max(Narm))
  dr <- matrix(NA,nrow=NS, ncol=max(Narm))
  study<-unique(sim_dat$study) 
  
  data_CB <- list('Narm'=N, 'Nstudy'=NS, 
                  'Ndrug'=NT, 'study'= s, 'drug'=t, 
                  'y'=y, 'n'=n ,'Omega'=diag(rep(0.2,times=4)))
  inits_CB <- list(list(mu=rep(0,max(NS)), d=c(NA,rep(1,max(t)-1))),
                   list(mu=rep(0,max(NS)), d=c(NA,rep(0,max(t)-1))))
  
  para_CB<-c( "d", "tau", "best1", 'best2', 'best3')
  fit_CB<-jags(data=data_CB, inits=inits_CB, para_CB,
               n.iter=5000, n.burnin = 2000, n.chains = 2, n.thin = 1,
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
  
  # return power results
  return(rbind(reject_correct_gemtc, reject_correct_LA, reject_correct_CB, reject_correct_AB))
}


result_smoke_all = result_smoke_all/S

result_smoke_all = matrix(result_smoke_all, nrow = 4)

colnames(result_smoke_all) <- c("Power AB","Power AC","Power AD")
result_smoke_all = result_smoke_all %>% as.data.frame() 
result_smoke_all$Model = c("GEMTC", "LA", "CB", "AB")
result_smoke_all

save(result_smoke_all, file = "result_smoke_all2.RData")

stopCluster(cl)



