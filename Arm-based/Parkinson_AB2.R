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

# Name
# 1: Placebo
# 2: Pramipexole
# 3: Ropinirole
# 4: Bromocriptine
# 5: Cabergoline

load("parkinson_ab.RData")
parkinson_ab$treatment = as.numeric(parkinson_ab$treatment)

# normal/identity:forcontinuous(meandifference)data. 
# Requiredcolumns:[mean,std.err]or[mean,std.dev,sampleSize]. 
# Result:relativemeandifference.

parkinson_net <- mtc.network(parkinson_ab)
cons.model <- mtc.model(parkinson_net, type="consistency", 
                        likelihood="normal", link="identity", linearModel="random",
                        hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                        re.prior.sd = 10)
cons.out <- mtc.run(cons.model, n.adapt=20000, n.iter=50000, thin=1)
prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
prob <- round(prob, digits=3)
sucra <- gemtc::sucra(prob)
correct_rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))


estimates = summary(gemtc::relative.effect(cons.out,"1", c("2","3","4","5")))
estimates$measure
estimates$summaries$statistics

d12 = estimates$summaries$statistics[1,1]
d13 = estimates$summaries$statistics[2,1]
d14 = estimates$summaries$statistics[3,1]
d15 = estimates$summaries$statistics[4,1]
tau = estimates$summaries$statistics[5,1]


########
## simulation
library(foreach)
library(doParallel)

# Register the parallel backend
cl <- makeCluster(detectCores())

# Register the parallel backend
registerDoParallel(cl)

S = 1000


n_study = n_distinct(parkinson_ab$study)

mu_1 = mean(parkinson_ab[which(parkinson_ab$treatment == 1), ]$mean)
# mu_2 = mu_1 + d12 
# mu_3 = mu_1 + d13 
# mu_4 = mu_1 + d14 

dk = c(0, d12,d13,d14,d15)
trt = c(1,2,3,4,5)



result_pk_all = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  dk = c(0, -1.8669717, -0.5365794, -0.5338494, -0.8342570)
  tau = 0.3550
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = parkinson_ab[parkinson_ab$study==i,]
    
    alpha_iB = rnorm(n = 1, mean=mu_1, sd=0.1)
    
    if (study_data$treatment[1] == "1"){
      for (k in 1:nrow(study_data)){
        if (k == 1) {
          delta_ik = alpha_iB
          sd_ik = study_data$std.dev[k]
          n_ik = study_data$sampleSize[k]
          y_ik = c(y_ik, rnorm(n = 1, mean = delta_ik, sd = sqrt(sd_ik^2/n_ik)))
        } else {
          delta_ik = alpha_iB + rnorm(n = 1, mean = dk[which(trt == study_data$treatment[k])], sd = tau)
          sd_ik = study_data$std.dev[k]
          n_ik = study_data$sampleSize[k]
          y_ik = c(y_ik, rnorm(n = 1, mean = delta_ik, sd = sqrt(sd_ik^2/n_ik)))
        }
      }
    } else{
      for (k in 1:nrow(study_data)){
        delta_ik = alpha_iB + rnorm(n = 1, mean = dk[which(trt == study_data$treatment[k])], sd = tau)
        sd_ik = study_data$std.dev[k]
        n_ik = study_data$sampleSize[k]
        y_ik = c(y_ik, rnorm(n = 1, mean = delta_ik, sd = sqrt(sd_ik^2/n_ik)))
      }
    }
  }
  
  sim_dat = parkinson_ab
  sim_dat$mean = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- gemtc::mtc.model(network, type="consistency", 
                                 likelihood="normal", link="identity", linearModel="random",
                                 hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                                 re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=5000, n.iter=2000, thin=1)
  
  prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
  prob <- round(prob, digits=3)
  sucra <- gemtc::sucra(prob)
  rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
  # A_last = as.numeric(rank_order[4] == "A")
  Correct_order = as.numeric(identical(rank_order,correct_rank_order))
  
  res = summary(gemtc::relative.effect(cons.out,"1", c("2","3","4","5")))
  
  d12_lower_95 = res$summaries$quantiles[1,1]
  d12_upper_95 = res$summaries$quantiles[1,5]
  d12_reject = 1-as.numeric(0 >= d12_lower_95 & 0 <= d12_upper_95)
  
  d13_lower_95 = res$summaries$quantiles[2,1]
  d13_upper_95 = res$summaries$quantiles[2,5]
  d13_reject = 1-as.numeric(0 >= d13_lower_95 & 0 <= d13_upper_95)
  
  d14_lower_95 = res$summaries$quantiles[3,1]
  d14_upper_95 = res$summaries$quantiles[3,5]
  d14_reject = 1-as.numeric(0 >= d14_lower_95 & 0 <= d14_upper_95)
  
  d15_lower_95 = res$summaries$quantiles[4,1]
  d15_upper_95 = res$summaries$quantiles[4,5]
  d15_reject = 1-as.numeric(0 >= d15_lower_95 & 0 <= d15_upper_95)
  
  reject_correct_gemtc = c(d12_reject, d13_reject, d14_reject, d15_reject)
  
  
  ### Data pre for jags, AB models
  dk = c(0, -1.7634840, -0.1174405, -0.5640385, -1.5103555)
  tau_AB = c(0.5211938, 0.5622814, 0.5999670, 0.7519857, 0.5761361)
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = parkinson_ab[parkinson_ab$study==i,]
    
    alpha_iB = rnorm(n = 1, mean=mu_1, sd=0.1)
    
    if (study_data$treatment[1] == "1"){
      for (k in 1:nrow(study_data)){
        if (k == 1) {
          delta_ik = alpha_iB
          sd_ik = study_data$std.dev[k]
          n_ik = study_data$sampleSize[k]
          y_ik = c(y_ik, rnorm(n = 1, mean = delta_ik, sd = sqrt(sd_ik^2/n_ik)))
        } else {
          delta_ik = alpha_iB + rnorm(n = 1, mean = dk[which(trt == study_data$treatment[k])], 
                                      sd = tau_AB[which(trt == study_data$treatment[k])])
          sd_ik = study_data$std.dev[k]
          n_ik = study_data$sampleSize[k]
          y_ik = c(y_ik, rnorm(n = 1, mean = delta_ik, sd = sqrt(sd_ik^2/n_ik)))
        }
      }
    } else{
      for (k in 1:nrow(study_data)){
        delta_ik = alpha_iB + rnorm(n = 1, mean = dk[which(trt == study_data$treatment[k])], 
                                    sd = tau_AB[which(trt == study_data$treatment[k])])
        sd_ik = study_data$std.dev[k]
        n_ik = study_data$sampleSize[k]
        y_ik = c(y_ik, rnorm(n = 1, mean = delta_ik, sd = sqrt(sd_ik^2/n_ik)))
      }
    }
  }
  
  sim_dat = parkinson_ab
  sim_dat$mean = y_ik
  
  NS = 7
  NT = 5
  N = nrow(sim_dat)
  s = sim_dat$study
  t = sim_dat$treatment
  y = sim_dat$mean
  sigma = sim_dat$std.dev
  n = sim_dat$sampleSize
  drug_list<- c(1,2,3,4,5)
  Narm <- as.numeric(table(sim_dat$study))
  
  # AB models
  data_AB <- list('Narm'=N, 'Nstudy'=NS, 
                  'Ndrug'=NT, 'study'= s, 'drug'=t, "sigma" = sigma, 
                  'y'=y, 'n'=n ,'Omega'=diag(rep(0.2,times=5)),
                  'zero.AB' = (rep(0, times=5)))
  inits_AB<- list(list(mu=rep(0,5)),
                  list(mu=rep(0,5)))
  para_AB<-c( "lor", "tau")
  fit_AB<-jags(data=data_AB, inits=inits_AB, para_AB,
               n.iter=5000, n.burnin = 2000, n.chains = 2, n.thin = 1,
               DIC=TRUE, model.file=ABWish_C)
   
  # #saving treatment effect output
  AB_trt_results<-data.frame(fit_AB$BUGSoutput$summary[,c(1, 3, 7)])
  AB_trt_results <- tibble::rownames_to_column(AB_trt_results, "drug_list")
  AB_trt_results<-AB_trt_results%>%
    filter(drug_list %in% c("lor[1]", "lor[2]", "lor[3]", "lor[4]", "lor[5]"))
  

  ABresults <- AB_trt_results%>%
    mutate(LL = as.numeric(X2.5.),
           UL = as.numeric(X97.5.),
           mean = as.numeric(mean))

  d12_lower_95 = ABresults[2,3]
  d12_upper_95 = ABresults[2,4]
  d12_reject = 1-as.numeric(0 >= d12_lower_95 & 0 <= d12_upper_95)

  d13_lower_95 = ABresults[3,3]
  d13_upper_95 = ABresults[3,4]
  d13_reject = 1-as.numeric(0 >= d13_lower_95 & 0 <= d13_upper_95)

  d14_lower_95 = ABresults[4,3]
  d14_upper_95 = ABresults[4,4]
  d14_reject = 1-as.numeric(0 >= d14_lower_95 & 0 <= d14_upper_95)
  
  d15_lower_95 = ABresults[5,3]
  d15_upper_95 = ABresults[5,4]
  d15_reject = 1-as.numeric(0 >= d15_lower_95 & 0 <= d15_upper_95)

  reject_correct_AB = c(d12_reject, d13_reject, d14_reject, d15_reject)
  

  # LA model
  dk = c(0, -1.8342801, -0.4891893, -0.5174223, -0.8127853)
  tau = 0.3975971
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = parkinson_ab[parkinson_ab$study==i,]
    
    alpha_iB = rnorm(n = 1, mean=mu_1, sd=0.1)
    
    if (study_data$treatment[1] == "1"){
      for (k in 1:nrow(study_data)){
        if (k == 1) {
          delta_ik = alpha_iB
          sd_ik = study_data$std.dev[k]
          n_ik = study_data$sampleSize[k]
          y_ik = c(y_ik, rnorm(n = 1, mean = delta_ik, sd = sqrt(sd_ik^2/n_ik)))
        } else {
          delta_ik = alpha_iB + rnorm(n = 1, mean = dk[which(trt == study_data$treatment[k])], sd = tau)
          sd_ik = study_data$std.dev[k]
          n_ik = study_data$sampleSize[k]
          y_ik = c(y_ik, rnorm(n = 1, mean = delta_ik, sd = sqrt(sd_ik^2/n_ik)))
        }
      }
    } else{
      for (k in 1:nrow(study_data)){
        delta_ik = alpha_iB + rnorm(n = 1, mean = dk[which(trt == study_data$treatment[k])], sd = tau)
        sd_ik = study_data$std.dev[k]
        n_ik = study_data$sampleSize[k]
        y_ik = c(y_ik, rnorm(n = 1, mean = delta_ik, sd = sqrt(sd_ik^2/n_ik)))
      }
    }
  }
  
  sim_dat = parkinson_ab
  sim_dat$mean = y_ik
  
  NS = 7
  NT = 5
  N = nrow(sim_dat)
  s = sim_dat$study
  t = sim_dat$treatment
  y = sim_dat$mean
  sigma = sim_dat$std.dev
  n = sim_dat$sampleSize
  drug_list<- c(1,2,3,4,5)
  Narm <- as.numeric(table(sim_dat$study))
  n.obs <- matrix(NA,nrow=NS, ncol=max(Narm))
  n.eve <- matrix(NA,nrow=NS, ncol=max(Narm))
  n.sd <- matrix(NA,nrow=NS, ncol=max(Narm))
  dr <- matrix(NA,nrow=NS, ncol=max(Narm))
  
  study<-unique(sim_dat$study)
  for (i in 1:NS){
    n.obs[i,1:Narm[i]] <- sim_dat$sampleSize[sim_dat$study==study[i]]
    n.eve[i,1:Narm[i]] <- sim_dat$mean[sim_dat$study==study[i]]
    n.sd[i,1:Narm[i]] <- sim_dat$std.dev[sim_dat$study==study[i]]
    dr[i,1:Narm[i]] <- match(sim_dat$treatment[sim_dat$study==study[i]],drug_list)
  }
  
  
  ##putting data into list form
  data_LA <- list('Narm'=Narm, 'Nstudy'=NS,'Ndrug'=NT, 'drug'=dr,'y'=n.eve,'n'=n.obs, "sigma" = n.sd) 
  init_LA <- list(list(mu=rep(0,max(NS)), d=c(NA,rep(0,max(t)-1))),
                  list(mu=rep(0,max(NS)), d=c(NA,rep(0,max(t)-1))))
  para_LA <- c('d','tau')
  fit_LA <- jags(data=data_LA, inits=init_LA, para_LA,
                 n.iter=5000, n.burnin = 2000, n.chains = 2, n.thin = 1,
                 DIC=TRUE, model.file=LARE_C)
  
  #saving treatment effect output
  LA_trt_results<-data.frame(fit_LA$BUGSoutput$summary[,c(1, 3, 7)])
  LA_trt_results <- tibble::rownames_to_column(LA_trt_results, "drug_list")
  LA_trt_results<-LA_trt_results%>%
    filter(drug_list %in% c("d[1]", "d[2]", "d[3]", "d[4]", "d[5]"))
  
  LAresults <- LA_trt_results%>%
    mutate(LL = as.numeric(X2.5.),
           UL = as.numeric(X97.5.),
           mean = as.numeric(mean))
  
  d12_lower_95 = LAresults[2,3]
  d12_upper_95 = LAresults[2,4]
  d12_reject = 1-as.numeric(0 >= d12_lower_95 & 0 <= d12_upper_95)
  
  d13_lower_95 = LAresults[3,3]
  d13_upper_95 = LAresults[3,4]
  d13_reject = 1-as.numeric(0 >= d13_lower_95 & 0 <= d13_upper_95)
  
  d14_lower_95 = LAresults[4,3]
  d14_upper_95 = LAresults[4,4]
  d14_reject = 1-as.numeric(0 >= d14_lower_95 & 0 <= d14_upper_95)
  
  d15_lower_95 = LAresults[5,3]
  d15_upper_95 = LAresults[5,4]
  d15_reject = 1-as.numeric(0 >= d15_lower_95 & 0 <= d15_upper_95)
  
  reject_correct_LA = c(d12_reject, d13_reject, d14_reject, d15_reject)
  
  
  # CB model
  dk = c(0, -1.8688443, -0.5181954, -0.5616805, -0.8566140)
  tau_CB = c(0.5729980, 0.5003060, 0.5094048, 0.4793562, 0.5167020)
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = parkinson_ab[parkinson_ab$study==i,]
    
    alpha_iB = rnorm(n = 1, mean=mu_1, sd=0.1)
    
    if (study_data$treatment[1] == "1"){
      for (k in 1:nrow(study_data)){
        if (k == 1) {
          delta_ik = alpha_iB
          sd_ik = study_data$std.dev[k]
          n_ik = study_data$sampleSize[k]
          y_ik = c(y_ik, rnorm(n = 1, mean = delta_ik, sd = sqrt(sd_ik^2/n_ik)))
        } else {
          delta_ik = alpha_iB + rnorm(n = 1, mean = dk[which(trt == study_data$treatment[k])], 
                                      sd = tau_CB[which(trt == study_data$treatment[k])])
          sd_ik = study_data$std.dev[k]
          n_ik = study_data$sampleSize[k]
          y_ik = c(y_ik, rnorm(n = 1, mean = delta_ik, sd = sqrt(sd_ik^2/n_ik)))
        }
      }
    } else{
      for (k in 1:nrow(study_data)){
        delta_ik = alpha_iB + rnorm(n = 1, mean = dk[which(trt == study_data$treatment[k])], 
                                    sd = tau_CB[which(trt == study_data$treatment[k])])
        sd_ik = study_data$std.dev[k]
        n_ik = study_data$sampleSize[k]
        y_ik = c(y_ik, rnorm(n = 1, mean = delta_ik, sd = sqrt(sd_ik^2/n_ik)))
      }
    }
  }
  
  sim_dat = parkinson_ab
  sim_dat$mean = y_ik
  
  NS = 7
  NT = 5
  N = nrow(sim_dat)
  s = sim_dat$study
  t = sim_dat$treatment
  y = sim_dat$mean
  sigma = sim_dat$std.dev
  n = sim_dat$sampleSize
  drug_list<- c(1,2,3,4,5)
  Narm <- as.numeric(table(sim_dat$study))
  
  data_CB <- list('Narm'=N, 'Nstudy'=NS, 
                  'Ndrug'=NT, 'study'= s, 'drug'=t, "sigma" = sigma,
                  'y'=y, 'n'=n ,'Omega'=diag(rep(0.2,times=5)))
  inits_CB <- list(list(mu=rep(0,max(NS)), d=c(NA,rep(1,max(t)-1))),
                   list(mu=rep(0,max(NS)), d=c(NA,rep(0,max(t)-1))))
  
  para_CB<-c( "d", "tau")
  fit_CB<-jags(data=data_CB, inits=inits_CB, para_CB,
               n.iter=5000, n.burnin = 2000, n.chains = 2, n.thin = 1,
               DIC=TRUE, model.file=CBWish_C)
  
  #saving treatment effect output
  CB_trt_results<-data.frame(fit_CB$BUGSoutput$summary[,c(1, 3, 7)])
  CB_trt_results <- tibble::rownames_to_column(CB_trt_results, "drug_list")
  CB_trt_results<-CB_trt_results%>%
    filter(drug_list %in% c("d[1]", "d[2]", "d[3]", "d[4]", "d[5]"))
  
  CBresults <- CB_trt_results%>%
    mutate(LL = as.numeric(X2.5.),
           UL = as.numeric(X97.5.),
           mean = as.numeric(mean))
  
  d12_lower_95 = CBresults[2,3]
  d12_upper_95 = CBresults[2,4]
  d12_reject = 1-as.numeric(0 >= d12_lower_95 & 0 <= d12_upper_95)
  
  d13_lower_95 = CBresults[3,3]
  d13_upper_95 = CBresults[3,4]
  d13_reject = 1-as.numeric(0 >= d13_lower_95 & 0 <= d13_upper_95)
  
  d14_lower_95 = CBresults[4,3]
  d14_upper_95 = CBresults[4,4]
  d14_reject = 1-as.numeric(0 >= d14_lower_95 & 0 <= d14_upper_95)
  
  d15_lower_95 = CBresults[5,3]
  d15_upper_95 = CBresults[5,4]
  d15_reject = 1-as.numeric(0 >= d15_lower_95 & 0 <= d15_upper_95)
  
  reject_correct_CB =  c(d12_reject, d13_reject, d14_reject, d15_reject)
  
  # return power results
  return(rbind(reject_correct_gemtc, reject_correct_LA, reject_correct_CB, reject_correct_AB))
}


result_pk_all = result_pk_all/S

result_pk_all = matrix(result_pk_all, nrow = 4)

colnames(result_pk_all) <- c("Power D12","Power D13","Power D14","Power D15")
result_pk_all = result_pk_all %>% as.data.frame() 
result_pk_all$Model = c("GEMTC", "LA", "CB", "AB")
result_pk_all

save(result_pk_all, file = "result_pk_all2.RData")

stopCluster(cl)



