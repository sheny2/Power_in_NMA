source("Model2.R")
library(R2jags)
library(doMC)
library(doParallel)

S = 500

k_ab = 6
k_ac = 6
k_bc = 6

pi_a = 0.5

OR_ab = 1.2
OR_ac = 1.8

tau = 0.1


N_cores = detectCores() - 1

# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)



result_four_model = foreach::foreach (i = 1:S, .combine = "rbind") %dopar% {
  
  source("Model2.R")
  library(R2jags)
  library(tidyverse)

  # simulate k_ab studies (indirect)
  dat_ab = data.frame(matrix(ncol = 4, nrow = 0))
  if (k_ab != 0){
    for (j in 1:k_ab)
    {
      n_j = runif(n = 1, min = 100, max = 500)
      n_aj = n_bj = round(n_j / 2)
      
      log_OR_ab_j = rnorm(n = 1, mean = log(OR_ab),  sd = tau)
      
      pi_aj = runif(n = 1, min = 1/2 * pi_a, max = 3/2 * pi_a)
      pi_bj = pi_aj * exp(log_OR_ab_j) / (1 - pi_aj + pi_aj * exp(log_OR_ab_j) )
      
      e_aj = rbinom(n = 1, size = n_aj, prob = pi_aj)
      e_bj = rbinom(n = 1, size = n_bj, prob = pi_bj)
      
      dat_ab = rbind(dat_ab, rbind(cbind(j, "A", n_aj, e_aj), cbind(j, "B", n_bj, e_bj)))
    }
  }
  
  
  # simulate k_ac studies (indirect)
  dat_ac = data.frame(matrix(ncol = 4, nrow = 0))
  if (k_ac != 0){
    for (j in 1:k_ac)
    {
      n_j = runif(n = 1, min = 100, max = 500)
      n_aj = n_cj = round(n_j / 2)
      
      log_OR_ac_j = rnorm(n = 1, mean = log(OR_ac),  sd = tau)
      
      pi_aj = runif(n = 1, min = 1/2 * pi_a, max = 3/2 * pi_a)
      pi_cj = pi_aj * exp(log_OR_ac_j) / (1 - pi_aj + pi_aj * exp(log_OR_ac_j) )
      
      e_aj = rbinom(n = 1, size = n_aj, prob = pi_aj)
      e_cj = rbinom(n = 1, size = n_cj, prob = pi_cj)
      
      dat_ac = rbind(dat_ac, rbind(cbind(j+k_ab, "A", n_aj, e_aj), cbind(j+k_ab, "C", n_cj, e_cj)))
    }
  }
  
  
  # simulate k_bc studies (direct)
  OR_bc = exp(log(OR_ac) - log(OR_ab))
  pi_b = (OR_ab * pi_a / (1 - pi_a))/(1 + OR_ab * pi_a / (1-pi_a))
  
  dat_bc = data.frame(matrix(ncol = 4, nrow = 0))
  if (k_bc != 0){
    for (j in 1:k_bc)
    {
      n_j = runif(n = 1, min = 100, max = 500)
      n_bj = n_cj = round(n_j / 2)
      
      log_OR_bc_j = rnorm(n = 1, mean = log(OR_bc),  sd = tau)
      
      pi_bj = runif(n = 1, min = 1/2 * pi_b, max = 3/2 * pi_b)
      pi_cj = pi_bj * exp(log_OR_bc_j) / (1 - pi_bj + pi_bj * exp(log_OR_bc_j) )
      
      e_bj = rbinom(n = 1, size = n_bj, prob = pi_bj)
      e_cj = rbinom(n = 1, size = n_cj, prob = pi_cj)
      
      dat_bc = rbind(dat_bc, rbind(cbind(j+k_ab+k_ac, "B", n_bj, e_bj), cbind(j+k_ab+k_ac, "C", n_cj, e_cj)))
    }
  }
  
  
  colnames(dat_ab) = colnames(dat_ac) = colnames(dat_bc) = c("study", "treatment", "sampleSize", "responders")
  all_data = rbind(dat_ab, dat_ac, dat_bc)
  
  
  
  all_data[,c('sampleSize', 'responders')] <-
    lapply(all_data[,c('sampleSize', 'responders')], as.numeric)
  
  
  all_data_pc = tibble(s.id = all_data$study,
                       t.id = all_data$treatment, 
                       r = all_data$responders, 
                       n = all_data$sampleSize) 
  
  
  
  # data pre for jags
  NS = n_distinct(all_data$study)
  NT = n_distinct(all_data$treatment)
  N = nrow(all_data)
  s = all_data$study
  # t = all_data$treatment
  t = as.integer(factor(all_data$treatment, levels = c("A", "B", "C"), labels = c(1,2,3)))
  y = all_data$responders
  n = all_data$sampleSize
  drug_list <- c("A", "B", "C")
  Narm <- as.numeric(table(all_data$study))
  n.obs <- matrix(NA,nrow=NS, ncol=max(Narm))
  n.eve <- matrix(NA,nrow=NS, ncol=max(Narm))
  dr <- matrix(NA,nrow=NS, ncol=max(Narm))
  study<-unique(all_data$study)
  for (i in 1:NS){
    n.obs[i,1:Narm[i]] <- all_data$sampleSize[all_data$study==study[i]]
    n.eve[i,1:Narm[i]] <- all_data$responders[all_data$study==study[i]]
    dr[i,1:Narm[i]] <- match(all_data$treatment[all_data$study==study[i]],drug_list)
  }
  
  
  ### Running Arm Based  Model  ########
  data_AB <- list('Narm'=N, 'Nstudy'=NS,
                  'Ndrug'=NT, 'study'= s, 'drug'=t,
                  'y'=y, 'n'=n ,'Omega'=diag(rep(0.2,times=3)),
                  'zero.AB' = (rep(0, times=3)))
  inits_AB<- list(list(mu=rep(0,3)),
                  list(mu=rep(0,3)))
  para_AB<-c( "or", "lor", "tau")
  fit_AB_IW <-jags(data=data_AB, inits=inits_AB, para_AB,
               n.iter= 5000, n.burnin = 2000, n.chains = 2, n.thin = 1,
               DIC=TRUE, model.file=ABWish.het.cor)
  
  
  # output data
  AB_het_IW_DIC = fit_AB_IW$BUGSoutput$DIC
  AB_het_IW<-data.frame(fit_AB_IW$BUGSoutput$summary[,c(1, 3, 7)])
  AB_het_IW <- tibble::rownames_to_column(AB_het_IW, "drug_list")
  AB_het_IW_core <-AB_het_IW%>%
    filter(drug_list %in% c("or[1]", "or[2]", "or[3]", "tau[1]", "tau[2]", "tau[3]"))
  # AB_het_IW
  
  
  
  
  ### Running Arm Based  Model 2 ########
  data_AB <- list('Narm'=N, 'Nstudy'=NS,
                  'Ndrug'=NT, 'study'= s, 'drug'=t,
                  'y'=y, 'n'=n ,
                  'zero.AB' = (rep(0, times=3)))
  inits_AB<- list(list(mu=rep(0,3)),
                  list(mu=rep(0,3)))
  para_AB<-c( "lor", "or", "rho", "sigma")
  fit_AB_het.eqcor <-jags(data=data_AB, inits=inits_AB, para_AB,
               n.iter= 5000, n.burnin = 2000, n.chains = 2, n.thin = 1,
               DIC=TRUE, model.file=NMA.het.eqcor)
  
  AB_het_eqcor_DIC = fit_AB_het.eqcor$BUGSoutput$DIC
  AB_het_eqcor<-data.frame(fit_AB_het.eqcor$BUGSoutput$summary[,c(1, 3, 7)])
  AB_het_eqcor <- tibble::rownames_to_column(AB_het_eqcor, "drug_list")
  AB_het_eqcor_core<-AB_het_eqcor%>%
    filter(drug_list %in% c("or[1]", "or[2]", "or[3]", "rho", "sigma[1]", "sigma[2]", "sigma[3]"))
  
  
  
  
  ### Running Arm Based  Model 3 ########
  data_AB <- list('Narm'=N, 'Nstudy'=NS,
                  'Ndrug'=NT, 'study'= s, 'drug'=t,
                  'y'=y, 'n'=n ,
                  'zero.AB' = (rep(0, times=3)))
  inits_AB<- list(list(mu=rep(0,3)),
                  list(mu=rep(0,3)))
  para_AB<-c( "lor", "or", "rho", "sigma")
  fit_AB_homo_eqcor<-jags(data=data_AB, inits=inits_AB, para_AB,
               n.iter= 5000, n.burnin = 2000, n.chains = 2, n.thin = 1,
               DIC=TRUE, model.file=NMA.homo.eqcor)
  
  AB_homo_eqcor_DIC = fit_AB_homo_eqcor$BUGSoutput$DIC
  AB_homo_eqcor<-data.frame(fit_AB_homo_eqcor$BUGSoutput$summary[,c(1, 3, 7)])
  AB_homo_eqcor <- tibble::rownames_to_column(AB_homo_eqcor, "drug_list")
  AB_homo_eqcor_core<-AB_homo_eqcor%>%
    filter(drug_list %in% c("or[1]", "or[2]", "or[3]", "rho", "sigma"))
  
  
  
  
  ### Running LA Model ########
  
  
  ##putting data into list form
  data_LA <- list('Narm'=Narm, 'Nstudy'=NS,'Ndrug'=NT, 'drug'=dr,'y'=n.eve,'n'=n.obs) 
  init_LA <- list(list(mu=rep(0,max(NS)), d=c(NA,rep(0,max(t)-1))),
                  list(mu=rep(0,max(NS)), d=c(NA,rep(0,max(t)-1))))
  para_LA <- c('d', "or", 'tau')
  fit_LA <- jags(data=data_LA, inits=init_LA, para_LA,
                 n.iter= 5000, n.burnin = 2000, n.chains = 2, n.thin = 1,
                 DIC=TRUE, model.file=LARE)
  
  #saving treatment effect output
  LA_results_DIC = fit_LA$BUGSoutput$DIC
  LA_results<-data.frame(fit_LA$BUGSoutput$summary[,c(1, 3, 7)])
  LA_results <- tibble::rownames_to_column(LA_results, "drug_list")
  LA_results_core <-LA_results%>%
    filter(drug_list %in% c("d[1]", "d[2]", "d[3]", "tau"))
  
  
  # save
  
  rbind(
  cbind(Model = "AB_het_IW", AB_het_IW, DIC = AB_het_IW_DIC),
  cbind(Model = "AB_het_eqcor", AB_het_eqcor, DIC = AB_het_eqcor_DIC),
  cbind(Model = "AB_homo_eqcor", AB_homo_eqcor, DIC = AB_homo_eqcor_DIC),
  cbind(Model = "LA_results", LA_results, DIC = LA_results_DIC)
  )

}


stopCluster(cl)



save(result_four_model, file = "result_four_model_tau_0.1.RData")


