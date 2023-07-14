library(doMC)
library(gemtc)
source("Models.R")

k_ab = k_ac = 6
k_bc = 6

pi_a = 0.5

OR_ab = 1.2
OR_ac = 1.8

tau = 0.2
S = 10

power.sim.AB_full <- function(S = 500, k_ab = 0, k_ac = 0, k_bc = 0, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.4, tau = 0.01, verbose = F){
  
  result = foreach::foreach (i = 1:S, .combine = "+") %dopar% {
  
    library(R2jags)
    library(tidyverse)
    source("Models.R")
    
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
    

    # data pre for jags
    NS = n_distinct(all_data$study)
    NT = n_distinct(all_data$treatment)
    N = nrow(all_data)
    s = all_data$study
    # t = all_data$treatment
    t = as.integer(factor(all_data$treatment, levels = c("B","A","C"), labels = c(1,2,3)))
    y = all_data$responders
    n = all_data$sampleSize
    drug_list <- c("B","A","C")
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
    para_AB<-c( "lor", "tau", "best1", "best2", "best3")
    fit_AB<-jags(data=data_AB, inits=inits_AB, para_AB,
                 n.iter = 5000, n.burnin = 2000, n.chains = 2, n.thin = 1,
                 DIC=TRUE, model.file=ABWish)

    # output data
    AB_trt_results<-data.frame(fit_AB$BUGSoutput$summary[,c(1, 3, 7)])
    AB_trt_results <- tibble::rownames_to_column(AB_trt_results, "drug_list")
    AB_trt_results<-AB_trt_results%>%
      filter(drug_list %in% c("lor[1]", "lor[2]", "lor[3]"))

    AB_rank_prob<-data.frame(fit_AB$BUGSoutput$summary[,c(1, 3, 7)])
    AB_rank_prob <- tibble::rownames_to_column(AB_rank_prob, "best")

    AB_rank_prob <- data.frame(cbind(AB_rank_prob[AB_rank_prob$best%in%c("best1[1]", "best1[2]", "best1[3]"),2],
                                   AB_rank_prob[AB_rank_prob$best%in%c("best2[1]", "best2[2]", "best2[3]"),2],
                                   AB_rank_prob[AB_rank_prob$best%in%c("best3[1]", "best3[2]", "best3[3]"),2]))
    rownames(AB_rank_prob)<-c("trt_1", "trt_2", "trt_3")


    lower_95 = AB_trt_results$X2.5.[3]
    upper_95 = AB_trt_results$X97.5.[3]


    # 1B 2A 3C
    sucra <- gemtc::sucra(AB_rank_prob)
    rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
    true_rank = c("trt_2", "trt_1", "trt_3")


    reject_null = as.numeric(0 > upper_95 || 0 < lower_95)
    rank_correct = as.numeric(identical(rank_order, true_rank))
    c(reject_null, rank_correct)
    
    # # LA model
    # 
    # ##putting data into list form
    # data_LA <- list('Narm'=Narm, 'Nstudy'=NS,'Ndrug'=NT, 'drug'=dr,'y'=n.eve,'n'=n.obs) 
    # init_LA <- list(list(mu=rep(0,max(NS)), d=c(NA,rep(0,max(t)-1))),
    #                 list(mu=rep(0,max(NS)), d=c(NA,rep(0,max(t)-1))))
    # para_LA <- c('d','tau','best1', 'best2', 'best3')
    # fit_LA <- jags(data=data_LA, inits=init_LA, para_LA,
    #                n.iter=5000, n.burnin = 2000, n.chains = 2, n.thin = 1,
    #                DIC=TRUE, model.file=LARE)
    # #output data 
    # # fit_LA$BUGSoutput$summary[,c(1, 3, 7)]
    # 
    # #saving treatment effect output
    # LA_trt_results<-data.frame(fit_LA$BUGSoutput$summary[,c(1, 3, 7)])
    # LA_trt_results <- tibble::rownames_to_column(LA_trt_results, "drug_list")
    # LA_trt_results<-LA_trt_results%>%
    #   filter(drug_list %in% c("d[1]", "d[2]", "d[3]"))
    # 
    # #saving rank probability of treatment to be the best treatment (number 1)
    # LA_rank_prob<-data.frame(fit_LA$BUGSoutput$summary[,c(1, 3, 7)])
    # LA_rank_prob <- tibble::rownames_to_column(LA_rank_prob, "best")
    # 
    # LA_rank_prob<-data.frame(cbind(LA_rank_prob[LA_rank_prob$best%in%c("best1[1]", "best1[2]", "best1[3]"),2],
    #                                LA_rank_prob[LA_rank_prob$best%in%c("best2[1]", "best2[2]", "best2[3]"),2],
    #                                LA_rank_prob[LA_rank_prob$best%in%c("best3[1]", "best3[2]", "best3[3]"),2]))
    # rownames(LA_rank_prob)<-c("trt_1", "trt_2", "trt_3")
    # 
    # 
    
    
    # #### GEMTC
    # 
    # network_abc <- gemtc::mtc.network(data.ab=all_data)
    # 
    # # fit NMA model
    # cons.model <- gemtc::mtc.model(network_abc, type="consistency", likelihood="binom", link="logit", linearModel="random")
    # 
    # cons.out = NULL
    # 
    # cons.out <- gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
    # 
    # # Rank order
    # prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
    # prob <- round(prob, digits=3)
    # sucra <- gemtc::sucra(prob)
    # rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
    # 
    # 
    # if (k_ab == 0 & k_ac == 0) {
    #   res = summary(gemtc::relative.effect(cons.out,"B",c("B","C")))
    #   lower_95 = res$summaries$quantiles[2,1]
    #   upper_95 = res$summaries$quantiles[2,5]
    #   true_rank = c("C","B")
    # } else {
    #   res = summary(gemtc::relative.effect(cons.out,"B",c("A","B","C")))
    #   lower_95 = res$summaries$quantiles[3,1]
    #   upper_95 = res$summaries$quantiles[3,5]
    #   true_rank = c("C","B","A")
    # }
    # 
    # reject_null = as.numeric(0 > upper_95 || 0 < lower_95)
    # rank_correct = as.numeric(identical(rank_order, true_rank))
    # c(reject_null, rank_correct)
  
  }
  
  power = result[1] / S 
  rank_correct_prob = result[2] / S
  
  
  if(verbose == T){
    print(paste("Power of k_ab =", k_ab, "k_ac =", k_ac, "pi_a =", pi_a, "OR_bc =", OR_bc, "tau =", tau, "is", power))
  }
  
  return(paste(power, rank_correct_prob))
}




k_ab = k_ac = 0
k_bc = 6

pi_a = 0.5

OR_ab = 1.2
OR_ac = 1.8

tau = 0.2
S = 10

power.sim.AB_direct <- function(S = 500, k_ab = 0, k_ac = 0, k_bc = 0, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.4, tau = 0.01, verbose = F){
  
  result = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
    
    library(R2jags)
    library(tidyverse)
    source("Models.R")
    
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
    
    
    # data pre for jags
    NS = n_distinct(all_data$study)
    NT = n_distinct(all_data$treatment)
    N = nrow(all_data)
    s = all_data$study
    # t = all_data$treatment
    t = as.integer(factor(all_data$treatment, levels = c("B","C"), labels = c(1,2)))
    y = all_data$responders
    n = all_data$sampleSize
    drug_list <- c("B","C")
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
                    'y'=y, 'n'=n ,'Omega'=diag(rep(0.2,times=2)),
                    'zero.AB' = (rep(0, times=2)))
    inits_AB<- list(list(mu=rep(0,2)),
                    list(mu=rep(0,2)))
    para_AB<-c( "lor", "tau", "best1", "best2")
    fit_AB<-jags(data=data_AB, inits=inits_AB, para_AB,
                 n.iter = 5000, n.burnin = 2000, n.chains = 2, n.thin = 1,
                 DIC=TRUE, model.file=ABWish)
    
    #output data
    
    AB_trt_results<-data.frame(fit_AB$BUGSoutput$summary[,c(1, 3, 7)])
    AB_trt_results <- tibble::rownames_to_column(AB_trt_results, "drug_list")
    AB_trt_results<-AB_trt_results%>%
      filter(drug_list %in% c("lor[1]", "lor[2]"))
    
    AB_rank_prob<-data.frame(fit_AB$BUGSoutput$summary[,c(1, 3, 7)])
    AB_rank_prob <- tibble::rownames_to_column(AB_rank_prob, "best")
    
    AB_rank_prob <- data.frame(cbind(AB_rank_prob[AB_rank_prob$best%in%c("best1[1]", "best1[2]"),2],
                                     AB_rank_prob[AB_rank_prob$best%in%c("best2[1]", "best2[2]"),2]))
    rownames(AB_rank_prob)<-c("trt_1", "trt_2")
    
    
    lower_95 = AB_trt_results$X2.5.[2]
    upper_95 = AB_trt_results$X97.5.[2]
    
    
    sucra <- gemtc::sucra(AB_rank_prob)
    rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
    true_rank = c("trt_1", "trt_2")
    
    
    reject_null = as.numeric(0 > upper_95 || 0 < lower_95)
    rank_correct = as.numeric(identical(rank_order, true_rank))
    c(reject_null, rank_correct)
    
    # 
    # 
    # #### GEMTC
    # 
    # network_abc <- gemtc::mtc.network(data.ab=all_data)
    # 
    # # fit NMA model
    # cons.model <- gemtc::mtc.model(network_abc, type="consistency", likelihood="binom", link="logit", linearModel="random")
    # 
    # cons.out = NULL
    # 
    # cons.out <- gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
    # 
    # # Rank order
    # prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
    # prob <- round(prob, digits=3)
    # sucra <- gemtc::sucra(prob)
    # rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
    # 
    # 
    # if (k_ab == 0 & k_ac == 0) {
    #   res = summary(gemtc::relative.effect(cons.out,"B",c("B","C")))
    #   lower_95 = res$summaries$quantiles[2,1]
    #   upper_95 = res$summaries$quantiles[2,5]
    #   true_rank = c("C","B")
    # } else {
    #   res = summary(gemtc::relative.effect(cons.out,"B",c("A","B","C")))
    #   lower_95 = res$summaries$quantiles[3,1]
    #   upper_95 = res$summaries$quantiles[3,5]
    #   true_rank = c("C","B","A")
    # }
    # 
    # reject_null = as.numeric(0 > upper_95 || 0 < lower_95)
    # rank_correct = as.numeric(identical(rank_order, true_rank))
    # c(reject_null, rank_correct)
  }
  
  power = result[1] / S 
  rank_correct_prob = result[2] / S
  
  
  if(verbose == T){
    print(paste("Power of k_ab =", k_ab, "k_ac =", k_ac, "pi_a =", pi_a, "OR_bc =", OR_bc, "tau =", tau, "is", power))
  }
  
  return(paste(power, rank_correct_prob))
}



# test
# power.sim.AB_direct(S = 20, k_ab = 0, k_ac = 0, k_bc = 3, pi_a = 0.1, OR_ab = 1.2, OR_ac = 1.6, tau = 0.1)
