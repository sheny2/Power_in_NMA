library(doMC)
# source("Model2.R")

k_ab = k_ac = 6
k_bc = 6

pi_a = 0.5

OR_ab = 1.2
OR_ac = 1.8

tau = 0.2
S = 5

power.sim.AB_full <- function(S = 500, k_ab = 0, k_ac = 0, k_bc = 0, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.4, tau = 0.01, verbose = F){
  
  result = foreach::foreach (i = 1:S, .combine = "+") %dopar% {
    
    library(R2jags)
    library(tidyverse)
    source("Model2.R")
    probability_to_log_odds <- function(probability) {
      if (probability <= 0 | probability >= 1) {
        stop("Probability must be between 0 and 1.")
      }
      log_odds <- log(probability / (1 - probability))
      return(log_odds)
    }
    
    
    inverse_logit <- function(x) {
      odds <- exp(x)  
      probability <- odds / (1 + odds) 
      return(probability)  
    }
    
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
    
    
    
    
    
    # Generate data using AB
    
    n_study = n_distinct(all_data$study)
    trt = unique(all_data$treatment)
    
    pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab ) 
    pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )
    
    mean_vec = c(probability_to_log_odds(pi_a), probability_to_log_odds(pi_b), probability_to_log_odds(pi_c))
    
    COV_mat = matrix(nrow = 3, ncol = 3)
    sigma = tau
    Ndrug = 3
    
    rho = runif(1, -1/(Ndrug - 1), 0.9999)
    diag <- (1 + (Ndrug - 2)*rho)/(1 + (Ndrug - 2)*rho - (Ndrug - 1)*rho^2)
    offdiag <- (-rho/(1 + (Ndrug - 2)*rho - (Ndrug - 1)*rho^2))
    
    for(j in 1:Ndrug){
      for(k in 1:Ndrug){ 
        COV_mat[j,k] <- 1/sigma^2*ifelse(j == k, diag, offdiag)
      }
    }
    
    COV_mat_T = solve(COV_mat)
    
    
    y_ik = c()
    for (i in 1:n_study){
      study_data = all_data[all_data$study==i,]    
      delta_sim <- MASS::mvrnorm(1, mean_vec, COV_mat_T)
      Arm_prob = inverse_logit(delta_sim)
      
      for (k in 1:nrow(study_data)){
        p_ik = Arm_prob[which(trt == study_data$treatment[k])]
        y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
      }
    }
    
    all_data$responders = y_ik
    
    
    
    # data pre for jags
    NS = n_distinct(all_data$study)
    NT = n_distinct(all_data$treatment)
    N = nrow(all_data)
    s = all_data$study
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
    
    data_AB <- list('Narm'=N, 'Nstudy'=NS, 
                    'Ndrug'=NT, 'study'= s, 'drug'=t, 
                    'y'=y, 'n'=n ,
                    'zero.AB' = (rep(0, times=3)))
    inits_AB<- list(list(mu=rep(0,3)),
                    list(mu=rep(0,3)))
    para_AB<-c( "lor","rho", "sigma", "or")
    fit_AB_homo.eqcor<-jags(data=data_AB, inits=inits_AB, para_AB,
                            n.iter=3000, n.burnin = 1000, n.chains = 2, n.thin = 1,
                            DIC=TRUE, model.file=NMA.homo.eqcor)
    
    
    AB_trt_results<-data.frame(fit_AB_homo.eqcor$BUGSoutput$summary[,c(1, 3, 7)])
    AB_trt_results <- tibble::rownames_to_column(AB_trt_results, "drug_list")
    AB_trt_results<-AB_trt_results%>%
      filter(drug_list %in% c("lor[1]", "lor[2]", "lor[3]", "or[1]", "or[2]", "or[3]"))
    
    point_estimate = AB_trt_results$mean[6]
    lower_95 = AB_trt_results$X2.5.[3]
    upper_95 = AB_trt_results$X97.5.[3]
    
    
    # all_data_pc = tibble(s.id = all_data$study,
    #                      t.id = all_data$treatment, 
    #                      r = all_data$responders, 
    #                      n = all_data$sampleSize) 
    # # het_cor
    # AB_Result_het_cor = nma.ab.bin(s.id, t.id, r, n, data = all_data_pc, param= "LOR",
    #                                model = "het_cor", n.adapt = 2000, n.iter = 5000, n.chains = 2)
    # AB_Result_het_cor_PointEst = as.numeric(sub("\\s*\\(.*", "", AB_Result_het_cor$LogOddsRatio$Mean_SD[3,2]))
    # CI95 = gsub(".*\\(([^,]+),\\s+([^,]+)\\).*", "\\1 \\2", AB_Result_het_cor$LogOddsRatio$Median_CI[3,2])
    # AB_Result_het_cor_lower_95 <- as.numeric(unlist(strsplit(CI95, " "))[1])
    # AB_Result_het_cor_upper_95 <- as.numeric(unlist(strsplit(CI95, " "))[2])
    # 
    # reject_null_het_cor = as.numeric(0 > AB_Result_het_cor_upper_95 || 0 < AB_Result_het_cor_lower_95)
    # bias_het_cor = abs(AB_Result_het_cor_PointEst - log(OR_bc))
    
    reject_null_het_cor = as.numeric(0 > upper_95 || 0 < lower_95)
    bias_het_cor = abs(point_estimate - log(OR_bc))
    c(reject_null_het_cor, bias_het_cor)
  }
  
  
  power = result[1] / S
  avg_bias_abs = result[2] / S
  
  return(paste(power, avg_bias_abs))
}


# 
# power.sim.AB_full(S = 5, k_ab = 2, k_ac = 2, k_bc = 2, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.4, tau = 0.01, verbose = F)
# power.sim.AB_full(S = 5, k_ab = 6, k_ac = 6, k_bc = 6, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.8, tau = 0.02)



k_ab = k_ac = 0
k_bc = 3

pi_a = 0.5

OR_ab = 1.2
OR_ac = 1.8

tau = 0.2
S = 5

power.sim.AB_direct <- function(S = 500, k_ab = 0, k_ac = 0, k_bc = 0, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.4, tau = 0.01, verbose = F){
  
  result = foreach::foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
    
    library(R2jags)
    library(tidyverse)
    source("Model2.R")
    probability_to_log_odds <- function(probability) {
      if (probability <= 0 | probability >= 1) {
        stop("Probability must be between 0 and 1.")
      }
      log_odds <- log(probability / (1 - probability))
      return(log_odds)
    }
    
    
    inverse_logit <- function(x) {
      odds <- exp(x)  
      probability <- odds / (1 + odds) 
      return(probability)  
    }
    
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
    
    # all_data_pc = tibble(s.id = all_data$study,
    #                      t.id = all_data$treatment, 
    #                      r = all_data$responders, 
    #                      n = all_data$sampleSize) 
    
    
    # Generate data using AB
    
    n_study = n_distinct(all_data$study)
    trt = unique(all_data$treatment)
    
    pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab ) 
    pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )
    
    mean_vec = c(probability_to_log_odds(pi_b), probability_to_log_odds(pi_c))
    
    COV_mat = matrix(nrow = 2, ncol = 2)
    sigma = tau
    Ndrug = 2
    
    rho = runif(1, -1/(Ndrug - 1), 0.9999)
    diag <- (1 + (Ndrug - 2)*rho)/(1 + (Ndrug - 2)*rho - (Ndrug - 1)*rho^2)
    offdiag <- (-rho/(1 + (Ndrug - 2)*rho - (Ndrug - 1)*rho^2))
    
    for(j in 1:Ndrug){
      for(k in 1:Ndrug){ 
        COV_mat[j,k] <- 1/sigma^2*ifelse(j == k, diag, offdiag)
      }
    }
    
    COV_mat_T = solve(COV_mat)
    
    
    y_ik = c()
    for (i in 1:n_study){
      study_data = all_data[all_data$study==i,]    
      delta_sim <- MASS::mvrnorm(1, mean_vec, COV_mat_T)
      Arm_prob = inverse_logit(delta_sim)
      
      for (k in 1:nrow(study_data)){
        p_ik = Arm_prob[which(trt == study_data$treatment[k])]
        y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
      }
    }
    
    all_data$responders = y_ik
    
    
    
    # data pre for jags
    NS = n_distinct(all_data$study)
    NT = n_distinct(all_data$treatment)
    N = nrow(all_data)
    s = all_data$study
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
    
    data_AB <- list('Narm'=N, 'Nstudy'=NS, 
                    'Ndrug'=NT, 'study'= s, 'drug'=t, 
                    'y'=y, 'n'=n ,
                    'zero.AB' = (rep(0, times=2)))
    inits_AB<- list(list(mu=rep(0,2)),
                    list(mu=rep(0,2)))
    para_AB<-c( "lor","rho", "sigma", "or")
    fit_AB_homo.eqcor<-jags(data=data_AB, inits=inits_AB, para_AB,
                            n.iter=3000, n.burnin = 1000, n.chains = 2, n.thin = 1,
                            DIC=TRUE, model.file=NMA.homo.eqcor)
    
    
    AB_trt_results<-data.frame(fit_AB_homo.eqcor$BUGSoutput$summary[,c(1, 3, 7)])
    AB_trt_results <- tibble::rownames_to_column(AB_trt_results, "drug_list")
    AB_trt_results<-AB_trt_results%>%
      filter(drug_list %in% c("lor[1]", "lor[2]", "or[1]", "or[2]"))
    
    point_estimate = AB_trt_results$mean[4]
    lower_95 = AB_trt_results$X2.5.[2]
    upper_95 = AB_trt_results$X97.5.[2]
    
    
    # all_data_pc = tibble(s.id = all_data$study,
    #                      t.id = all_data$treatment, 
    #                      r = all_data$responders, 
    #                      n = all_data$sampleSize) 
    # # het_cor
    # AB_Result_het_cor = nma.ab.bin(s.id, t.id, r, n, data = all_data_pc, param= "LOR",
    #                                model = "het_cor", n.adapt = 2000, n.iter = 5000, n.chains = 2)
    # AB_Result_het_cor_PointEst = as.numeric(sub("\\s*\\(.*", "", AB_Result_het_cor$LogOddsRatio$Mean_SD[3,2]))
    # CI95 = gsub(".*\\(([^,]+),\\s+([^,]+)\\).*", "\\1 \\2", AB_Result_het_cor$LogOddsRatio$Median_CI[3,2])
    # AB_Result_het_cor_lower_95 <- as.numeric(unlist(strsplit(CI95, " "))[1])
    # AB_Result_het_cor_upper_95 <- as.numeric(unlist(strsplit(CI95, " "))[2])
    # 
    # reject_null_het_cor = as.numeric(0 > AB_Result_het_cor_upper_95 || 0 < AB_Result_het_cor_lower_95)
    # bias_het_cor = abs(AB_Result_het_cor_PointEst - log(OR_bc))
    
    reject_null_het_cor = as.numeric(0 > upper_95 || 0 < lower_95)
    bias_het_cor = abs(point_estimate - log(OR_bc))
    c(reject_null_het_cor, bias_het_cor)
  }
  
  
  power = result[1] / S
  avg_bias_abs = result[2] / S
  
  return(paste(power, avg_bias_abs))
}


# test
# power.sim.AB_direct(S = 3, k_ab = 0, k_ac = 0, k_bc = 3, pi_a = 0.1, OR_ab = 1.2, OR_ac = 1.6, tau = 0.1)
