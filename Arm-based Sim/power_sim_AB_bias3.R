library(doMC)
library(gemtc)
# source("Models.R")

k_ab = k_ac = 6
k_bc = 6

pi_a = 0.5

OR_ab = 1.2
OR_ac = 1.8

tau = 0.2
S = 5

power.sim.AB_full <- function(S = 500, k_ab = 0, k_ac = 0, k_bc = 0, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.4, tau = 0.01, verbose = F){
  
  result = foreach::foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
    
    library(R2jags)
    library(tidyverse)
    library(pcnetmeta)
    
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
    
    # het_cor
    AB_Result_het_cor = nma.ab.bin(s.id, t.id, r, n, data = all_data_pc, param= "LOR",
                                   model = "het_cor", n.adapt = 2000, n.iter = 5000, n.chains = 2)
    AB_Result_het_cor_PointEst = as.numeric(sub("\\s*\\(.*", "", AB_Result_het_cor$LogOddsRatio$Mean_SD[3,2]))
    CI95 = gsub(".*\\(([^,]+),\\s+([^,]+)\\).*", "\\1 \\2", AB_Result_het_cor$LogOddsRatio$Median_CI[3,2])
    AB_Result_het_cor_lower_95 <- as.numeric(unlist(strsplit(CI95, " "))[1])
    AB_Result_het_cor_upper_95 <- as.numeric(unlist(strsplit(CI95, " "))[2])
    
    reject_null_het_cor = as.numeric(0 > AB_Result_het_cor_upper_95 || 0 < AB_Result_het_cor_lower_95)
    bias_het_cor = abs(AB_Result_het_cor_PointEst - log(OR_bc))
    
    # het_eqcor
    AB_Result_het_eqcor = nma.ab.bin(s.id, t.id, r, n, data = all_data_pc, param= "LOR",
                                     model = "het_eqcor", n.adapt = 2000, n.iter = 5000, n.chains = 2)
    AB_Result_het_eqcor_PointEst = as.numeric(sub("\\s*\\(.*", "", AB_Result_het_eqcor$LogOddsRatio$Mean_SD[3,2]))
    CI95 = gsub(".*\\(([^,]+),\\s+([^,]+)\\).*", "\\1 \\2", AB_Result_het_eqcor$LogOddsRatio$Median_CI[3,2])
    AB_Result_het_eqcor_lower_95 <- as.numeric(unlist(strsplit(CI95, " "))[1])
    AB_Result_het_eqcor_upper_95 <- as.numeric(unlist(strsplit(CI95, " "))[2])
    
    reject_null_het_eqcor = as.numeric(0 > AB_Result_het_eqcor_upper_95 || 0 < AB_Result_het_eqcor_lower_95)
    bias_het_eqcor = abs(AB_Result_het_eqcor_PointEst - log(OR_bc))
    
    # hom_eqcor
    AB_Result_hom_eqcor = nma.ab.bin(s.id, t.id, r, n, data = all_data_pc, param= "LOR",
                                     model = "hom_eqcor", n.adapt = 2000, n.iter = 5000, n.chains = 2)
    AB_Result_hom_eqcor_PointEst = as.numeric(sub("\\s*\\(.*", "", AB_Result_hom_eqcor$LogOddsRatio$Mean_SD[3,2]))
    CI95 = gsub(".*\\(([^,]+),\\s+([^,]+)\\).*", "\\1 \\2", AB_Result_hom_eqcor$LogOddsRatio$Median_CI[3,2])
    AB_Result_hom_eqcor_lower_95 <- as.numeric(unlist(strsplit(CI95, " "))[1])
    AB_Result_hom_eqcor_upper_95 <- as.numeric(unlist(strsplit(CI95, " "))[2])
    
    reject_null_hom_eqcor = as.numeric(0 > AB_Result_hom_eqcor_upper_95 || 0 < AB_Result_hom_eqcor_lower_95)
    bias_hom_eqcor = abs(AB_Result_hom_eqcor_PointEst - log(OR_bc))
    
    
    c(reject_null_het_cor, reject_null_het_eqcor, reject_null_hom_eqcor,
      bias_het_cor, bias_het_eqcor, bias_hom_eqcor)
  }
  

  power_het_cor = result[1] / S
  power_het_eqcor = result[2] / S
  power_hom_eqcor = result[3] / S
  avg_bias_het_cor = result[4] / S
  avg_bias_het_eqcor = result[5] / S
  avg_bias_hom_eqcor = result[6] / S
  
  return(paste(power_het_cor, power_het_eqcor, power_hom_eqcor,
               avg_bias_het_cor, avg_bias_het_eqcor, avg_bias_hom_eqcor))
}


# 
# power.sim.AB_full(S = 5, k_ab = 2, k_ac = 2, k_bc = 2, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.4, tau = 0.01, verbose = F)
# power.sim.AB_full(S = 5, k_ab = 6, k_ac = 6, k_bc = 6, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.8, tau = 0.02, verbose = F)



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
    library(pcnetmeta)
    
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
    
    
    # het_cor
    AB_Result_het_cor = nma.ab.bin(s.id, t.id, r, n, data = all_data_pc, param= "LOR",
                                   model = "het_cor", n.adapt = 1000, n.iter = 3000, n.chains = 2)
    AB_Result_het_cor_PointEst = as.numeric(sub("\\s*\\(.*", "", AB_Result_het_cor$LogOddsRatio$Mean_SD[2,1]))
    CI95 = gsub(".*\\(([^,]+),\\s+([^,]+)\\).*", "\\1 \\2", AB_Result_het_cor$LogOddsRatio$Median_CI[2,1])
    AB_Result_het_cor_lower_95 <- as.numeric(unlist(strsplit(CI95, " "))[1])
    AB_Result_het_cor_upper_95 <- as.numeric(unlist(strsplit(CI95, " "))[2])
    
    reject_null_het_cor = as.numeric(0 > AB_Result_het_cor_upper_95 || 0 < AB_Result_het_cor_lower_95)
    bias_het_cor = abs(AB_Result_het_cor_PointEst - log(OR_bc))
    
    # het_eqcor
    AB_Result_het_eqcor = nma.ab.bin(s.id, t.id, r, n, data = all_data_pc, param= "LOR",
                                     model = "het_eqcor", n.adapt = 1000, n.iter = 3000, n.chains = 2)
    AB_Result_het_eqcor_PointEst = as.numeric(sub("\\s*\\(.*", "", AB_Result_het_eqcor$LogOddsRatio$Mean_SD[2,1]))
    CI95 = gsub(".*\\(([^,]+),\\s+([^,]+)\\).*", "\\1 \\2", AB_Result_het_eqcor$LogOddsRatio$Median_CI[2,1])
    AB_Result_het_eqcor_lower_95 <- as.numeric(unlist(strsplit(CI95, " "))[1])
    AB_Result_het_eqcor_upper_95 <- as.numeric(unlist(strsplit(CI95, " "))[2])
    
    reject_null_het_eqcor = as.numeric(0 > AB_Result_het_eqcor_upper_95 || 0 < AB_Result_het_eqcor_lower_95)
    bias_het_eqcor = abs(AB_Result_het_eqcor_PointEst - log(OR_bc))
    
    # hom_eqcor
    AB_Result_hom_eqcor = nma.ab.bin(s.id, t.id, r, n, data = all_data_pc, param= "LOR",
                                     model = "hom_eqcor", n.adapt = 1000, n.iter = 3000, n.chains = 2)
    AB_Result_hom_eqcor_PointEst = as.numeric(sub("\\s*\\(.*", "", AB_Result_hom_eqcor$LogOddsRatio$Mean_SD[2,1]))
    CI95 = gsub(".*\\(([^,]+),\\s+([^,]+)\\).*", "\\1 \\2", AB_Result_hom_eqcor$LogOddsRatio$Median_CI[2,1])
    AB_Result_hom_eqcor_lower_95 <- as.numeric(unlist(strsplit(CI95, " "))[1])
    AB_Result_hom_eqcor_upper_95 <- as.numeric(unlist(strsplit(CI95, " "))[2])
    
    reject_null_hom_eqcor = as.numeric(0 > AB_Result_hom_eqcor_upper_95 || 0 < AB_Result_hom_eqcor_lower_95)
    bias_hom_eqcor = abs(AB_Result_hom_eqcor_PointEst - log(OR_bc))
    
    
    c(reject_null_het_cor, reject_null_het_eqcor, reject_null_hom_eqcor,
      bias_het_cor, bias_het_eqcor, bias_hom_eqcor)
  }
  
  
  power_het_cor = result[1] / S
  power_het_eqcor = result[2] / S
  power_hom_eqcor = result[3] / S
  avg_bias_het_cor = result[4] / S
  avg_bias_het_eqcor = result[5] / S
  avg_bias_hom_eqcor = result[6] / S
  
  return(paste(power_het_cor, power_het_eqcor, power_hom_eqcor,
               avg_bias_het_cor, avg_bias_het_eqcor, avg_bias_hom_eqcor))
}



# test
# power.sim.AB_direct(S = 3, k_ab = 0, k_ac = 0, k_bc = 3, pi_a = 0.1, OR_ab = 1.2, OR_ac = 1.6, tau = 0.1)
