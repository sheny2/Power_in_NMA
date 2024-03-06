library(metafor)
library(tidyverse)

# Construct Table output
set.seed(12345)

S = 1000

k_ab_all <- c(5, 10, 25, 100)
k_ac <- 5

tau_all <- c(0.001, 0.2, 0.4)

OR_ab <- 1.2  # the true relative effect of B versus A
OR_ac <- 1.4  # the true relative effect of C versus A

pi_a <- 0.3   # the true average event rate in the common comparator group A

power_result = matrix(NA, nrow = length(k_ab_all), ncol = length(tau_all))

for (k in 1:length(k_ab_all)){
  for (t in 1:length(tau_all)){
    k_ab = k_ab_all[k]
    tau = tau_all[t]
    
    reject_correctly = 0
    
    for (i in 1:S) {
      # simulate k_ab studies
      dat_ab = data.frame()
      for (j in 1:k_ab)
      {
        n_j = runif(n = 1, min = 20, max = 500)
        n_aj = n_bj = round(n_j / 2)
        
        log_OR_ab_j = rnorm(n = 1, mean = log(OR_ab),  sd = tau)
        
        pi_aj = runif(n = 1, min = 1/2 * pi_a, max = 3/2 * pi_a)
        pi_bj = pi_aj * exp(log_OR_ab_j) / (1 - pi_aj + pi_aj * exp(log_OR_ab_j) )
        
        e_aj = rbinom(n = 1, size = n_aj, prob = pi_aj)
        e_bj = rbinom(n = 1, size = n_bj, prob = pi_bj)
        
        dat_ab = rbind(dat_ab, cbind(e_aj, n_aj, e_bj, n_bj))
      }
      
      # fit random effect MA
      result_ab = rma(measure = "OR", ai = e_aj, bi = n_aj - e_aj, ci = e_bj, di = n_bj - e_bj, 
                      data = dat_ab, method = "DL")
      est_log_OR_ab = result_ab$beta[,1]
      se_log_OR_ab = result_ab$se
      
      # simulate k_ac studies
      dat_ac = data.frame()
      for (j in 1:k_ac)
      {
        n_j = runif(n = 1, min = 20, max = 500)
        n_aj = n_cj = round(n_j / 2)
        
        log_OR_ac_j = rnorm(n = 1, mean = log(OR_ac),  sd = tau)
        
        pi_aj = runif(n = 1, min = 1/2 * pi_a, max = 3/2 * pi_a)
        pi_cj = pi_aj * exp(log_OR_ac_j) / (1 - pi_aj + pi_aj * exp(log_OR_ac_j) )
        
        e_aj = rbinom(n = 1, size = n_aj, prob = pi_aj)
        e_cj = rbinom(n = 1, size = n_cj, prob = pi_cj)  
        
        dat_ac = rbind(dat_ac, cbind(e_aj, n_aj, e_cj, n_cj))
      }
      
      # fit random effect MA
      result_ac = rma(measure = "OR", ai = e_aj, bi = n_aj - e_aj, ci = e_cj, di = n_cj - e_cj, 
                      data = dat_ac, method = "DL")
      est_log_OR_ac = result_ac$beta[,1]
      se_log_OR_ac = result_ac$se
      
      
      # Bucher's Method: indirect estimate of log(OR_bc)
      est_log_OR_bc = est_log_OR_ab - est_log_OR_ac
      se_log_OR_bc = sqrt(se_log_OR_ab^2 + se_log_OR_ac^2)
      
      # Compute 95% CI
      upper_95 = exp(est_log_OR_bc + 1.96 * se_log_OR_bc)
      lower_95 = exp(est_log_OR_bc - 1.96 * se_log_OR_bc)
      
      if (1 > upper_95 || 1 < lower_95)
      { reject_correctly = reject_correctly + 1 }
    }
    
    power_result[k,t] = reject_correctly / S
  }
  
}

colnames(power_result) = c("tau = 0.001", "tau = 0.2", "tau = 0.4")
rownames(power_result) = c("k_ab = 5", "k_ab = 10", "k_ab = 25", "k_ab = 100")
power_result  # k_ac = 5, pi_a = 0.3, OR_ab <- 1.2, OR_ac <- 1.4

