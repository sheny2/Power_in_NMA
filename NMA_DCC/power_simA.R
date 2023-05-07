library(metafor)
library(tidyverse)
library(doMC)

# power simulation function

## Using indirect evidence (AB, AC) to estimate missing direct comparison (BC)

# all meta-analysis was fit in DerSimonian and Laird ("DL") method

power.sim_indirect <- function(S = 5000, k_ab, k_ac, pi_a, OR_ab, OR_ac, tau, verbose = F){

  reject_correctly = foreach (i = 1:S, .combine = "+") %dopar% {
    # simulate k_ab studies
    dat_ab = data.frame()
    for (j in 1:k_ab)
    {
      n_j = runif(n = 1, min = 100, max = 500)
      n_aj = n_bj = round(n_j / 2)
      
      log_OR_ab_j = rnorm(n = 1, mean = log(OR_ab),  sd = tau)
      
      pi_aj = runif(n = 1, min = 1/2 * pi_a, max = 3/2 * pi_a)
      pi_bj = pi_aj * exp(log_OR_ab_j) / (1 - pi_aj + pi_aj * exp(log_OR_ab_j) )
      
      e_aj = rbinom(n = 1, size = n_aj, prob = pi_aj)
      e_bj = rbinom(n = 1, size = n_bj, prob = pi_bj)
      
      dat_ab = rbind(dat_ab, cbind(e_aj, n_aj, e_bj, n_bj))
    }
    
    # fit random effect MA
    result_ab = metafor::rma(measure = "OR", ai = e_aj, bi = n_aj - e_aj, ci = e_bj, di = n_bj - e_bj, 
                             data = dat_ab, method = "DL")
    est_log_OR_ab = result_ab$beta[,1]
    se_log_OR_ab = result_ab$se
    
    # simulate k_ac studies
    dat_ac = data.frame()
    for (j in 1:k_ac)
    {
      n_j = runif(n = 1, min = 100, max = 500)
      n_aj = n_cj = round(n_j / 2)
      
      log_OR_ac_j = rnorm(n = 1, mean = log(OR_ac),  sd = tau)
      
      pi_aj = runif(n = 1, min = 1/2 * pi_a, max = 3/2 * pi_a)
      pi_cj = pi_aj * exp(log_OR_ac_j) / (1 - pi_aj + pi_aj * exp(log_OR_ac_j) )
      
      e_aj = rbinom(n = 1, size = n_aj, prob = pi_aj)
      e_cj = rbinom(n = 1, size = n_cj, prob = pi_cj)  
      
      dat_ac = rbind(dat_ac, cbind(e_aj, n_aj, e_cj, n_cj))
    }
    
    # fit random effect MA
    result_ac = metafor::rma(measure = "OR", ai = e_aj, bi = n_aj - e_aj, ci = e_cj, di = n_cj - e_cj, 
                    data = dat_ac, method = "DL")
    est_log_OR_ac = result_ac$beta[,1]
    se_log_OR_ac = result_ac$se
    
    
    # Bucher's Method: indirect estimate of log(OR_bc)
    est_log_OR_bc = est_log_OR_ab - est_log_OR_ac
    se_log_OR_bc = sqrt(se_log_OR_ab^2 + se_log_OR_ac^2)
    
    # Compute 95% CI
    upper_95 = est_log_OR_bc + 1.96 * se_log_OR_bc
    lower_95 = est_log_OR_bc - 1.96 * se_log_OR_bc
    
    as.numeric(0 > upper_95 || 0 < lower_95)
  }

  
  power = reject_correctly / S
  
  if(verbose == T){
    print(paste("Power of k_ab =", k_ab, "k_ac =", k_ac, "pi_a =", pi_a, "OR_bc =", OR_bc, "tau =", tau, "is", power))
  }
  
  return(power)
}


# test
# power.sim_indirect(S = 1000, k_ab = 10, k_ac = 10, pi_a = 0.1, OR_ab = 1.2, OR_ac = 1.4, tau = 0.1)
# power.sim_indirect(S = 1000, k_ab = 6, k_ac = 6, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.4, tau = 0.4)