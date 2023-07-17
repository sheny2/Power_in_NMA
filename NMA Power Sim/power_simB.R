library(metafor)
library(tidyverse)
library(doMC)

## Have direct evidence (BC)

# Assume consistency and similarities


# pi_b <- (OR_ab) * (pi_a) / (1-pi_a) / (1 + (OR_ab) * (pi_a) / (1-pi_a))
# pi_c <- (OR_ac) * (pi_a) / (1-pi_a) / (1 + (OR_ac) * (pi_a) / (1-pi_a))
# OR_bc = exp(log(OR_ac) - log(OR_ab))


power.sim_direct <- function(S = 5000, k_bc, pi_a, OR_ab, OR_ac, tau, verbose = F){
  
  reject_correctly = foreach (i = 1:S, .combine = "+") %dopar% {
    
    # simulate k_bc studies (direct)
    OR_bc = exp(log(OR_ac) - log(OR_ab))
    pi_b = (OR_ab * pi_a / (1 - pi_a))/(1 + OR_ab * pi_a / (1-pi_a))
    
    dat_bc = data.frame()
    for (j in 1:k_bc)
    {
      n_j = runif(n = 1, min = 100, max = 500)
      n_bj = n_cj = round(n_j / 2)
      
      log_OR_bc_j = rnorm(n = 1, mean = log(OR_bc),  sd = tau)
      
      pi_bj = runif(n = 1, min = 1/2 * pi_b, max = 3/2 * pi_b)
      pi_cj = pi_bj * exp(log_OR_bc_j) / (1 - pi_bj + pi_bj * exp(log_OR_bc_j) )
      
      e_bj = rbinom(n = 1, size = n_bj, prob = pi_bj)
      e_cj = rbinom(n = 1, size = n_cj, prob = pi_cj)
      
      dat_bc = rbind(dat_bc, cbind(e_bj, n_bj, e_cj, n_cj))
    }
    
    # fit random effect MA for bc
    result_bc = metafor::rma(measure = "OR", ai = e_bj, bi = n_bj - e_bj, ci = e_cj, di = n_cj - e_cj,
                    data = dat_bc, method = "DL")
    dir_est_log_OR_bc = result_bc$beta[,1]
    dir_se_log_OR_bc = result_bc$se
    
    # Compute 95% CI
    lower_95 = dir_est_log_OR_bc - 1.96 * dir_se_log_OR_bc
    upper_95 = dir_est_log_OR_bc + 1.96 * dir_se_log_OR_bc
    
    as.numeric(0 > upper_95 || 0 < lower_95)
  }
  
  power = reject_correctly / S
  
  if(verbose == T){
    print(paste("Power of k_bc =", k_bc, "pi_a =", pi_a, "OR_bc =", OR_bc, "tau =", tau, "is", power))
  }
  
  return(power)
}


# test
# power.sim_direct(S = 5000, k_bc = 10, pi_a = 0.1, OR_ab = 1.2, OR_ac = 1.4, tau = 0.1)
