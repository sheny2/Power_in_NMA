library(metafor)
library(tidyverse)


# power simulation function

set.seed(8848)

## Using indirect evidence (AB, AC) to estimate missing direct comparison (BC)

# all meta-analysis was fit in DerSimonian and Laird ("DL") method

power.sim_indirect <- function(S = 5000, k_ab, k_ac, pi_a, OR_ab, OR_ac, tau, verbose = F){
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
    result_ab = metafor::rma(measure = "OR", ai = e_aj, bi = n_aj - e_aj, ci = e_bj, di = n_bj - e_bj, 
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
  
  power = reject_correctly / S
  
  if(verbose == T){
  print(paste("Power of k_ab =", k_ab, "k_ac =", k_ac, "pi_a =", pi_a, "OR_bc =", "tau =", tau, "is", power))
  }
  
  return(power)
}


# test
power.sim_indirect(S = 5000, k_ab = 10, k_ac = 10, pi_a = 0.1, OR_ab = 1.2, OR_ac = 1.4, tau = 0.2)





## Have direct evidence (BC)

# Assume consistency and similarities


# pi_b <- (OR_ab) * (pi_a) / (1-pi_a) / (1 + (OR_ab) * (pi_a) / (1-pi_a))
# pi_c <- (OR_ac) * (pi_a) / (1-pi_a) / (1 + (OR_ac) * (pi_a) / (1-pi_a))
# OR_bc = exp(log(OR_ac) - log(OR_ab))


power.sim_direct <- function(S = 5000, k_bc, pi_a, OR_ab, OR_ac, tau, verbose = F){
  
  reject_correctly = 0

  for (i in 1:S) {
    
    # simulate k_bc studies (direct)
    OR_bc = exp(log(OR_ac) - log(OR_ab))
    pi_b = (OR_ab * pi_a / (1 - pi_a))/(1 + OR_ab * pi_a / (1-pi_a))

    dat_bc = data.frame()
    for (j in 1:k_bc)
    {
      n_j = runif(n = 1, min = 20, max = 500)
      n_bj = n_cj = round(n_j / 2)

      log_OR_bc_j = rnorm(n = 1, mean = log(OR_bc),  sd = tau)

      pi_bj = runif(n = 1, min = 1/2 * pi_b, max = 3/2 * pi_b)
      pi_cj = pi_bj * exp(log_OR_bc_j) / (1 - pi_bj + pi_bj * exp(log_OR_bc_j) )

      e_bj = rbinom(n = 1, size = n_bj, prob = pi_bj)
      e_cj = rbinom(n = 1, size = n_cj, prob = pi_cj)

      dat_bc = rbind(dat_bc, cbind(e_bj, n_bj, e_cj, n_cj))
    }

    # fit random effect MA for bc
    result_bc = rma(measure = "OR", ai = e_bj, bi = n_bj - e_bj, ci = e_cj, di = n_cj - e_cj,
                    data = dat_bc, method = "DL")
    dir_est_log_OR_bc = result_bc$beta[,1]
    dir_se_log_OR_bc = result_bc$se
    
    # Compute 95% CI
    lower_95 = exp(dir_est_log_OR_bc - 1.96 * dir_se_log_OR_bc)
    upper_95 = exp(dir_est_log_OR_bc + 1.96 * dir_se_log_OR_bc)
    
    if (1 > upper_95 || 1 < lower_95)
    { reject_correctly = reject_correctly + 1 }
    
  }
  
  power = reject_correctly / S
  
  if(verbose == T){
  print(paste("Power of k_bc =", k_bc, "pi_a =", pi_a, "OR_bc =", OR_bc, "tau =", tau, "is", power))
  }
  
  return(power)
}


power.sim_direct(S = 5000, k_bc = 10, pi_a = 0.1, OR_ab = 1.2, OR_ac = 1.4, tau = 0.2)
