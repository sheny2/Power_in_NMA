library(doMC)
library(gemtc)

k_ab = k_ac = k_bc = 2

pi_a = 0.2

OR_ab = 1.2
OR_ac = 1.4

tau = 0.1

power.sim.NMA <- function(S = 5000, k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, tau, verbose = F){
  
  reject_correctly = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
    
    # simulate k_ab studies (indirect)
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
      
      dat_ab = rbind(dat_ab, rbind(cbind(j, "A", n_aj, e_aj), cbind(j, "B", n_bj, e_bj)))
    }
    
    
    # simulate k_ac studies (indirect)
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
      
      dat_ac = rbind(dat_ac, rbind(cbind(j+k_ab, "A", n_aj, e_aj), cbind(j+k_ab, "C", n_cj, e_cj)))
    }
    
    
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
      
      dat_bc = rbind(dat_bc, rbind(cbind(j+k_ab+k_ac, "B", n_bj, e_bj), cbind(j+k_ab+k_ac, "C", n_cj, e_cj)))
    }
    
    
    colnames(dat_ab) = colnames(dat_ac) = colnames(dat_bc) = c("study", "treatment", "sampleSize", "responders")
    all_data = rbind(dat_ab, dat_ac, dat_bc)
    
    all_data[,c('sampleSize', 'responders')] <-
      lapply(all_data[,c('sampleSize', 'responders')], as.numeric)
    
    network_abc <- gemtc::mtc.network(data.ab=all_data)
    
    # fit NMA model
    cons.model <- gemtc::mtc.model(network_abc, type="consistency", likelihood="binom", link="logit", linearModel="random")
    
    cons.out = NULL
    
    cons.out <- 
      tryCatch({
        gemtc::mtc.run(cons.model, n.adapt=500, n.iter=2000, thin=1)
      }, error = function(e) {
        cat("Error:", conditionMessage(e), "\n")
      })
    
    
    if (is.null(cons.out)) 
    {print(paste("A model didn't converge for k_ab =", k_ab, "k_ac =", k_ac, "pi_a =", pi_a, "OR_bc =", OR_bc, "tau =", tau))
      0}
    else 
    {res = summary(gemtc::relative.effect(cons.out,"B",c("A","B","C")))
    lower_95 = res$summaries$quantiles[3,1]
    upper_95 = res$summaries$quantiles[3,5]
    
    as.numeric(0 > upper_95 || 0 < lower_95)}
  }
  
  power = reject_correctly / S
  
  if(verbose == T){
    print(paste("Power of k_ab =", k_ab, "k_ac =", k_ac, "pi_a =", pi_a, "OR_bc =", OR_bc, "tau =", tau, "is", power))
  }
  return(power)
}

# test
# power.sim.NMA(S = 1000, k_ab = 10, k_ac = 10, k_bc = 10, pi_a = 0.2, OR_ab = 1.2, OR_ac = 1.4, tau = 0.1)
