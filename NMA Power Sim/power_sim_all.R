library(doMC)
library(gemtc)

k_ab = k_ac = 6
k_bc = 6

pi_a = 0.5

OR_ab = 1.2
OR_ac = 1.8

tau = 0.2
S = 10

power.sim.NMA1 <- function(S = 500, k_ab = 0, k_ac = 0, k_bc = 0, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.4, tau = 0.01, verbose = F){
  
  result = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
    
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
    
    network_abc <- gemtc::mtc.network(data.ab=all_data)
    
    # fit NMA model
    cons.model <- gemtc::mtc.model(network_abc, type="consistency", likelihood="binom", link="logit", linearModel="random")
    
    cons.out = NULL
    
    cons.out <- gemtc::mtc.run(cons.model, n.adapt=500, n.iter=3000, thin=1)
    
    # Rank order
    prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
    prob <- round(prob, digits=3)
    sucra <- gemtc::sucra(prob)
    rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
    
    
    if (k_ab == 0 & k_ac == 0) {
      res = summary(gemtc::relative.effect(cons.out,"B",c("B","C")))
      lower_95 = res$summaries$quantiles[2,1]
      upper_95 = res$summaries$quantiles[2,5]
      true_rank = c("C","B")
    } else {
      res = summary(gemtc::relative.effect(cons.out,"B",c("A","B","C")))
      lower_95 = res$summaries$quantiles[3,1]
      upper_95 = res$summaries$quantiles[3,5]
      true_rank = c("C","B","A")
    }
    
    reject_null = as.numeric(0 > upper_95 || 0 < lower_95)
    rank_correct = as.numeric(identical(rank_order, true_rank))
    c(reject_null, rank_correct)
  }
  
  power = result[1] / S 
  rank_correct_prob = result[2] / S
  
  
  if(verbose == T){
    print(paste("Power of k_ab =", k_ab, "k_ac =", k_ac, "pi_a =", pi_a, "OR_bc =", OR_bc, "tau =", tau, "is", power))
  }
  
  return(paste(power,rank_correct_prob))
}






library(doMC)
library(gemtc)

k_ab = k_ac = 6
k_bc = 6

pi_a = 0.5

OR_ab = 1.2
OR_ac = 1.8

tau = 0.2
S = 10

power.sim.NMA2 <- function(S = 500, k_ab = 0, k_ac = 0, k_bc = 0, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.4, tau = 0.01, verbose = F){
  
  result = foreach (i = 1:S, .combine = rbind, .errorhandling='remove') %dopar% {
    
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
    
    network_abc <- gemtc::mtc.network(data.ab=all_data)
    
    # fit NMA model
    cons.model <- gemtc::mtc.model(network_abc, type="consistency", likelihood="binom", link="logit", linearModel="random")
    
    cons.out = NULL
    
    cons.out <- gemtc::mtc.run(cons.model, n.adapt=500, n.iter=5000, thin=1)
    
    # Rank order
    prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
    prob <- round(prob, digits=3)
    sucra <- gemtc::sucra(prob)
    rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
    
    
    if (k_ab == 0 & k_ac == 0) {
      res = summary(gemtc::relative.effect(cons.out,"B",c("B","C")))
      lower_95 = res$summaries$quantiles[2,1]
      upper_95 = res$summaries$quantiles[2,5]
      point_est = res$summaries$statistics[2,1]
      true_rank = c("C","B")
    } else {
      res = summary(gemtc::relative.effect(cons.out,"B",c("A","B","C")))
      lower_95 = res$summaries$quantiles[3,1]
      upper_95 = res$summaries$quantiles[3,5]
      point_est = res$summaries$statistics[3,1]
      true_rank = c("C","B","A")
    }
    
    reject_null = as.numeric(0 > upper_95 || 0 < lower_95)
    rank_correct = as.numeric(identical(rank_order, true_rank))
    
    point_true = log(OR_bc)
    
    c(reject_null, rank_correct, point_est, point_true)
  }
  # 
  # power = sum(result[,1]) / S 
  # rank_correct_prob = sum(result[,2]) / S
  # 
  # bias = sum(result[,3]-result[,4]) / S
  
  
  if(verbose == T){
    print(paste("Power of k_ab =", k_ab, "k_ac =", k_ac, "pi_a =", pi_a, "OR_bc =", OR_bc, "tau =", tau, "is", power))
  }
  
  return(result)
}




library(doMC)
library(gemtc)

k_ab = k_ac = 6
k_bc = 6

pi_a = 0.5

OR_ab = 1.2
OR_ac = 1.8

tau = 0.2
S = 10

power.sim.NMA3 <- function(k_ab = 0, k_ac = 0, k_bc = 0, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.4, tau = 0.01){

    # simulate k_ab studies (indirect)
    dat_ab = data.frame(matrix(ncol = 4, nrow = 0))
    if (as.numeric(k_ab) != 0){
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
    
    network_abc <- gemtc::mtc.network(data.ab=all_data)
    
    # fit NMA model
    cons.model <- gemtc::mtc.model(network_abc, type="consistency", likelihood="binom", link="logit", linearModel="random")
    
    cons.out = NULL
    
    cons.out <- gemtc::mtc.run(cons.model, n.adapt=500, n.iter=5000, thin=1)
    
    # Rank order
    prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
    prob <- round(prob, digits=3)
    sucra <- gemtc::sucra(prob)
    rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
    
    
    if (k_ab == 0 & k_ac == 0) {
      res = summary(gemtc::relative.effect(cons.out,"B",c("B","C")))
      lower_95 = res$summaries$quantiles[2,1]
      upper_95 = res$summaries$quantiles[2,5]
      point_est = res$summaries$statistics[2,1]
      true_rank = c("C","B")
    } else {
      res = summary(gemtc::relative.effect(cons.out,"B",c("A","B","C")))
      lower_95 = res$summaries$quantiles[3,1]
      upper_95 = res$summaries$quantiles[3,5]
      point_est = res$summaries$statistics[3,1]
      true_rank = c("C","B","A")
    }
    
    reject_null = as.numeric(0 > upper_95 || 0 < lower_95)
    rank_correct = as.numeric(identical(rank_order, true_rank))
    
    point_true = log(OR_bc)
    
  
  return(paste(reject_null, rank_correct, point_est, point_true))
}



power.sim.NMA3 <- function(k_ab = 0, k_ac = 0, k_bc = 0, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.4, tau = 0.01){
  
  result = foreach (i = 1, .combine = rbind, .errorhandling='remove') %dopar% {
    # simulate k_ab studies (indirect)
    dat_ab = data.frame(matrix(ncol = 4, nrow = 0))
    if (as.numeric(k_ab) != 0){
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
    
    network_abc <- gemtc::mtc.network(data.ab=all_data)
    
    # fit NMA model
    cons.model <- gemtc::mtc.model(network_abc, type="consistency", likelihood="binom", link="logit", linearModel="random")
    
    cons.out = NULL
    
    
    cons.out <- gemtc::mtc.run(cons.model, n.adapt=500, n.iter=5000,thin=1)
    
    
    # Rank order
    prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
    prob <- round(prob, digits=3)
    sucra <- gemtc::sucra(prob)
    rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
    
    
    if (k_ab == 0 & k_ac == 0) {
      res = summary(gemtc::relative.effect(cons.out,"B",c("B","C")))
      lower_95 = res$summaries$quantiles[2,1]
      upper_95 = res$summaries$quantiles[2,5]
      point_est = res$summaries$statistics[2,1]
      true_rank = c("C","B")
    } else {
      res = summary(gemtc::relative.effect(cons.out,"B",c("A","B","C")))
      lower_95 = res$summaries$quantiles[3,1]
      upper_95 = res$summaries$quantiles[3,5]
      point_est = res$summaries$statistics[3,1]
      true_rank = c("C","B","A")
    }
    
    reject_null = as.numeric(0 > upper_95 || 0 < lower_95)
    rank_correct = as.numeric(identical(rank_order, true_rank))
    
    point_true = log(OR_bc)
    c(reject_null, rank_correct, point_est, point_true)
    
  }
  
  if (is.null(result))
    {result = c(0,0,0,0)}
  
  return(paste(result[1], result[2], result[3], result[4]))
}


power.sim.NMA3 <- function(k_ab = 0, k_ac = 0, k_bc = 0, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.4, tau = 0.01){
  
  # simulate k_ab studies (indirect)
  dat_ab = data.frame(matrix(ncol = 4, nrow = 0))
  if (as.numeric(k_ab) != 0){
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
  
  network_abc <- gemtc::mtc.network(data.ab=all_data)
  
  # fit NMA model
  cons.model <- gemtc::mtc.model(network_abc, type="consistency", likelihood="binom", link="logit", linearModel="random")
  
  cons.out = NULL
  
  cons.out <- gemtc::mtc.run(cons.model, n.adapt=500, n.iter=5000, thin=1)
  
  # Rank order
  prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
  prob <- round(prob, digits=3)
  sucra <- gemtc::sucra(prob)
  rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
  
  
  if (k_ab == 0 & k_ac == 0) {
    res = summary(gemtc::relative.effect(cons.out,"B",c("B","C")))
    lower_95 = res$summaries$quantiles[2,1]
    upper_95 = res$summaries$quantiles[2,5]
    point_est = res$summaries$statistics[2,1]
    true_rank = c("C","B")
  } else {
    res = summary(gemtc::relative.effect(cons.out,"B",c("A","B","C")))
    lower_95 = res$summaries$quantiles[3,1]
    upper_95 = res$summaries$quantiles[3,5]
    point_est = res$summaries$statistics[3,1]
    true_rank = c("C","B","A")
  }
  
  reject_null = as.numeric(0 > upper_95 || 0 < lower_95)
  rank_correct = as.numeric(identical(rank_order, true_rank))
  
  point_true = log(OR_bc)
  
  
  return(paste(reject_null, rank_correct, point_est, point_true))
}



power.sim.NMA4 <- function(S = 3, k_ab = 0, k_ac = 0, k_bc = 0, pi_a = 0.5, OR_ab = 1.2, OR_ac = 1.4, tau = 0.01){
  
  result = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
    # simulate k_ab studies (indirect)
    dat_ab = data.frame(matrix(ncol = 4, nrow = 0))
    if (as.numeric(k_ab) != 0){
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
    
    network_abc <- gemtc::mtc.network(data.ab=all_data)
    
    # fit NMA model
    cons.model <- gemtc::mtc.model(network_abc, type="consistency", likelihood="binom", link="logit", linearModel="random")
    
    cons.out = NULL
    
    
    cons.out <- gemtc::mtc.run(cons.model, n.adapt=500, n.iter=2000,thin=1)
    
    
    # Rank order
    prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
    prob <- round(prob, digits=3)
    sucra <- gemtc::sucra(prob)
    rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
    
    
    if (k_ab == 0 & k_ac == 0) {
      res = summary(gemtc::relative.effect(cons.out,"B",c("B","C")))
      lower_95 = res$summaries$quantiles[2,1]
      upper_95 = res$summaries$quantiles[2,5]
      point_est = res$summaries$statistics[2,1]
      true_rank = c("C","B")
    } else {
      res = summary(gemtc::relative.effect(cons.out,"B",c("A","B","C")))
      lower_95 = res$summaries$quantiles[3,1]
      upper_95 = res$summaries$quantiles[3,5]
      point_est = res$summaries$statistics[3,1]
      true_rank = c("C","B","A")
    }
    
    reject_null = as.numeric(0 > upper_95 || 0 < lower_95)
    rank_correct = as.numeric(identical(rank_order, true_rank))
    
    point_true = log(OR_bc)
    c(reject_null, rank_correct, 
      point_est - point_true, 
      abs(point_est - point_true))
    
  }
  
  # if (is.null(result))
  #   {result = c(0,0,0,0)}
  
  
  power = result[1] / S 
  rank_correct_prob = result[2] / S
  avg_bias = result[3] / S
  avg_bias_abs = result[4] / S
  
  return(paste(power,rank_correct_prob, avg_bias, avg_bias_abs))
}


# test
# power.sim.NMA1(S = 100, k_ab = 6, k_ac = 6, k_bc = 3, pi_a = 0.1, OR_ab = 1.2, OR_ac = 1.4, tau = 0.1)


# test
# power.sim.NMA2(S = 100, k_ab = 6, k_ac = 6, k_bc = 3, pi_a = 0.1, OR_ab = 1.2, OR_ac = 1.4, tau = 0.1)

# power.sim.NMA4(k_ab = 6, k_ac = 6, k_bc = 3, pi_a = 0.1, OR_ab = 1.2, OR_ac = 1.4, tau = 0.1)

