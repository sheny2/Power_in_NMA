library(foreach)
library(doParallel)
# Initialize parallel backend
cl <- makeCluster(detectCores())

# Register the parallel backend
registerDoParallel(cl)

S = 1000

library(tidyverse)
library(gemtc)

# source("ifplot.fun.R")
# source("BayesDiagnos.fun.R")
# source("ranko_sucra.R")

set.seed(2023)

library(pcnetmeta)
data(diabetes)
diabetes_ab = data.frame(study = factor(diabetes$s.id), treatment = factor(diabetes$t.id), 
                         sampleSize = diabetes$n, responders = diabetes$r)

diabetes_ab$study = as.numeric(diabetes_ab$study)
trt = unique(diabetes_ab$treatment)

# Helper function

calculate_pi <- function(odds_ratio, p_A) {
  odds_A <- p_A / (1 - p_A)
  odds_B <- odds_A * odds_ratio
  p_B <- odds_B / (1 + odds_B)
  
  return(p_B)
}


logit <- function(p) {
  return(log(p/(1-p)))
}


inverse_logit <- function(x) {
  odds <- exp(x)  
  probability <- odds / (1 + odds) 
  return(probability)  
}




# Generate Data (Lu & Ades model)

pi_a = sum(diabetes_ab[diabetes_ab$treatment=="1",]$responders) / 
  sum(diabetes_ab[diabetes_ab$treatment=="1",]$sampleSize)

# True LOR
Log_OR_dat = c(0, -0.2934780,-0.0712749,-0.2427150, -0.4125047, -0.4914017)
tau = 0.1345282

OR_ab = exp(Log_OR_dat[2])
OR_ac = exp(Log_OR_dat[3])
OR_ad = exp(Log_OR_dat[4])
OR_ae = exp(Log_OR_dat[5])
OR_af = exp(Log_OR_dat[6])

pi_b =  calculate_pi(OR_ab, pi_a)
pi_c =  calculate_pi(OR_ac, pi_a)
pi_d =  calculate_pi(OR_ad, pi_a)
pi_e =  calculate_pi(OR_ae, pi_a)
pi_f =  calculate_pi(OR_af, pi_a)


result_diabetes = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  
  # LA model
  Log_OR_dat = c(0, -0.2934780,-0.0712749,-0.2427150, -0.4125047, -0.4914017)
  tau = 0.1345282
  n_study = n_distinct(diabetes_ab$study)
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = diabetes_ab[diabetes_ab$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    if (study_data$treatment[1] == "1") {
      for (k in 1:nrow(study_data)){
        if (k == 1) {
          logit_p_ik = alpha_iB 
          p_ik = inverse_logit(logit_p_ik)
          y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
        } else {
          d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
          delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
          logit_p_ik = alpha_iB + delta_iBk
          p_ik = inverse_logit(logit_p_ik)
          y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
        }
      }
    } else {for (k in 1:nrow(study_data)){
      d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
      delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
      logit_p_ik = alpha_iB + delta_iBk
      p_ik = inverse_logit(logit_p_ik)
      y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
    }
    }
  }
  
  sim_dat = diabetes_ab
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res56 = summary(gemtc::relative.effect(cons.out,"5", "6"))
  res23 = summary(gemtc::relative.effect(cons.out,"2", "3")) 
  res14 = summary(gemtc::relative.effect(cons.out,"1", "4")) 
  res34 = summary(gemtc::relative.effect(cons.out,"3", "4")) 
  
  lower_56 = res56$summaries$quantiles[1,1]
  upper_56 = res56$summaries$quantiles[1,5]
  reject_56 = 1-as.numeric(0 >= lower_56 & 0 <= upper_56)
  
  lower_23 = res23$summaries$quantiles[1,1]
  upper_23 = res23$summaries$quantiles[1,5]
  reject_23 = 1-as.numeric(0 >= lower_23 & 0 <= upper_23)
  
  lower_14 = res14$summaries$quantiles[1,1]
  upper_14 = res14$summaries$quantiles[1,5]
  reject_14 = 1-as.numeric(0 >= lower_14 & 0 <= upper_14)
  
  lower_34 = res34$summaries$quantiles[1,1]
  upper_34 = res34$summaries$quantiles[1,5]
  reject_34 = 1-as.numeric(0 >= lower_34 & 0 <= upper_34)
  
  reject_correct = c(reject_56, reject_23, reject_14, reject_34)
  
  return(c(reject_correct))
}

result_diabetes = result_diabetes/S

result_diabetes = matrix(result_diabetes, nrow = 1)
colnames(result_diabetes) <- c("Power.5-6","Power.2-3","Power.1-4","Power.3-4")

save(result_diabetes, file = "result_diabetes_4power_original.RData")






# add 56
diabetes_ab_56 = rbind(diabetes_ab, c(23, "5", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_e)) )
diabetes_ab_56 = rbind(diabetes_ab_56, c(23, "6", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_f)) )
diabetes_ab_56 = rbind(diabetes_ab_56, c(24, "5", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_e)) )
diabetes_ab_56 = rbind(diabetes_ab_56, c(24, "6", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_f)) )
diabetes_ab_56 = rbind(diabetes_ab_56, c(25, "5", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_e)) )
diabetes_ab_56 = rbind(diabetes_ab_56, c(25, "6", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_f)) )
diabetes_ab_56 = rbind(diabetes_ab_56, c(26, "5", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_e)) )
diabetes_ab_56 = rbind(diabetes_ab_56, c(26, "6", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_f)) )
diabetes_ab_56 = rbind(diabetes_ab_56, c(27, "5", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_e)) )
diabetes_ab_56 = rbind(diabetes_ab_56, c(27, "6", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_f)) )

diabetes_ab_56$sampleSize = as.numeric(diabetes_ab_56$sampleSize)
diabetes_ab_56$responders = as.numeric(diabetes_ab_56$responders)
n_study = n_distinct(diabetes_ab_56$study)

result_diabetes = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  
  # LA model
  Log_OR_dat = c(0, -0.2934780,-0.0712749,-0.2427150, -0.4125047, -0.4914017)
  tau = 0.1345282
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = diabetes_ab_56[diabetes_ab_56$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    if (study_data$treatment[1] == "1") {
      for (k in 1:nrow(study_data)){
        if (k == 1) {
          logit_p_ik = alpha_iB 
          p_ik = inverse_logit(logit_p_ik)
          y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
        } else {
          d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
          delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
          logit_p_ik = alpha_iB + delta_iBk
          p_ik = inverse_logit(logit_p_ik)
          y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
        }
      }
    } else {for (k in 1:nrow(study_data)){
      d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
      delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
      logit_p_ik = alpha_iB + delta_iBk
      p_ik = inverse_logit(logit_p_ik)
      y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
    }
    }
  }
  
  sim_dat = diabetes_ab_56
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res56 = summary(gemtc::relative.effect(cons.out,"5", "6"))
  res23 = summary(gemtc::relative.effect(cons.out,"2", "3")) 
  res14 = summary(gemtc::relative.effect(cons.out,"1", "4")) 
  res34 = summary(gemtc::relative.effect(cons.out,"3", "4")) 
  
  lower_56 = res56$summaries$quantiles[1,1]
  upper_56 = res56$summaries$quantiles[1,5]
  reject_56 = 1-as.numeric(0 >= lower_56 & 0 <= upper_56)
  
  lower_23 = res23$summaries$quantiles[1,1]
  upper_23 = res23$summaries$quantiles[1,5]
  reject_23 = 1-as.numeric(0 >= lower_23 & 0 <= upper_23)
  
  lower_14 = res14$summaries$quantiles[1,1]
  upper_14 = res14$summaries$quantiles[1,5]
  reject_14 = 1-as.numeric(0 >= lower_14 & 0 <= upper_14)
  
  lower_34 = res34$summaries$quantiles[1,1]
  upper_34 = res34$summaries$quantiles[1,5]
  reject_34 = 1-as.numeric(0 >= lower_34 & 0 <= upper_34)
  
  reject_correct = c(reject_56, reject_23, reject_14, reject_34)
  
  return(c(reject_correct))
}

result_diabetes = result_diabetes/S

result_diabetes = matrix(result_diabetes, nrow = 1)
colnames(result_diabetes) <- c("Power.5-6","Power.2-3","Power.1-4","Power.3-4")

save(result_diabetes, file = "result_diabetes_4power_add56five.RData")

# stopCluster(cl)




# add 23
diabetes_ab_23 = rbind(diabetes_ab, c(23, "2", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_b)) )
diabetes_ab_23 = rbind(diabetes_ab_23, c(23, "3", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_c)) )
diabetes_ab_23 = rbind(diabetes_ab_23, c(24, "2", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_b)) )
diabetes_ab_23 = rbind(diabetes_ab_23, c(24, "3", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_c)) )
diabetes_ab_23 = rbind(diabetes_ab_23, c(25, "2", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_b)) )
diabetes_ab_23 = rbind(diabetes_ab_23, c(25, "3", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_c)) )
diabetes_ab_23 = rbind(diabetes_ab_23, c(26, "2", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_b)) )
diabetes_ab_23 = rbind(diabetes_ab_23, c(26, "3", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_c)) )
diabetes_ab_23 = rbind(diabetes_ab_23, c(27, "2", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_b)) )
diabetes_ab_23 = rbind(diabetes_ab_23, c(27, "3", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_c)) )

diabetes_ab_23$sampleSize = as.numeric(diabetes_ab_23$sampleSize)
diabetes_ab_23$responders = as.numeric(diabetes_ab_23$responders)
n_study = n_distinct(diabetes_ab_23$study)

result_diabetes = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  
  # LA model
  Log_OR_dat = c(0, -0.2934780,-0.0712749,-0.2427150, -0.4125047, -0.4914017)
  tau = 0.1345282
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = diabetes_ab_23[diabetes_ab_23$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    if (study_data$treatment[1] == "1") {
      for (k in 1:nrow(study_data)){
        if (k == 1) {
          logit_p_ik = alpha_iB 
          p_ik = inverse_logit(logit_p_ik)
          y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
        } else {
          d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
          delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
          logit_p_ik = alpha_iB + delta_iBk
          p_ik = inverse_logit(logit_p_ik)
          y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
        }
      }
    } else {for (k in 1:nrow(study_data)){
      d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
      delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
      logit_p_ik = alpha_iB + delta_iBk
      p_ik = inverse_logit(logit_p_ik)
      y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
    }
    }
  }
  
  sim_dat = diabetes_ab_23
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res56 = summary(gemtc::relative.effect(cons.out,"5", "6"))
  res23 = summary(gemtc::relative.effect(cons.out,"2", "3")) 
  res14 = summary(gemtc::relative.effect(cons.out,"1", "4")) 
  res34 = summary(gemtc::relative.effect(cons.out,"3", "4")) 
  
  lower_56 = res56$summaries$quantiles[1,1]
  upper_56 = res56$summaries$quantiles[1,5]
  reject_56 = 1-as.numeric(0 >= lower_56 & 0 <= upper_56)
  
  lower_23 = res23$summaries$quantiles[1,1]
  upper_23 = res23$summaries$quantiles[1,5]
  reject_23 = 1-as.numeric(0 >= lower_23 & 0 <= upper_23)
  
  lower_14 = res14$summaries$quantiles[1,1]
  upper_14 = res14$summaries$quantiles[1,5]
  reject_14 = 1-as.numeric(0 >= lower_14 & 0 <= upper_14)
  
  lower_34 = res34$summaries$quantiles[1,1]
  upper_34 = res34$summaries$quantiles[1,5]
  reject_34 = 1-as.numeric(0 >= lower_34 & 0 <= upper_34)
  
  reject_correct = c(reject_56, reject_23, reject_14, reject_34)
  
  return(c(reject_correct))
}

result_diabetes = result_diabetes/S

result_diabetes = matrix(result_diabetes, nrow = 1)
colnames(result_diabetes) <- c("Power.5-6","Power.2-3","Power.1-4","Power.3-4")

save(result_diabetes, file = "result_diabetes_4power_add23five.RData")

# stopCluster(cl)



# add 14
diabetes_ab_14 = rbind(diabetes_ab, c(23, "1", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_a)) )
diabetes_ab_14 = rbind(diabetes_ab_14, c(23, "4", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_d)) )
diabetes_ab_14 = rbind(diabetes_ab_14, c(24, "1", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_a)) )
diabetes_ab_14 = rbind(diabetes_ab_14, c(24, "4", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_d)) )
diabetes_ab_14 = rbind(diabetes_ab_14, c(25, "1", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_a)) )
diabetes_ab_14 = rbind(diabetes_ab_14, c(25, "4", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_d)) )
diabetes_ab_14 = rbind(diabetes_ab_14, c(26, "1", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_a)) )
diabetes_ab_14 = rbind(diabetes_ab_14, c(26, "4", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_d)) )
diabetes_ab_14 = rbind(diabetes_ab_14, c(27, "1", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_a)) )
diabetes_ab_14 = rbind(diabetes_ab_14, c(27, "4", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_d)) )


diabetes_ab_14$sampleSize = as.numeric(diabetes_ab_14$sampleSize)
diabetes_ab_14$responders = as.numeric(diabetes_ab_14$responders)
n_study = n_distinct(diabetes_ab_14$study)

result_diabetes = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  
  # LA model
  Log_OR_dat = c(0, -0.2934780,-0.0712749,-0.2427150, -0.4125047, -0.4914017)
  tau = 0.1345282
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = diabetes_ab_14[diabetes_ab_14$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    if (study_data$treatment[1] == "1") {
      for (k in 1:nrow(study_data)){
        if (k == 1) {
          logit_p_ik = alpha_iB 
          p_ik = inverse_logit(logit_p_ik)
          y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
        } else {
          d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
          delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
          logit_p_ik = alpha_iB + delta_iBk
          p_ik = inverse_logit(logit_p_ik)
          y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
        }
      }
    } else {for (k in 1:nrow(study_data)){
      d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
      delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
      logit_p_ik = alpha_iB + delta_iBk
      p_ik = inverse_logit(logit_p_ik)
      y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
    }
    }
  }
  
  sim_dat = diabetes_ab_14
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res56 = summary(gemtc::relative.effect(cons.out,"5", "6"))
  res23 = summary(gemtc::relative.effect(cons.out,"2", "3")) 
  res14 = summary(gemtc::relative.effect(cons.out,"1", "4")) 
  res34 = summary(gemtc::relative.effect(cons.out,"3", "4")) 
  
  lower_56 = res56$summaries$quantiles[1,1]
  upper_56 = res56$summaries$quantiles[1,5]
  reject_56 = 1-as.numeric(0 >= lower_56 & 0 <= upper_56)
  
  lower_23 = res23$summaries$quantiles[1,1]
  upper_23 = res23$summaries$quantiles[1,5]
  reject_23 = 1-as.numeric(0 >= lower_23 & 0 <= upper_23)
  
  lower_14 = res14$summaries$quantiles[1,1]
  upper_14 = res14$summaries$quantiles[1,5]
  reject_14 = 1-as.numeric(0 >= lower_14 & 0 <= upper_14)
  
  lower_34 = res34$summaries$quantiles[1,1]
  upper_34 = res34$summaries$quantiles[1,5]
  reject_34 = 1-as.numeric(0 >= lower_34 & 0 <= upper_34)
  
  reject_correct = c(reject_56, reject_23, reject_14, reject_34)
  
  return(c(reject_correct))
}

result_diabetes = result_diabetes/S

result_diabetes = matrix(result_diabetes, nrow = 1)
colnames(result_diabetes) <- c("Power.5-6","Power.2-3","Power.1-4","Power.3-4")

save(result_diabetes, file = "result_diabetes_4power_add14five.RData")

# stopCluster(cl)




# add 34
diabetes_ab_34 = rbind(diabetes_ab, c(23, "3", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_c)) )
diabetes_ab_34 = rbind(diabetes_ab_34, c(23, "4", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_d)) )
diabetes_ab_34 = rbind(diabetes_ab_34, c(24, "3", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_c)) )
diabetes_ab_34 = rbind(diabetes_ab_34, c(24, "4", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_d)) )
diabetes_ab_34 = rbind(diabetes_ab_34, c(25, "3", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_c)) )
diabetes_ab_34 = rbind(diabetes_ab_34, c(25, "4", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_d)) )
diabetes_ab_34 = rbind(diabetes_ab_34, c(26, "3", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_c)) )
diabetes_ab_34 = rbind(diabetes_ab_34, c(26, "4", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_d)) )
diabetes_ab_34 = rbind(diabetes_ab_34, c(27, "3", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_c)) )
diabetes_ab_34 = rbind(diabetes_ab_34, c(27, "4", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_d)) )

diabetes_ab_34$sampleSize = as.numeric(diabetes_ab_34$sampleSize)
diabetes_ab_34$responders = as.numeric(diabetes_ab_34$responders)
n_study = n_distinct(diabetes_ab_34$study)

result_diabetes = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  
  # LA model
  Log_OR_dat = c(0, -0.2934780,-0.0712749,-0.2427150, -0.4125047, -0.4914017)
  tau = 0.1345282
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = diabetes_ab_34[diabetes_ab_34$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    if (study_data$treatment[1] == "1") {
      for (k in 1:nrow(study_data)){
        if (k == 1) {
          logit_p_ik = alpha_iB 
          p_ik = inverse_logit(logit_p_ik)
          y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
        } else {
          d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
          delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
          logit_p_ik = alpha_iB + delta_iBk
          p_ik = inverse_logit(logit_p_ik)
          y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
        }
      }
    } else {for (k in 1:nrow(study_data)){
      d_k = Log_OR_dat[which(trt == study_data$treatment[k])]
      delta_iBk = rnorm(n=1, mean=d_k, sd=tau)
      logit_p_ik = alpha_iB + delta_iBk
      p_ik = inverse_logit(logit_p_ik)
      y_ik = c(y_ik, rbinom(n = 1, size = study_data$sampleSize[k], prob = p_ik))
    }
    }
  }
  
  sim_dat = diabetes_ab_34
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res56 = summary(gemtc::relative.effect(cons.out,"5", "6"))
  res23 = summary(gemtc::relative.effect(cons.out,"2", "3")) 
  res14 = summary(gemtc::relative.effect(cons.out,"1", "4")) 
  res34 = summary(gemtc::relative.effect(cons.out,"3", "4")) 
  
  lower_56 = res56$summaries$quantiles[1,1]
  upper_56 = res56$summaries$quantiles[1,5]
  reject_56 = 1-as.numeric(0 >= lower_56 & 0 <= upper_56)
  
  lower_23 = res23$summaries$quantiles[1,1]
  upper_23 = res23$summaries$quantiles[1,5]
  reject_23 = 1-as.numeric(0 >= lower_23 & 0 <= upper_23)
  
  lower_14 = res14$summaries$quantiles[1,1]
  upper_14 = res14$summaries$quantiles[1,5]
  reject_14 = 1-as.numeric(0 >= lower_14 & 0 <= upper_14)
  
  lower_34 = res34$summaries$quantiles[1,1]
  upper_34 = res34$summaries$quantiles[1,5]
  reject_34 = 1-as.numeric(0 >= lower_34 & 0 <= upper_34)
  
  reject_correct = c(reject_56, reject_23, reject_14, reject_34)
  
  return(c(reject_correct))
}

result_diabetes = result_diabetes/S

result_diabetes = matrix(result_diabetes, nrow = 1)
colnames(result_diabetes) <- c("Power.5-6","Power.2-3","Power.1-4","Power.3-4")

save(result_diabetes, file = "result_diabetes_4power_add34five.RData")
stopCluster(cl)



