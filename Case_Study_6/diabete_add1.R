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
  
  res = summary(gemtc::relative.effect(cons.out,"1", c("2","3","4","5","6")))
  
  AB_lower_95 = res$summaries$quantiles[1,1]
  AB_upper_95 = res$summaries$quantiles[1,5]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)
  
  AC_lower_95 = res$summaries$quantiles[2,1]
  AC_upper_95 = res$summaries$quantiles[2,5]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)
  
  AD_lower_95 = res$summaries$quantiles[3,1]
  AD_upper_95 = res$summaries$quantiles[3,5]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)
  
  AE_lower_95 = res$summaries$quantiles[4,1]
  AE_upper_95 = res$summaries$quantiles[4,5]
  AE_reject = 1-as.numeric(0 >= AE_lower_95 & 0 <= AE_upper_95)
  
  AF_lower_95 = res$summaries$quantiles[5,1]
  AF_upper_95 = res$summaries$quantiles[5,5]
  AF_reject = 1-as.numeric(0 >= AF_lower_95 & 0 <= AF_upper_95)
  
  
  reject_correct = c(AB_reject, AC_reject, AD_reject, AE_reject, AF_reject)
  
  return(c(reject_correct))
}

# first three powers / next three Type I errors / last one A being ranked last
result_diabetes = result_diabetes/S

result_diabetes = matrix(result_diabetes, nrow = 1)
colnames(result_diabetes) <- c("AB_power", "AC_power", "AD_power", "AE_power", "AF_power")

save(result_diabetes, file = "result_diabetes_original.RData")






# add AB
diabetes_ab_AB = rbind(diabetes_ab, c(23, "1", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_a)) )
diabetes_ab_AB = rbind(diabetes_ab_AB, c(23, "2", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_b)) )
diabetes_ab_AB$sampleSize = as.numeric(diabetes_ab_AB$sampleSize)
diabetes_ab_AB$responders = as.numeric(diabetes_ab_AB$responders)
n_study = n_distinct(diabetes_ab_AB$study)

result_diabetes = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  
  # LA model
  Log_OR_dat = c(0, -0.2934780,-0.0712749,-0.2427150, -0.4125047, -0.4914017)
  tau = 0.1345282
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = diabetes_ab_AB[diabetes_ab_AB$study==i,]
    
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
  
  sim_dat = diabetes_ab_AB
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res = summary(gemtc::relative.effect(cons.out,"1", c("2","3","4","5","6")))
  
  AB_lower_95 = res$summaries$quantiles[1,1]
  AB_upper_95 = res$summaries$quantiles[1,5]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)
  
  AC_lower_95 = res$summaries$quantiles[2,1]
  AC_upper_95 = res$summaries$quantiles[2,5]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)
  
  AD_lower_95 = res$summaries$quantiles[3,1]
  AD_upper_95 = res$summaries$quantiles[3,5]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)
  
  AE_lower_95 = res$summaries$quantiles[4,1]
  AE_upper_95 = res$summaries$quantiles[4,5]
  AE_reject = 1-as.numeric(0 >= AE_lower_95 & 0 <= AE_upper_95)
  
  AF_lower_95 = res$summaries$quantiles[5,1]
  AF_upper_95 = res$summaries$quantiles[5,5]
  AF_reject = 1-as.numeric(0 >= AF_lower_95 & 0 <= AF_upper_95)
  
  
  reject_correct = c(AB_reject, AC_reject, AD_reject, AE_reject, AF_reject)
  
  return(c(reject_correct))
}

# first three powers / next three Type I errors / last one A being ranked last
result_diabetes = result_diabetes/S

result_diabetes = matrix(result_diabetes, nrow = 1)
colnames(result_diabetes) <- c("AB_power", "AC_power", "AD_power", "AE_power", "AF_power")

save(result_diabetes, file = "result_diabetes_AB.RData")

# stopCluster(cl)




# add AC
diabetes_ab_AC = rbind(diabetes_ab, c(23, "1", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_a)) )
diabetes_ab_AC = rbind(diabetes_ab_AC, c(23, "3", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_c)) )
diabetes_ab_AC$sampleSize = as.numeric(diabetes_ab_AC$sampleSize)
diabetes_ab_AC$responders = as.numeric(diabetes_ab_AC$responders)
n_study = n_distinct(diabetes_ab_AC$study)

result_diabetes = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  
  # LA model
  Log_OR_dat = c(0, -0.2934780,-0.0712749,-0.2427150, -0.4125047, -0.4914017)
  tau = 0.1345282
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = diabetes_ab_AC[diabetes_ab_AC$study==i,]
    
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
  
  sim_dat = diabetes_ab_AC
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res = summary(gemtc::relative.effect(cons.out,"1", c("2","3","4","5","6")))
  
  AB_lower_95 = res$summaries$quantiles[1,1]
  AB_upper_95 = res$summaries$quantiles[1,5]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)
  
  AC_lower_95 = res$summaries$quantiles[2,1]
  AC_upper_95 = res$summaries$quantiles[2,5]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)
  
  AD_lower_95 = res$summaries$quantiles[3,1]
  AD_upper_95 = res$summaries$quantiles[3,5]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)
  
  AE_lower_95 = res$summaries$quantiles[4,1]
  AE_upper_95 = res$summaries$quantiles[4,5]
  AE_reject = 1-as.numeric(0 >= AE_lower_95 & 0 <= AE_upper_95)
  
  AF_lower_95 = res$summaries$quantiles[5,1]
  AF_upper_95 = res$summaries$quantiles[5,5]
  AF_reject = 1-as.numeric(0 >= AF_lower_95 & 0 <= AF_upper_95)
  
  
  reject_correct = c(AB_reject, AC_reject, AD_reject, AE_reject, AF_reject)
  
  return(c(reject_correct))
}

# first three powers / next three Type I errors / last one A being ranked last
result_diabetes = result_diabetes/S

result_diabetes = matrix(result_diabetes, nrow = 1)
colnames(result_diabetes) <- c("AB_power", "AC_power", "AD_power", "AE_power", "AF_power")

save(result_diabetes, file = "result_diabetes_AC.RData")

# stopCluster(cl)



# add AC
diabetes_ab_AD = rbind(diabetes_ab, c(23, "1", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_a)) )
diabetes_ab_AD = rbind(diabetes_ab_AD, c(23, "4", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_d)) )
diabetes_ab_AD$sampleSize = as.numeric(diabetes_ab_AD$sampleSize)
diabetes_ab_AD$responders = as.numeric(diabetes_ab_AD$responders)
n_study = n_distinct(diabetes_ab_AD$study)

result_diabetes = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  
  # LA model
  Log_OR_dat = c(0, -0.2934780,-0.0712749,-0.2427150, -0.4125047, -0.4914017)
  tau = 0.1345282
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = diabetes_ab_AD[diabetes_ab_AD$study==i,]
    
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
  
  sim_dat = diabetes_ab_AD
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res = summary(gemtc::relative.effect(cons.out,"1", c("2","3","4","5","6")))
  
  AB_lower_95 = res$summaries$quantiles[1,1]
  AB_upper_95 = res$summaries$quantiles[1,5]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)
  
  AC_lower_95 = res$summaries$quantiles[2,1]
  AC_upper_95 = res$summaries$quantiles[2,5]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)
  
  AD_lower_95 = res$summaries$quantiles[3,1]
  AD_upper_95 = res$summaries$quantiles[3,5]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)
  
  AE_lower_95 = res$summaries$quantiles[4,1]
  AE_upper_95 = res$summaries$quantiles[4,5]
  AE_reject = 1-as.numeric(0 >= AE_lower_95 & 0 <= AE_upper_95)
  
  AF_lower_95 = res$summaries$quantiles[5,1]
  AF_upper_95 = res$summaries$quantiles[5,5]
  AF_reject = 1-as.numeric(0 >= AF_lower_95 & 0 <= AF_upper_95)
  
  
  reject_correct = c(AB_reject, AC_reject, AD_reject, AE_reject, AF_reject)
  
  return(c(reject_correct))
}

# first three powers / next three Type I errors / last one A being ranked last
result_diabetes = result_diabetes/S

result_diabetes = matrix(result_diabetes, nrow = 1)
colnames(result_diabetes) <- c("AB_power", "AC_power", "AD_power", "AE_power", "AF_power")

save(result_diabetes, file = "result_diabetes_AD.RData")

# stopCluster(cl)




# add AE
diabetes_ab_AE = rbind(diabetes_ab, c(23, "1", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_a)) )
diabetes_ab_AE = rbind(diabetes_ab_AE, c(23, "5", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_e)) )
diabetes_ab_AE$sampleSize = as.numeric(diabetes_ab_AE$sampleSize)
diabetes_ab_AE$responders = as.numeric(diabetes_ab_AE$responders)
n_study = n_distinct(diabetes_ab_AE$study)

result_diabetes = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  
  # LA model
  Log_OR_dat = c(0, -0.2934780,-0.0712749,-0.2427150, -0.4125047, -0.4914017)
  tau = 0.1345282
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = diabetes_ab_AE[diabetes_ab_AE$study==i,]
    
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
  
  sim_dat = diabetes_ab_AE
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res = summary(gemtc::relative.effect(cons.out,"1", c("2","3","4","5","6")))
  
  AB_lower_95 = res$summaries$quantiles[1,1]
  AB_upper_95 = res$summaries$quantiles[1,5]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)
  
  AC_lower_95 = res$summaries$quantiles[2,1]
  AC_upper_95 = res$summaries$quantiles[2,5]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)
  
  AD_lower_95 = res$summaries$quantiles[3,1]
  AD_upper_95 = res$summaries$quantiles[3,5]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)
  
  AE_lower_95 = res$summaries$quantiles[4,1]
  AE_upper_95 = res$summaries$quantiles[4,5]
  AE_reject = 1-as.numeric(0 >= AE_lower_95 & 0 <= AE_upper_95)
  
  AF_lower_95 = res$summaries$quantiles[5,1]
  AF_upper_95 = res$summaries$quantiles[5,5]
  AF_reject = 1-as.numeric(0 >= AF_lower_95 & 0 <= AF_upper_95)
  
  
  reject_correct = c(AB_reject, AC_reject, AD_reject, AE_reject, AF_reject)
  
  return(c(reject_correct))
}

# first three powers / next three Type I errors / last one A being ranked last
result_diabetes = result_diabetes/S

result_diabetes = matrix(result_diabetes, nrow = 1)
colnames(result_diabetes) <- c("AB_power", "AC_power", "AD_power", "AE_power", "AF_power")

save(result_diabetes, file = "result_diabetes_AE.RData")

# stopCluster(cl)





# add AF
diabetes_ab_AF = rbind(diabetes_ab, c(23, "1", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_a)) )
diabetes_ab_AF = rbind(diabetes_ab_AF, c(23, "6", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_f)) )
diabetes_ab_AF$sampleSize = as.numeric(diabetes_ab_AF$sampleSize)
diabetes_ab_AF$responders = as.numeric(diabetes_ab_AF$responders)
n_study = n_distinct(diabetes_ab_AF$study)

result_diabetes = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  
  # LA model
  Log_OR_dat = c(0, -0.2934780,-0.0712749,-0.2427150, -0.4125047, -0.4914017)
  tau = 0.1345282
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = diabetes_ab_AF[diabetes_ab_AF$study==i,]
    
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
  
  sim_dat = diabetes_ab_AF
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res = summary(gemtc::relative.effect(cons.out,"1", c("2","3","4","5","6")))
  
  AB_lower_95 = res$summaries$quantiles[1,1]
  AB_upper_95 = res$summaries$quantiles[1,5]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)
  
  AC_lower_95 = res$summaries$quantiles[2,1]
  AC_upper_95 = res$summaries$quantiles[2,5]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)
  
  AD_lower_95 = res$summaries$quantiles[3,1]
  AD_upper_95 = res$summaries$quantiles[3,5]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)
  
  AE_lower_95 = res$summaries$quantiles[4,1]
  AE_upper_95 = res$summaries$quantiles[4,5]
  AE_reject = 1-as.numeric(0 >= AE_lower_95 & 0 <= AE_upper_95)
  
  AF_lower_95 = res$summaries$quantiles[5,1]
  AF_upper_95 = res$summaries$quantiles[5,5]
  AF_reject = 1-as.numeric(0 >= AF_lower_95 & 0 <= AF_upper_95)
  
  
  reject_correct = c(AB_reject, AC_reject, AD_reject, AE_reject, AF_reject)
  
  return(c(reject_correct))
}

# first three powers / next three Type I errors / last one A being ranked last
result_diabetes = result_diabetes/S

result_diabetes = matrix(result_diabetes, nrow = 1)
colnames(result_diabetes) <- c("AB_power", "AC_power", "AD_power", "AE_power", "AF_power")

save(result_diabetes, file = "result_diabetes_AF.RData")





# add EF
diabetes_ab_EF = rbind(diabetes_ab, c(23, "5", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_e)) )
diabetes_ab_EF = rbind(diabetes_ab_EF, c(23, "6", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_f)) )
diabetes_ab_EF$sampleSize = as.numeric(diabetes_ab_EF$sampleSize)
diabetes_ab_EF$responders = as.numeric(diabetes_ab_EF$responders)
n_study = n_distinct(diabetes_ab_EF$study)

result_diabetes = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  
  # LA model
  Log_OR_dat = c(0, -0.2934780,-0.0712749,-0.2427150, -0.4125047, -0.4914017)
  tau = 0.1345282
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = diabetes_ab_EF[diabetes_ab_EF$study==i,]
    
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
  
  sim_dat = diabetes_ab_EF
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res = summary(gemtc::relative.effect(cons.out,"1", c("2","3","4","5","6")))
  
  AB_lower_95 = res$summaries$quantiles[1,1]
  AB_upper_95 = res$summaries$quantiles[1,5]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)
  
  AC_lower_95 = res$summaries$quantiles[2,1]
  AC_upper_95 = res$summaries$quantiles[2,5]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)
  
  AD_lower_95 = res$summaries$quantiles[3,1]
  AD_upper_95 = res$summaries$quantiles[3,5]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)
  
  AE_lower_95 = res$summaries$quantiles[4,1]
  AE_upper_95 = res$summaries$quantiles[4,5]
  AE_reject = 1-as.numeric(0 >= AE_lower_95 & 0 <= AE_upper_95)
  
  AF_lower_95 = res$summaries$quantiles[5,1]
  AF_upper_95 = res$summaries$quantiles[5,5]
  AF_reject = 1-as.numeric(0 >= AF_lower_95 & 0 <= AF_upper_95)
  
  
  reject_correct = c(AB_reject, AC_reject, AD_reject, AE_reject, AF_reject)
  
  return(c(reject_correct))
}

# first three powers / next three Type I errors / last one A being ranked last
result_diabetes = result_diabetes/S

result_diabetes = matrix(result_diabetes, nrow = 1)
colnames(result_diabetes) <- c("AB_power", "AC_power", "AD_power", "AE_power", "AF_power")

save(result_diabetes, file = "result_diabetes_EF.RData")






# add EF
diabetes_ab_EF2 = rbind(diabetes_ab, c(23, "5", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_e * 0.5)) )
diabetes_ab_EF2 = rbind(diabetes_ab_EF2, c(23, "6", round(mean(diabetes_ab$sampleSize)), round(mean(diabetes_ab$sampleSize) * pi_f * 3)) )
diabetes_ab_EF2$sampleSize = as.numeric(diabetes_ab_EF2$sampleSize)
diabetes_ab_EF2$responders = as.numeric(diabetes_ab_EF2$responders)
n_study = n_distinct(diabetes_ab_EF2$study)

result_diabetes = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  
  # LA model
  Log_OR_dat = c(0, -0.2934780,-0.0712749,-0.2427150, -0.4125047, -0.4914017)
  tau = 0.1345282
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = diabetes_ab_EF2[diabetes_ab_EF2$study==i,]
    
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
  
  sim_dat = diabetes_ab_EF2
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res = summary(gemtc::relative.effect(cons.out,"1", c("2","3","4","5","6")))
  
  AB_lower_95 = res$summaries$quantiles[1,1]
  AB_upper_95 = res$summaries$quantiles[1,5]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)
  
  AC_lower_95 = res$summaries$quantiles[2,1]
  AC_upper_95 = res$summaries$quantiles[2,5]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)
  
  AD_lower_95 = res$summaries$quantiles[3,1]
  AD_upper_95 = res$summaries$quantiles[3,5]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)
  
  AE_lower_95 = res$summaries$quantiles[4,1]
  AE_upper_95 = res$summaries$quantiles[4,5]
  AE_reject = 1-as.numeric(0 >= AE_lower_95 & 0 <= AE_upper_95)
  
  AF_lower_95 = res$summaries$quantiles[5,1]
  AF_upper_95 = res$summaries$quantiles[5,5]
  AF_reject = 1-as.numeric(0 >= AF_lower_95 & 0 <= AF_upper_95)
  
  
  reject_correct = c(AB_reject, AC_reject, AD_reject, AE_reject, AF_reject)
  
  return(c(reject_correct))
}

# first three powers / next three Type I errors / last one A being ranked last
result_diabetes = result_diabetes/S

result_diabetes = matrix(result_diabetes, nrow = 1)
colnames(result_diabetes) <- c("AB_power", "AC_power", "AD_power", "AE_power", "AF_power")

save(result_diabetes, file = "result_diabetes_EF2.RData")

stopCluster(cl)





