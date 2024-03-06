library(foreach)
library(doParallel)
# Initialize parallel backend
cl <- makeCluster(22)

# Register the parallel backend
registerDoParallel(cl)

S = 1000

library(tidyverse)
library(gemtc)

# source("ifplot.fun.R")
# source("BayesDiagnos.fun.R")
# source("ranko_sucra.R")

set.seed(2023)

######## Import data
load("smokingcessation_ab.RData")

smokingcessation_ab$study = as.numeric(smokingcessation_ab$study)

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

pi_a = sum(smokingcessation_ab[smokingcessation_ab$treatment=="A",]$responders) / 
  sum(smokingcessation_ab[smokingcessation_ab$treatment=="A",]$sampleSize)

# True LOR
Log_OR_dat = c(0, 0.4944, 0.8378, 1.0983)
tau = 0.8458

OR_ab = exp(Log_OR_dat[2])
OR_ac = exp(Log_OR_dat[3])
OR_ad = exp(Log_OR_dat[4])

pi_b =  calculate_pi(OR_ab, pi_a)
pi_c =  calculate_pi(OR_ac, pi_a)
pi_d =  calculate_pi(OR_ad, pi_a)

trt = c("A", "B", "C", "D")


result_smoke = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)

  Log_OR_dat = c(0, 0.4944, 0.8378, 1.0983)
  tau = 0.8458
  n_study = n_distinct(smokingcessation_ab$study)
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = smokingcessation_ab[smokingcessation_ab$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    if (study_data$treatment[1] == "A") {
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
  
  sim_dat = smokingcessation_ab
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res = summary(gemtc::relative.effect(cons.out,"A", c("B","C","D")))
  
  AB_lower_95 = res$summaries$quantiles[1,1]
  AB_upper_95 = res$summaries$quantiles[1,5]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)
  
  AC_lower_95 = res$summaries$quantiles[2,1]
  AC_upper_95 = res$summaries$quantiles[2,5]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)
  
  AD_lower_95 = res$summaries$quantiles[3,1]
  AD_upper_95 = res$summaries$quantiles[3,5]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)
  
  reject_correct = c(AB_reject, AC_reject, AD_reject)
  
  return(c(reject_correct))
}

result_smoke = result_smoke/S

result_smoke = matrix(result_smoke, nrow = 1)
colnames(result_smoke) <- c("AB_power", "AC_power", "AD_power")

save(result_smoke, file = "result_smoke_original.RData")

# stopCluster(cl)





# AB
smokingcessation_ab_AB = rbind(smokingcessation_ab, c(25, "A", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_a)) )
smokingcessation_ab_AB = rbind(smokingcessation_ab_AB, c(25, "B", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_b)) )
smokingcessation_ab_AB = rbind(smokingcessation_ab_AB, c(26, "A", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_a)) )
smokingcessation_ab_AB = rbind(smokingcessation_ab_AB, c(26, "B", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_b)) )
smokingcessation_ab_AB = rbind(smokingcessation_ab_AB, c(27, "A", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_a)) )
smokingcessation_ab_AB = rbind(smokingcessation_ab_AB, c(27, "B", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_b)) )
smokingcessation_ab_AB = rbind(smokingcessation_ab_AB, c(28, "A", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_a)) )
smokingcessation_ab_AB = rbind(smokingcessation_ab_AB, c(28, "B", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_b)) )
smokingcessation_ab_AB = rbind(smokingcessation_ab_AB, c(29, "A", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_a)) )
smokingcessation_ab_AB = rbind(smokingcessation_ab_AB, c(29, "B", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_b)) )

smokingcessation_ab_AB$sampleSize = as.numeric(smokingcessation_ab_AB$sampleSize)
smokingcessation_ab_AB$responders = as.numeric(smokingcessation_ab_AB$responders)
n_study = n_distinct(smokingcessation_ab_AB$study)

result_smoke = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  Log_OR_dat = c(0, 0.4944, 0.8378, 1.0983)
  tau = 0.8458
  n_study = n_distinct(smokingcessation_ab_AB$study)
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = smokingcessation_ab_AB[smokingcessation_ab_AB$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    if (study_data$treatment[1] == "A") {
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
  
  sim_dat = smokingcessation_ab_AB
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res = summary(gemtc::relative.effect(cons.out,"A", c("B","C","D")))
  
  AB_lower_95 = res$summaries$quantiles[1,1]
  AB_upper_95 = res$summaries$quantiles[1,5]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)
  
  AC_lower_95 = res$summaries$quantiles[2,1]
  AC_upper_95 = res$summaries$quantiles[2,5]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)
  
  AD_lower_95 = res$summaries$quantiles[3,1]
  AD_upper_95 = res$summaries$quantiles[3,5]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)
  
  reject_correct = c(AB_reject, AC_reject, AD_reject)
  
  return(c(reject_correct))
}

result_smoke = result_smoke/S

result_smoke = matrix(result_smoke, nrow = 1)
colnames(result_smoke) <- c("AB_power", "AC_power", "AD_power")

save(result_smoke, file = "result_smoke_AB_add5.RData")




# AC
smokingcessation_ab_AC = rbind(smokingcessation_ab, c(25, "A", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_a)) )
smokingcessation_ab_AC = rbind(smokingcessation_ab_AC, c(25, "C", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_c)) )
smokingcessation_ab_AC = rbind(smokingcessation_ab_AC, c(26, "A", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_a)) )
smokingcessation_ab_AC = rbind(smokingcessation_ab_AC, c(26, "C", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_c)) )
smokingcessation_ab_AC = rbind(smokingcessation_ab_AC, c(27, "A", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_a)) )
smokingcessation_ab_AC = rbind(smokingcessation_ab_AC, c(27, "C", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_c)) )
smokingcessation_ab_AC = rbind(smokingcessation_ab_AC, c(28, "A", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_a)) )
smokingcessation_ab_AC = rbind(smokingcessation_ab_AC, c(28, "C", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_c)) )
smokingcessation_ab_AC = rbind(smokingcessation_ab_AC, c(29, "A", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_a)) )
smokingcessation_ab_AC = rbind(smokingcessation_ab_AC, c(29, "C", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_c)) )

smokingcessation_ab_AC$sampleSize = as.numeric(smokingcessation_ab_AC$sampleSize)
smokingcessation_ab_AC$responders = as.numeric(smokingcessation_ab_AC$responders)
n_study = n_distinct(smokingcessation_ab_AC$study)

result_smoke = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  Log_OR_dat = c(0, 0.4944, 0.8378, 1.0983)
  tau = 0.8458
  n_study = n_distinct(smokingcessation_ab_AC$study)
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = smokingcessation_ab_AC[smokingcessation_ab_AC$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    if (study_data$treatment[1] == "A") {
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
  
  sim_dat = smokingcessation_ab_AC
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res = summary(gemtc::relative.effect(cons.out,"A", c("B","C","D")))
  
  AB_lower_95 = res$summaries$quantiles[1,1]
  AB_upper_95 = res$summaries$quantiles[1,5]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)
  
  AC_lower_95 = res$summaries$quantiles[2,1]
  AC_upper_95 = res$summaries$quantiles[2,5]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)
  
  AD_lower_95 = res$summaries$quantiles[3,1]
  AD_upper_95 = res$summaries$quantiles[3,5]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)
  
  reject_correct = c(AB_reject, AC_reject, AD_reject)
  
  return(c(reject_correct))
}

result_smoke = result_smoke/S

result_smoke = matrix(result_smoke, nrow = 1)
colnames(result_smoke) <- c("AB_power", "AC_power", "AD_power")

save(result_smoke, file = "result_smoke_AC_add5.RData")






# AD
smokingcessation_ab_AD = rbind(smokingcessation_ab, c(25, "A", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_a)) )
smokingcessation_ab_AD = rbind(smokingcessation_ab_AD, c(25, "D", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_d)) )
smokingcessation_ab_AD = rbind(smokingcessation_ab_AD, c(26, "A", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_a)) )
smokingcessation_ab_AD = rbind(smokingcessation_ab_AD, c(26, "D", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_d)) )
smokingcessation_ab_AD = rbind(smokingcessation_ab_AD, c(27, "A", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_a)) )
smokingcessation_ab_AD = rbind(smokingcessation_ab_AD, c(27, "D", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_d)) )
smokingcessation_ab_AD = rbind(smokingcessation_ab_AD, c(28, "A", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_a)) )
smokingcessation_ab_AD = rbind(smokingcessation_ab_AD, c(28, "D", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_d)) )
smokingcessation_ab_AD = rbind(smokingcessation_ab_AD, c(29, "A", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_a)) )
smokingcessation_ab_AD = rbind(smokingcessation_ab_AD, c(29, "D", round(mean(smokingcessation_ab$sampleSize)), round(mean(smokingcessation_ab$sampleSize) * pi_d)) )
smokingcessation_ab_AD$sampleSize = as.numeric(smokingcessation_ab_AD$sampleSize)
smokingcessation_ab_AD$responders = as.numeric(smokingcessation_ab_AD$responders)
n_study = n_distinct(smokingcessation_ab_AD$study)

result_smoke = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  library(R2jags)
  library(tidyverse)
  library(gemtc)
  
  Log_OR_dat = c(0, 0.4944, 0.8378, 1.0983)
  tau = 0.8458
  n_study = n_distinct(smokingcessation_ab_AD$study)
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = smokingcessation_ab_AD[smokingcessation_ab_AD$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
    if (study_data$treatment[1] == "A") {
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
  
  sim_dat = smokingcessation_ab_AD
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random",
                          hy.prior =  mtc.hy.prior(type="std.dev", distr="dunif", 0.01, 10),
                          re.prior.sd = 10)
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=2000, n.iter=5000, thin=1)
  
  res = summary(gemtc::relative.effect(cons.out,"A", c("B","C","D")))
  
  AB_lower_95 = res$summaries$quantiles[1,1]
  AB_upper_95 = res$summaries$quantiles[1,5]
  AB_reject = 1-as.numeric(0 >= AB_lower_95 & 0 <= AB_upper_95)
  
  AC_lower_95 = res$summaries$quantiles[2,1]
  AC_upper_95 = res$summaries$quantiles[2,5]
  AC_reject = 1-as.numeric(0 >= AC_lower_95 & 0 <= AC_upper_95)
  
  AD_lower_95 = res$summaries$quantiles[3,1]
  AD_upper_95 = res$summaries$quantiles[3,5]
  AD_reject = 1-as.numeric(0 >= AD_lower_95 & 0 <= AD_upper_95)
  
  reject_correct = c(AB_reject, AC_reject, AD_reject)
  
  return(c(reject_correct))
}

result_smoke = result_smoke/S

result_smoke = matrix(result_smoke, nrow = 1)
colnames(result_smoke) <- c("AB_power", "AC_power", "AD_power")

save(result_smoke, file = "result_smoke_AD_add5.RData")


stopCluster(cl)






