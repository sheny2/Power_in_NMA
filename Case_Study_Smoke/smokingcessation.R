# use data(smokingcessation) 

library(netmeta)
data(smokingcessation)
colnames(smokingcessation) = c("event", "n", "event", "n", "event", "n", "treat1", "treat2", "treat3")

smokingcessation_ab <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(smokingcessation_ab) <- c("study_id", "treatments", "event", "n")

for (i in 1:24) {
  if (smokingcessation[i,9] != ""){
    treatments = as.character(smokingcessation[i,7:9])
    event_n = rbind(smokingcessation[i,1:2],smokingcessation[i,3:4],smokingcessation[i,5:6])
    study_id = rep(i, length(treatments))
    smokingcessation_ab = rbind(smokingcessation_ab, cbind(study_id, treatments, event_n))
    
  } else{ 
    treatments = as.character(smokingcessation[i,7:8])
    event_n = rbind(smokingcessation[i,1:2],smokingcessation[i,3:4])
    study_id = rep(i, length(treatments))
    smokingcessation_ab = rbind(smokingcessation_ab, cbind(study_id, treatments, event_n))
  }
}

smokingcessation_ab = data.frame(cbind(smokingcessation_ab[,1],smokingcessation_ab[,2],
                                   smokingcessation_ab[,4],smokingcessation_ab[,3]))

colnames(smokingcessation_ab) = c("study","treatment", "sampleSize","responders")
smokingcessation_ab[,c('sampleSize', 'responders')] <-
  lapply(smokingcessation_ab[,c('sampleSize', 'responders')], as.numeric)



library(tidyverse)
library(gemtc)
data(smoking)
set.seed(8848)

trt = c("A","B","C","D")
network <- mtc.network(data.ab=smokingcessation_ab)
plot(network)
plot(smoking)
# summary(network)


cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random")
cons.out <- mtc.run(cons.model, n.adapt=20000, n.iter=50000, thin=1)
estimates = summary(cons.out)
estimates

prob <- rank.probability(cons.out, preferredDirection=-1)
prob <- round(prob, digits=3)
prob




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


# simulation

library(foreach)
library(doParallel)
# Initialize parallel backend
cl <- makeCluster(12)

# Register the parallel backend
registerDoParallel(cl)


S = 1000


Log_OR_dat = c(log(1), estimates$summaries$statistics[,1][1:3])
OR_dat = exp(Log_OR_dat)


OR_ab = exp(estimates$summaries$statistics[,1][1])
OR_ac = exp(estimates$summaries$statistics[,1][2])
OR_ad = exp(estimates$summaries$statistics[,1][3])

pi_a = sum(smokingcessation_ab[smokingcessation_ab$treatment=="A",]$responders) / 
  sum(smokingcessation_ab[smokingcessation_ab$treatment=="A",]$sampleSize)

n_study = n_distinct(smokingcessation_ab$study)

tau = estimates$summaries$statistics[,1][4]

pi_b =  calculate_pi(OR_ab, pi_a)
pi_c =  calculate_pi(OR_ac, pi_a)
pi_d =  calculate_pi(OR_ad, pi_a)


result_smoke = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  
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
  cons.model <- gemtc::mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random")
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=500, n.iter=2000, thin=1)
  
  prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
  prob <- round(prob, digits=3)
  sucra <- gemtc::sucra(prob)
  rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
  A_last = as.numeric(rank_order[4] == "A")
  Correct_order = as.numeric(identical(rank_order,c("D","C","B", "A")))
  
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
  
  # OR_lower_95 = exp(lower_95)
  # OR_upper_95 = exp(upper_95)
  # 1-as.numeric(1 > OR_lower_95 & 1 < OR_upper_95)
  
  res = summary(gemtc::relative.effect(cons.out,c("B","B","C"), c("C","D","D")))
  
  BC_lower_95 = res$summaries$quantiles[1,1]
  BC_upper_95 = res$summaries$quantiles[1,5]
  BC_reject = 1-as.numeric(0 >= BC_lower_95 & 0 <= BC_upper_95)
  
  BD_lower_95 = res$summaries$quantiles[2,1]
  BD_upper_95 = res$summaries$quantiles[2,5]
  BD_reject = 1-as.numeric(0 >= BD_lower_95 & 0 <= BD_upper_95)
  
  CD_lower_95 = res$summaries$quantiles[3,1]
  CD_upper_95 = res$summaries$quantiles[3,5]
  CD_reject = 1-as.numeric(0 >= CD_lower_95 & 0 <= CD_upper_95)
  
  reject_correct = c(AB_reject, AC_reject, AD_reject)
  reject_incorrect = c(BC_reject, BD_reject, CD_reject)
  
  return(c(reject_correct, reject_incorrect, A_last, Correct_order))
}


result_smoke = result_smoke/S

result_smoke = matrix(result_smoke, nrow = 1)
colnames(result_smoke) <- c("Power AB","Power AC","Power AD",
                            "Power BC","Power BD","Power CD",
                            "Rank A last", "Al Correct order")

# Compare with relative effect
# estimates

save(result_smoke, file = "result_smoke.RData")




# Eff

mapping <- c('A' = 1, 'B' = 2, 'C' = 3, 'D' = 4)

smokingcessation_eff = smokingcessation_ab
smokingcessation_eff$treatment <- mapping[smokingcessation_eff$treatment]


smoke_eff = data.frame(sid = smokingcessation_ab$study, tid = smokingcessation_eff$treatment, 
                       r = smokingcessation_ab$responders, n = smokingcessation_ab$sampleSize)


source("functions.R")

dir.binary(smoke_eff)
eff.binary(smoke_eff)



