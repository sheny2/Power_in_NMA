library(foreach)
library(doParallel)
# Initialize parallel backend
cl <- makeCluster(20)

# Register the parallel backend
registerDoParallel(cl)

S = 5000

library(tidyverse)
library(gemtc)

source("ifplot.fun.R")
source("BayesDiagnos.fun.R")
source("ranko_sucra.R")

set.seed(2023)

######## Import data
bleed <- read.table("Data_Bleeding.csv", sep=",", header=TRUE)
mace <- read.table("Data_MACE.csv", sep=",", header=TRUE)


######## Row 7 and 8 need to be combined
bleed[7,3:7] = bleed[7,3:7] + bleed[8,3:7]
bleed <- bleed[-8,]
mace[7,3:10] = mace[7,3:10] + mace[8,3:10]
mace <- mace[-8,]


Bleed_data <- bleed[,c("study", "treatment", "n", "TIMImajor")]
Bleed_data = Bleed_data[1:(nrow(Bleed_data)-2), ]

MACE_data <- mace[,c("study", "treatment", "n", "MACE")]        
MACE_data = MACE_data[1:(nrow(MACE_data)-2), ]


## Note: A is reference 

trt <- c("A", "B", "C","D")
# trt <- c("VKA + DAPT", "VKA + P2Y12", "NOAC + DAPT", "NOAC + P2Y12")

trts <- read.table(textConnection('id description
                                  A "VKA + DAPT"
                                  B "VKA + P2Y12"
                                  C "NOAC + DAPT"
                                  D "NOAC + P2Y12"'), header=TRUE)


Bleed_data$treatment <- trt[Bleed_data$treatment]
colnames(Bleed_data) = c("study","treatment", "sampleSize","responders")
MACE_data$treatment <- trt[MACE_data$treatment]
colnames(MACE_data) = c("study","treatment", "sampleSize","responders")


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


### Major Bleeding
# pi_a = sum(bleed[bleed$treatment==1,]$TIMImajor) / sum(bleed[bleed$treatment==1,]$n)
pi_a = sum(Bleed_data[Bleed_data$treatment=="A",]$responders) / sum(Bleed_data[Bleed_data$treatment=="A",]$sampleSize)

tau = 0.24

# True OR (from paper)
OR_ab = 0.58
OR_ac = 0.70
OR_ad = 0.49
Log_OR_dat = c(1, log(c(OR_ab,OR_ac,OR_ad)))

pi_b =  calculate_pi(OR_ab, pi_a)
pi_c =  calculate_pi(OR_ac, pi_a)
pi_d =  calculate_pi(OR_ad, pi_a)


# add AB
Bleed_data_AB = rbind(Bleed_data, c(5, "A", round(mean(Bleed_data$sampleSize)), round(mean(Bleed_data$sampleSize) * pi_a)) )
Bleed_data_AB = rbind(Bleed_data_AB, c(5, "B", round(mean(Bleed_data$sampleSize)), round(mean(Bleed_data$sampleSize) * pi_b)) )
Bleed_data_AB$sampleSize = as.numeric(Bleed_data_AB$sampleSize)
Bleed_data_AB$responders = as.numeric(Bleed_data_AB$responders)
n_study = n_distinct(Bleed_data_AB$study)

result_bleeding = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = Bleed_data_AB[Bleed_data_AB$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
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
  }
  
  sim_dat = Bleed_data_AB
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat, treatments=trts)
  cons.model <- gemtc::mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random")
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=500, n.iter=2000, thin=1)
  
  prob <- gemtc::rank.probability(cons.out, preferredDirection = -1)
  prob <- round(prob, digits=3)
  sucra <- gemtc::sucra(prob)
  rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
  A_last = as.numeric(rank_order[4] == "A")
  
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
  
  return(c(reject_correct, reject_incorrect, A_last))
}

# first three powers / next three Type I errors / last one A being ranked last
result_rate_bleeding = result_bleeding/S

result_rate_bleeding = matrix(result_rate_bleeding, nrow = 1)
colnames(result_rate_bleeding) <- c("Power AB","Power AC","Power AD",
                                    "Type I err BC","Type I err BD","Type I err CD",
                                    "Rank A last")

save(result_rate_bleeding, file = "result_rate_bleeding_AB.RData")

# load("result_rate_bleeding.RData")



# add AC
Bleed_data_AC = rbind(Bleed_data, c(5, "A", round(mean(Bleed_data$sampleSize)), round(mean(Bleed_data$sampleSize) * pi_a)) )
Bleed_data_AC = rbind(Bleed_data_AC, c(5, "C", round(mean(Bleed_data$sampleSize)), round(mean(Bleed_data$sampleSize) * pi_c)) )
Bleed_data_AC$sampleSize = as.numeric(Bleed_data_AC$sampleSize)
Bleed_data_AC$responders = as.numeric(Bleed_data_AC$responders)
n_study = n_distinct(Bleed_data_AC$study)

result_bleeding = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = Bleed_data_AC[Bleed_data_AC$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
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
  }
  
  sim_dat = Bleed_data_AC
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat, treatments=trts)
  cons.model <- gemtc::mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random")
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=500, n.iter=2000, thin=1)
  
  prob <- gemtc::rank.probability(cons.out, preferredDirection = -1)
  prob <- round(prob, digits=3)
  sucra <- gemtc::sucra(prob)
  rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
  A_last = as.numeric(rank_order[4] == "A")
  
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
  
  return(c(reject_correct, reject_incorrect, A_last))
}

# first three powers / next three Type I errors / last one A being ranked last
result_rate_bleeding = result_bleeding/S

result_rate_bleeding = matrix(result_rate_bleeding, nrow = 1)
colnames(result_rate_bleeding) <- c("Power AB","Power AC","Power AD",
                                    "Type I err BC","Type I err BD","Type I err CD",
                                    "Rank A last")

save(result_rate_bleeding, file = "result_rate_bleeding_AC.RData")




# add AD
Bleed_data_AD = rbind(Bleed_data, c(5, "A", round(mean(Bleed_data$sampleSize)), round(mean(Bleed_data$sampleSize) * pi_a)) )
Bleed_data_AD = rbind(Bleed_data_AD, c(5, "D", round(mean(Bleed_data$sampleSize)), round(mean(Bleed_data$sampleSize) * pi_d)) )
Bleed_data_AD$sampleSize = as.numeric(Bleed_data_AD$sampleSize)
Bleed_data_AD$responders = as.numeric(Bleed_data_AD$responders)
n_study = n_distinct(Bleed_data_AD$study)

result_bleeding = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = Bleed_data_AD[Bleed_data_AD$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
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
  }
  
  sim_dat = Bleed_data_AD
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat, treatments=trts)
  cons.model <- gemtc::mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random")
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=500, n.iter=2000, thin=1)
  
  prob <- gemtc::rank.probability(cons.out, preferredDirection = -1)
  prob <- round(prob, digits=3)
  sucra <- gemtc::sucra(prob)
  rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
  A_last = as.numeric(rank_order[4] == "A")
  
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
  
  return(c(reject_correct, reject_incorrect, A_last))
}

# first three powers / next three Type I errors / last one A being ranked last
result_rate_bleeding = result_bleeding/S

result_rate_bleeding = matrix(result_rate_bleeding, nrow = 1)
colnames(result_rate_bleeding) <- c("Power AB","Power AC","Power AD",
                                    "Type I err BC","Type I err BD","Type I err CD",
                                    "Rank A last")

save(result_rate_bleeding, file = "result_rate_bleeding_AD.RData")





# add CD
Bleed_data_CD = rbind(Bleed_data, c(5, "C", round(mean(Bleed_data$sampleSize)), round(mean(Bleed_data$sampleSize) * pi_c)) )
Bleed_data_CD = rbind(Bleed_data_CD, c(5, "D", round(mean(Bleed_data$sampleSize)), round(mean(Bleed_data$sampleSize) * pi_d)) )
Bleed_data_CD$sampleSize = as.numeric(Bleed_data_CD$sampleSize)
Bleed_data_CD$responders = as.numeric(Bleed_data_CD$responders)
n_study = n_distinct(Bleed_data_CD$study)

result_bleeding = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = Bleed_data_CD[Bleed_data_CD$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
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
  }
  
  sim_dat = Bleed_data_CD
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat, treatments=trts)
  cons.model <- gemtc::mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random")
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=500, n.iter=2000, thin=1)
  
  prob <- gemtc::rank.probability(cons.out, preferredDirection = -1)
  prob <- round(prob, digits=3)
  sucra <- gemtc::sucra(prob)
  rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
  A_last = as.numeric(rank_order[4] == "A")
  
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
  
  return(c(reject_correct, reject_incorrect, A_last))
}

# first three powers / next three Type I errors / last one A being ranked last
result_rate_bleeding = result_bleeding/S

result_rate_bleeding = matrix(result_rate_bleeding, nrow = 1)
colnames(result_rate_bleeding) <- c("Power AB","Power AC","Power AD",
                                    "Type I err BC","Type I err BD","Type I err CD",
                                    "Rank A last")

save(result_rate_bleeding, file = "result_rate_bleeding_CD.RData")




# add BC
Bleed_data_BC = rbind(Bleed_data, c(5, "B", round(mean(Bleed_data$sampleSize)), round(mean(Bleed_data$sampleSize) * pi_b)) )
Bleed_data_BC = rbind(Bleed_data_BC, c(5, "C", round(mean(Bleed_data$sampleSize)), round(mean(Bleed_data$sampleSize) * pi_c)) )
Bleed_data_BC$sampleSize = as.numeric(Bleed_data_BC$sampleSize)
Bleed_data_BC$responders = as.numeric(Bleed_data_BC$responders)
n_study = n_distinct(Bleed_data_BC$study)

result_bleeding = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = Bleed_data_BC[Bleed_data_BC$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
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
  }
  
  sim_dat = Bleed_data_BC
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat, treatments=trts)
  cons.model <- gemtc::mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random")
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=500, n.iter=2000, thin=1)
  
  prob <- gemtc::rank.probability(cons.out, preferredDirection = -1)
  prob <- round(prob, digits=3)
  sucra <- gemtc::sucra(prob)
  rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
  A_last = as.numeric(rank_order[4] == "A")
  
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
  
  return(c(reject_correct, reject_incorrect, A_last))
}

# first three powers / next three Type I errors / last one A being ranked last
result_rate_bleeding = result_bleeding/S

result_rate_bleeding = matrix(result_rate_bleeding, nrow = 1)
colnames(result_rate_bleeding) <- c("Power AB","Power AC","Power AD",
                                    "Type I err BC","Type I err BD","Type I err CD",
                                    "Rank A last")

save(result_rate_bleeding, file = "result_rate_bleeding_BC.RData")





# add BC
Bleed_data_BD = rbind(Bleed_data, c(5, "B", round(mean(Bleed_data$sampleSize)), round(mean(Bleed_data$sampleSize) * pi_b)) )
Bleed_data_BD = rbind(Bleed_data_BD, c(5, "D", round(mean(Bleed_data$sampleSize)), round(mean(Bleed_data$sampleSize) * pi_d)) )
Bleed_data_BD$sampleSize = as.numeric(Bleed_data_BD$sampleSize)
Bleed_data_BD$responders = as.numeric(Bleed_data_BD$responders)
n_study = n_distinct(Bleed_data_BD$study)

result_bleeding = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = Bleed_data_BD[Bleed_data_BD$study==i,]
    
    alpha_iB = rnorm(n=1, mean=logit(pi_a), sd=0.1)
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
  }
  
  sim_dat = Bleed_data_BD
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat, treatments=trts)
  cons.model <- gemtc::mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random")
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=500, n.iter=2000, thin=1)
  
  prob <- gemtc::rank.probability(cons.out, preferredDirection = -1)
  prob <- round(prob, digits=3)
  sucra <- gemtc::sucra(prob)
  rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
  A_last = as.numeric(rank_order[4] == "A")
  
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
  
  return(c(reject_correct, reject_incorrect, A_last))
}

# first three powers / next three Type I errors / last one A being ranked last
result_rate_bleeding = result_bleeding/S

result_rate_bleeding = matrix(result_rate_bleeding, nrow = 1)
colnames(result_rate_bleeding) <- c("Power AB","Power AC","Power AD",
                                    "Type I err BC","Type I err BD","Type I err CD",
                                    "Rank A last")

save(result_rate_bleeding, file = "result_rate_bleeding_BD.RData")