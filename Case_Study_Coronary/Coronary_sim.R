library(foreach)
library(doParallel)
# Initialize parallel backend
cl <- makeCluster(12)

# Register the parallel backend
registerDoParallel(cl)

S = 50


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

trt <- c("A","B", "C","D")
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

network <- mtc.network(data.ab=Bleed_data, treatments=trts)
plot(network)
summary(network)


cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random")
cons.out <- mtc.run(cons.model, n.adapt=20000, n.iter=50000, thin=1)
summary(cons.out)
gemtc::forest(cons.out)

prob <- rank.probability(cons.out, preferredDirection=-1)
prob <- round(prob, digits=3)



# Helper function

calculate_pi <- function(odds_ratio, p_A) {
  odds_A <- p_A / (1 - p_A)
  odds_B <- odds_A * odds_ratio
  p_B <- odds_B / (1 + odds_B)
  
  return(p_B)
}

# from prob to log odds
logit <- function(p) {
  return(log(p/(1-p)))
}

# from log odds to prob
inverse_logit <- function(x) {
  odds <- exp(x)  
  probability <- odds / (1 + odds) 
  return(probability)  
}




# Generate Data (Lu & Ades model)


### Major Bleeding
# pi_a = sum(bleed[bleed$treatment==1,]$TIMImajor) / sum(bleed[bleed$treatment==1,]$n)
pi_a = sum(Bleed_data[Bleed_data$treatment=="A",]$responders) / sum(Bleed_data[Bleed_data$treatment=="A",]$sampleSize)

n_study = n_distinct(Bleed_data$study)

tau = 0.24

# True OR (from paper)
OR_ab = 0.58
OR_ac = 0.70
OR_ad = 0.49
Log_OR_dat = log(c(1, OR_ab,OR_ac,OR_ad))

pi_b =  calculate_pi(OR_ab, pi_a)
pi_c =  calculate_pi(OR_ac, pi_a)
pi_d =  calculate_pi(OR_ad, pi_a)


result_bleeding = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = Bleed_data[Bleed_data$study==i,]
    
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
  
  sim_dat = Bleed_data
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
  
  return(c(reject_correct, reject_incorrect, A_last))
}

# first three powers / next three Type I errors / last one A being ranked last
result_rate_bleeding = result_bleeding/S

result_rate_bleeding = matrix(result_rate_bleeding, nrow = 1)
colnames(result_rate_bleeding) <- c("Power AB","Power AC","Power AD",
                                    "Type I err BC","Type I err BD","Type I err CD",
                                    "Rank A last")

save(result_rate_bleeding, file = "result_rate_bleeding.RData")

# load("result_rate_bleeding.RData")







### MACE

pi_a = sum(mace[mace$treatment==1,]$MACE) / sum(mace[mace$treatment==1,]$n)

tau = 0.23

OR_ab = 0.96
OR_ac = 0.94
OR_ad = 1.02
Log_OR_dat = log(c(1, OR_ab,OR_ac,OR_ad))

pi_b =  calculate_pi(OR_ab, pi_a)
pi_c =  calculate_pi(OR_ac, pi_a)
pi_d =  calculate_pi(OR_ad, pi_a)



order_MACE = foreach (i = 1:10000, .combine = "c", .errorhandling='remove') %dopar% {
  
  y_ik = c()
  for (i in 1:4){
    study_data = MACE_data[MACE_data$study==i,]
    
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
  
  sim_dat = MACE_data
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat, treatments=trts)
  cons.model <- gemtc::mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random")
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=500, n.iter=2000, thin=1)
  # summary(cons.out)
  
  prob <- gemtc::rank.probability(cons.out, preferredDirection = -1)
  prob <- round(prob, digits=3)
  sucra <- gemtc::sucra(prob)
  rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
  rank_order = paste(rank_order,collapse=' ')
  
  return(rank_order)
}

# order_MACE
table(order_MACE)
prop.table(table(order_MACE))

save(order_MACE, file = "order_MACE.RData")



result_MACE = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  
  y_ik = c()
  for (i in 1:4){
    study_data = MACE_data[MACE_data$study==i,]
    
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
  
  sim_dat = MACE_data
  sim_dat$responders = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat, treatments=trts)
  cons.model <- gemtc::mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random")
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=500, n.iter=2000, thin=1)
  # summary(cons.out)
  
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
  
  type_1_err = c(AB_reject, AC_reject, AD_reject, BC_reject, BD_reject, CD_reject)
  
  return(type_1_err)
}

result_rate_MACE = result_MACE/S
result_rate_MACE = matrix(result_rate_MACE, nrow = 1)
colnames(result_rate_MACE) <- c("Type I err AB","Type I err AC","Type I err AD",
                                "Type I err BC","Type I err BD","Type I err CD")

save(result_rate_MACE, file = "result_rate_MACE.RData")




