library(gemtc)
library(tidyverse)
load("parkinson_ab.RData")

set.seed(8848)

# Name
# 1: Placebo
# 2: Pramipexole
# 3: Ropinirole
# 4: Bromocriptine
# 5: Cabergoline

parkinson_ab

# normal/identity:forcontinuous(meandifference)data. 
# Requiredcolumns:[mean,std.err]or[mean,std.dev,sampleSize]. 
# Result:relativemeandifference.

parkinson_net <- mtc.network(parkinson_ab)
cons.model <- mtc.model(parkinson_net, type="consistency", 
                        likelihood="normal", link="identity", linearModel="random")
cons.out <- mtc.run(cons.model, n.adapt=20000, n.iter=50000, thin=1)
prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
prob <- round(prob, digits=3)
sucra <- gemtc::sucra(prob)
correct_rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))


estimates = summary(gemtc::relative.effect(cons.out,"1", c("2","3","4","5")))
estimates$measure
estimates$summaries$statistics

d12 = estimates$summaries$statistics[1,1]
d13 = estimates$summaries$statistics[2,1]
d14 = estimates$summaries$statistics[3,1]
d15 = estimates$summaries$statistics[4,1]
tau = estimates$summaries$statistics[5,1]



# simulation

library(foreach)
library(doParallel)
# Initialize parallel backend
cl <- makeCluster(8)

# Register the parallel backend
registerDoParallel(cl)


S = 5000


n_study = n_distinct(parkinson_ab$study)

mu_1 = mean(parkinson_ab[which(parkinson_ab$treatment == 1), ]$mean)
# mu_2 = mu_1 + d12 
# mu_3 = mu_1 + d13 
# mu_4 = mu_1 + d14 

dk = c(0, d12,d13,d14,d15)
trt = c("1","2","3","4","5")


result_pk = foreach (i = 1:S, .combine = "+", .errorhandling='remove') %dopar% {
  
  y_ik = c()
  for (i in 1:n_study){
    study_data = parkinson_ab[parkinson_ab$study==i,]
    
    alpha_iB = rnorm(n = 1, mean=mu_1, sd=0.1)
    
    if (study_data$treatment[1] == "1"){
      for (k in 1:nrow(study_data)){
        if (k == 1) {
          delta_ik = alpha_iB
          sd_ik = study_data$std.dev[k]
          n_ik = study_data$sampleSize[k]
          y_ik = c(y_ik, rnorm(n = 1, mean = delta_ik, sd = sqrt(sd_ik^2/n_ik)))
        } else {
          delta_ik = alpha_iB + rnorm(n = 1, mean = dk[which(trt == study_data$treatment[k])], sd = tau)
          sd_ik = study_data$std.dev[k]
          n_ik = study_data$sampleSize[k]
          y_ik = c(y_ik, rnorm(n = 1, mean = delta_ik, sd = sqrt(sd_ik^2/n_ik)))
        }
      }
    } else{
      for (k in 1:nrow(study_data)){
        delta_ik = alpha_iB + rnorm(n = 1, mean = dk[which(trt == study_data$treatment[k])], sd = tau)
        sd_ik = study_data$std.dev[k]
        n_ik = study_data$sampleSize[k]
        y_ik = c(y_ik, rnorm(n = 1, mean = delta_ik, sd = sqrt(sd_ik^2/n_ik)))
      }
    }
  }
    
  sim_dat = parkinson_ab
  sim_dat$mean = y_ik
  
  network <- gemtc::mtc.network(data.ab=sim_dat)
  cons.model <- gemtc::mtc.model(network, type="consistency", 
                                 likelihood="normal", link="identity", linearModel="random")
  cons.out <-gemtc::mtc.run(cons.model, n.adapt=500, n.iter=2000, thin=1)
  
  prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
  prob <- round(prob, digits=3)
  sucra <- gemtc::sucra(prob)
  rank_order = rownames(as.matrix(sort(sucra, decreasing = T)))
  # A_last = as.numeric(rank_order[4] == "A")
  Correct_order = as.numeric(identical(rank_order,correct_rank_order))
  
  res = summary(gemtc::relative.effect(cons.out,"1", c("2","3","4","5")))
    
  d12_lower_95 = res$summaries$quantiles[1,1]
  d12_upper_95 = res$summaries$quantiles[1,5]
  d12_reject = 1-as.numeric(0 >= d12_lower_95 & 0 <= d12_upper_95)
  
  d13_lower_95 = res$summaries$quantiles[2,1]
  d13_upper_95 = res$summaries$quantiles[2,5]
  d13_reject = 1-as.numeric(0 >= d13_lower_95 & 0 <= d13_upper_95)
  
  d14_lower_95 = res$summaries$quantiles[3,1]
  d14_upper_95 = res$summaries$quantiles[3,5]
  d14_reject = 1-as.numeric(0 >= d14_lower_95 & 0 <= d14_upper_95)
  
  d15_lower_95 = res$summaries$quantiles[4,1]
  d15_upper_95 = res$summaries$quantiles[4,5]
  d15_reject = 1-as.numeric(0 >= d15_lower_95 & 0 <= d15_upper_95)

  res_2 = summary(gemtc::relative.effect(cons.out,"2", c("3","4","5")))
  
  d23_lower_95 = res_2$summaries$quantiles[1,1]
  d23_upper_95 = res_2$summaries$quantiles[1,5]
  d23_reject = 1-as.numeric(0 >= d23_lower_95 & 0 <= d23_upper_95)
  
  d24_lower_95 = res_2$summaries$quantiles[2,1]
  d24_upper_95 = res_2$summaries$quantiles[2,5]
  d24_reject = 1-as.numeric(0 >= d24_lower_95 & 0 <= d24_upper_95)
  
  d25_lower_95 = res_2$summaries$quantiles[3,1]
  d25_upper_95 = res_2$summaries$quantiles[3,5]
  d25_reject = 1-as.numeric(0 >= d25_lower_95 & 0 <= d25_upper_95)

  
  res_34 = summary(gemtc::relative.effect(cons.out,c("3","3","4"), c("4","5","5")))
  
  d34_lower_95 = res_34$summaries$quantiles[1,1]
  d34_upper_95 = res_34$summaries$quantiles[1,5]
  d34_reject = 1-as.numeric(0 >= d34_lower_95 & 0 <= d34_upper_95)
  
  d35_lower_95 = res_34$summaries$quantiles[2,1]
  d35_upper_95 = res_34$summaries$quantiles[2,5]
  d35_reject = 1-as.numeric(0 >= d35_lower_95 & 0 <= d35_upper_95)
  
  d45_lower_95 = res_34$summaries$quantiles[3,1]
  d45_upper_95 = res_34$summaries$quantiles[3,5]
  d45_reject = 1-as.numeric(0 >= d45_lower_95 & 0 <= d45_upper_95)
  
  
  reject_correct = c(d12_reject, d13_reject, d14_reject, d15_reject, 
                     d23_reject, d24_reject, d25_reject, 
                     d34_reject, d35_reject, d45_reject, Correct_order)
  
  return(c(reject_correct))
}


result_pk = result_pk/S

result_pk = matrix(result_pk, nrow = 1)
colnames(result_pk) <- c("Power D12","Power D13","Power D14", "Power D15", 
                         "Power D23","Power D24","Power D25", 
                         "Power D34","Power D35","Power D45", 
                         "Correct_order")


result_pk

save(result_pk, file = "result_pk.RData")













