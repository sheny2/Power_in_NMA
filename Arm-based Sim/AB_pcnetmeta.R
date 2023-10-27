library(pcnetmeta)


# 
# 
# data("smoke")
# # increase n.iter to reach convergence
# set.seed(1234)
# nma.out <- nma.ab.bin(s.id, t.id, r, n, data = smoke,
#                       trtname = c("NC", "SH", "IC", "GC"), param= "AR",
#                       model = "het_cor", n.adapt = 1000, n.iter = 100, n.chains = 1)
# absolute.plot(nma.out, save = FALSE)
# #absolute.plot(nma.out)
# absolute.plot(nma.out, alphabetic = FALSE, save = FALSE)
# nma.out



################################# 


k_ab = k_ac = 6
k_bc = 6

pi_a = 0.5

OR_ab = 1.2
OR_ac = 1.8

tau = 0.2
S = 10


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



# data pre for jags
NS = n_distinct(all_data$study)
NT = n_distinct(all_data$treatment)
N = nrow(all_data)
s = all_data$study
# t = all_data$treatment
t = as.integer(factor(all_data$treatment, levels = c("B","A","C"), labels = c(1,2,3)))
y = all_data$responders
n = all_data$sampleSize
drug_list <- c("B","A","C")
Narm <- as.numeric(table(all_data$study))
n.obs <- matrix(NA,nrow=NS, ncol=max(Narm))
n.eve <- matrix(NA,nrow=NS, ncol=max(Narm))
dr <- matrix(NA,nrow=NS, ncol=max(Narm))
study<-unique(all_data$study)

for (i in 1:NS){
  n.obs[i,1:Narm[i]] <- all_data$sampleSize[all_data$study==study[i]]
  n.eve[i,1:Narm[i]] <- all_data$responders[all_data$study==study[i]]
  dr[i,1:Narm[i]] <- match(all_data$treatment[all_data$study==study[i]],drug_list)
}


### Running Arm Based  Model  ########
data_AB <- list('Narm'=N, 'Nstudy'=NS,
                'Ndrug'=NT, 'study'= s, 'drug'=t,
                'y'=y, 'n'=n ,'Omega'=diag(rep(0.2,times=3)),
                'zero.AB' = (rep(0, times=3)))
inits_AB<- list(list(mu=rep(0,3)),
                list(mu=rep(0,3)))
para_AB<-c( "lor", "tau", "best1", "best2", "best3")
fit_AB<-jags(data=data_AB, inits=inits_AB, para_AB,
             n.iter = 5000, n.burnin = 2000, n.chains = 2, n.thin = 1,
             DIC=TRUE, model.file=ABWish)


# output data
AB_trt_results<-data.frame(fit_AB$BUGSoutput$summary[,c(1, 3, 7)])
AB_trt_results <- tibble::rownames_to_column(AB_trt_results, "drug_list")
AB_trt_results<-AB_trt_results%>%
  filter(drug_list %in% c("lor[1]", "lor[2]", "lor[3]"))

AB_rank_prob<-data.frame(fit_AB$BUGSoutput$summary[,c(1, 3, 7)])
AB_rank_prob <- tibble::rownames_to_column(AB_rank_prob, "best")

AB_rank_prob <- data.frame(cbind(AB_rank_prob[AB_rank_prob$best%in%c("best1[1]", "best1[2]", "best1[3]"),2],
                                 AB_rank_prob[AB_rank_prob$best%in%c("best2[1]", "best2[2]", "best2[3]"),2],
                                 AB_rank_prob[AB_rank_prob$best%in%c("best3[1]", "best3[2]", "best3[3]"),2]))
rownames(AB_rank_prob)<-c("trt_1", "trt_2", "trt_3")

AB_trt_results[,1:2]


all_data_pc = tibble(s.id = all_data$study,
                    t.id = all_data$treatment, 
                    r = all_data$responders, 
                    n = all_data$sampleSize) 


par(mfrow = c(1,3))
AB_Result_het_cor = nma.ab.bin(s.id, t.id, r, n, data = all_data_pc, param= "LOR",
           model = "het_cor", n.adapt = 2000, n.iter = 5000, n.chains = 2)
AB_Result_het_cor$LogOddsRatio$Mean_SD[1,2]
AB_Result_het_cor$LogOddsRatio$Mean_SD[3,2]
contrast.plot(AB_Result_het_cor, save = FALSE, reference = "B")



AB_Result_het_eqcor = nma.ab.bin(s.id, t.id, r, n, data = all_data_pc, param= "LOR",
                          model = "het_eqcor", n.adapt = 2000, n.iter = 5000, n.chains = 2)

AB_Result_het_eqcor$LogOddsRatio$Mean_SD[1,2]
AB_Result_het_eqcor$LogOddsRatio$Mean_SD[3,2]
contrast.plot(AB_Result_het_eqcor, save = FALSE, reference = "B")


AB_Result_hom_eqcor = nma.ab.bin(s.id, t.id, r, n, data = all_data_pc, param= "LOR",
                          model = "hom_eqcor", n.adapt = 2000, n.iter = 5000, n.chains = 2)

AB_Result_hom_eqcor$LogOddsRatio$Mean_SD[1,2]
AB_Result_hom_eqcor$LogOddsRatio$Mean_SD[3,2]
contrast.plot(AB_Result_hom_eqcor, save = FALSE, reference = "B")


par(mfrow = c(1,3))
contrast.plot(AB_Result_het_cor, save = FALSE, reference = "B")
contrast.plot(AB_Result_het_eqcor, save = FALSE, reference = "B")
contrast.plot(AB_Result_hom_eqcor, save = FALSE, reference = "B")