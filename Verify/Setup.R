library(tidyverse)
library(foreach)
library(doParallel)
library(magrittr)
library(R2jags)

set.seed(123456)

N_cores = detectCores() 
N_sim = 500


load("df_indirect_blank.RData")
load("df_direct_blank.RData")
load("df_BNMA_blank.RData")
load("df_indirect_bias_additional_blank.RData")

df_direct_blank = df_direct_blank[,1:8]
df_indirect_blank = df_indirect_blank[,1:9]
df_indirect_bias_additional_blank = df_indirect_bias_additional_blank[,1:9]



df_BNMA_blank = df_BNMA_blank[1:9,1:10]


NMA_Net4 = do.call(rbind, replicate(6, df_BNMA_blank, simplify = FALSE))
NMA_Net5 = do.call(rbind, replicate(6, df_BNMA_blank, simplify = FALSE))
NMA_Net6 = do.call(rbind, replicate(6, df_BNMA_blank, simplify = FALSE))


NMA_Net4$k_ab = 2
NMA_Net4$k_ac = c(rep(1,9),rep(2,9),rep(3,9),rep(6,9),rep(9,9),rep(12,9))
NMA_Net4$k_bc = 2


NMA_Net5$k_ab = c(rep(1,9),rep(2,9),rep(3,9),rep(6,9),rep(9,9),rep(12,9))
NMA_Net5$k_ac = c(rep(1,9),rep(2,9),rep(3,9),rep(6,9),rep(9,9),rep(12,9))
NMA_Net5$k_bc = 2


NMA_Net6$k_ab = 2
NMA_Net6$k_ac = 2
NMA_Net6$k_bc = c(rep(1,9),rep(2,9),rep(3,9),rep(6,9),rep(9,9),rep(12,9))


