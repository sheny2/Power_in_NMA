library(tidyverse)
library(foreach)
library(doParallel)
library(magrittr)
library(R2jags)
library(gemtc)

set.seed(123456)

N_cores = detectCores() 
N_sim = 500

source("Power_sim_fun.R")

cl <- makeCluster(N_cores)

registerDoParallel(cl)

power_mutate <- function(df) {
  df %>% mutate(power = power.sim.NMA(N_sim, k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, tau))
}

NMA_Net4_res <- NMA_Net4 %>% 
  mutate(tau = readr::parse_number(tau)) %>% 
  group_by(k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, OR_bc, tau) %>%
  do(power_mutate(.)) %>% 
  separate(power, c("power", "rank_correct", "avg_bias", "avg_bias_abs"), " ", convert = TRUE)


stopCluster(cl)

saveRDS(NMA_Net4_res, file = "NMA_Net4_res.rds")






cl <- makeCluster(N_cores)

registerDoParallel(cl)

power_mutate <- function(df) {
  df %>% mutate(power = power.sim.NMA(N_sim, k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, tau))
}

NMA_Net5_res <- NMA_Net5 %>% 
  mutate(tau = readr::parse_number(tau)) %>% 
  group_by(k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, OR_bc, tau) %>%
  do(power_mutate(.)) %>% 
  separate(power, c("power", "rank_correct", "avg_bias", "avg_bias_abs"), " ", convert = TRUE)


stopCluster(cl)

saveRDS(NMA_Net5_res, file = "NMA_Net5_res.rds")







cl <- makeCluster(N_cores)

registerDoParallel(cl)

power_mutate <- function(df) {
  df %>% mutate(power = power.sim.NMA(N_sim, k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, tau))
}

NMA_Net6_res <- NMA_Net6 %>% 
  mutate(tau = readr::parse_number(tau)) %>% 
  group_by(k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, OR_bc, tau) %>%
  do(power_mutate(.)) %>% 
  separate(power, c("power", "rank_correct", "avg_bias", "avg_bias_abs"), " ", convert = TRUE)


stopCluster(cl)

saveRDS(NMA_Net6_res, file = "NMA_Net6_res.rds")