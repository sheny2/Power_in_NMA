# Simulation to compare (several) big studies and (lots of) small studies


library(tidyverse)
library(dplyr)
library(foreach)
library(doParallel)
library(magrittr)
library(gemtc)

set.seed(123456)

N_cores = detectCores() - 1
N_sim = 1000


source("simulation_big_small.R")


pi_a <- c(0.05, 0.5)

tau <- c(0.001, 0.2, 0.4)

OR_ab <- c(1.2)
OR_ac <- c(1.4, 1.6, 1.8)

k_ab = c(2,5,10,25)
k_bc = 0


# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = bigsmall.sim(N_sim, k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, tau, total = 3000))
}

# Apply the function in parallel using foreach
indirect_result_big_small <- expand_grid(pi_a, OR_ab, OR_ac, tau, k_ab) %>% 
  mutate(k_bc = k_bc, k_ac = k_ab) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
  mutate(pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )) %>% 
  mutate(pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )) %>% 
  group_by(k_ab, k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.)) %>% separate(power, c("power", "rank_correct"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(indirect_result_big_small, file = "indirect_result_big_small.RData")




pi_a <- c(0.05, 0.5)

tau <- c(0.001, 0.2, 0.4)

OR_ab <- c(1.2)
OR_ac <- c(1.4, 1.6, 1.8)

k_ab = 0
k_ac = 0
k_bc = c(3,6,10,15,25)


# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = bigsmall.sim(N_sim, k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, tau,total = 3000))
}

# Apply the function in parallel using foreach
direct_result_big_small <- expand_grid(pi_a, OR_ab, OR_ac, tau, k_bc) %>% 
  mutate(k_ab = k_ab, k_ac = k_ac) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
  mutate(pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )) %>% 
  mutate(pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )) %>% 
  group_by(k_ab, k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.)) %>% separate(power, c("power", "rank_correct"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(direct_result_big_small, file = "direct_result_big_small.RData")





# overall

pi_a <- c(0.05, 0.5)

tau <- c(0.001, 0.2, 0.4)

OR_ab <- c(1.2)
OR_ac <- c(1.4, 1.6, 1.8)

k_ab = 0
k_ac = 0
k_bc = c(2,4,6,8)


# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = bigsmall.sim(N_sim, k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, tau,total = 3000))
}

# Apply the function in parallel using foreach
overall_result_big_small <- expand_grid(pi_a, OR_ab, OR_ac, tau, k_bc) %>% 
  mutate(k_ab = k_bc, k_ac = k_bc) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
  mutate(pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )) %>% 
  mutate(pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )) %>% 
  group_by(k_ab, k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.)) %>% separate(power, c("power", "rank_correct"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(overall_result_big_small, file = "overall_result_big_small.RData")
