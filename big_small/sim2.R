# Simulation to compare (several) big studies and (lots of) small studies


library(tidyverse)
library(dplyr)
library(foreach)
library(doParallel)
library(magrittr)
library(gemtc)

set.seed(123456)

N_cores = 20
N_sim = 5000

source("simulation_big.R")
source("simulation_small.R")



# big 4 * 600
# small 12 * 200



### In the setting of only indirect evidence

# 4 big studies approx 600 each
# two studies of AB, two studies of BC

# 12 big studies approx 600 each
# six studies of AB, six studies of BC


# big 
pi_a <- c(0.05, 0.5)

tau <- c(0.001, 0.2, 0.4)

OR_ab <- c(1.2)
OR_ac <- c(1.4, 1.6, 1.8)

k_ab = 2
k_ac = 2
k_bc = 0


# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = big.sim(N_sim, k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
indirect_result_big <- expand_grid(pi_a, OR_ab, OR_ac, tau) %>% 
  mutate(k_ab = k_ab, k_ac = k_ab, k_bc = k_ab) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
  mutate(pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )) %>% 
  mutate(pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )) %>% 
  group_by(k_ab, k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.)) %>% separate(power, c("power", "rank_correct"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(indirect_result_big, file = "indirect_result_big.RData")




# small 
pi_a <- c(0.05, 0.5)

tau <- c(0.001, 0.2, 0.4)

OR_ab <- c(1.2)
OR_ac <- c(1.4, 1.6, 1.8)

k_ab = 6
k_ac = 6
k_bc = 0


# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = small.sim(N_sim, k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
indirect_result_small <- expand_grid(pi_a, OR_ab, OR_ac, tau) %>% 
  mutate(k_ab = k_ab, k_ac = k_ab, k_bc = k_ab) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
  mutate(pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )) %>% 
  mutate(pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )) %>% 
  group_by(k_ab, k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.)) %>% separate(power, c("power", "rank_correct"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(indirect_result_small, file = "indirect_result_small.RData")

