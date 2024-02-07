library(metafor)
library(tidyverse)
library(foreach)
library(doParallel)
library(magrittr)
library(gemtc)

set.seed(123456)

N_cores = detectCores() 
N_sim = 400


load("df_indirect_blank.RData")
load("df_direct_blank.RData")
load("df_BNMA_blank.RData")
load("df_indirect_bias_additional_blank.RData")


source("power_sim_AB.R")


### Indirect Evidence Only

# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = power.sim.AB_full(N_sim, k_ab, k_ac, k_bc = 0, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df_indirect_AB <- df_indirect_blank %>% 
  mutate(tau = readr::parse_number(tau)) %>% 
  group_by(k_ab, k_ac, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.)) %>% separate(power, c("power", "avg_bias_abs"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(df_indirect_AB, file = "df_indirect_AB.RData")




# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = power.sim.AB_full(N_sim, k_ab, k_ac, k_bc = 0, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df_indirect_additional_AB <- df_indirect_bias_additional_blank %>% 
  mutate(tau = readr::parse_number(tau)) %>% 
  group_by(k_ab, k_ac, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.)) %>% separate(power, c("power", "avg_bias_abs"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(df_indirect_additional_AB, file = "df_indirect_additional_AB.RData")





### Direct Evidence Only


# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = power.sim.AB_direct(N_sim, k_ab = 0, k_ac = 0, k_bc, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df_direct_AB <- df_direct_blank %>% 
  mutate(tau = readr::parse_number(tau)) %>% 
  group_by(k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.)) %>% separate(power, c("power", "avg_bias_abs"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(df_direct_AB, file = "df_direct_AB.RData")





### BNMA Evidence 

# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = power.sim.AB_full(N_sim, k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df_BNMA_AB <- df_BNMA_blank %>%
  # mutate(tau = readr::parse_number(tau)) %>% 
  group_by(k_ab, k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.)) %>% separate(power, c("power", "avg_bias_abs"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(df_BNMA_AB, file = "df_BNMA_AB.RData")



