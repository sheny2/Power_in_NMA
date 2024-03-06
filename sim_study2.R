library(tidyverse)
library(dplyr)
library(foreach)
library(doParallel)
library(magrittr)
library(gemtc)



pi_a <- c(0.1, 0.3)

tau <- c(0.001, 0.2, 0.4)

OR_ab <- c(1.2)
OR_ac <- c(1.4, 1.6, 1.8)

k_ab <- c(1:10*2)
# k_ac <- c(1:10*2)
# k_bc <- c(1:10*2)


# expand_grid(pi_a, OR_ab, OR_ac, tau, k_ab) %>% 
#   mutate(k_ac = k_ab, k_bc = k_ab) %>% 
#   mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
#   rowwise() %>% 
#   mutate(power = power.sim2.foreach(1000, k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, tau))



set.seed(12536)

# Initialize parallel backend
cl <- makeCluster(7)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = power.sim2.foreach(10, k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df <- expand_grid(pi_a, OR_ab, OR_ac, tau, k_ab) %>% 
  mutate(k_ac = k_ab, k_bc = k_ab) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
  group_by(k_ab, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.))

# Stop the parallel backend
stopCluster(cl)


