library(metafor)
library(tidyverse)
library(foreach)
library(doParallel)
library(magrittr)
library(gemtc)

set.seed(123456)

N_cores = 24
N_sim = 5000


source("power_simA.R")
source("power_simB.R")
source("power_simC.R")

pi_a <- c(0.1, 0.3, 0.5)
# pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )
# pi_c = pi_a * OR_bc / (1 - pi_a + pi_a * OR_bc )

tau <- c(0.001, 0.2, 0.4)

OR_ab <- c(1.2)
OR_ac <- c(1.4, 1.6, 1.8, 2.0)
# OR_bc = exp(log(OR_ac) - log(OR_ab))

k_ab = c(1:10*2) 


### Indirect Evidence Only

# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = power.sim_indirect(N_sim, k_ab, k_ac, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df_indirect <- expand_grid(pi_a, OR_ab, OR_ac, tau, k_ab) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>%
  mutate(pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )) %>% 
  mutate(pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )) %>% 
  mutate(k_ac = k_ab) %>% 
  group_by(k_ab, k_ac, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.))

# Stop the parallel backend
stopCluster(cl)

save(df_indirect, file = "df_indirect.RData")



### Direct Evidence Only

k_bc <- c(1:10*2) 

# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = power.sim_direct(N_sim, k_bc, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df_direct <- expand_grid(k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
  mutate(pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )) %>% 
  mutate(pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )) %>% 
  group_by(k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.))

# Stop the parallel backend
stopCluster(cl)

save(df_direct, file = "df_direct.RData")




### BNMA Evidence 

# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = power.sim.NMA(N_sim, k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df_BNMA <- expand_grid(pi_a, OR_ab, OR_ac, tau, k_ab) %>% 
  mutate(k_ac = k_ab, k_bc = k_ab) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
  mutate(pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )) %>% 
  mutate(pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )) %>% 
  group_by(k_ab, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.))

# Stop the parallel backend
stopCluster(cl)

save(df_BNMA, file = "df_BNMA.RData")
