library(tidyverse)
library(foreach)
library(doParallel)
library(magrittr)
library(gemtc)

set.seed(123456)

N_cores = 20
N_sim = 500

source("power_sim_all.R")

# pi_a <- c(0.2, 0.4, 0.6)
pi_a <- c(0.1, 0.5)
# pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )
# pi_c = pi_a * OR_bc / (1 - pi_a + pi_a * OR_bc )

tau <- c(0.001, 0.2, 0.4)

OR_ab <- c(1.2)
OR_ac <- c(1.4, 1.6, 1.8)
# OR_bc = exp(log(OR_ac) - log(OR_ab))

k_ab = c(6,12) 
DR_INDR = c(1,2,3,6)



# expand_grid(pi_a, OR_ab, OR_ac, tau, k_ab, DR_INDR) %>%
#   mutate(k_ac = k_ab, k_bc = k_ab / DR_INDR) %>%
#   select(-DR_INDR)%>%
#   mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>%
#   mutate(pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )) %>%
#   mutate(pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )) %>% view()




### Indirect Evidence Only

# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>% mutate(power = power.sim.NMA3(k_ab, k_ac, k_bc = 0, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df_indirect_bias <- expand_grid(S = 1:N_sim, pi_a, OR_ab, OR_ac, tau, k_ab) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>%
  mutate(pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )) %>% 
  mutate(pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )) %>% 
  mutate(k_ac = k_ab) %>% 
  group_by(S, k_ab, k_ac, pi_a, OR_ab, OR_ac, tau) %>%
  do(power_mutate(.))%>% separate(power, c("power", "rank_correct", "point_est", "point_true"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(df_indirect_bias, file = "df_indirect_BNMA_bias.RData")



### Direct Evidence Only

k_bc <- c(1,2,3,6,
          2,4,6,12) %>% unique()

# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = power.sim.NMA3(k_ab = 0, k_ac = 0, k_bc, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df_direct_bias <- expand_grid(S = 1:N_sim, k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
  mutate(pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )) %>% 
  mutate(pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )) %>% 
  group_by(S, k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.))%>% separate(power, c("power", "rank_correct", "point_est", "point_true"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(df_direct_bias, file = "df_direct_BNMA_bias.RData")






### BNMA Evidence 

# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = power.sim.NMA3(k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df_BNMA_bias <- expand_grid(S = 1:N_sim, pi_a, OR_ab, OR_ac, tau, k_ab, DR_INDR) %>% 
  mutate(k_ac = k_ab, k_bc = k_ab / DR_INDR) %>% 
  select(-DR_INDR)%>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
  mutate(pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )) %>% 
  mutate(pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )) %>% 
  group_by(S, k_ab, k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.)) %>% separate(power, c("power", "rank_correct", "point_est", "point_true"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(df_BNMA_bias, file = "df_overall_BNMA_bias.RData")
