library(metafor)
library(tidyverse)
library(foreach)
library(doParallel)
library(magrittr)
library(gemtc)

set.seed(123456)

N_cores = detectCores() 
N_sim = 50

# source("power_sim_AB_bias.R")
# source("power_sim_AB_bias.R")
source("power_sim_AB_bias3.R")







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
df_indirect_new <- expand_grid(pi_a, OR_ab, OR_ac, tau, k_ab) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>%
  mutate(pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )) %>% 
  mutate(pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )) %>% 
  mutate(k_ac = k_ab) %>% 
  group_by(k_ab, k_ac, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.))%>% separate(power, c("power_het_cor","power_het_eqcor", "power_hom_eqcor",
                                           "avg_bias_het_cor", "avg_bias_het_eqcor", "avg_bias_hom_eqcor"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(df_indirect_new, file = "df_indirect_BNMA_AB_bias3.RData")





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
df_direct_new <- expand_grid(k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
  mutate(pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )) %>% 
  mutate(pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )) %>% 
  group_by(k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.))%>% separate(power, c("power_het_cor","power_het_eqcor", "power_hom_eqcor",
                                           "avg_bias_het_cor", "avg_bias_het_eqcor", "avg_bias_hom_eqcor"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(df_direct_new, file = "df_direct_BNMA_AB_bias3.RData")







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
df_BNMA_new <- expand_grid(pi_a, OR_ab, OR_ac, tau, k_ab, DR_INDR) %>% 
  mutate(k_ac = k_ab, k_bc = k_ab / DR_INDR) %>% 
  select(-DR_INDR)%>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
  mutate(pi_b = pi_a * OR_ab / (1 - pi_a + pi_a * OR_ab )) %>% 
  mutate(pi_c = pi_a * OR_ac / (1 - pi_a + pi_a * OR_ac )) %>% 
  group_by(k_ab, k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.)) %>% separate(power, c("power_het_cor","power_het_eqcor", "power_hom_eqcor",
                                            "avg_bias_het_cor", "avg_bias_het_eqcor", "avg_bias_hom_eqcor"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(df_BNMA_new, file = "df_overall_BNMA_AB_bias3_100.RData")



