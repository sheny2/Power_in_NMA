library(tidyverse)
library(foreach)
library(doParallel)
library(magrittr)
library(gemtc)

set.seed(123456)

N_cores = detectCores()
N_sim = 1000

source("power_sim_all.R")

load("df_indirect_BNMA_bias2.RData")

df_indirect_bias_add = df_indirect_bias[,1:9]
df_indirect_bias_add$k_ac = ifelse(df_indirect_bias_add$k_ab == 6, 12, 6)

### Indirect Evidence Only

# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>% mutate(power = power.sim.NMA4(N_sim, k_ab, k_ac, k_bc = 0, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df_indirect_bias_additional <- df_indirect_bias_add %>% 
  group_by(k_ab, k_ac, pi_a, OR_ab, OR_ac, tau) %>%
  do(power_mutate(.))%>% separate(power, c("power", "rank_correct", "avg_bias", "avg_bias_abs"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(df_indirect_bias_additional, file = "df_indirect_BNMA_bias_2_additional.RData")


