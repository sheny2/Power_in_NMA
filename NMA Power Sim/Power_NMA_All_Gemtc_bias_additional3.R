library(tidyverse)
library(foreach)
library(doParallel)
library(magrittr)
library(gemtc)

set.seed(123456)

N_cores = detectCores()
N_sim = 700

source("power_sim_all.R")

load("df_overall_BNMA_bias2.RData")

prepare0 = df_BNMA_bias[19:36,-(11:14)]
prepare0$k_ab = prepare0$k_ac = 1

prepare1 = df_BNMA_bias[19:36,-(11:14)]
prepare1$k_ab = prepare1$k_ac = 2

prepare2 = df_BNMA_bias[19:36,-(11:14)]
prepare2$k_ab = prepare2$k_ac = 3

prepare3 = df_BNMA_bias[19:36,-(11:14)]
prepare3$k_ab = prepare3$k_ac = 4


df_BNMA_bias_pre = rbind(prepare0, prepare1, prepare2, prepare3)



# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>% mutate(power = power.sim.NMA4(N_sim, k_ab, k_ac, k_bc = 0, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df_BNMA_bias_additional <- df_BNMA_bias_pre %>% 
  group_by(k_ab, k_ac, pi_a, OR_ab, OR_ac, tau) %>%
  do(power_mutate(.))%>% separate(power, c("power", "rank_correct", "avg_bias", "avg_bias_abs"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(df_BNMA_bias_additional, file = "df_overall_BNMA_bias2_additional.RData")







N_cores = detectCores()
N_sim = 500

df_BNMA_additional = df_BNMA_blank %>% filter(k_ab == 6)
df_BNMA_additional$k_ab = 2
df_BNMA_additional$k_ac = 2

# Initialize parallel backend
cl <- makeCluster(N_cores)

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel

power_mutate <- function(df) {
  df %>% mutate(power = power.sim.NMA4(N_sim, k_ab, k_ac, k_bc, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df_BNMA_bias_additional_Net5 <- df_BNMA_additional %>%
  mutate(tau = readr::parse_number(tau)) %>% 
  group_by(k_ab, k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.))%>% separate(power, c("power", "rank_correct", "avg_bias", "avg_bias_abs"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(df_BNMA_bias_additional_Net5, file = "df_BNMA_bias_additional_Net5.RData")





