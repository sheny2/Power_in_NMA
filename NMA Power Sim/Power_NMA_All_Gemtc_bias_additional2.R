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

df_indirect_bias_example = df_indirect_bias %>% filter(k_ab == 12)

pre0 = df_indirect_bias_example[1:9]
pre0$k_ab = 2; pre0$k_ac = 2

pre1 = df_indirect_bias_example[1:9]
pre1$k_ab = 3; pre1$k_ac = 3

pre2 = df_indirect_bias_example[1:9]
pre2$k_ab = 4; pre2$k_ac = 4

df_indirect_bias_newpre = rbind(pre0, pre1, pre2)


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
df_indirect_bias_additional2 <- df_indirect_bias_newpre %>% 
  group_by(k_ab, k_ac, pi_a, OR_ab, OR_ac, tau) %>%
  do(power_mutate(.))%>% separate(power, c("power", "rank_correct", "avg_bias", "avg_bias_abs"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(df_indirect_bias_additional2, file = "df_indirect_BNMA_bias_2_additional2.RData")







