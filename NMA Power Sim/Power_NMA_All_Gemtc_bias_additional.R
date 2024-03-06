library(tidyverse)
library(foreach)
library(doParallel)
library(magrittr)
library(gemtc)

set.seed(123456)

N_cores = detectCores()
N_sim = 500

source("power_sim_all.R")

load("df_indirect_BNMA_bias2.RData")

df_indirect_bias_add1 = df_indirect_bias[1:18,1:9]
df_indirect_bias_add1$k_ab = 2
df_indirect_bias_add1$k_ac = 2

df_indirect_bias_add2 = df_indirect_bias[1:18,1:9]
df_indirect_bias_add2$k_ab = 3
df_indirect_bias_add2$k_ac = 2

df_indirect_bias_add3 = df_indirect_bias[1:18,1:9]
df_indirect_bias_add3$k_ab = 6
df_indirect_bias_add3$k_ac = 2

df_indirect_bias_add4 = df_indirect_bias[1:18,1:9]
df_indirect_bias_add4$k_ab = 12
df_indirect_bias_add4$k_ac = 2

df_indirect_bias_add5 = df_indirect_bias[1:18,1:9]
df_indirect_bias_add5$k_ab = 2
df_indirect_bias_add5$k_ac = 3

df_indirect_bias_add6 = df_indirect_bias[1:18,1:9]
df_indirect_bias_add6$k_ab = 2
df_indirect_bias_add6$k_ac = 6

df_indirect_bias_add7 = df_indirect_bias[1:18,1:9]
df_indirect_bias_add7$k_ab = 2
df_indirect_bias_add7$k_ac = 12

df_indirect_bias_add8 = df_indirect_bias[1:18,1:9]
df_indirect_bias_add8$k_ab = 12
df_indirect_bias_add8$k_ac = 12


df_indirect_bias_add8 = df_indirect_bias[1:18,1:9]
df_indirect_bias_add8$k_ab = 2
df_indirect_bias_add8$k_ac = 1


df_indirect_bias_add_new = rbind(df_indirect_bias_add1, df_indirect_bias_add2, df_indirect_bias_add3, 
                                 df_indirect_bias_add4, df_indirect_bias_add5, df_indirect_bias_add6,
                                 df_indirect_bias_add7, df_indirect_bias_add8)



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
df_indirect_BNMA_bias_2_additional_new <- df_indirect_bias_add_new %>% 
  group_by(k_ab, k_ac, pi_a, OR_ab, OR_ac, tau) %>%
  do(power_mutate(.))%>% separate(power, c("power", "rank_correct", "avg_bias", "avg_bias_abs"), " ", convert = TRUE)

# Stop the parallel backend
stopCluster(cl)

save(df_indirect_BNMA_bias_2_additional_new, file = "df_indirect_BNMA_bias_2_additional_new.RData")





