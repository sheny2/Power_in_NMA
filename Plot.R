library(metafor)
library(tidyverse)

load("sim_result_1000.RData")
df = sim_result_1000

a <- df %>% filter(pi_a == 0.1, k_ac == 1) %>% 
  ggplot(aes(x = k_ab, y = power, color = factor(OR_bc))) + geom_line(size = 1.5) + facet_wrap(~tau) + 
  labs(title = "pi_a = 0.1, k_ac = 1", color = "True OR_bc") + theme_bw()

b <- df %>% filter(pi_a == 0.3, k_ac == 1) %>% 
  ggplot(aes(x = k_ab, y = power, color = factor(OR_bc))) + geom_line(size = 1.5) + facet_wrap(~tau) + 
  labs(title = "pi_a = 0.3, k_ac = 1", color = "True OR_bc") + theme_bw()


c <- df %>% filter(pi_a == 0.1, k_ac == 5) %>% 
  ggplot(aes(x = k_ab, y = power, color = factor(OR_bc))) + geom_line(size = 1.5) + facet_wrap(~tau) + 
  labs(title = "pi_a = 0.1, k_ac = 5", color = "True OR_bc") + theme_bw()

d <- df %>% filter(pi_a == 0.3, k_ac == 5) %>% 
  ggplot(aes(x = k_ab, y = power, color = factor(OR_bc))) + geom_line(size = 1.5) + facet_wrap(~tau) + 
  labs(title = "pi_a = 0.3, k_ac = 5", color = "True OR_bc") + theme_bw()


gridExtra::grid.arrange(a,b,c,d)
