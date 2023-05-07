library(metafor)
library(tidyverse)

# Simulation goal:  The performance of the indirect inference on OR_bc

# OR_bc = exp(log(OR_ac) - log(OR_ab))
# random-effects meta-analysis is based on the DerSimonian and Laird method (DL)



# Simulation Setup

## data generation

# full
k_ab <- c(5, 10, 25, 100)
k_ac <- c(1, 5)

tau <- c(0.001, 0.2, 0.4)


# one-row simulation setup
k_ab <- 25    # the number of trials pertaining to the B versus A
k_ac <- 1     # the number of trials pertaining to the C versus A

pi_a <- 0.1   # the true average event rate in the common comparator group A

OR_ab <- 1.2  # the true relative effect of B versus A
OR_ac <- 1.4  # the true relative effect of C versus A

tau <- 0.4    # the between-study standard deviation


# power simulation function

set.seed(8848)

power.sim <- function(S = 5000, k_ab, k_ac, pi_a, OR_ab, OR_ac, tau){
  reject_correctly = 0
  
  for (i in 1:S) {
    # simulate k_ab studies
    dat_ab = data.frame()
    for (j in 1:k_ab)
    {
      n_j = runif(n = 1, min = 20, max = 500)
      n_aj = n_bj = round(n_j / 2)
      
      log_OR_ab_j = rnorm(n = 1, mean = log(OR_ab),  sd = tau)
      
      pi_aj = runif(n = 1, min = 1/2 * pi_a, max = 3/2 * pi_a)
      pi_bj = pi_aj * exp(log_OR_ab_j) / (1 - pi_aj + pi_aj * exp(log_OR_ab_j) )
      
      e_aj = rbinom(n = 1, size = n_aj, prob = pi_aj)
      e_bj = rbinom(n = 1, size = n_bj, prob = pi_bj)
      
      dat_ab = rbind(dat_ab, cbind(e_aj, n_aj, e_bj, n_bj))
    }
    
    # fit random effect MA
    result_ab = rma(measure = "OR", ai = e_aj, bi = n_aj - e_aj, ci = e_bj, di = n_bj - e_bj, 
                    data = dat_ab, method = "DL")
    est_log_OR_ab = result_ab$beta[,1]
    se_log_OR_ab = result_ab$se
    
    # simulate k_ac studies
    dat_ac = data.frame()
    for (j in 1:k_ac)
    {
      n_j = runif(n = 1, min = 20, max = 500)
      n_aj = n_cj = round(n_j / 2)
      
      log_OR_ac_j = rnorm(n = 1, mean = log(OR_ac),  sd = tau)
      
      pi_aj = runif(n = 1, min = 1/2 * pi_a, max = 3/2 * pi_a)
      pi_cj = pi_aj * exp(log_OR_ac_j) / (1 - pi_aj + pi_aj * exp(log_OR_ac_j) )
      
      e_aj = rbinom(n = 1, size = n_aj, prob = pi_aj)
      e_cj = rbinom(n = 1, size = n_cj, prob = pi_cj)  
      
      dat_ac = rbind(dat_ac, cbind(e_aj, n_aj, e_cj, n_cj))
    }
    
    # fit random effect MA
    result_ac = rma(measure = "OR", ai = e_aj, bi = n_aj - e_aj, ci = e_cj, di = n_cj - e_cj, 
                    data = dat_ac, method = "DL")
    est_log_OR_ac = result_ac$beta[,1]
    se_log_OR_ac = result_ac$se
    
    
    # Bucher's Method: indirect estimate of log(OR_bc)
    est_log_OR_bc = est_log_OR_ab - est_log_OR_ac
    se_log_OR_bc = sqrt(se_log_OR_ab^2 + se_log_OR_ac^2)
    
    # Compute 95% CI
    upper_95 = exp(est_log_OR_bc + 1.96 * se_log_OR_bc)
    lower_95 = exp(est_log_OR_bc - 1.96 * se_log_OR_bc)
    
    if (1 > upper_95 || 1 < lower_95)
      { reject_correctly = reject_correctly + 1 }
    
  }
  
  power = reject_correctly / S
  # print(paste("Power of k_ab =", k_ab, "k_ac =", k_ac, "pi_a =", pi_a, "OR_ab =", OR_ab, "OR_ac =", OR_ac, "tau =", tau, "is", power))
  
  return(power)

}



# test
power.sim(S = 1000, k_ab = 20, k_ac = 20, pi_a = 0.1, OR_ab = 1.2, OR_ac = 1.4, tau = 0.2)


########## Simulation

# possible values

k_ab <- c(5, 10, 25, 100)
k_ac <- c(1, 5)

pi_a <- c(0.1, 0.3)
# pi_a <- c(0.1, 0.3, 0.5, 0.7)

tau <- c(0.001, 0.2, 0.4)

OR_ab <- c(1.2)
OR_ac <- c(1.4, 1.6, 1.8)

# OR_ab <- c(1.8)
# OR_ac <- c(1.6, 1.4, 1.2)

k_bc <- c(2, 5, 10, 25, 50, 100)



# conduct simulations and store results
sim_result = expand_grid(k_ab, k_ac, pi_a, OR_ab, OR_ac, tau) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
  rowwise() %>% 
  mutate(power = power.sim(100, k_ab, k_ac, pi_a, OR_ab, OR_ac, tau))

a <- sim_result %>% filter(pi_a == 0.1, k_ac == 1) %>% 
  ggplot(aes(x = k_ab, y = power, color = factor(OR_bc))) + geom_line(size = 1.5) + facet_wrap(~tau) + 
  labs(title = "pi_a = 0.1, k_ac = 1", color = "True OR_bc")

b <- sim_result %>% filter(pi_a == 0.3, k_ac == 1) %>% 
  ggplot(aes(x = k_ab, y = power, color = factor(OR_bc))) + geom_line(size = 1.5) + facet_wrap(~tau) + 
  labs(title = "pi_a = 0.3, k_ac = 1", color = "True OR_bc")


c <- sim_result %>% filter(pi_a == 0.1, k_ac == 5) %>% 
  ggplot(aes(x = k_ab, y = power, color = factor(OR_bc))) + geom_line(size = 1.5) + facet_wrap(~tau) + 
  labs(title = "pi_a = 0.1, k_ac = 5", color = "True OR_bc")

d <- sim_result %>% filter(pi_a == 0.3, k_ac == 5) %>% 
  ggplot(aes(x = k_ab, y = power, color = factor(OR_bc))) + geom_line(size = 1.5) + facet_wrap(~tau) + 
  labs(title = "pi_a = 0.3, k_ac = 5", color = "True OR_bc")

gridExtra::grid.arrange(a,b,c,d)


# ignore sketch
# pi_b <- (OR_ab) * (pi_a) / (1-pi_a) / (1 + (OR_ab) * (pi_a) / (1-pi_a))
# pi_c <- (OR_ac) * (pi_a) / (1-pi_a) / (1 + (OR_ac) * (pi_a) / (1-pi_a))

# escalc(measure = "OR", ai = e_aj, bi = n_aj - e_aj, ci = e_bj, di = n_bj - e_bj, data = dat_ab, append = TRUE)



###### test parallel

# runtime
system.time(expand_grid(k_ab, k_ac, pi_a, OR_ab, OR_ac, tau) %>% 
              mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
              rowwise() %>% 
              mutate(power = power.sim(1000, k_ab, k_ac, pi_a, OR_ab, OR_ac, tau)))


library(dplyr)
library(foreach)
library(doParallel)
library(magrittr)

set.seed(12536)

# Initialize parallel backend
cl <- makeCluster(detectCores())

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = power.sim(500, k_ab, k_ac, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df <- expand_grid(k_ab, k_ac, pi_a, OR_ab, OR_ac, tau) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
  group_by(k_ab, k_ac, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.))

# Stop the parallel backend
stopCluster(cl)



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

# Save data
# sim_result_1000_1lower = df
# save(sim_result_1000_1lower, file = "sim_result_1000_1lower.RData")




# only direct

library(dplyr)
library(foreach)
library(doParallel)
library(magrittr)

set.seed(12536)

# Initialize parallel backend
cl <- makeCluster(detectCores())

# Register the parallel backend
registerDoParallel(cl)

# Define a function to be applied in parallel
power_mutate <- function(df) {
  df %>%  mutate(power = power.sim2(1000, k_ab = 1, k_ac = 1, k_bc, pi_a, OR_ab, OR_ac, tau))
}

# Apply the function in parallel using foreach
df <- expand_grid(k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  mutate(OR_bc = round(exp(log(OR_ac) - log(OR_ab)), 2)) %>% 
  group_by(k_bc, pi_a, OR_ab, OR_ac, tau) %>% 
  do(power_mutate(.))

# Stop the parallel backend
stopCluster(cl)




