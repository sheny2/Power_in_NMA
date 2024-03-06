

# one-row simulation setup
k_ab <- 25    # the number of trials pertaining to the B versus A
k_ac <- 1     # the number of trials pertaining to the C versus A

pi_a <- 0.1   # the true average event rate in the common comparator group A

OR_ab <- 1.2  # the true relative effect of B versus A
OR_ac <- 1.4  # the true relative effect of C versus A

tau <- 0.4    # the between-study standard deviation


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


# all meta-analysis was fit using odds ratio in DerSimonian and Laird ("DL") method

rma1 = rma(measure = "OR", ai = e_aj, bi = n_aj - e_aj, ci = e_bj, di = n_bj - e_bj, 
    data = dat_ab, method = "DL")

summary(rma1)

rma1$beta
rma1$se







