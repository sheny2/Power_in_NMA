library(dmetar)
library(tidyverse)

# Define assumptions
theta <- 0.2
K <- 10
n1 <- 25
n2 <- 25
# Calculate pooled effect standard error
sigma <- sqrt(((n1+n2)/(n1*n2)+(theta^2/(2*n1+n2)))/K)
# Calculate z
z = theta/sigma
# Calculate the power
1 - pnorm(1.96-z) + pnorm(-1.96-z)


power.analysis(d = 0.2,
               k = 10,
               n1 = 25,
               n2 = 25,
               p = 0.05)


power.analysis(d = 0.2,
               k = 10,
               n1 = 10,
               n2 = 25,
               p = 0.05,
               heterogeneity = "moderate")




success <- 100:400

plot(success, dbinom(success, size=500, prob=.5),type='h')


plot(success, dbinom(success, size=500, prob=.6),type='h')

qbinom(0.00001,
       size = 500,
       prob = 0.5,
       lower.tail = F)





