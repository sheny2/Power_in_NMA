library(tidyverse)
library(INLA)

library(nmaINLA)
data("Smokdat", package = "nmaINLA")
head(Smokdat)


SmokdatINLA <- create_INLA_dat(dat = Smokdat,
                               armVars = c('treatment' = 't', 'responders' = 'r',
                                           'sampleSize' = 'n'),
                               nArmsVar = 'na',
                               design = 'des')
head(SmokdatINLA)

plot_nma(s.id = study, t.id = treatment, data = SmokdatINLA)

fit.consistency <- nma_inla(SmokdatINLA, likelihood = "binomial",
                            fixed.par = c(0, 1000), tau.prior = "uniform",
                            tau.par = c(0, 5), type = "consistency")



d12.inla <- inla.smarginal(marginal = fit.consistency$marginals.fixed$d12)
plot(d12.inla, type = "l", xlab = expression(paste(d[12])), ylab = " ")

log.prec.het <- fit.consistency$internal.marginals.hyperpar$`Log precision for het`
tau2.inla <- inla.tmarginal(function(x) 1/exp(x), log.prec.het, n = 20000)
plot(tau2.inla, type = "l", xlab = expression(paste(tau)), ylab = " ")

fit.jackson <- nma_inla(SmokdatINLA, likelihood = "binomial",
                        fixed.par = c(0, 1000), tau.prior = "uniform",
                        tau.par = c(0, 5), kappa.prior = "uniform",
                        kappa.par = c(0, 5), type = "jackson")



data("Strokedat", package = "nmaINLA")
# deleting 13th study
Strokedat.mreg <- Strokedat[-c(13),]
# centering the covariate
Strokedat.mreg$age <- Strokedat.mreg$age - mean(Strokedat.mreg$age)
# data preparation for INLA
StrokedatINLA.mreg <- create_INLA_dat(dat = Strokedat.mreg,
                                      armVars = c('treatment' = 't','responders' = 'r',
                                                  'sampleSize' = 'n'),
                                      nArmsVar = 'na',
                                      design = 'des',
                                      covariate = 'age')


fit.Stroke.CONS.MREG.INLA <- nma_inla(StrokedatINLA.mreg, likelihood = "binomial",
                                      fixed.par = c(0, 1000), tau.prior = "uniform",
                                      tau.par = c(0, 2), type = 'consistency',
                                      mreg = TRUE)






library(gemtc)
bleed <- read.table("../Data_Bleeding.csv", sep=",", header=TRUE)
mace <- read.table("../Data_MACE.csv", sep=",", header=TRUE)

bleed[7,3:7] = bleed[7,3:7] + bleed[8,3:7]
bleed <- bleed[-8,]
mace[7,3:10] = mace[7,3:10] + mace[8,3:10]
mace <- mace[-8,]

trt <- c("A", "B", "C", "D")
trts <- read.table(textConnection('
                                  id description
                                  A "VKA + P2Y12 + Asprin"
                                  B "VKA + P2Y12"
                                  C "NOAC + P2Y12 + Asprin"
                                  D "NOAC + P2Y12"'), header=TRUE)


outcomename = "TIMImajor"
data <- bleed[,c("study", "treatment", "n", outcomename)]
data$treatment <- trt[data$treatment]
colnames(data) = c("study","treatment", "sampleSize","responders")
network <- mtc.network(data.ab=data, treatments=trts)
plot(network)
summary(network)

## Run a random effects NMA model under consistency wity binary outcome
cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random")
cons.out <- mtc.run(cons.model, n.adapt=20000, n.iter=50000, thin=1)
summary(cons.out)
## Rank probability
prob <- rank.probability(cons.out, preferredDirection=-1)
prob <- round(prob, digits=3)



datINLA = create_INLA_dat(dat = data %>% group_by(study) %>% mutate(na = n()) %>% ungroup(),
                          nArmsVar = 'na')
datINLA$treatment = data$treatment
datINLA$responders = data$responders
datINLA$sampleSize = data$sampleSize

# datINLA = data %>% group_by(study) %>% mutate(na = n()) %>% group_by(study) %>% mutate(g = 1:mean(na)) %>% mutate(g = ifelse(treatment == "A", NA, g-1)) %>% ungroup()

plot_nma(s.id = study, t.id = treatment, data = datINLA)

fit.consistency <- nma_inla(datINLA, likelihood = "binomial",type = "consistency")


# data %>% group_by(study) %>% mutate(na = n()) %>% group_by(study) %>% mutate(g = 1:mean(na)) %>% mutate(g = ifelse(treatment == "A", NA, g-1)) %>% ungroup()
