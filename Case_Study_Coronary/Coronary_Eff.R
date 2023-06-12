# library(foreach)
# library(doParallel)
# # Initialize parallel backend
# cl <- makeCluster(12)
# 
# # Register the parallel backend
# registerDoParallel(cl)



library(tidyverse)
# library(gemtc)


set.seed(2023)

######## Import data
bleed <- read.table("Data_Bleeding.csv", sep=",", header=TRUE)
mace <- read.table("Data_MACE.csv", sep=",", header=TRUE)


######## Row 7 and 8 need to be combined
bleed[7,3:7] = bleed[7,3:7] + bleed[8,3:7]
bleed <- bleed[-8,]
mace[7,3:10] = mace[7,3:10] + mace[8,3:10]
mace <- mace[-8,]


Bleed_data <- bleed[,c("study", "treatment", "n", "TIMImajor")]
Bleed_data = Bleed_data[1:(nrow(Bleed_data)-2), ]

MACE_data <- mace[,c("study", "treatment", "n", "MACE")]        
MACE_data = MACE_data[1:(nrow(MACE_data)-2), ]


## Note: A is reference

trt <- c("A","B", "C","D")
# trt <- c("VKA + DAPT", "VKA + P2Y12", "NOAC + DAPT", "NOAC + P2Y12")

trts <- read.table(textConnection('id description
                                  A "VKA + DAPT"
                                  B "VKA + P2Y12"
                                  C "NOAC + DAPT"
                                  D "NOAC + P2Y12"'), header=TRUE)



# Bleed_data$treatment <- trt[Bleed_data$treatment]
colnames(Bleed_data) = c("study","treatment", "sampleSize","responders")
# MACE_data$treatment <- trt[MACE_data$treatment]
colnames(MACE_data) = c("study","treatment", "sampleSize","responders")



Bleed_data_eff = data.frame(sid = Bleed_data$study, tid = Bleed_data$treatment, r = Bleed_data$responders, n = Bleed_data$sampleSize)
MACE_data_eff = data.frame(sid = MACE_data$study, tid = MACE_data$treatment, r = MACE_data$responders, n = MACE_data$sampleSize)

# Bleed_data_eff
# MACE_data_eff


library("igraph")
source("functions.R")


dir.binary(Bleed_data_eff)
eff.binary(Bleed_data_eff)

dir.binary(MACE_data_eff)
eff.binary(MACE_data_eff)


# stargazer::stargazer(eff.binary(Bleed_data_eff))

