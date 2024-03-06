# normal/identity: for continuous (mean difference) data.
# Required columns: [mean, std.err] or [mean, std.dev, sampleSize].


# data <- read.table(textConnection('
#   study  treatment  mean   std.dev  sampleSize
#   01     A          -1.12  0.6      15
#   01     B          -1.55  0.5      16
#   02     A          -0.8   0.7      33
#   02     B          -1.1   0.5      31'), header=TRUE)
# network <- mtc.network(data)


library(gemtc)
data(parkinson)
parkinson_rn = parkinson
colnames(parkinson_rn) = c("study", "treatment", "mean", "std.dev", "sampleSize", 
                           "treatment", "mean", "std.dev", "sampleSize",
                           "treatment", "mean", "std.dev", "sampleSize")

parkinson_ab <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(parkinson_ab) <- c("study", "treatment", "mean", "std.dev", "sampleSize")

for (i in 1:7) {
  if (!is.na(parkinson_rn[i,10])){
    treatment = as.character(parkinson_rn[i,c(2,6,10)])
    event_n = rbind(parkinson_rn[i,3:5],parkinson_rn[i,7:9],parkinson_rn[i,11:13])
    study = rep(i, length(treatment))
    parkinson_ab = rbind(parkinson_ab, cbind(study, treatment, event_n))
  } else{ 
    treatment = as.character(parkinson_rn[i,c(2,6)])
    event_n = rbind(parkinson_rn[i,3:5],parkinson_rn[i,7:9])
    study = rep(i, length(treatment))
    parkinson_ab = rbind(parkinson_ab, cbind(study, treatment, event_n))
  }
}

parkinson_net <- mtc.network(parkinson_ab)
plot(parkinson_net)


cons.model <- mtc.model(parkinson_net, type="consistency", 
                        likelihood="normal", link="identity", linearModel="random")
cons.out <- mtc.run(cons.model, n.adapt=20000, n.iter=50000, thin=1)
estimates <- summary(cons.out)

estimates$summaries$statistics[,1]



