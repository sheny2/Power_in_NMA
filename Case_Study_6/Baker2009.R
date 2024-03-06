library(pcnetmeta)
data(diabetes)
diabetes


diabetes_ab = data.frame(study = factor(diabetes$s.id), treatment = factor(diabetes$t.id), 
                         sampleSize = diabetes$n, responders = diabetes$r)

network <- mtc.network(diabetes_ab)
plot(network)




nma.networkplot(study, treatment, data = smokingcessation_ab, 
                title = "Smoking Sessation Treatments", node.col = "orange", edge.col = "gray", adjust.thick = 10,
                trtname = c("A: No intervention", "B: Self help", "C: Individual counseling", "D: Group counseling"))

nma.networkplot(study, treatment, data = smokingcessation_ab, title = "Smoking Sessation", node.col = "orange")


