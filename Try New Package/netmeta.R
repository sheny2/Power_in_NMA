library(netmeta)


# R package netmeta (Balduzzi et al., 2023) provides frequentist methods for network meta-analysis


data(smokingcessation)

smokingcessation

p1 <- pairwise(list(treat1, treat2, treat3),
               event = list(event1, event2, event3), n = list(n1, n2, n3),
               data = smokingcessation, sm = "OR")

net1 <- netmeta(p1, common = FALSE)


data(Baker2009)
Baker2009



data(Senn2013)
# Only consider first five studies (to reduce runtime of example)
#
studies <- unique(Senn2013$studlab)
Senn2013.5 <- subset(Senn2013, studlab %in% studies[1:5])
# Conduct network meta-analysis with placebo as reference treatment
#
net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data = Senn2013.5, sm = "MD", reference = "plac")
# Decomposition of Cochran's Q
#
decomp.design(net1)




# Artificial dataset
#
t1 <- c("A + B", "A + C", "A"    , "A"    , "D", "D", "E")
t2 <- c("C"    , "B"    , "B + C", "A + D", "E", "F", "F")
#
mean    <- c(4.1, 2.05, 0, 0, 0.1, 0.1, 0.05)
se.mean <- rep(0.1, 7)
#
study <- paste("study", c(1:4, 5, 5, 5))
#
dat <- data.frame(mean, se.mean, t1, t2, study,
                  stringsAsFactors = FALSE)
#
trts <- c("A", "A + B", "A + C", "A + D",
          "B", "B + C", "C", "D", "E", "F")
#
comps <- LETTERS[1:6]
# Use netconnection() to display network information
#
netconnection(t1, t2, study)
dc1 <- discomb(mean, se.mean, t1, t2, study, seq = trts)
dc1
forest(dc1, ref = "F")
# Define C matrix manually (which will produce the same results) #
C<-rbind(c(1,0,0,0,0,0), #A
         c(1,1,0,0,0,0), #A+B 
         c(1,0,1,0,0,0), #A+C
         c(1,0,0,1,0,0), #A+D 
         c(0,1,0,0,0,0), #B 
         c(0,1,1,0,0,0), #B+C 
         c(0,0,1,0,0,0), #C 
         c(0,0,0,1,0,0), #D 
         c(0,0,0,0,1,0), #E 
         c(0,0,0,0,0,1)) #F
         #
colnames(C) <- comps
rownames(C) <- trts
#
dc2 <- discomb(mean, se.mean, t1, t2, study, seq = trts,
              C.matrix = C)
#
# Compare C matrices
#
all.equal(dc1$C.matrix, dc2$C.matrix)






data(Linde2016)
# Only consider studies including Face-to-face PST (to reduce
# runtime of example)
#
face <- subset(Linde2016, id %in% c(16, 24, 49, 118))
# Conduct random effects network meta-analysis
#
net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
                data = face, ref = "placebo", sm = "OR", common = FALSE)
# Additive model for treatment components (with placebo as inactive
# treatment)
#
nc1 <- netcomb(net1, inactive = "placebo")

# Some complex interventions
#
ints <- c("F + TCA", "F + Plac", "SSRI + Plac + TCA")
netcomplex(nc1, ints)
#
forest(netcomplex(nc1, ints))
forest(netcomplex(nc1, ints), nchar.comps = 4)
# Component effects
#
forest(netcomplex(nc1, nc1$comps))





data(Linde2016)
# Only consider studies including Face-to-face PST (to reduce
# runtime of example)
#
face <- subset(Linde2016, id %in% c(16, 24, 49, 118))
# Conduct random effects network meta-analysis
#
net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
                data = face, ref = "placebo", sm = "OR", common = FALSE)
# Additive model for treatment components (with placebo as inactive
# treatment)
#
nc1 <- netcomb(net1, inactive = "placebo")
# Some comparisons
#
t1 <- c("F + TCA", "F + Plac", "SSRI + Plac + TCA")
t2 <- c("UC", "Plac", "UC")
#
netcomparison(nc1, t1, t2)
#
forest(netcomparison(nc1, t1, t2))
forest(netcomparison(nc1, t1, t2), nchar.comps = 4)
forest(netcomparison(nc1, c("F", "TCA"), "UC"), nchar.comps = 4)



data(Franchini2012)
Franchini2012



data(Gurusamy2011)
# Only consider three studies (to reduce runtime of example)
#
studies <- c("Findlay 2001", "Garcia-Huete 1997", "Dalmau 2000")
three <- subset(Gurusamy2011, study %in% studies)
# Transform data from long arm-based format to contrast-based
# format. Argument 'sm' has to be used for odds ratio as summary
# measure; by default the risk ratio is used in the metabin
# function called internally.
#
p1 <- pairwise(treatment, death, n, studlab = study,
               data = three, sm = "OR")
# Conduct Mantel-Haenszel network meta-analysis
#
netmetabin(p1, ref = "cont")
## Not run:
p2 <- pairwise(treatment, death, n, studlab = study,
               data = Gurusamy2011, sm = "OR")
# Conduct Mantel-Haenszel network meta-analysis
netmetabin(p2, ref = "cont")
## End(Not run)




# use netmeta
data(smokingcessation)
# Transform data from arm-based format to contrast-based format
#
p1 <- pairwise(list(treat1, treat2, treat3),
               event = list(event1, event2, event3), n = list(n1, n2, n3),
               data = smokingcessation, sm = "OR")
# Conduct random effects network meta-analysis
#
net1 <- netmeta(p1, common = FALSE)
net1

## Not run:
data(Senn2013)
# Conduct common effects network meta-analysis
#
net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data = Senn2013, sm = "MD", random = FALSE)
net2
net2$Q.decomp
# Comparison with reference group
#
print(net2, reference = "plac")
# Conduct random effects network meta-analysis
#
net3 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data = Senn2013, sm = "MD", common = FALSE)
net3
# Change printing order of treatments with placebo last and use
# long treatment names
#
trts <- c("acar", "benf", "metf", "migl", "piog",
          "rosi", "sita", "sulf", "vild", "plac")
net4 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
                data = Senn2013, sm = "MD", common = FALSE,
                seq = trts, reference = "Placebo")
print(net4, digits = 2)
## End(Not run)
