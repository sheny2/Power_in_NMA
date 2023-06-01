# setwd("E:/Dropbox/Research/NMAEN/Code") # specify your working directory

library("igraph")
source("functions.R")

################################################
## Illustrative Example: Trikalinos et al., 2009
################################################
tri <- read.csv("Trikalinos.csv")
tri.dir <- dir.binary(tri)

# effective number of studies
N12 <- tri.dir$studies[1,2]
N13 <- tri.dir$studies[1,3]
N23 <- tri.dir$studies[2,3]
N34 <- tri.dir$studies[3,4]

E12 <- N12 + 1/(1/N13 + 1/N23)
E12
E13 <- N13 + 1/(1/N12 + 1/N23)
E13
E23 <- N23 + 1/(1/N12 + 1/N13)
E23
E34 <- N34
E34
E14 <- 1/(1/E13 + 1/N34)
E14
E14 <- sum(solve(matrix(c(1/N13 + 1/N34, 1/N34, 1/N34, 1/N12 + 1/N23 + 1/N34), 2, 2)))
E14
E24 <- 1/(1/E23 + 1/N34)
E24
E24 <- sum(solve(matrix(c(1/N23 + 1/N34, 1/N34, 1/N34, 1/N12 + 1/N13 + 1/N34), 2, 2)))
E24

# effective sample size
n12 <- tri.dir$sz[1,2]
n13 <- tri.dir$sz[1,3]
n23 <- tri.dir$sz[2,3]
n34 <- tri.dir$sz[3,4]

ESS12 <- n12 + 1/(1/n13 + 1/n23)
ESS12
ESS13 <- n13 + 1/(1/n12 + 1/n23)
ESS13
ESS23 <- n23 + 1/(1/n12 + 1/n13)
ESS23
ESS34 <- n34
ESS34
ESS14 <- sum(solve(matrix(c(1/n13 + 1/n34, 1/n34, 1/n34, 1/n12 + 1/n23 + 1/n34), 2, 2)))
ESS14
ESS24 <- sum(solve(matrix(c(1/n23 + 1/n34, 1/n34, 1/n34, 1/n12 + 1/n13 + 1/n34), 2, 2)))
ESS24

# effective precision
prec12 <- tri.dir$prec[1,2]
prec13 <- tri.dir$prec[1,3]
prec23 <- tri.dir$prec[2,3]
prec34 <- tri.dir$prec[3,4]

EP12 <- prec12 + 1/(1/prec13 + 1/prec23)
EP12
EP13 <- prec13 + 1/(1/prec12 + 1/prec23)
EP13
EP23 <- prec23 + 1/(1/prec12 + 1/prec13)
EP23
EP34 <- prec34
EP34
EP14 <- sum(solve(matrix(c(1/prec13 + 1/prec34, 1/prec34, 1/prec34, 1/prec12 + 1/prec23 + 1/prec34), 2, 2)))
EP14
EP24 <- sum(solve(matrix(c(1/prec23 + 1/prec34, 1/prec34, 1/prec34, 1/prec12 + 1/prec13 + 1/prec34), 2, 2)))
EP24

# alternatively, using the eff.binary() function
tri.eff <- eff.binary(tri)
tri.eff

# plot network
networkplot(tri.dir$studies, title = "(a) Number of studies of direct evidence", adjust.title = 0.5)
networkplot(tri.dir$sz, title = "(b) Sample size of direct evidence", adjust.title = 0.5)
networkplot(tri.dir$prec, title = "(c) Precision of direct evidence", adjust.title = 0.5)
wt.min.studies <- min(tri.dir$studies[tri.dir$studies > 0], na.rm = TRUE)
networkplot(tri.eff$eff.studies, wt.min = wt.min.studies, title = "(d) Effective number of studies of overall evidence", adjust.title = 0.5)
wt.min.sz <- min(tri.dir$sz[tri.dir$sz > 0], na.rm = TRUE)
networkplot(tri.eff$eff.sz, wt.min = wt.min.sz, title = "(e) Effective sample size of overall evidence", adjust.title = 0.5)
wt.min.prec <- min(tri.dir$prec[tri.dir$prec > 0], na.rm = TRUE)
networkplot(tri.eff$eff.prec, wt.min = wt.min.prec, title = "(f) Effective precision of overall evidence", adjust.title = 0.5)

###############
##More Examples
###############

wel <- read.csv("Welton.csv")
wel.dir <- dir.binary(wel)
wel.dir
wel.eff <- eff.binary(wel)
wel.eff

wel.K <- dim(wel.dir$studies)[1]
wel.compare.studies <- matrix(NA, wel.K, wel.K)
wel.compare.sz <- matrix(NA, wel.K, wel.K)
wel.compare.prec <- matrix(NA, wel.K, wel.K)
wel.compare.studies[upper.tri(wel.compare.studies)] <- wel.dir$studies[upper.tri(wel.dir$studies)]
wel.compare.studies[lower.tri(wel.compare.studies)] <- round(wel.eff$eff.studies[lower.tri(wel.eff$eff.studies)], 1)
wel.compare.studies
wel.compare.sz[upper.tri(wel.compare.sz)] <- wel.dir$sz[upper.tri(wel.dir$sz)]
wel.compare.sz[lower.tri(wel.compare.sz)] <- round(wel.eff$eff.sz[lower.tri(wel.eff$eff.sz)])
wel.compare.sz
wel.compare.prec[upper.tri(wel.compare.prec)] <- round(wel.dir$prec[upper.tri(wel.dir$prec)], 1)
wel.compare.prec[lower.tri(wel.compare.prec)] <- round(wel.eff$eff.prec[lower.tri(wel.eff$eff.prec)], 1)
wel.compare.prec

networkplot(wel.dir$studies, title = "(a) Number of studies of direct evidence", adjust.title = 0.5)
networkplot(wel.dir$sz, title = "(b) Sample size of direct evidence", adjust.title = 0.5)
networkplot(wel.dir$prec, title = "(c) Precision of direct evidence", adjust.title = 0.5)
wt.min.studies <- min(wel.dir$studies[wel.dir$studies > 0], na.rm = TRUE)
networkplot(wel.eff$eff.studies, wt.min = wt.min.studies, title = "(d) Effective number of studies of overall evidence", adjust.title = 0.5)
wt.min.sz <- min(wel.dir$sz[wel.dir$sz > 0], na.rm = TRUE)
networkplot(wel.eff$eff.sz, wt.min = wt.min.sz, title = "(e) Effective sample size of overall evidence", adjust.title = 0.5)
wt.min.prec <- min(wel.dir$prec[wel.dir$prec > 0], na.rm = TRUE)
networkplot(wel.eff$eff.prec, wt.min = wt.min.prec, title = "(f) Effective precision of overall evidence", adjust.title = 0.5)

ell <- read.csv("Elliott.csv")
ell.dir <- dir.binary(ell)
ell.dir
ell.eff <- eff.binary(ell)
ell.eff

ell.K <- dim(ell.dir$studies)[1]
ell.compare.studies <- matrix(NA, ell.K, ell.K)
ell.compare.sz <- matrix(NA, ell.K, ell.K)
ell.compare.prec <- matrix(NA, ell.K, ell.K)
ell.compare.studies[upper.tri(ell.compare.studies)] <- ell.dir$studies[upper.tri(ell.dir$studies)]
ell.compare.studies[lower.tri(ell.compare.studies)] <- round(ell.eff$eff.studies[lower.tri(ell.eff$eff.studies)], 1)
ell.compare.studies
ell.compare.sz[upper.tri(ell.compare.sz)] <- ell.dir$sz[upper.tri(ell.dir$sz)]
ell.compare.sz[lower.tri(ell.compare.sz)] <- round(ell.eff$eff.sz[lower.tri(ell.eff$eff.sz)])
ell.compare.sz
ell.compare.prec[upper.tri(ell.compare.prec)] <- round(ell.dir$prec[upper.tri(ell.dir$prec)])
ell.compare.prec[lower.tri(ell.compare.prec)] <- round(ell.eff$eff.prec[lower.tri(ell.eff$eff.prec)])
ell.compare.prec

networkplot(ell.dir$studies, title = "(a) Number of studies of direct evidence", adjust.title = 0.5)
networkplot(ell.dir$sz, title = "(b) Sample size of direct evidence", adjust.title = 0.5)
networkplot(ell.dir$prec, title = "(c) Precision of direct evidence", adjust.title = 0.5)
wt.min.studies <- min(ell.dir$studies[ell.dir$studies > 0], na.rm = TRUE)
networkplot(ell.eff$eff.studies, wt.min = wt.min.studies, title = "(d) Effective number of studies of overall evidence", adjust.title = 0.5)
wt.min.sz <- min(ell.dir$sz[ell.dir$sz > 0], na.rm = TRUE)
networkplot(ell.eff$eff.sz, wt.min = wt.min.sz, title = "(e) Effective sample size of overall evidence", adjust.title = 0.5)
wt.min.prec <- min(ell.dir$prec[ell.dir$prec > 0], na.rm = TRUE)
networkplot(ell.eff$eff.prec, wt.min = wt.min.prec, title = "(f) Effective precision of overall evidence", adjust.title = 0.5)

pic <- read.csv("Picard.csv")
pic.dir <- dir.binary(pic)
pic.dir
pic.eff <- eff.binary(pic)
pic.eff
#$eff.studies
#          [,1]      [,2]      [,3]      [,4]      [,5]      [,6]     [,7]      [,8]
#[1,]        NA 13.984789 16.481821 22.704897 11.285224 16.322689 6.624252 11.120446
#[2,] 13.984789        NA 10.236616 12.240570  8.466909 11.120001 4.840021  9.461512
#[3,] 16.481821 10.236616        NA 11.638738  8.590710  9.811129 5.804974  7.517411
#[4,] 22.704897 12.240570 11.638738        NA  8.943549 11.863220 5.357450  9.365078
#[5,] 11.285224  8.466909  8.590710  8.943549        NA  8.565217 4.881378  6.838898
#[6,] 16.322689 11.120001  9.811129 11.863220  8.565217        NA 4.956015  8.143302
#[7,]  6.624252  4.840021  5.804974  5.357450  4.881378  4.956015       NA  4.324364
#[8,] 11.120446  9.461512  7.517411  9.365078  6.838898  8.143302 4.324364        NA
#
#$eff.sz
#          [,1]     [,2]      [,3]      [,4]     [,5]      [,6]     [,7]     [,8]
#[1,]        NA 935.3156 1381.3158 2243.2158 668.5828 1209.9846 710.0824 776.8170
#[2,]  935.3156       NA  684.9081  848.0952 474.3579  765.3295 433.6763 664.8671
#[3,] 1381.3158 684.9081        NA  963.6215 527.5077  759.3789 591.3004 543.1345
#[4,] 2243.2158 848.0952  963.6215        NA 559.4914  933.6471 558.3448 699.8743
#[5,]  668.5828 474.3579  527.5077  559.4914       NA  541.8999 440.5678 412.1604
#[6,] 1209.9846 765.3295  759.3789  933.6471 541.8999        NA 479.4515 578.0132
#[7,]  710.0824 433.6763  591.3004  558.3448 440.5678  479.4515       NA 386.8009
#[8,]  776.8170 664.8671  543.1345  699.8743 412.1604  578.0132 386.8009       NA
#
#$eff.prec
#         [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
#[1,]       NA 41.29283 64.70412 91.15781 24.86405 54.63320 28.50595 29.86656
#[2,] 41.29283       NA 31.55661 37.40149 19.03413 34.80261 17.89094 28.98488
#[3,] 64.70412 31.55661       NA 42.27226 20.77004 35.71434 23.98579 22.55049
#[4,] 91.15781 37.40149 42.27226       NA 21.22582 41.54452 22.28380 27.24971
#[5,] 24.86405 19.03413 20.77004 21.22582       NA 20.69090 15.55038 15.67023
#[6,] 54.63320 34.80261 35.71434 41.54452 20.69090       NA 19.70512 24.06783
#[7,] 28.50595 17.89094 23.98579 22.28380 15.55038 19.70512       NA 15.09140
#[8,] 29.86656 28.98488 22.55049 27.24971 15.67023 24.06783 15.09140       NA

pic.K <- dim(pic.dir$studies)[1]
pic.compare.studies <- matrix(NA, pic.K, pic.K)
pic.compare.sz <- matrix(NA, pic.K, pic.K)
pic.compare.prec <- matrix(NA, pic.K, pic.K)
pic.compare.studies[upper.tri(pic.compare.studies)] <- pic.dir$studies[upper.tri(pic.dir$studies)]
pic.compare.studies[lower.tri(pic.compare.studies)] <- round(pic.eff$eff.studies[lower.tri(pic.eff$eff.studies)], 1)
pic.compare.studies
pic.compare.sz[upper.tri(pic.compare.sz)] <- pic.dir$sz[upper.tri(pic.dir$sz)]
pic.compare.sz[lower.tri(pic.compare.sz)] <- round(pic.eff$eff.sz[lower.tri(pic.eff$eff.sz)])
pic.compare.sz
pic.compare.prec[upper.tri(pic.compare.prec)] <- round(pic.dir$prec[upper.tri(pic.dir$prec)], 1)
pic.compare.prec[lower.tri(pic.compare.prec)] <- round(pic.eff$eff.prec[lower.tri(pic.eff$eff.prec)], 1)
pic.compare.prec

networkplot(pic.dir$studies, title = "(a) Number of studies of direct evidence", adjust.title = 0.5)
networkplot(pic.dir$sz, title = "(b) Sample size of direct evidence", adjust.title = 0.5)
networkplot(pic.dir$prec, title = "(c) Precision of direct evidence", adjust.title = 0.5)
wt.min.studies <- min(pic.dir$studies[pic.dir$studies > 0], na.rm = TRUE)
networkplot(pic.eff$eff.studies, wt.min = wt.min.studies, title = "(d) Effective number of studies of overall evidence", adjust.title = 0.5)
wt.min.sz <- min(pic.dir$sz[pic.dir$sz > 0], na.rm = TRUE)
networkplot(pic.eff$eff.sz, wt.min = wt.min.sz, title = "(e) Effective sample size of overall evidence", adjust.title = 0.5)
wt.min.prec <- min(pic.dir$prec[pic.dir$prec > 0], na.rm = TRUE)
networkplot(pic.eff$eff.prec, wt.min = wt.min.prec, title = "(f) Effective precision of overall evidence", adjust.title = 0.5)
