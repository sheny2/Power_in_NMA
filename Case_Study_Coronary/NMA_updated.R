######## NMA AF ACS (updated)
######## ENTRUST-AF is added
######## Last update by Hwanhee Hong: Sep 29, 2019
######## Data were checked again; some changes below
######## PIONEER: MI counts were updated
######## PIONEER: All cause death used safety outcome population
######## AUGUSTUS: The 4 safety outcomes were updated
######## 
######## Set your working directory
setwd("/Users/hh190/Dropbox/2.Hwanhee_Duke/02.Collaboration/1.RenatoLopes/9.UpdatedNMA")
set.seed(2918)

library(gemtc)
library(fields)
source("ifplot.fun.R")
source("BayesDiagnos.fun.R")
source("ranko_sucra.R")

######## Import data
bleed <- read.table("Data_Bleeding.csv", sep=",", header=TRUE)
mace <- read.table("Data_MACE.csv", sep=",", header=TRUE)

######## Row 7 and 8 need to be combined
bleed[7,3:7] = bleed[7,3:7] + bleed[8,3:7]
bleed <- bleed[-8,]
mace[7,3:10] = mace[7,3:10] + mace[8,3:10]
mace <- mace[-8,]

######## Function: NMA under the consistency assumption  
fun.NMA.consistency <- function(data, trt, trts, NT, outcomename, figlab, color) {
  
  ## Structure data
  if(outcomename=="AllDeath") {
    data <- data[,c("study", "treatment", "n_AllDeath", outcomename)]
  } else {
    data <- data[,c("study", "treatment", "n", outcomename)]
  }
  data$treatment <- trt[data$treatment]
  colnames(data) = c("study","treatment", "sampleSize","responders")
  network <- mtc.network(data.ab=data, treatments=trts)
  plot(network)
  summary(network)
  
  ## Run a random effects NMA model under consistency wity binary outcome
  cons.model <- mtc.model(network, type="consistency", likelihood="binom", link="logit", linearModel="random")
  cons.out <- mtc.run(cons.model, n.adapt=20000, n.iter=50000, thin=1)
  ## Rank probability
  prob <- rank.probability(cons.out, preferredDirection=-1)
  prob <- round(prob, digits=3)
  ## Results
  sink(paste("Output/",outcomename,"_res.txt",sep=""), append=FALSE)
  summary(cons.out)
  cat("#### d.G.A = A - G", "\n")
  cat("#### relative.effect shows column - row format for a league table", "\n")
  summary(relative.effect(cons.out,"A",trt))
  summary(relative.effect(cons.out,"B",trt))
  summary(relative.effect(cons.out,"C",trt))
  summary(relative.effect(cons.out,"D",trt))
  cat("#### Rank probability", "\n")
  prob
  sink()
  
  ## Diagnostic test
  pdf(paste("Output/",outcomename,"_diag.pdf",sep=""))
  gelman.plot(cons.out)
  plot(cons.out)
  dev.off()
  
  ## League table
  # Make league dataset
  league <- data.frame()
  ind <- c("A","B","C","D")
  for(i in 1:4){
    a <- summary(relative.effect(cons.out,ind[i],trt))$summaries
    temp <- as.data.frame(a$statistics[1:4,])
    temp$upper <- a$quantiles[1:4,5]
    temp$lower <- a$quantiles[1:4,1]
    league <- rbind(league, temp[,c("Mean","upper","lower")])
  }
  league <- round(league,digits = 2)
  # Make table
  head(league)
  league.table <- NULL
  for (i in 1:4) {
    res <- round(exp(league$Mean[(4*(i-1)+1):(4*i)]), digits = 2)
    ci <- paste("(",round(exp(league$lower[(4*(i-1)+1):(4*i)]),digits = 2),";",round(exp(league$upper[(4*(i-1)+1):(4*i)]),digits = 2),")", sep="")
    league.table <- rbind(league.table,res,ci)
  }
  write.table(league.table, paste("Output/",outcomename,"_league.csv",sep=""), sep=",", quote = FALSE, row.names = F)
  
  #### Forest plot
  y <- c(3.5,2.5,1.5,0.5)
  x = round(exp(league$Mean[1:4]),digits = 3)
  
  pdf(paste("Output/",outcomename,"_forest.pdf",sep=""), width=4, height=2.8)
  par(mfrow=c(1,1), mar=c(4,7,1,0), mgp=c(3, 0.1, 0))
  plot(NULL, xlim=c(0,3), ylim=c(0,4), axes=FALSE, xlab="", ylab="")
  lines(c(1,1),c(-0.3,4), col=1,cex.axis = 0.65)
  points(x, y, cex=0.7, pch=16, col=color)
  #VKA+P2Y12+Aspirin: reference, no 95% CI
  #VKA+P2Y12
  lines(c(exp(league$lower[2]),exp(league$upper[2])),c(2.5,2.5),lty=1, col=color)
  #NOAC+P2Y12+Aspirin
  lines(c(exp(league$lower[3]),exp(league$upper[3])),c(1.5,1.5),lty=1, col=color)
  #NOAC+P2Y12
  lines(c(exp(league$lower[4]),exp(league$upper[4])),c(0.5,0.5),lty=1, col=color)
  #
  axis(2, at=y, 
       lab=c(expression("VKA + DAPT"),
             expression("VKA + P2Y"[12]*" inhibitor"),
             expression("NOAC + DAPT"),
             expression("NOAC + P2Y"[12]*" inhibitor")),
       tck=0, las=2, cex.axis = 0.65, col.ticks="white")
  axis(2, at=3.2, lab="(Reference)", tck=0, las=2, cex.axis = 0.65, col.ticks="white")
  axis(1, at=c(0,1,2,3), lab=c(0,1,2,3), cex.axis=0.8, tck=-0.03)
  mtext("Favors non-reference strategy", side=1, line=1, at=0.8, cex=0.6, adj=1)
  mtext("Favors reference", side=1, line=1, at=2, cex=0.6, adj=0)
  mtext(paste("Odds Ratio for ",figlab,sep=""), side=1, line=2, at=1, cex=0.7)
  dev.off()

  ## SUCRA
  cumprob.data <- apply(t(prob), 2, cumsum)
  # Calculate SUCRA
  sucra.res <- c()
  for(i in 1:NT) {
    sucra <- (sum(cumprob.data[1:(NT-1),i])/(NT-1))*100
    sucra.res <- cbind(sucra.res,sucra) 
  }
  colnames(sucra.res) <- trts[,2]
  sucra.res <- round(sucra.res, digits=1)
  write.table(sucra.res, paste("Output/",outcomename,"_sucra.csv",sep=""), sep=",", quote = FALSE, row.names = FALSE)
  
  ## Save results
  RESULTS <- NULL
  RESULTS[[1]] <- summary(cons.out)
  RESULTS[[2]] <- league
  RESULTS[[3]] <- prob
  RESULTS[[4]] <- sucra.res
  names(RESULTS) <- c("Summary","League table","Rank probability","SUCRA")
  
  return(RESULTS)
}

######## Run models

## Treatment coding
# 1	VKA + P2Y12 + Asprin
# 2	VKA + P2Y12
# 3	NOAC + P2Y12 + Asprin
# 4	NOAC + P2Y12
trt <- c("A", "B", "C", "D")
trts <- read.table(textConnection('
                                  id description
                                  A "VKA + P2Y12 + Asprin"
                                  B "VKA + P2Y12"
                                  C "NOAC + P2Y12 + Asprin"
                                  D "NOAC + P2Y12"'), header=TRUE)
NT <- 4 ## number of treatments

#### First, check the JAMA NMA for MI; MI counts were 17/21/14 and now 19/17/21
JAMA.MI <- fun.NMA.consistency(data=mace[1:11,], trt=trt, trts=trts, NT=NT, outcomename="MI", figlab="Myocardial Infarction", color=4)
JAMA.AllDeath <- fun.NMA.consistency(data=mace[1:11,], trt=trt, trts=trts, NT=NT, outcomename="AllDeath", figlab="All-Cause Death", color=4)

#### Bleeding outcomes (safety)
TIMImajor.res <- fun.NMA.consistency(data=bleed, trt=trt, trts=trts, NT=NT, outcomename="TIMImajor", figlab="TIMI Major Bleeding", color=2)
TIMImajorminor.res <- fun.NMA.consistency(data=bleed, trt=trt, trts=trts, NT=NT, outcomename="TIMImajorminor", figlab="TIMI Major or Minor Bleeding", color=2)
TrialPrimary.res <- fun.NMA.consistency(data=bleed, trt=trt, trts=trts, NT=NT, outcomename="TrialPrimary", figlab="Trial-Defined Primary Safety Outcome", color=2)
ICH.res <- fun.NMA.consistency(data=bleed, trt=trt, trts=trts, NT=NT, outcomename="ICH", figlab="Intracranial Hemorrhage", color=2)

#### MACE outcomes (efficacy)
MACE.res <- fun.NMA.consistency(data=mace, trt=trt, trts=trts, NT=NT, outcomename="MACE", figlab="Trial-Defined Primary MACE Outcome", color=4)
AllDeath.res <- fun.NMA.consistency(data=mace, trt=trt, trts=trts, NT=NT, outcomename="AllDeath", figlab="All-Cause Death", color=4)
CardioDeath.res <- fun.NMA.consistency(data=mace, trt=trt, trts=trts, NT=NT, outcomename="CardioDeath", figlab="Cardiovascular Death", color=4)
MI.res <- fun.NMA.consistency(data=mace, trt=trt, trts=trts, NT=NT, outcomename="MI", figlab="Myocardial Infarction", color=4)
Stroke.res <- fun.NMA.consistency(data=mace, trt=trt, trts=trts, NT=NT, outcomename="Stroke", figlab="Stroke", color=4)
ST.res <- fun.NMA.consistency(data=mace, trt=trt, trts=trts, NT=NT, outcomename="ST", figlab="Stent Thrombosis", color=4)

save.image(file="NMA_updated.RData")

######## Table: SUCRA

tab.SUCRA <- rbind(TIMImajor.res$SUCRA, TIMImajorminor.res$SUCRA, TrialPrimary.res$SUCRA, ICH.res$SUCRA,
                   MACE.res$SUCRA, AllDeath.res$SUCRA, CardioDeath.res$SUCRA, MI.res$SUCRA, Stroke.res$SUCRA, ST.res$SUCRA)
row.names(tab.SUCRA) <- c("TIMImajor", "TIMImajorminor", "TrialPrimary", "ICH",
                          "MACE", "AllDeath", "CardioDeath", "MI", "Stroke", "ST")
write.table(tab.SUCRA, "Output/SUCRA_table.csv", sep=",", quote = FALSE, row.names = T)

######## Plot: OR for TIMI major bleeding and MACE

league <- TIMImajor.res[[2]]
league_eff <- MACE.res[[2]] 
x1 = round(exp(league$Mean[1:4]),digits = 3)
x2 = round(exp(league_eff$Mean[1:4]),digits = 3)

pdf("Output/ORsBleedxMACE.pdf", width=4, height=4)
par(mfrow=c(1,1), mar=c(3,4,2,1), mgp=c(1, 0.2, 0))
par(cex.axis=0.4, cex.lab=0.4)
plot(x1,x2,pch=c(1,17,15,18), xlim = c(0,1.6), ylim = c(0,1.6), xlab = "Odds Ratio for TIMI Major Bleeding", ylab = "Odds Ratio for MACE", axes = FALSE, cex=0.5)

#VKA+P2Y12
lines(c(exp(league$lower[2]),exp(league$upper[2])), c(x2[2],x2[2]), col=2, lty=1)
lines(c(x1[2],x1[2]),c(exp(league_eff$lower[2]),exp(league_eff$upper[2])),lty=1, col=4)

#NOAC+P2Y12+Aspirin
lines(c(exp(league$lower[3]),exp(league$upper[3])), c(x2[3],x2[3]), col=2, lty=1)
lines(c(x1[3],x1[3]),c(exp(league_eff$lower[3]),exp(league_eff$upper[3])),lty=1, col=4)

#NOAC+P2Y12
lines(c(exp(league$lower[4]),exp(league$upper[4])), c(x2[4],x2[4]), col=2, lty=1)
lines(c(x1[4],x1[4]),c(exp(league_eff$lower[4]),exp(league_eff$upper[4])),lty=1, col=4)

lines(c(-0.05,1.6),c(1,1),lty=2, lwd=0.6)
lines(c(1,1),c(-0.05,1.6),lty=2, lwd=0.6)

points(x1[2:3],x2[2:3],pch=c(17,15),cex=0.6)
points(x1[4],x2[4],pch=18, cex=0.8)

axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6), 
     lab=c("0.0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6"),
     tck=-0.01, las=2)
axis(1, at=c(0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6), 
     lab=c("0.0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6"),
     tck=-0.01, las=1)

legend("topright", c(expression("VKA + DAPT"),
                     expression("VKA + P2Y"[12]*" inhibitor"),
                     expression("NOAC + DAPT"),
                     expression("NOAC + P2Y"[12]*" inhibitor"), "95% CI for TIMI major bleeding", "95% CI for MACE"), 
       pch = c(1,17,15,18,NA,NA), lty = c(NA,NA,NA,NA,1,1),col = c(1,1,1,1,2,4),cex = 0.3, lwd = 0.6)
dev.off()















