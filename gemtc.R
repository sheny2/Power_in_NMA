set.seed(2918)


library(gemtc)
library(fields)

source("ifplot.fun.R")
source("BayesDiagnos.fun.R")
source("ranko_sucra.R")


bleed <- read.table("Data_Bleeding.csv", sep=",", header=TRUE)

######## Row 7 and 8 need to be combined
bleed[7,3:7] = bleed[7,3:7] + bleed[8,3:7]
bleed <- bleed[-8,]

######## Function: NMA under the consistency assumption  
fun.NMA.consistency <- function(data, trt, trts, NT, outcomename, figlab, color) {
  
  data = bleed
  outcomename = "TIMImajor" 
  
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
  
  
  # 
  # sucra_plot_value<- function(nma_model, filter){
  #   prob <- rank.probability(cons.out, preferredDirection=-1)
  #   prob = round(prob, digits=3)
  #   sucra<-sucra(prob)
  #   
  #   ### SUCRA
  #   cumprob.data <- apply(t(prob), 2, cumsum)
  #   n.treat = dim(cumprob.data)[2]
  #   rank<-seq(1:n.treat)
  #   cum_df1<-as.data.frame(cbind(cumprob.data,rank))
  #   
  #   cum_df<-cum_df1%>%
  #     pivot_longer(!rank, names_to = "trt")
  #   
  #   p<-(ggplot(cum_df,aes(rank,value))+
  #         geom_line()+
  #         theme_bw() + 
  #         theme(panel.border = element_blank(), 
  #               # panel.grid.major = element_blank(),
  #               #panel.grid.minor = element_blank(),
  #               axis.line = element_line(colour = "black"))+
  #         geom_area( fill = "darkblue"))+
  #     facet_wrap(~trt)+
  #     labs(title=paste0("SUCRA Plot for Filter: ",filter))
  #   print(p)
  #   return(sucra) }
    
  
  return(RESULTS)
}


trt <- c("A", "B", "C", "D")
trts <- read.table(textConnection('
                                  id description
                                  A "VKA + P2Y12 + Asprin"
                                  B "VKA + P2Y12"
                                  C "NOAC + P2Y12 + Asprin"
                                  D "NOAC + P2Y12"'), header=TRUE)
NT <- 4 ## number of treatments

TIMImajor.res <- fun.NMA.consistency(data=bleed, trt=trt, trts=trts, NT=NT, outcomename="TIMImajor", figlab="TIMI Major Bleeding", color=2)

TIMImajorminor.res <- fun.NMA.consistency(data=bleed, trt=trt, trts=trts, NT=NT, outcomename="TIMImajorminor", figlab="TIMI Major or Minor Bleeding", color=2)






