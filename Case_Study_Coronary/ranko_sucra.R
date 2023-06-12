#### Functions for plotting rankogram and sucra gram

### Input: 
# saf/eff: probability tables for safety and efficacy results
# name: used in file names
# trt.name: vector of treatment names

ranko <- function(saf, eff,name,lab.name,dir){
  n.treat = max(dim(saf)[2],dim(eff)[2])
  ind = c('A','B','C','D','E','F','G','H')
  
  for(i in 1:n.treat){
    pdf(paste(dir,name,"_Rank_",ind[i],".pdf",sep = ""),width = 3.2, height = 3.2)
    par(mar=c(4,4,2,1)+0.1)
    par(cex.axis=0.6, cex.lab=0.6)
    if(!is.null(saf))
      dry.pla <- saf[,ind[i]]
    if(!is.null(eff))
      imp.pla <- eff[,ind[i]]
    
    if(!is.null(saf)){
      plot(c(1:n.treat),dry.pla,type="l",xlab = "",ylab="",ylim=c(0,1),
           axes=FALSE,col="red", lwd=2)
      if(!is.null(eff))
        points(c(1:n.treat),imp.pla,type="l",lty=1,col="blue", lwd=2)
    }
    else
      plot(c(1:n.treat),imp.pla,type="l",xlab = "",ylab="",ylim=c(0,1),
           axes=FALSE,col="blue", lwd=2)
    
    axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0), 
         lab=rep("",6),
         tck=-0.01, las=2)
    axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0), 
         lab=c("0.0","0.2","0.4","0.6","0.8","1.0"),
         las=2, line = -0.2, lwd=0)
    title(ylab="Rank probability", line=1.8, cex.lab=0.6)
    
    axis(1, at=c(1,2,3,4,5,6,7,8), 
         lab=rep("",8),
         tck=-0.01, las=1)
    axis(1, at=c(1,2,3,4,5,6,7,8), 
         lab=c(expression("1"^"st"),expression("2"^"nd"),expression("3"^"rd"),expression("4"^"th"),expression("5"^"th"),expression("6"^"th"),expression("7"^"th"),expression("8"^"th")),
         las=1, line=-0.8, lwd = 0)
    mtext(lab.name[[i]], side = 1, at=n.treat/2+0.5, line = 1.2, cex=0.6)
    
    if(i==1){
      if(!is.null(saf)&!is.null(eff))
        legend(1,1,c("Major bleeding","MACE"),lty=c(1,1),bty = "n",cex=0.5,col = c("red","blue"))
      else if(!is.null(saf))
        legend(1,1,"Major bleeding",lty=1,bty = "n",cex=0.5,col = "red")
      else
        legend(1,1,"MACE",lty=1,bty = "n",cex=0.5,col = "blue")
    }
    #mtext(trt.name[i], NORTH<-3, at=1, line=0.25, cex=1.2)
    dev.off()
    
  }
  
}


### Input: 
# saf/eff: cumulative probability tables for safety and efficacy results
# name: used in file names
# trt.name: vector of treatment names

sucra.fun <- function(saf, eff,name,lab.name,dir){
  n.treat = max(dim(saf)[2],dim(eff)[2])
  ind = c('A','B','C','D','E','F','G','H')
  x = c(1,seq(from=1.5,by=1,length.out = n.treat-1),n.treat)
  
  for(i in 1:n.treat){
    pdf(paste(dir,name,"_SUCRA_",ind[i],".pdf",sep = ""),width = 3.2, height = 3.2)
    par(mar=c(4,4,2,1)+0.1)
    par(cex.axis=0.6, cex.lab=0.6)
    if(!is.null(saf))
      dry.pla <- c(saf[1,ind[i]],saf[1:(n.treat-1),ind[i]],saf[n.treat-1,ind[i]])
    if(!is.null(eff))
      imp.pla <- c(eff[1,ind[i]],eff[1:(n.treat-1),ind[i]],eff[n.treat-1,ind[i]])
    
    if(!is.null(saf)){
      plot(x,dry.pla,type="l",xlab="",ylab="",ylim=c(0,1),
         axes=FALSE,col="red", lwd=2)
      if(!is.null(eff))
        points(x,imp.pla,type="l",lty=1,col="blue", lwd=2)
    }
    else
      plot(x,imp.pla,type="l",xlab="",ylab="",ylim=c(0,1),
           axes=FALSE,col="blue", lwd=2)
    
    axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0), 
         lab=rep("",6),
         tck=-0.01, las=2)
    axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0), 
         lab=c("0.0","0.2","0.4","0.6","0.8","1.0"),
         las=2, line = -0.2, lwd=0)
    title(ylab="Cumulative rank probability", line=1.8, cex.lab=0.6)
    
    axis(1, at=c(1,2,3,4,5,6,7,8), 
         lab=rep("",8),
         tck=-0.01, las=1)
    axis(1, at=c(1,2,3,4,5,6,7,8), 
         lab=c(expression("1"^"st"),expression("2"^"nd"),expression("3"^"rd"),expression("4"^"th"),expression("5"^"th"),expression("6"^"th"),expression("7"^"th"),expression("8"^"th")),
         las=1, line=-0.8, lwd = 0)
    mtext(lab.name[[i]], side = 1, at=n.treat/2+0.5, line = 1.2, cex=0.6)
    
    if(i==1){
      if(!is.null(saf)&!is.null(eff))
        legend(1,1,c("Major bleeding","MACE"),lty=c(1,1),bty = "n",cex=0.5,col = c("red","blue"))
      else if(!is.null(saf))
        legend(1,1,"Major bleeding",lty=1,bty = "n",cex=0.5,col = "red")
      else
        legend(1,1,"MACE",lty=1,bty = "n",cex=0.5,col = "blue")
    }
    #mtext(trt.name[i], NORTH<-3, at=1, line=0.25, cex=1.2)
    dev.off()
    
  }
}

### Input: 
# saf/eff: sucra values
# name: used in file names
# trt.name: vector of treatment names

rank_plot <- function(saf, eff,name,lab.name,dir){
  n = dim(saf)[2]
  ind = c('A','B','C','D','E','F','G','H')
  lab = c(expression('8'^'th'),expression('7'^'th'),expression('6'^'th'),expression('5'^'th'),expression('4'^'th'),expression('3'^'rd'),expression('2'^'nd'),expression('1'^'st'))
  
    pdf(paste(dir,name,"_rank_plot",".pdf",sep = ""))
    par(mar=c(4,4,2,1)+0.1)
    
    x = rank(saf)
    y = rank(eff)

    
    plot(x,y, xlab="Ranking for major bleeding",ylab = "Ranking for MACE",ylim=c(0.8,n+0.2),xlim=c(0.8,n+0.2),xaxt='n',yaxt='n')
    for (i in 1:n)
      text(x[i],y[i], labels = c(lab.name[[i]]), cex=0.7,pos=1)
    
    axis(2, at=c(1:n), 
         lab=lab[8-n+1:n],
         tck=-0.01, las=2)
    
    axis(1, at=c(1:n), 
         lab=lab[8-n+1:n],
         tck=-0.01, las=1)
    dev.off()
  
}

