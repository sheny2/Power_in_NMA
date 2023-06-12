"ifplot.fun" = function(id,theta,var.theta,treatA,treatB, plot = T, common.tau=F,method="DL", network.tau.sq=NULL){
 	options(warn = -1)
	library(combinat)
	library(metafor)
	Dataset<-data.frame(theta,var.theta,id,treatA,treatB)
	no.rows <- length(theta)
	if(is.numeric(treatA)){treatA<-letters[treatA]}
	if(is.numeric(treatB)){treatB<-letters[treatB]}
	if(data.class(treatA) !="character" | data.class(treatB) !="character") {
		treatA <- as.character(treatA)
		treatB <- as.character(treatB)
	}
	inlex.order <- as.numeric(apply(Dataset[, 4:5], 1, order)[1,  ] == 1)
	treat.A <- treatA
	treat.A[inlex.order == 0] <- treatB[inlex.order == 0]
	treat.B <- treatB
	treat.B[inlex.order == 0] <- treatA[inlex.order == 0]
	treatA <- treat.A
	treatB <- treat.B
	theta1 <- theta
	theta1[inlex.order == 0] <-  - theta1[inlex.order == 0]
	theta <- theta1
	comp <- c()
	for(i in 1:no.rows) {
		comp <- append(comp, paste(sort(c(treatA[i], treatB[i])),sep = "", collapse = ""))
	}
	Dataset <- cbind.data.frame(theta1, var.theta, id, comp)
	treatments <- unique(c(treatA, treatB))
	av.treat <- unique(comp)
	cat("\n","\n  *-----  Evaluating the consistency of the network ------*","\n")
	cat("\n  Nr of treatments: ", length(treatments))
	pos.l <- .possible.loops(treatments)[1]
	cat("\n  Nr of all possible first order loops (triangles): ", length(pos.l[[1]]))
	availability <- apply(pos.l[[1]], 2, match, av.treat)
	availability <- !is.na(apply(availability, 2, sum))
	pos.l <- pos.l[[1]][, availability]
	dim <- dim(pos.l)
	if(!is.null(dim)) {
		pos.l <- lapply(apply(pos.l, 2, list), unlist)
		no.loops <- length(pos.l)
		loops <- length(pos.l)
		cat("\n  Nr of available first order loops: ", length(pos.l), "\n", "\n")
	}
	if(is.null(dim)) {
		int<-as.integer(availability)
		for(i in 1:length(availability)){
			if(int[i]==1){
			name<-names(availability[i])
			}
		}		
		no.loops <- 1
		loops <- 1
		cat("\n  Nr of available first order loops: ", 1, "\n", "\n")
	}
	pos.l.q <- .possible.loops.q(treatments)
	cat("\n  Nr of all possible second order loops (quadrilaterals): ", length(pos.l.q[[1]]))
	tab<-data.frame()
	if(length(treatments)==4){
		pos.l.q <- c(apply(combn(2,treatments[1:4],sort), 2,paste,collapse="",sep=","))
		tab<-as.data.frame(av.treat[!is.na(match(av.treat,pos.l.q))])
		names(tab)<-paste(treatments[1:4],collapse="")
		name.q<-names(tab)
		tab1<-as.matrix(tab)
		data.comp<-data.frame()
		for(t in 1:4){
			split.comp<-unlist(strsplit(tab1[t,1],""))
			data.comp[t,1]<-split.comp[1]
			data.comp[t,2]<-split.comp[2]
		}
		mat<-as.matrix(data.comp)
		if(max(table(mat))==3){
			availability.q <-0
			tab<-data.frame()
			no.loops.q<-0
		}else{
		availability.q<-1
		no.loops.q <- 1
		}
		cat("\n  Nr of available second order loops (quadrilaterals): ", availability.q, "\n", "\n")
	}
	if(length(treatments)!=4){
		availability.q <-data.frame( apply(pos.l.q[[1]], 2, match, av.treat))		
		availability.final<-data.frame()
		pos.l.q<-data.frame(pos.l.q)
		z<-0
		for(i in 1:(dim(availability.q)[2])){		
			table<-data.frame()
			k<-0
			for(j in 1:6){	
				if(!is.na(availability.q[j,i])){
					k<-k+1
					table[k,1]<-pos.l.q[j,i]
				}
			}
			if(k>=4){
				l<-c()
				l<-c(substring(table[,1],1,1:1),substring(table[,1],2,2:2))
				m<-c(table(l))
				if(m[1]==2 & m[2]==2 & m[3]==2 & m[4]==2){
					availability.final[1,i]<-"TRUE"
					z<-z+1
					for(p in 1:dim(table)[1]){
						tab[p,z]<-table[p,1]
					}
				}else{
					availability.final[1,i]<-"FALSE"
				}
			}else{
				availability.final[1,i]<-"FALSE"
			}
		}
		names(availability.final)<-names(availability.q)
		no.loops.q <- dim(tab)[2]
		k<-1
		name.q<-c()
		for(i in 1:dim(availability.final)[2]){
			if(availability.final[i]=="TRUE"){
				name.q[k]<-names(availability.final[i])
				k<-k+1
			}
		}
		names(tab)<-name.q	
		cat("\n  Nr of available second order loops (quadrilaterals): ", dim(tab)[2], "\n", "\n")
	}
    EVALUATE.LOOP <- c()
    LOOPS <- c()
	tau.squared1<-c()
	tau.squared2<-c()
	tau.squared3<-c()
	tau.squared1.q<-c()
	tau.squared2.q<-c()
	tau.squared3.q<-c()
	tau.squared4.q<-c()
	if(!is.null(dim)){
    	results<-c()
		for(i in 1:no.loops) {
			cat("\n", i, ": Evaluation of the loop", names(pos.l[i]))
			loop <- c()
			loop.dimension <- length(pos.l[[i]])
			table.comp <- table(match(comp, pos.l[[i]]))
			names(table.comp) <- pos.l[[i]]
			cat("\n", "Direct comparisons in the loop:", "\n")
			print(table.comp)
			duplicat <- table(id[!is.na(match(comp, pos.l[[i]]))]) > 2
			if(sum(duplicat) > 0) {
				duplication <- names(duplicat[duplicat])
				DatasetOnuse <- Dataset
				leaveout <- c()
				for(k in 1:length(duplication)) {
					most.rep.comp <- names(table.comp[table.comp ==max(table.comp)])
					cat("\n The study with id=", duplication[k], sep = "")
					cat(" has more than two treatments in this loop")
					if(length(most.rep.comp) > 1) {
						most.rep.comp <- as.character(most.rep.comp[1])
					}
					leaveout <- append(leaveout, c(1:no.rows)[id == duplication[k] & comp ==most.rep.comp])
					cat(" (out is row nr ", leaveout[k], " for comparison ", as.character(comp[leaveout[k]]), ")", sep= "")
			}
			DatasetOnuse <- DatasetOnuse[ - leaveout,  ]				
			}else {
				DatasetOnuse <- Dataset
			}
			results1<-c()
			tau.squared<-c()
			if(common.tau==F){	
				for(j in 1:loop.dimension) {
					dta <- DatasetOnuse[DatasetOnuse$comp == pos.l[[i]][j], 1:2]
					Dnew <-Dataset[Dataset$comp == pos.l[[i]][j], 1:2]
					cat("\n   Meta-analysis for the", pos.l[[i]][j], "comparison")
					if(dim(dta)[1] == 0){
						cat("\n   The loop", names(pos.l[i]), "is described by only ONE MULTI-ARM trial.")
						cat("\n   There is no inconsistency for the loop", names(pos.l[i]))
						metaanalysis<-c(Value = Dnew[, 1], se = sqrt(Dnew[, 2]))
						tau.squared[j]<-0
					}
					if(dim(dta)[1] == 1) {
						metaanalysis <- c(Value = dta[, 1], se = sqrt(dta[, 2]))
						tau.squared[j]<-0
					}else {
						if(dim(dta)[1] > 1) {
							o<-rma(dta[,1],dta[,2],data=dta,method=method)
							metaanalysis <- c(o$b, o$se)
							tau.squared[j]<-o$tau
						}
					}
					tau.squared1[i]<-tau.squared[1]
					tau.squared2[i]<-tau.squared[2]
					tau.squared3[i]<-tau.squared[3]
					cat("\n   mean(se)=", round(metaanalysis[1], 3))
					cat("(", round(metaanalysis[2], 3), ")", sep = "")
					loop <- rbind(loop, metaanalysis)
					results1<- rbind(results1,c(pos.l[[i]][j],round(metaanalysis[1], 3),round(metaanalysis[2], 3)))
				}
				
			}
			if(is.null(network.tau.sq) & common.tau==T){
				data.metareg<-c()
				for(j in 1:loop.dimension) {
					dta <- DatasetOnuse[DatasetOnuse$comp == pos.l[[i]][j], 1:4]
					if(dim(dta)[1] != 0){
						D <- DatasetOnuse[DatasetOnuse$comp == pos.l[[i]][j], 1:4]
						data.metareg<-rbind(data.metareg, D)		
					}else{
						D <- Dataset[Dataset$comp == pos.l[[i]][j], 1:4]
						data.metareg<-rbind(data.metareg, D)	
					}
				}
				available.comps<-c()
				available.comps<-unique(data.metareg$comp)
				available.comps<-matrix(available.comps,nrow=length(available.comps),byrow=T)
				tr<-1:nrow(available.comps)
				available.comps<-as.data.frame(available.comps)
				available.comps[,2]<-tr
				comp2<-match(data.metareg$comp,available.comps[,1])	
				data.metareg$comp.a<-as.numeric(comp2==1)
				data.metareg$comp.b<-as.numeric(comp2==2)
				data.metareg$comp.c<-as.numeric(comp2==3)
				if(nrow(data.metareg)==3){
					tau.squared<-0
				}else{
					o<-rma(data.metareg$theta1, data.metareg$var.theta, mods=cbind(comp.a,comp.b,comp.c),data=data.metareg,intercept=FALSE,method=method)
					tau.squared<-o$tau
				}
				for(j in 1:loop.dimension) {
					dta <- DatasetOnuse[DatasetOnuse$comp == pos.l[[i]][j], 1:4]
					Dnew <-Dataset[Dataset$comp == pos.l[[i]][j], 1:2]
					cat("\n   Meta-analysis for the", pos.l[[i]][j], "comparison")	
					if(dim(dta)[1] == 0){
						cat("\n   The loop", names(pos.l[i]), "is described by only ONE MULTI-ARM trial.")
						cat("\n   There is no inconsistency for the loop", names(pos.l[i]))
						o<-new.meta(Dnew[, 1], sqrt(Dnew[, 2]),tau.squared)
					}else{
						o<-new.meta(dta[,1],sqrt(dta[,2]),tau.squared)
					}
					metaanalysis <- c(o$TE.random, o$seTE.random)
					cat("\n   mean(se)=", round(metaanalysis[1], 3))
					cat("(", round(metaanalysis[2], 3), ")", sep = "")
					loop <- rbind(loop, metaanalysis)
					results1<- rbind(results1,c(pos.l[[i]][j],round(metaanalysis[1], 3),round(metaanalysis[2], 3)))	
				}
				tau.squared1[i]<-tau.squared
				tau.squared2[i]<-tau.squared
				tau.squared3[i]<-tau.squared
		}
		if(!is.null(network.tau.sq) & common.tau==T){
			for(j in 1:loop.dimension) {
				dta <- DatasetOnuse[DatasetOnuse$comp == pos.l[[i]][j], 1:4]
				Dnew <-Dataset[Dataset$comp == pos.l[[i]][j], 1:2]
				cat("\n   Meta-analysis for the", pos.l[[i]][j], "comparison")		
				if(dim(dta)[1] == 0){
					o<-new.meta(Dnew[, 1], sqrt(Dnew[, 2]),network.tau.sq)
					cat("\n   The loop", names(pos.l[i]), "is described by only ONE MULTI-ARM trial.")
					cat("\n   There is no inconsistency for the loop", names(pos.l[i]))
				}else{
					o<-new.meta(dta[,1],sqrt(dta[,2]),network.tau.sq)
				}
				metaanalysis <- c(o$TE.random, o$seTE.random)
				cat("\n   mean(se)=", round(metaanalysis[1], 3))
				cat("(", round(metaanalysis[2], 3), ")", sep = "")
				loop <- rbind(loop, metaanalysis)
				results1<- rbind(results1,c(pos.l[[i]][j],round(metaanalysis[1], 3),round(metaanalysis[2], 3)))	
			}
			tau.squared1[i]<-network.tau.sq
			tau.squared2[i]<-network.tau.sq
			tau.squared3[i]<-network.tau.sq
		}
			dimnames(loop)[1] <- list(pos.l[[i]])
			loop[, 2] <- loop[, 2]^2
			leftLoop <- c(sum(loop[1:(loop.dimension - 1), 1]), sum(loop[1:(loop.dimension - 1), 2]))
			cat("\n   Indirect summary estimate for the", names(table.comp)[3], "comparison")
			cat("\n   Mean(se)=", round(leftLoop[1], 3))
			cat("(", round(sqrt(leftLoop[2]), 3), ")", sep = "")
			EvaluateLoop <- c(leftLoop[1] - loop[loop.dimension, 1],leftLoop[2] + loop[loop.dimension,2])
			cat("\n", "\n   Inconsistency within the loop")
			cat(":  Mean(se)=", round(EvaluateLoop[1], 3))
			cat("(", round(sqrt(EvaluateLoop[2]), 3), ")", sep = "", "\n","\n")
			LOOPS <- append(LOOPS, list(loop))
			EVALUATE.LOOP <- rbind(EVALUATE.LOOP, EvaluateLoop)	
			results<- rbind(results,c(names(pos.l[i]),results1,as.numeric(table.comp),round(tau.squared1[i],4),round(tau.squared2[i],4),round(tau.squared3[i],4)))	
		}
		names(LOOPS) <- names(pos.l)
		dimnames(EVALUATE.LOOP)[1] <- list(names(pos.l))
	}
    if(dim(tab)[2]!=0){	
	EVALUATE.LOOP.q <- c()
	LOOPS.q <- c()
	results.q<-c()
		for(i in 1:no.loops.q) {
			cat("\n", i, ": Evaluation of the loop", name.q[i])
			loop.q <- c()
			loop.dimension.q <- dim(tab)[1]
			table.comp <- table(match(comp, tab[,i]))
			names(table.comp) <- tab[,i]
			cat("\n", "Direct comparisons in the loop:", "\n")
			print(table.comp)
			DatasetOnuse <- Dataset
			results1.q<-c()
			tau.squared.q<-c()
			tab1<-as.matrix(tab)
			if(common.tau==F){	
				for(j in 1:loop.dimension.q) {
					cat("\n   Meta-analysis for the", tab1[j,i], "comparison")
					dta.q <- DatasetOnuse[DatasetOnuse$comp == tab1[j,i], 1:2]
					if(dim(dta.q)[1] == 1) {
						metaanalysis.q <- c(Value = dta.q[, 1], se = sqrt(dta.q[, 2]))
						tau.squared.q[j]<-0
					}else {
						o.q<-rma(dta.q[,1],dta.q[,2],data=dta.q,method=method)
						metaanalysis.q <- c(o.q$b, o.q$se)
						tau.squared.q[j]<-o.q$tau
					}
					cat("\n   mean(se)=", round(metaanalysis.q[1], 3))
					cat("(", round(metaanalysis.q[2], 3), ")", sep = "")
					loop.q <- rbind(loop.q, metaanalysis.q)
					results1.q<- rbind(results1.q,c(tab1[j,i],round(metaanalysis.q[1], 3),round(metaanalysis.q[2], 3)))	
				}
				tau.squared1.q[i]<-tau.squared.q[1]
				tau.squared2.q[i]<-tau.squared.q[2]
				tau.squared3.q[i]<-tau.squared.q[3]
				tau.squared4.q[i]<-tau.squared.q[4]
			}
			if(is.null(network.tau.sq) & common.tau==T){
				data.metareg.q<-c()
				for(j in 1:loop.dimension.q) {
					D <- Dataset[Dataset$comp == tab1[j,i], 1:4]
					data.metareg.q<-rbind(data.metareg.q, D)	
				}
				available.comps.q<-c()
				available.comps.q<-unique(data.metareg.q$comp)
				available.comps.q<-matrix(available.comps.q,nrow=length(available.comps.q),byrow=T)
				tr.q<-1:nrow(available.comps.q)
				available.comps.q<-as.data.frame(available.comps.q)
				available.comps.q[,2]<-tr.q
				comp2.q<-match(data.metareg.q$comp,available.comps.q[,1])	
				data.metareg.q$comp.a<-as.numeric(comp2.q==1)
				data.metareg.q$comp.b<-as.numeric(comp2.q==2)
				data.metareg.q$comp.c<-as.numeric(comp2.q==3)
				if(nrow(data.metareg.q)==4){
					tau.squared.q<-0
				}else{
					o.q<-rma(data.metareg.q$theta1, data.metareg.q$var.theta, mods=cbind(comp.a,comp.b,comp.c),data=data.metareg.q,intercept=FALSE,method=method)
					tau.squared.q<-o.q$tau
				}
				results1.q<-c()
				for(j in 1:loop.dimension.q) {
					dta.q <- DatasetOnuse[DatasetOnuse$comp == tab1[j,i], 1:4]
					cat("\n   Meta-analysis for the", tab1[j,i], "comparison")		
					Dnew <-Dataset[Dataset$comp ==tab1[j,i], 1:2]
					o.q<-new.meta(dta.q[,1],sqrt(dta.q[,2]),tau.squared.q)
					metaanalysis.q <- c(o.q$TE.random, o.q$seTE.random)
					cat("\n   mean(se)=", round(metaanalysis.q[1], 3))
					cat("(", round(metaanalysis.q[2], 3), ")", sep = "")
					loop.q <- rbind(loop.q, metaanalysis.q)
					results1.q<- rbind(results1.q,c(tab1[j,i],round(metaanalysis.q[1], 3),round(metaanalysis.q[2], 3)))	
				}
				tau.squared1.q[i]<-tau.squared.q
				tau.squared2.q[i]<-tau.squared.q
				tau.squared3.q[i]<-tau.squared.q
				tau.squared4.q[i]<-tau.squared.q
			}
			if(!is.null(network.tau.sq) & common.tau==T){
				results1.q<-c()
				for(j in 1:loop.dimension.q) {
					dta.q <- DatasetOnuse[DatasetOnuse$comp == tab1[j,i], 1:4]
					cat("\n   Meta-analysis for the", tab1[j,i], "comparison")		
					Dnew <-Dataset[Dataset$comp ==tab1[j,i], 1:2]
					o.q<-new.meta(dta.q[,1],sqrt(dta.q[,2]),network.tau.sq)
					metaanalysis.q <- c(o.q$TE.random, o.q$seTE.random)
					cat("\n   mean(se)=", round(metaanalysis.q[1], 3))
					cat("(", round(metaanalysis.q[2], 3), ")", sep = "")
					loop.q <- rbind(loop.q, metaanalysis.q)
					results1.q<- rbind(results1.q,c(tab1[j,i],round(metaanalysis.q[1], 3),round(metaanalysis.q[2], 3)))	
				}
				tau.squared1.q[i]<-network.tau.sq
				tau.squared2.q[i]<-network.tau.sq
				tau.squared3.q[i]<-network.tau.sq
				tau.squared4.q[i]<-network.tau.sq
			}
			dimnames(loop.q)[1] <- list(tab1[,i])
			loop.q[, 2] <- loop.q[, 2]^2
			data.comp<-data.frame()
			for(t in 1:4){
				split.comp<-unlist(strsplit(tab1[t,1],""))
				data.comp[t,1]<-split.comp[1]
				data.comp[t,2]<-split.comp[2]
			}
			data.comp[,3]<-tab1
			indirect<-c()
			for(t in 1:3){
				if((data.comp[4,1]==data.comp[t,1])||(data.comp[4,2]==data.comp[t,2])){
					indirect[t]<-as.numeric(results1.q[t,2])
				}
				if((data.comp[4,1]==data.comp[t,2])||(data.comp[4,2]==data.comp[t,1])){
					p<- as.numeric(results1.q[t,2])
					indirect[t]<- -p
				}else{
					if(t==1){
						if((data.comp[t,1]==data.comp[t+1,1] & data.comp[4,2]==data.comp[t+1,2])||(data.comp[t,1]==data.comp[t+1,2] & data.comp[4,2]==data.comp[t+1,1])){
							p<- as.numeric(results1.q[t,2])
							indirect[t]<- -p					
						}
						if((data.comp[t,1]==data.comp[t+1,1] & data.comp[4,1]==data.comp[t+1,2])||(data.comp[t,1]==data.comp[t+1,2] & data.comp[4,1]==data.comp[t+1,1])){
							indirect[t]<-as.numeric(results1.q[t,2])				
						}		
					}
					if(t>1){
						if((data.comp[t,1]==data.comp[t-1,1] & data.comp[4,2]==data.comp[t-1,2])||(data.comp[t,1]==data.comp[t-1,2] & data.comp[4,2]==data.comp[t-1,1])||(data.comp[t,1]==data.comp[t-1,1] & data.comp[4,1]==data.comp[t-1,2])||(data.comp[t,1]==data.comp[t-1,2] & data.comp[4,1]==data.comp[t-1,1])){
							p<- as.numeric(results1.q[t,2])
							indirect[t]<- -p					
						}
					}
					if(t>2){
						if((data.comp[t,1]==data.comp[t-2,1] & data.comp[4,2]==data.comp[t-2,2])||(data.comp[t,1]==data.comp[t-2,2] & data.comp[4,2]==data.comp[t-2,1])||(data.comp[t,1]==data.comp[t-2,1] & data.comp[4,1]==data.comp[t-2,2])||(data.comp[t,1]==data.comp[t-2,2] & data.comp[4,1]==data.comp[t-2,1])){
							p<- as.numeric(results1.q[t,2])
							indirect[t]<- -p					
						}
					}
				}
			}
			leftLoop.q<-c(sum(indirect),sum(loop.q[1:(loop.dimension.q - 1), 2]))
			cat("\n   Indirect summary estimate for the", names(table.comp)[4], "comparison")
			cat("\n   Mean(se)=", round(leftLoop.q[1], 3))
			cat("(", round(sqrt(leftLoop.q[2]), 3), ")", sep = "")
			EvaluateLoop.q <- c(leftLoop.q[1] - loop.q[loop.dimension.q, 1], leftLoop.q[2] + loop.q[loop.dimension.q, 2])
			cat("\n", "\n   Inconsistency within the loop")
			cat(":  Mean(se)=", round(EvaluateLoop.q[1], 3))
			cat("(", round(sqrt(EvaluateLoop.q[2]), 3), ")", sep = "", "\n","\n")
			LOOPS.q <- append(LOOPS.q, list(loop.q))
			EVALUATE.LOOP.q <- rbind(EVALUATE.LOOP.q, EvaluateLoop.q)
			results.q<- rbind(results.q,c(name.q[i],results1.q,as.numeric(table.comp),round(tau.squared1.q[i],4),round(tau.squared2.q[i],4),round(tau.squared3.q[i],4),round(tau.squared4.q[i],4)))	
		}
		names(LOOPS.q) <- name.q
		dimnames(EVALUATE.LOOP.q)[1] <- list(name.q)
	}
 if(is.null(dim)){
    results<-c()
    loop <- c()
	cat("\n", 1, ": Evaluation of the loop", name)
	loop.dimension <- length(pos.l)
	table.comp <- table(match(comp, pos.l))
	names(table.comp) <- pos.l
	cat("\n", "Direct comparisons in the loop:", "\n")
	print(table.comp)
	duplicat <- table(id[!is.na(match(comp, pos.l))]) > 2
	if(sum(duplicat) > 0) {
		duplication <- names(duplicat[duplicat])
		DatasetOnuse <- Dataset
		leaveout <- c()
		for(k in 1:length(duplication)) {
			most.rep.comp <- names(table.comp[table.comp == max(table.comp)])
			cat("\n The study with id=", duplication[k], sep = "")
			cat(" has more than two treatments in this loop")
			if(length(most.rep.comp) > 1) {
				most.rep.comp <- as.character(most.rep.comp[1])
			}
			leaveout <- append(leaveout, c(1:no.rows)[id == duplication[k] & comp == most.rep.comp])
			cat(" (out is row nr ", leaveout[k], " for comparison ", as.character(comp[leaveout[k]]), ")", sep = "")
		}
		DatasetOnuse <- DatasetOnuse[ - leaveout,  ]
	}else{
		DatasetOnuse <- Dataset
	}
	results1<-c()
	tau.squared	<-c()
	if(common.tau==F){	
		for(j in 1:loop.dimension) {
			dta <- DatasetOnuse[DatasetOnuse$comp == pos.l[j], 1:2]
			cat("\n   Meta-analysis for the", pos.l[j], "comparison")
			if(dim(dta)[1] == 0){
				cat("\n   The loop", name, "is described by only ONE MULTI-ARM trial.")
				cat("\n   There is no inconsistency for the loop", name)
				Dnew <-Dataset[Dataset$comp == pos.l[j], 1:2]
				metaanalysis<-c(Value = Dnew[, 1], se = sqrt(Dnew[, 2]))
				tau.squared[j]<-0
			}
			if(dim(dta)[1] == 1) {
				metaanalysis <- c(Value = dta[, 1], se = sqrt(dta[, 2]))
				tau.squared[j]<-0
			}
			else {
				if(dim(dta)[1] > 1) {
					o<-rma(dta[,1],dta[,2],data=dta,method=method)
					metaanalysis <- c(o$b, o$se)
					tau.squared[j]<-o$tau
				}
			}
			cat("\n   mean(se)=", round(metaanalysis[1], 3))
			cat("(", round(metaanalysis[2], 3), ")", sep = "")
			loop <- rbind(loop, metaanalysis)
			results1<- rbind(results1,c(pos.l[j],round(metaanalysis[1], 3),round(metaanalysis[2], 3)))	
		}
		tau.squared1<-tau.squared[1]
		tau.squared2<-tau.squared[2]
		tau.squared3<-tau.squared[3]
	}
	if(is.null(network.tau.sq) & common.tau==T){
		data.metareg<-c()
		for(j in 1:loop.dimension) {
			dta <- DatasetOnuse[DatasetOnuse$comp == pos.l[j], 1:4]
			if(dim(dta)[1] != 0){
				D <- DatasetOnuse[DatasetOnuse$comp == pos.l[j], 1:4]
				data.metareg<-rbind(data.metareg, D)		
			}else{
				D <- Dataset[Dataset$comp == pos.l[j], 1:4]
				data.metareg<-rbind(data.metareg, D)	
			}
		}
		available.comps<-c()
		available.comps<-unique(data.metareg$comp)
		available.comps<-matrix(available.comps,nrow=length(available.comps),byrow=T)
		tr<-1:nrow(available.comps)
		available.comps<-as.data.frame(available.comps)
		available.comps[,2]<-tr
		comp2<-match(data.metareg$comp,available.comps[,1])	
		data.metareg$comp.a<-as.numeric(comp2==1)
		data.metareg$comp.b<-as.numeric(comp2==2)
		data.metareg$comp.c<-as.numeric(comp2==3)
		if(nrow(data.metareg)==3){
			tau.squared<-0
		}else{
			o<-rma(data.metareg$theta1, data.metareg$var.theta, mods=cbind(comp.a,comp.b,comp.c),data=data.metareg,intercept=FALSE,method=method)
			tau.squared<-o$tau
		}
		for(j in 1:loop.dimension) {
			dta <- DatasetOnuse[DatasetOnuse$comp == pos.l[j], 1:2]
			cat("\n   Meta-analysis for the", pos.l[j], "comparison")	
			if(dim(dta)[1] == 0){
				cat("\n   The loop", name, "is described by only ONE MULTI-ARM trial.")
				cat("\n   There is no inconsistency for the loop", name)
				Dnew <-Dataset[Dataset$comp == pos.l[j], 1:2]
				o<-new.meta(Dnew[, 1], sqrt(Dnew[, 2]),tau.squared)
			}else{
				o<-new.meta(dta[,1],sqrt(dta[,2]),tau.squared)
			}
			metaanalysis <- c(o$TE.random, o$seTE.random)
			cat("\n   mean(se)=", round(metaanalysis[1], 3))
			cat("(", round(metaanalysis[2], 3), ")", sep = "")
			loop <- rbind(loop, metaanalysis)
			results1<- rbind(results1,c(pos.l[j],round(metaanalysis[1], 3),round(metaanalysis[2], 3)))	
		}
		tau.squared1<-tau.squared
		tau.squared2<-tau.squared
		tau.squared3<-tau.squared
	}
	if(!is.null(network.tau.sq) & common.tau==T){
		for(j in 1:loop.dimension) {
			dta <- DatasetOnuse[DatasetOnuse$comp == pos.l[j], 1:2]
			cat("\n   Meta-analysis for the", pos.l[j], "comparison")	
			if(dim(dta)[1] == 0){
				cat("\n   The loop", name, "is described only by ONE MULTI-ARM trial.")
				cat("\n   There is no inconsistency for the loop", name)
				Dnew <-Dataset[Dataset$comp == pos.l[j], 1:2]
				o<-new.meta(Dnew[, 1], sqrt(Dnew[, 2]),network.tau.sq)
			}else{
				o<-new.meta(dta[,1],sqrt(dta[,2]),network.tau.sq)
			}
			metaanalysis <- c(o$TE.random, o$seTE.random)
			cat("\n   mean(se)=", round(metaanalysis[1], 3))
			cat("(", round(metaanalysis[2], 3), ")", sep = "")
			loop <- rbind(loop, metaanalysis)
			results1<- rbind(results1,c(pos.l[j],round(metaanalysis[1], 3),round(metaanalysis[2], 3)))	
		}
		tau.squared1<-network.tau.sq
		tau.squared2<-network.tau.sq
		tau.squared3<-network.tau.sq
	}
	loop[, 2] <- loop[, 2]^2
	leftLoop <- c(sum(loop[1:(loop.dimension - 1), 1]), sum(loop[1:(loop.dimension - 1), 2]))
	cat("\n   Indirect summary estimate for the", names(table.comp)[3], "comparison")
	cat("\n   Mean(se)=", round(leftLoop[1], 3))
	cat("(", round(sqrt(leftLoop[2]), 3), ")", sep = "")
	EvaluateLoop <- c(leftLoop[1] - loop[loop.dimension, 1], leftLoop[2] + loop[loop.dimension, 2])
	cat("\n", "\n   Inconsistency within the loop")
	cat(":  Mean(se)=", round(EvaluateLoop[1], 3))
	cat("(", round(sqrt(EvaluateLoop[2]), 3), ")", sep = "", "\n","\n")
	LOOPS <- append(LOOPS, list(loop))
	EVALUATE.LOOP <- rbind(EVALUATE.LOOP, EvaluateLoop)
	results<- rbind(results,c(name,results1,as.numeric(table.comp),round(tau.squared1,4),round(tau.squared2,4),round(tau.squared3,4)))	
 names(LOOPS) <- name
 dimnames(EVALUATE.LOOP)[1] <- list(name)
 }
 if(plot) {
	if(dim(tab)[2]!=0){
		par(mfrow=c(1,2))
		forest.default(abs(EVALUATE.LOOP[,1]),vi=EVALUATE.LOOP[,2], slab=names(LOOPS),xlab="Inconsistency for triangular loops")
		forest.default(abs(EVALUATE.LOOP.q[,1]),vi=EVALUATE.LOOP.q[,2], slab=names(LOOPS.q),xlab="Inconsistency for quadrilateral loops")
	}else{
		forest.default(abs(EVALUATE.LOOP[,1]),vi=EVALUATE.LOOP[,2], slab=names(LOOPS),xlab="Inconsistency for triangular loops")
	}
 }
 z.score<-round(EVALUATE.LOOP[,1]/sqrt(EVALUATE.LOOP[,2]),3)
 p.value<-round(2*pnorm(-abs(z.score)),3)
 out <- cbind.data.frame(names(LOOPS),round(abs(EVALUATE.LOOP[,1]),3),round(sqrt(EVALUATE.LOOP[,2]),3),z.score,p.value)
 names(out)<-c("Loop","IF","seIF","z.score","p.value")
 invisible()
 options(warn = 1)
 out<-as.data.frame(out)
 results<-as.data.frame(results)
 final<-data.frame(results,out[,2],out[,3],out[,4],out[,5])
 names(final)<-c("loop","comp1","comp2","comp3","mean1","mean2","mean3","se1","se2","se3","n1","n2","n3","tausq1","tausq2","tausq3","abs(IF)","seIF","z.score","p.value")
 if(dim(tab)[2]!=0){
	z.score.q<-round(EVALUATE.LOOP.q[,1]/sqrt(EVALUATE.LOOP.q[,2]),3)
	p.value.q<-round(2*pnorm(-abs(z.score.q)),3)
	out.q <- cbind.data.frame(names(LOOPS.q),round(abs(EVALUATE.LOOP.q[,1]),3),round(sqrt(EVALUATE.LOOP.q[,2]),3),z.score.q,p.value.q)
	names(out.q)<-c("Loop.q","IF.q","seIF.q","z.score(q)","p.value(q)")
	out.q<-as.data.frame(out.q)
	results.q<-as.data.frame(results.q)
	final.q<-data.frame(results.q,out.q[,2],out.q[,3],out.q[,4],out.q[,5])
	names(final.q)<-c("loop","comp1","comp2","comp3","comp4","mean1","mean2","mean3","mean4","se1","se2","se3","se4","n1","n2","n3","n4","tausq1","tausq2","tausq3","tausq4","abs(IF.q)","seIF.q","z.score(q)","p.value(q)")
 }
 if(dim(tab)[2]!=0){
	total<-list(triangles=final,quadrilaterals=final.q)
 }else{
	total<-final
 }
 return(total)
}
#########************************************ Sub-routines *****************************************************############
".possible.loops" = function(treatments){
	t <- length(treatments)
	t <- length(treatments)
	z <- list()
	u <- c()
	for(i in 1:t) {
		u <- cbind(u, (rbind(treatments[i], combn(2, x = treatments[- i]))))
	}
	u <- apply(u, 2, sort)
	price.in.vectors <- apply(u, 2, .needed.vectors)
	dimnames(price.in.vectors) <- list(NULL, apply(u, 2, paste, collapse = ""))
	price.in.vectors <- price.in.vectors[, unique(apply(u, 2, paste, collapse = ""))]
	z <- append(z, list(price.in.vectors))
	z
}
".possible.loops.q" = function(treatments){
	t <- length(treatments)
	z <- list()
	u <- c()
	if(t==4){
		u<-rbind(u,treatments[1:4])
		u <- apply(u, 2, sort)
		price.in.vectors <- c(apply(combn(2,treatments[1:4]), 2,paste,collapse="",sep=","))
		names(price.in.vectors)<-paste(treatments[1:4],collapse="")
		z <-price.in.vectors
	}
	
	if(t!=4){
		for(i in 1:t) {
			u <- cbind(u, (rbind(treatments[i], combn(3, x = treatments[- i]))))
		}
	u <- apply(u, 2, sort)
	price.in.vectors <- apply(u, 2, .needed.vectors.q)
	dimnames(price.in.vectors) <- list(NULL, apply(u, 2, paste, collapse = ""))
	price.in.vectors <- price.in.vectors[, unique(apply(u, 2, paste, collapse = ""))]
	z <- append(z, list(price.in.vectors))
	}
	z
}
"combn" = function(m, x, fun = NULL, simplify = TRUE, ...){
	if(length(m) > 1) {
		warning(paste("Argument m has", length(m), "elements: only the first used"))
		m <- m[1]
	}
	if(m < 0)
		stop("m < 0")
	if(m == 0)
		return(if(simplify) vector(mode(x), 0) else list())
	if(is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == x)
		x <- seq(x)
	n <- length(x)
	if(n < m)
		stop("n < m")
	e <- 0
	h <- m
	a <- 1:m
	nofun <- is.null(fun)
	count <- nCm(n, m, 0.1)
	out <- vector("list", count)
	out[[1]] <- if(nofun) x[a] else fun(x[a], ...)
	if(simplify) {
		dim.use <- NULL
		if(nofun) {
			if(count > 1)
				dim.use <- c(m, count)
		}
		else {
			out1 <- out[[1]]
			d <- dim(out1)
			if(count > 1) {
				if(length(d) > 1)
					dim.use <- c(d, count)
				else if(length(out1) > 1)
					dim.use <- c(length(out1), count)
			}
			else if(length(d) > 1)
				dim.use <- d
		}
	}
	i <- 2
	nmmp1 <- n - m + 1
	mp1 <- m + 1
	while(a[1] != nmmp1) {
		if(e < n - h) {
			h <- 1
			e <- a[m]
			j <- 1
		}else {
			h <- h + 1
			e <- a[mp1 - h]
			j <- 1:h
		}
		a[m - h + j] <- e + j
		out[[i]] <- if(nofun) x[a] else fun(x[a], ...)
		i <- i + 1
	}
	if(simplify) {
		if(is.null(dim.use))
			out <- unlist(out)
		else out <- array(unlist(out), dim.use)
	}
	out
}
"nCm" = function(n, m, tol = 1e-008){
	len <- max(length(n), length(m))
	out <- numeric(len)
	n <- rep(n, length = len)
	m <- rep(m, length = len)
	mint <- (trunc(m) == m)
	out[!mint] <- NA
	out[m == 0] <- 1
	whichm <- (mint & m > 0)
	whichn <- (n < 0)
	which <- (whichm & whichn)
	if(any(which)) {
		nnow <- n[which]
		mnow <- m[which]
		out[which] <- ((-1)^mnow) * Recall(mnow - nnow - 1, mnow)
	}
	whichn <- (n > 0)
	nint <- (trunc(n) == n)
	which <- (whichm & whichn & !nint & n < m)
	if(any(which)) {
		nnow <- n[which]
		mnow <- m[which]
		foo <- function(j, nn, mm)
		{
			n <- nn[j]
			m <- mm[j]
			iseq <- seq(n - m + 1, n)
			negs <- sum(iseq < 0)
			((-1)^negs) * exp(sum(log(abs(iseq))) - lgamma(m +
				1))
		}
		out[which] <- unlist(lapply(seq(along = nnow), foo, nn = nnow,
			mm = mnow))
	}
	which <- (whichm & whichn & n >= m)
	nnow <- n[which]
	mnow <- m[which]
	out[which] <- exp(lgamma(nnow + 1) - lgamma(mnow + 1) - lgamma(nnow -mnow + 1))
	nna <- !is.na(out)
	outnow <- out[nna]
	rout <- round(outnow)
	smalldif <- abs(rout - outnow) < tol
	outnow[smalldif] <- rout[smalldif]
	out[nna] <- outnow
	out
}
".needed.vectors" = function(loop){
	leng <- length(loop)
	need.vec <- apply(apply(rbind(loop, c(loop[-1], loop[1])), 2, sort),2, paste, collapse = "")
	need.vec
}
".needed.vectors.q" = function(loop){
	leng <- length(loop)
	need.vec <- apply(apply(combn(2,loop), 2, sort),2, paste, collapse = "")
	need.vec
}
new.meta<-function(theta,se.theta,tau.squared){
	w<-c()
	product<-c()
	var.theta<-c()
	for(i in 1:length(theta)){
		var.theta[i]<-(se.theta[i])^2
		w[i]<-1/(var.theta[i]+tau.squared)
		product[i]<-w[i]*theta[i]
	}
	TE.random<-sum(product)/sum(w)
	VarTE.random<-1/sum(w)
	seTE.random<-sqrt(VarTE.random)
	results<-list(TE.random=TE.random,seTE.random=seTE.random)
	return(results)
}