"BayesDiagnos" = function(binary = F, y, var, post.delta, r, n, post.p, post.D, Sid=NULL,ns.points=FALSE)
{
	# this function estimates diagnostics and model fit for Bayesian hierarchical models
	# y the (normal) observed effect size
	# var the variance of the effect size
	# post.delta matrix or vector of the mean posterior fitted values to y (from WinBUGS)
	# post.D matrix or vector of the mean posterior residual deviance (from WinBUGS)
	# n sample size
	# r events
	# post.p posterior probabilities(from WinBUGS)
	# BE CAREFUL with binary data! don't use the R database for r and n, better copy the two matrices from the winBUGS data file!
	if(binary) {
		if(is.matrix(Sid)) {
			Sid <- c(t(Sid))
		}
		Sid.zero<-Sid[which(r==0 | r==n)]
		if(!identical(Sid.zero, numeric(0))){
			for(i in 1: length(Sid)){
				for(j in 1:length(Sid.zero)){
				if(Sid[i]==Sid.zero[j]){
					n[i] <- n[i] + 1
					r[i] <- r[i] + 0.5		}
				}
			}
		}
		if(is.matrix(post.D)) {
			post.D <- c(t(post.D))
			post.D <- post.D[!is.na(post.D)]
		}
		if(is.matrix(post.p)) {
			post.p <- c(t(post.p))
			post.p <- post.p[!is.na(post.p)]
		}
		if(is.matrix(n)) {
			n <- c(t(n))
			n <- n[n > 1]
		}
		if(is.matrix(r)) {
			r <- c(t(r))
			r <- r[!is.na(c((r)))]
		}
		D.post <- -2 * (r * log((post.p * n)/r) + (n - r) * log(((1 - post.p) * n)/(n - r)))
	}else {
		D.post <- ((y - post.delta) * (y - post.delta))/var
	}
	pD <- post.D - D.post
	#leverage
	DIC <- sum(post.D) + sum(pD)
	x <- seq(-3, 3, 0.1)
	y1 <- 1 - x^2
	y2 <- 2 - x^2
	y3 <- 3 - x^2
	maxx<-round(max(sqrt(abs(post.D))),0)+0.5
	maxy<-round(max(pD),0)+0.5
	if(!is.null(Sid) & ns.points==TRUE){
		plot(sqrt(abs(post.D)), pD, pch = 1,xaxs = "i", yaxs = "i", main = "Fit of the model", xlab = "Residual Deviance (postDi)", ylab = "Leverage (pDi)",type="n",xlim=c(0,maxx),ylim=c(0,maxy))
		text(sqrt(abs(post.D)), pD,Sid)
	}
	if(ns.points==FALSE){
		plot(sqrt(abs(post.D)), pD, pch = 1,xaxs = "i", yaxs = "i", main = "Fit of the model", xlab = "Residual Deviance (postDi)", ylab = "Leverage (pDi)",xlim=c(0,maxx),ylim=c(0,maxy))
	}
	matlines(x, cbind(y1, y2, y3))
	cat(" ", "\n")
	cat("    *Model fit measures and diagnostics*",  "\n")
	cat("Residual deviance=", round(sum(post.D), 2), "\n")
	if(binary) {
		cat("Data points=", length(n), "\n")
	}else {
		cat("Data points=", length(y), "\n")
	}
	cat("(Note that total residual deviance should approximate the number of data points for a good fit)", "\n")
	cat("Effective number of parameters=", round(sum(pD), 2), "\n")
	cat("DIC=", round(DIC, 2), "\n")

}