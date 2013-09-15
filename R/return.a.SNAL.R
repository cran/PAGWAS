return.a.SNAL <- function(y,G,P,no.simulations,no.paths,no.snps){
	q=seq(1,1e-07,length.out=10)
	n=length(y)
	snp.parameters=c(10^{-1},10^{-2},10^{-3},10^{-4},10^{-5},10^{-6})
	cat("Different values for the hyper-parameter a are tested:","\n")
	cat(snp.parameters)
	cat("\n")
	results.frame=matrix(NA,ncol=length(snp.parameters),nrow=no.simulations)
	vec=c()

	object<-internal.function.simulated.data(y,G,P,no.paths,no.snps)
	y.pred=object$y.pred
	msq=object$msq
	effects.positions.paths=object$effects.positions.paths
	effects.positions.snps=object$effects.positions.snps	

	for(k in 1:no.simulations){
		set.seed(k)
    		y.sim <- rnorm(length(y.pred),y.pred,sqrt(msq))
		for(w in 1:length(snp.parameters)){
  			aa=snp.parameters[w]
  			ww=1
  			TPR.FPR=c()
  				for(ww in 1:length(q)){
    				s2=q[ww]
				gamma.star=SNAL(y.sim,G,P,a=aa,s2)
				gamma.star=gamma.star[[1]]
	    			vector.gamma=which(gamma.star!=0)
    				TPR.FPR=c(TPR.FPR,TPR.FPR(P,vector=vector.gamma,effects.positions.paths))
  			}
  		vector.TPR=seq(1,2*length(q),by=2)
  		vector.TNR=seq(2,2*length(q),by=2)
  		vector.TPR=TPR.FPR[vector.TPR]
  		vector.TNR=TPR.FPR[vector.TNR]
  		sens.wipf=c(0,vector.TPR,1)
  		spec.wipf=c(1,vector.TNR,0)
  		conv.wipf=roc.convex(sens=sens.wipf,spec=spec.wipf)
  		results.frame[k,w]=conv.wipf
		}
		vec=c(vec,which(results.frame[k,]==max(results.frame[k,])))
		}
	tt=table(vec)
	tt=which(tt==max(tt))
	aa.max=snp.parameters[tt]
	rm(tt,vec,conv.wipf,sens.wipf,spec.wipf,vector.TPR,vector.TNR,w,ww)
	rm(TPR.FPR,s2,k,object,y.pred,msq,effects.positions.paths,effects.positions.snps,n,snp.parameters,q)
	aa.max
	}

