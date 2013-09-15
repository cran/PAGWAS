return.bf.NBF<-function(y,G,P,a,b,s2,nu,no.simulations,no.paths,no.snps,BF.cutoffs){
	BF.cutoff<-match.arg(BF.cutoffs)
	if(missing(BF.cutoffs)){
		BF.cutoffs=seq(1e-07,1,length.out=10)
	}
	
	object<-internal.function.simulated.data(y,G,P,no.paths,no.snps)
	y.pred=object$y.pred
	msq=object$msq
	effects.positions.paths=object$effects.positions.paths
	effects.positions.snps=object$effects.positions.snps		
	
	results.frame=matrix(NA,nrow=no.simulations,ncol=2*length(q)+1)
	w=1
	for(w in 1:no.simulations){
  		set.seed(w)
  		y.sim <- rnorm(length(y.pred),y.pred,sqrt(msq))
  		results.frame[w,1]=var(y.sim)
  		i=1
  		BF=NBF(y.sim,G,P,a,b,s2,nu)
  		rm(i)
  		bf.TPR.FPR=lapply(q,bf.TPR.FPR,BF)
  		results.frame[w,2:ncol(results.frame)]=unlist(bf.TPR.FPR)
  		print(w)
	}
	mean.median=data.frame("TPR"=NA,"sd(TPR)"=NA,"FPR"=NA,"sd(FPR)"=NA)
	w=1
	for(w in 1:length(BF.cutoffs)){
  		mean.median[w,1]=mean(results.frame[,2*w])
  		mean.median[w,2]=sd(results.frame[,2*w])
  		mean.median[w,3]=1-mean(results.frame[,2*w+1])
  		mean.median[w,4]=sd(1-results.frame[,2*w+1])
	}
	rownames(mean.median)=BF.cutoffs
	mean.median
}


