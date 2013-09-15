return.s2.SNAL <- function(y,G,P,a,no.simulations,no.paths,no.snps,s2.cutoffs=c()){
	if(length(s2.cutoffs)==0){
		s2.cutoffs=seq(1,1e-07,length.out=10)
	}
	
	cat("Different values for the tuning parameter s2 are tested:","\n")
	cat(s2.cutoffs)
	cat("\n")
	
	aa=a	
	
	object<-internal.function.simulated.data(y,G,P,no.paths,no.snps)
	y.pred=object$y.pred
	msq=object$msq
	effects.positions.paths=object$effects.positions.paths
	effects.positions.snps=object$effects.positions.snps		
	
	results.frame=matrix(NA,nrow=no.simulations,ncol=2*length(s2.cutoffs)+1)
	for(w in 1:no.simulations){
  		set.seed(w)
  		y.sim <- rnorm(length(y.pred),y.pred,sqrt(msq))
  		ww=1
  		for(ww in 1:length(s2.cutoffs)){
    		s2=s2.cutoffs[ww]
    		gamma.star=SNAL(y.sim,G,P,a=aa,s2)
    		gamma.star=gamma.star[[1]]
    		vector.gamma=which(gamma.star!=0)
    		results.frame[w,((2*ww):(2*ww+1))]=TPR.FPR(P,vector=vector.gamma,effects.positions.paths)
  		}
	}
	rm(w,ww,y.sim,gamma.star,vector.gamma)
	mean.median=data.frame("TPR"=NA,"sd(TPR)"=NA,"FPR"=NA,"sd(FPR)"=NA)
	w=1
	for(w in 1:length(s2.cutoffs)){
  		mean.median[w,1]=mean(results.frame[,2*w])
  		mean.median[w,2]=sd(results.frame[,2*w])
  		mean.median[w,3]=1-mean(results.frame[,2*w+1])
  		mean.median[w,4]=sd(1-results.frame[,2*w+1])
	}
	rownames(mean.median)=s2.cutoffs
	
	rm(y.pred,object,aa,msq,effects.positions.paths,effects.positions.snps,results.frame)
	
	cat("The output matrix gives the true positive rate (TPR) and false positive rate (FPR) for the different tuning parameters. The user can choose s2 according to their objectives.")
	cat("\n")
    cat("\n")

	mean.median
}

