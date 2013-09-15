bf.area <- function(y,G,P,a,b,s2,nu,effects.positions.paths,bf.cutoffs){	   	
	i=1
	BF=NBF(y,G,P,a,b,s2=s2,nu=nu)
	M=ncol(P)
  	TPR.FPR=lapply(bf.cutoffs,bf.TPR.FPR,BF,effects.positions.paths,M)
  	vector.TPR=seq(1,2*length(bf.cutoffs),by=2)
  	vector.TNR=seq(2,2*length(bf.cutoffs),by=2)
  	vector.TPR=unlist(TPR.FPR)[vector.TPR]
  	vector.TNR=unlist(TPR.FPR)[vector.TNR]
  	sens.bf=c(0,vector.TPR,1)
  	spec.bf=c(1,vector.TNR,0)
  	conv.bf=roc.convex(sens=sens.bf,spec=spec.bf)
	conv.bf
}


