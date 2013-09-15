TPR.FPR <-function(P,vector,effects.positions.paths){
	ma=match(vector,effects.positions.paths,nomatch=0)
     	TP=length(which(ma!=0))
     	FP=length(which(ma==0))
     	FN=length(effects.positions.paths)-TP
     	TN=ncol(P)-length(vector)-FN
	c(TP/(TP+FN),TN/(FP+TN))
	}

