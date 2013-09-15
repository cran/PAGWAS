bf.TPR.FPR <- function(bf.cutoffs,BF,effects.positions.paths,M){
	threshold=bf.cutoffs
	vector=which(BF<=threshold)
	ma=match(vector,effects.positions.paths,nomatch=0)
     	TP=length(which(ma!=0))
     	FP=length(which(ma==0))
     	FN=length(effects.positions.paths)-TP
     	TN=M-length(vector)-FN
	c(TP/(TP+FN),TN/(FP+TN))
	}

