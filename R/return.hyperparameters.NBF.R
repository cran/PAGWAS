return.hyperparameters.NBF <- function(y,G,P,no.simulations,no.paths,no.snps,list.parameters){
	q=seq(1,1e-07,length.out=50)
	n=length(y)
	vec=c()
	
	cat(length(list.parameters)," different combinations of the four hyper-parameters are tested. Different combinations of the hyper-parameters are tested:","\n",sep="")
	print(list.parameters)
	cat("\n")
	
	object<-internal.function.simulated.data(y,G,P,no.paths,no.snps)
	y.pred=object$y.pred
	msq=object$msq
	effects.positions.paths=object$effects.positions.paths
	effects.positions.snps=object$effects.positions.snps	
	
	results.bf=matrix(NA,ncol=length(list.parameters),nrow=no.simulations)
	counter=1
	for(i in 1:no.simulations){
		set.seed(counter)
    		y.sim <- rnorm(length(y.pred),y.pred,sqrt(msq))
		for(w in 1:length(list.parameters)){
			par=list.parameters[[w]]
			aa=par[1]
			bb=par[2]
			s2=par[3]
			nu=par[4]
			results.bf[i,w]=bf.area(y.sim,G,P,aa,bb,s2,nu,effects.positions.paths,bf.cutoffs=q)
		}
	vec=c(vec,which(results.bf[i,]==max(results.bf[i,])))
	}
			
	tt=table(vec)
	tt=which(tt==max(tt))
	par=list.parameters[[tt]]
	A=par[1]
	B=par[2]
	s2=par[3]
	nu=par[4]
		
	output=list()
	output=list(A,B,s2,nu)
	names(output)=c("a","b","s2_0","nu_0")
	output
	}
