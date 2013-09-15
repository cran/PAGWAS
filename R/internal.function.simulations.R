internal.function.simulated.data <- function(y,Z,X,no.paths,no.snps){
	
	linear=lm(y~X)
	coef=summary(linear)$coefficients[,4]
	names(coef)=colnames(X)
	coef=sort(coef)
	top.paths=names(coef)[1:no.paths]

	frame=GWA.analysis(Z,y)
	frame=frame[order(frame[,3]),]
	top.snps=frame[1:no.snps,1]
	
	ma=match(top.snps,colnames(Z),nomatch=0)
	effects.positions.snps=ma
	rm(ma)

	ma=match(top.paths,colnames(X),nomatch=0)
	effects.positions.paths=ma
	rm(ma)
	
	linear<- lm(y~0+X[,effects.positions.paths] + Z[,effects.positions.snps])
	y.pred <- predict(linear)
	msq <- anova(linear)$"Mean Sq"[3]
	
	object=list("msq"=msq,"y.pred"=y.pred,"effects.positions.paths"=effects.positions.paths,"effects.positions.snps"=effects.positions.snps)
	object
}
