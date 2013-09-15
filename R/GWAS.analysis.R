GWA.analysis<-function(genotypes,response){

 frame<-data.frame("SNPname"=NA,"teststatistic"=NA,"pvalue"=NA)

 y=response
 X <- scale(genotypes,scale=FALSE)
 y <- y - mean(y)
 n=length(y)
 X[is.na(X)] <- 0
 yX <- as.vector(y %*% X)
 xx <- apply(X*X,2,sum)
 b1 <- yX/xx

 X.pred <- t(t(X) * b1)
 X.diff <- X.pred - y
 ssx <- apply(X.pred * X.pred,2,sum)
 ssres <- apply(X.diff * X.diff,2,sum)

 t.vals<-b1/(sqrt((1/xx)*ssres/(nrow(X)-1)))
 f.vals <- ssx/ ssres * (nrow(X)-1)
 p.vals <- pf(f.vals,1,nrow(X)-1,lower.tail=FALSE)

 frame[1:ncol(genotypes),1]<-colnames(genotypes)
 frame[1:ncol(genotypes),3]<-as.vector(p.vals[1:ncol(genotypes)])
 frame[1:ncol(genotypes),2]<-as.vector(sign(b1[1:ncol(genotypes)])*sqrt(f.vals[1:ncol(genotypes)]))

 frame
 }

