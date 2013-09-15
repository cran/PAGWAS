SNAL.calculation<-function(Y,Phi,s2){
	beta=c() 
	z=c()
	gamma=c()
	
	a.check=gamma
	b.check=c()
	
	#initialize chain
	beta=rep(0,ncol(Phi))
	z=rep(1,ncol(Phi))
	gamma=rep(0,ncol(Phi))

	beta.star=lars.calculation(x=Phi,y=Y,z=z,s2=s2[1])
	gamma[1:length(beta)]=z[1:length(beta)]^(-1/2)*abs(beta.star[1:length(beta)])
	Ggamma=diag(gamma)

	S=Phi%*%Ggamma%*%t(Phi)
	S=S+diag(length(Y))*s2[1]
	S.inv=matrix.inv.calculation(V=S)

	z.star=diag(t(Phi)%*%S.inv%*%Phi)
	beta=beta.star
	z=z.star
	counter=0
while(length(which((a.check-b.check)<=1e-6 & (a.check-b.check)>=-1e-6))!=length(a.check)){
    a.check=gamma
    beta.star=lars.calculation(x=Phi,y=Y,z=z,s2=s2[1])
    gamma[1:length(beta)]=z[1:length(beta)]^(-1/2)*abs(beta.star[1:length(beta)])
    Ggamma=diag(gamma)

    S=Phi%*%Ggamma%*%t(Phi)
    S=S+diag(length(Y))*s2[1]
    S.inv=matrix.inv.calculation(V=S)
    z.star=diag(t(Phi)%*%S.inv%*%Phi)
    beta=beta.star
    z=z.star
    b.check=gamma
    counter=counter+1  }

S.star=Phi%*%diag(gamma)%*%t(Phi)
S.star=S.star+diag(length(Y))
S.star.inv=matrix.inv.calculation(V=S.star)
ARD=as.vector(diag(gamma)%*%t(Phi)%*%S.star.inv%*%Y)
gamma.star=gamma

names(gamma.star)=colnames(Phi)
names(ARD)=colnames(Phi)
results=list("gamma"=gamma.star,"ARD"=ARD)
}