SNAL <- function(y,G,P,a,s2){
	Z=G
	X=G%*%P
	ZZ=Z%*%t(Z)
	V=(a*ZZ+diag(length(y)))
	W=matrix.invroot.calculation(V)
	Y=W%*%y
	Phi=W%*%X
	results=SNAL.calculation(Y,Phi,s2)
	results
}
