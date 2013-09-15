create.pathway.matrix <- function(genotypes,pathway.snps){
	unique.snps=unique(unlist(pathway.snps))
	P=matrix(0,nrow=length(unique.snps),ncol=length(pathway.snps))
	rownames(P)=unique.snps
	colnames(P)=names(pathway.snps)
	for(i in 1:length(pathway.snps)){
		m=match(pathway.snps[[i]],rownames(P))
		P[m,i]=1
		rm(m)
		}
	rm(i)
	
	L=nrow(genotypes)
	new=matrix(0,nrow=(L-nrow(P)),ncol=length(pathway.snps))
	
	colnames(new)=colnames(P)
	rownames(new)=colnames(genotypes)[!(colnames(genotypes) %in% rownames(P))]
	P=rbind(P,new)
	
	m=match(colnames(genotypes),rownames(P))
	P=P[m,]
	rm(new,m,L)
	P
	}

