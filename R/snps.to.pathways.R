snps.to.pathways <-function(pathway,gene.snps){
	pathway.snps.list=list()
		for(i in 1:length(pathway)){
			x=pathway[[i]]
			m=match(x,names(gene.snps))
			dat=c()
			if(length(m)>0){
				kk=1
				for(kk in 1:length(m)){
					dat=c(dat,gene.snps[[m[[kk]]]])
					}
				dat=unique(dat)
				}
			pathway.snps.list[[i]]=dat
			}
		names(pathway.snps.list)=names(pathway)
		rm(i)
		x=which(pathway.snps.list=="NULL")
		if(length(x)>0){
			pathway.snps.list=pathway.snps.list[-x]
			}
			pathway.snps.list
			}

