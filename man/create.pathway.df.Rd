\name{create.pathway.df}
\alias{create.pathway.df}
\title{
Creates a pathway data frame
}
\description{
Returns a data frame with L rows and M columns. L is the number of SNPs in the genotypes data frame and M is the number of tested pathways. 
}
\usage{
create.pathway.df(genotypes,snps.paths)
}

\arguments{
\item{genotypes}{
Genotype matrix, with L SNPs (columns) and N individuals (rows)}
  \item{snps.paths}{
	A list with entries the SNP members of each pathway. The size of the list is M}
}

\value{
A data frame with columns equal to the number of pathways in the pathway.snps list and rows equal to the number of tested SNPs
}

\seealso{
\code{\link{SNPs}}, \code{\link{genes}}, \code{\link{snps.to.pathways}}
\code{\link{snps.to.genes}}
}

\examples{
data(SNPs)
data(genes)
data(pathways)
data(genotypes)
snps.genes <- snps.to.genes(snp.info=SNPs,gene.info=genes, distance=0)
pathway.snps <- snps.to.pathways(pathways,snps.genes)
P <- create.pathway.df(genotypes=genotypes,snps.paths=pathway.snps)
}


