\name{create.pathway.matrix}
\alias{create.pathway.matrix}
\title{
Creates a pathway matrix
}
\description{
Creates a pathway matrix, with rows the SNPs assigned to each pathway
}
\usage{
create.pathway.matrix(genotypes,pathway.snps)
}

\arguments{
\item{genotypes}{
Genotype matrix, with L SNPs (columns) and N individuals (rows)}
  \item{pathway.snps}{
A list of the SNPs members of each pathway}
}

\value{
A matrix with columns equal to the number of pathways in the pathway.snps list and rows equal to the number of SNPs in the genotypes data-frame
}

\seealso{
\code{\link{SNPs}}, \code{\link{genes}}, \code{\link{snps.to.pathways}},
\code{\link{snps.to.genes}}
}

\examples{
\dontrun{
data(SNPs)
data(genes)
data(pathways)
data(genotypes)
snps.genes <- snps.to.genes(snp.info=SNPs,gene.info=genes, distance=50)
pathway.snps <- snps.to.pathways(pathways,snps.genes)
P <- create.pathway.matrix(genotypes,pathway.snps)
}
}

