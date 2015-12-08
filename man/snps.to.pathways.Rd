\name{snps.to.pathways}
\alias{snps.to.pathways}
\title{
Assigns SNPs to pathways
}
\description{
Assigns SNPs to pathways, using the pathway gene members and the SNPs assigned to each gene.}

\usage{
snps.to.pathways(pathways,gene.snps)
}

\arguments{
  \item{pathways}{A list of pathways with their gene members}
  \item{gene.snps}{A list of genes with the SNPs assigned to them according to their physical distance on the genome}
}

\value{
A list of the same size as the number of pathways in the pathway list. The names of the SNPs assigned to each pathway are returned. Empty pathways are also returned.
}


\seealso{
\code{\link{SNPs}}, \code{\link{genes}}, \code{\link{snps.to.genes}}
}

\examples{
data(SNPs)
data(genes)
data(pathways)
snps.genes <- snps.to.genes(snp.info=SNPs,gene.info=genes, distance=50)
pathway.snps <- snps.to.pathways(pathways,snps.genes)
}


