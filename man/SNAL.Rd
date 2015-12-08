\name{SNAL}
\alias{SNAL}
\title{Sparse Normal/Adaptive lasso method for finding associated pathways}

\description{Sparse Normal/Adaptive lasso method applied for finding the associated pathways. The iterative algorithm suggested by Wipf and Nagarajan (2008) is applied. A vector equal to the number of tested pathways is returned, the zero entries of the vector correspond to the pathways that are not associated. The posterior estimates of the beta coefficients are also returned as they are described by Wipf and Nagarajan (2008).
}

\usage{
SNAL(y, G, P, a, s2)
}

\arguments{
  \item{y}{
Response vector of length N}
  \item{G}{
Genotype matrix, with N rows and L columns (number of tested SNPs) }
  \item{P}{
Pathway matrix, with L columns and M columns (number of tested pathways) }
  \item{a}{
Hyper-parameter of the variance assumed for the integrated out SNP effects}
  \item{s2}{
Variance assumed for the response variable, the tuning parameter of adaptive lasso}
}


\value{
\item{gamma.star}{Estimates of gamma hyper-parameters}
\item{ARD}{Posterior estimates of beta coefficients}
}

\references{
Evangelou, M., Dudbridge, F., Wernisch, L. (2014). Two novel pathway analysis methods based on a hierarchical model. Bioinformatics, 30(5), 690 - 697.

Wipf, D. and Nagarajan, S. (2008). A new view of automatic relevance determination. Advances in Neural Information Processing Systems, 20 
}


\seealso{
\code{\link{SNAL.calculation}}
}

\examples{
\dontrun{
	data(genotypes)
	G=genotypes
	data(pathways)
	data(SNPs)
	data(genes)
	snps.genes=snps.to.genes(SNPs,genes,distance=0)
	snps.paths=snps.to.pathways(pathways,snps.genes)
	P=create.pathway.df(G,snps.paths)
	y=rnorm(nrow(G),mean=0,sd=10)
	SNAL(y,G,P,a,s2)}
}
