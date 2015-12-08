\name{NBF}
\alias{NBF}
\title{
Normal/Bayes factors method for finding associated pathways}

\description{
A vector of the computed Bayes factors for the tested pathways.}

\usage{
NBF(y, G, P, a, b, s2, nu)
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
  \item{b}{
Hyper-parameter of the variance assumed for the pathway effects}
  \item{s2}{
Hyper-parameter of the Inverse-Chi-squared distribution assumed for the variance of the response vector }
  \item{nu}{
Hyper-parameter of the Inverse-Chi-squared distribution assumed for the variance of the response vector }
}

\value{
A vector of the computed Bayes factors of the same length as the number of tested pathways 
}

\references{
Evangelou, M., Dudbridge, F., Wernisch, L. (2014). Two novel pathway analysis methods based on a hierarchical model. Bioinformatics, 30(5), 690 - 697.
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
	NBF(y,G,P,a,b,s2,nu)}
}
