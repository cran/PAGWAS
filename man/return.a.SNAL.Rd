\name{return.a.SNAL}
\alias{return.a.SNAL}
\title{
Returns the hyper-parameter a of SNAL 
}

\description{
This function returns the appropriate hyper-parameter a for the integrated out SNP effects of SNAL
}

\usage{
return.a.SNAL(y, G, P, no.simulations, no.paths, no.snps)
}

\arguments{
  \item{y}{Response vector of length N}
  \item{G}{Genotype matrix, with N rows and L columns (L is the number of tested SNPs) }
  \item{P}{
Pathway matrix, with L columns and M columns (M is the number of tested pathways) }
  \item{no.simulations}{Number of simulations to run}
  \item{no.paths}{Number of pathways assumed to be causal for the simulations}
  \item{no.snps}{Number of SNPs assumed to be causal for the simulations}
}

\value{
Returns the value of the hyper-parameter a that gives the highest the power of SNAL
}

\references{
Evangelou, M., Dudbridge, F., Wernisch, L. (2013). Two novel pathway analysis methods based on a hierarchical model. Bioinformatics (to appear)
}

\seealso{
\code{\link{SNAL}}, \code{\link{return.s2.SNAL}}
}

