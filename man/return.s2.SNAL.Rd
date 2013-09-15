\name{return.s2.SNAL}
\alias{return.s2.SNAL}
\title{
Returns the s2 tuning parameter of SNAL
}

\description{
This function returns a data-frame of true positive rates and false positive rates of SNAL for the s2 thresholds chosen by the user
}

\usage{
return.s2.SNAL(y, G, P, a, no.simulations, no.paths, no.snps, s2.cutoffs)
}

\arguments{
  \item{y}{Response vector of length N}
  \item{G}{Genotype matrix, with N rows and L columns (L is the number of tested SNPs) }
  \item{P}{
Pathway matrix, with L columns and M columns (M is the number of tested pathways) }
    \item{a}{ Hyper-parameter of the variance assumed for the integrated out SNP effects}
  \item{no.simulations}{Number of simulations to run
}
  \item{no.paths}{Number of pathways assumed to be causal for the simulations}
  \item{no.snps}{Number of SNPs assumed to be causal for the simulations
}
  \item{s2.cutoffs}{A vector of s2 thresholds to test
}
}

\value{
Returns a data-frame of true positive and false positive rates for the s2 thresholds given by the user
}

\references{
Evangelou, M., Dudbridge, F., Wernisch, L. (2013). Two novel pathway analysis methods based on a hierarchical model. Bioinformatics (to appear)
}


\seealso{
\code{\link{SNAL}},\code{\link{return.a.SNAL}}
}

\examples{
\dontrun{return.s2.SNAL(y,G,P,a=1e-2,no.simulations=100,
		no.paths=10,no.snps=20,s2.cutoffs=c(0.5,1,2))
}
}