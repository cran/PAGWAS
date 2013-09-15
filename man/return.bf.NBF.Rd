\name{return.bf.NBF}
\alias{return.bf.NBF}
\title{
Returns the Bayes factor threshold of NBF
}
\description{
This function returns a data-frame of true positive rates and false positive rates of NBF for the Bayes factors thresholds chosen by the user
}

\usage{
return.bf.NBF(y, G, P, a, b, s2, nu, 
	no.simulations, no.paths, no.snps, BF.cutoffs)
}

\arguments{
  \item{y}{
Response vector of length N}
  \item{G}{
Genotype matrix, with N rows and L columns (L is the number of tested SNPs) 
}
  \item{P}{
Pathway matrix, with L columns and M columns (M is the number of tested pathways)
}   
\item{a}{ 
Hyper-parameter a of the variance assumed for the integrated out SNP effects
}
  \item{b}{
Hyper-parameter b of the variance assumed for the pathway effects
}
\item{s2}{
Hyper-parameter s2_0 of the Inverse-Chi-squared distribution assumed for the variance of the response vector 
}
  \item{nu}{
Hyper-parameter nu_0 of the Inverse-Chi-squared distribution assumed for the variance of the response vector 
}   
\item{no.simulations}{
Number of simulations to run
}
\item{no.paths}{
Number of pathways assumed to be causal for the simulations
}
\item{no.snps}{
Number of SNPs assumed to be causal for the simulations
}
\item{BF.cutoffs}{
A vector of Bayes factors (BFs) thresholds to test
}
}


\value{
Returns a data-frame of true positive and false positive rates for the Bayes factors (BFs) thresholds given by the user
}

\references{
Evangelou, M., Dudbridge, F., Wernisch, L. (2013). Two novel pathway analysis methods based on a hierarchical model. Bioinformatics (to appear)
}


\seealso{
\code{\link{NBF}}, \code{\link{return.hyperparameters.NBF}}
}

\examples{
\dontrun{return.bf.NBF(y,G,P,a=1e-4,b=1e-2,s2=0.25,nu=200,
no.simulations=100,no.paths=10,no.snps=20,BF.cutoffs=c(0.5,0.75,0.9,1))}
}
