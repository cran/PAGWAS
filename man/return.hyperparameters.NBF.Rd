\name{return.hyperparameters.NBF}
\alias{return.hyperparameters.NBF}
\title{Returns the four hyper-parameters of NBF
}
\description{This function returns the four hyper-parameters a, b, s2_0, nu_0 of NBF}

\usage{return.hyperparameters.NBF(y, G, P, no.simulations, 
no.paths, no.snps, list.parameters)
}

\arguments{
  \item{y}{Response vector of length N}
  \item{G}{Genotype matrix, with N rows and L columns (L is the number of tested SNPs) }
  \item{P}{
Pathway matrix, with L columns and M columns (M is the number of tested pathways) }
   \item{no.simulations}{Number of simulations to run}
  \item{no.paths}{Number of pathways assumed to be causal for the simulations}
  \item{no.snps}{Number of SNPs assumed to be causal for the simulations}
  \item{list.parameters}{A list of various combinations of the four hyper-parameters, see \code{\link{list.of.parameters}}}
}

\value{
Returns a list with the four hyper-parameters of NBF
\item{a}{Hyper-parameter of the variance assumed for the integrated out SNP effects}
\item{b}{Hyper-parameter of the variance assumed for the pathway effects}
\item{s2_0}{Hyper-parameter of the Inverse-Chi-squared distribution assumed for the variance of the response vector }
\item{nu_0}{
Hyper-parameter of the Inverse-Chi-squared distribution assumed for the variance of the response vector 
}
}

\references{
Evangelou, M., Dudbridge, F., Wernisch, L. (2013). Two novel pathway analysis methods based on a hierarchical model. Bioinformatics (to appear)
}

\seealso{
\code{\link{NBF}}, \code{\link{return.bf.NBF}}
}


\examples{
\dontrun{hyper=return.hyperparameters.NBF(y,G,P,no.simulations=100,
no.paths=10,no.snps=20,
list.parameters=list(c(1e-3,1e-3,0.25,200),c(1e-2,1e-3,0.25,100)))
}
}