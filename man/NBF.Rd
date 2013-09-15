\name{NBF}
\alias{NBF}
\title{
Normal/Bayes factors (NBF) method for finding associated pathways
}

\description{
An implementation of the Normal/Bayes factors method proposed by Evangelou et al (2013) for finding associated pathways. A vector of the computed Bayes factors for each tested pathway is returned
}

\usage{
NBF(y, G, P, a, b, s2, nu)
}

\arguments{
  \item{y}{
Response vector of length N}
  \item{G}{
Genotype matrix, with N rows and L columns (L is the number of tested SNPs) }
  \item{P}{
Pathway matrix, with L columns and M columns (M is the number of tested pathways) }
  \item{a}{
Hyper-parameter a of the variance assumed for the integrated out SNP effects}
  \item{b}{
Hyper-parameter b of the variance assumed for the pathway effects}
  \item{s2}{
Hyper-parameter s2_0 of the Inverse-Chi-squared distribution assumed for the variance of the response vector }
  \item{nu}{
Hyper-parameter nu_0 of the Inverse-Chi-squared distribution assumed for the variance of the response vector }
}

\value{
A vector of the computed Bayes factors of the same size as the number of tested pathways 
}

\references{
Evangelou, M., Dudbridge, F., Wernisch, L. (2013). Two novel pathway analysis methods based on a hierarchical model. Bioinformatics (to appear)
}


\seealso{
\code{\link{return.hyperparameters.NBF}}, \code{\link{return.bf.NBF}}, \code{\link{SNAL}}
}

\examples{
\dontrun{NBF(y,G,P,a,b,s2,nu)}
}
