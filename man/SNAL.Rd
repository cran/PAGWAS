\name{SNAL}
\alias{SNAL}
\title{Sparse Normal/Adaptive lasso (SNAL) method for finding associated pathways}

\description{
An implementation of the Sparse Normal/Adaptive lasso method proposed by Evangelou et al (2013) for finding associated pathways. A list of two vectors, both of the same size as the number of tested pathways, is returned. The zero entries of the first vector, named gamma.star, correspond to the pathways that are not associated and the non-zero entries of the vector correspond to the associated pathways. The posterior estimates of the beta coefficients of the Bayesian hierarchical model are also returned as they described by Wipf and Nagarajan (2008)
}

\usage{
SNAL(y, G, P, a, s2)
}

\arguments{
  \item{y}{
Response vector of length N}
  \item{G}{
Genotype matrix, with N rows and L columns (L is the number of tested SNPs) }
  \item{P}{
Pathway matrix, with L columns and M columns (M is the number of tested pathways) }
  \item{a}{
Hyper-parameter of the variance assumed for the integrated out SNP effects}
  \item{s2}{
Variance assumed for the response variable that corresponds to the tuning parameter of the adaptive lasso problem}
}


\value{
\item{gamma.star}{Estimates of gamma hyper-parameters}
\item{ARD}{Posterior estimates of the beta coefficients}
}

\references{
Evangelou, M., Dudbridge, F., Wernisch, L. (2013). Two novel pathway analysis methods based on a hierarchical model. Bioinformatics (to appear)

Wipf, D. and Nagarajan, S. (2008). A new view of automatic relevance determination. Advances in Neural Information Processing Systems, 20 
}


\seealso{
\code{\link{pathways}}, \code{\link{genotypes}}, \code{\link{create.pathway.matrix}}, \code{\link{NBF}}
}

\examples{
\dontrun{
s=SNAL(y,G,P,a,s2)
associated_paths=which(s$gamma.star!=0)
}
}
