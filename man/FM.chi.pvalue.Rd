\name{FM.chi.pvalue}
\alias{FM.chi.pvalue}
\title{
Calculates the Fisher's method p-value for each tested pathway
}

\description{
Calculates the Fisher's method \emph{p}-value for a set of \emph{p}-values. It returns both the \emph{p}-value and the test statistic value of the Fisher's product method.
}

\usage{
FM.chi.pvalue(x)
}

\arguments{
  \item{x}{
A vector of \emph{p}-values. These \emph{p}-values can be either gene or SNP \emph{p}-values of a tested pathway}
}

\value{
\item{FMstatistic}{Fisher's product method test statistic}
\item{FMpvalue}{Fisher's method \emph{p}-value, computed using the exact distribution of the Fisher's method test statistic which is a \emph{Chi^2} distribution with degrees of freedom twice the size of vector x}
}

\references{
Evangelou M, Rendon A, Ouwehand WH, Wernisch L, Dudbridge F (2012) Comparison of Methods for Competitive Tests of Pathway Analysis. PLoS ONE 7(7): e41018. doi:10.1371/journal.pone.0041018}

\seealso{
\code{\link{pathways}}, \code{\link{snps.to.pathways}}
}

\examples{
FM.chi.pvalue(x=c(0.05,0.1))}


