\name{FM.chi.pvalue}
\alias{FM.chi.pvalue}
\title{
Calculates the Fisher's method p-value for each tested pathway
}

\description{
Calculates the Fisher's method \emph{p}-value for the tested pathway. It returns both the \emph{p}-value and the test statistic value of the Fisher's method for the tested pathway
}

\usage{
FM.chi.pvalue(x)
}

\arguments{
  \item{x}{
A vector of \emph{p}-values for either the genes or the SNPs of the tested pathway}
}

\value{
\item{FMstatistic}{Fisher's method test statistic}
\item{FMpvalue}{Fisher's method \emph{p}-value, computed using the exact distribution of the Fisher's method test statistic which is a Chi-squared distribution with degrees of freedom twice the size of the vector x}
}

\references{
Evangelou, M., Rendon, A., Ouhewand, W. H., Wernisch, L., Dudbridge, F. (2012). Comparison of methods for competitive tests of pathway analysis. PLoS ONE 7(7), e41018

Evangelou, M., Dudbridge, F., Wernisch, L. (2013). Two novel pathway analysis methods based on a hierarchical model. Bioinformatics (to appear)
}

\seealso{
\code{\link{pathways}}, \code{\link{snps.to.pathways}}
}

\examples{
FM.chi.pvalue(x=c(0.05,0.1))
}


