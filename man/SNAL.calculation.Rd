\name{SNAL.calculation}
\alias{SNAL.calculation}
\title{
Sparse Normal/Adaptive lasso method for finding associated variables. The SNAL method is applied to the linear regression Y= Phi beta + epsilon
}

\description{
For more details please read \code{\link{SNAL}.}
}

\usage{
SNAL.calculation(Y, Phi, s2)
}

\arguments{
  \item{Y}{Response vector of length N}
  \item{Phi}{Design matrix, with N rows and M columns (number of tested variables)}
  \item{s2}{
Variance assumed for the response variable, the tuning parameter of the adaptive lasso problem}
}


\value{
\item{gamma.star}{Estimates of gamma hyper-parameters}
\item{ARD}{Posterior estimates of beta coefficients}
}

\references{
Evangelou, M., Dudbridge, F., Wernisch, L. (2014). Two novel pathway analysis methods based on a hierarchical model. Bioinformatics, 30(5), 690 - 697

Wipf, D. and Nagarajan, S. (2008). A new view of automatic relevance determination. Advances in Neural Information Processing Systems, 20  }



\seealso{
\code{\link{SNAL}}
}

\examples{
\dontrun{SNAL.calculation(Y,Phi,s2=0.5)}
}
