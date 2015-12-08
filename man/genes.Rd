\name{genes}
\alias{genes}
\docType{data}
\title{
A data frame of 20 artificial genes with their chromosomes and positions on the genome
}

\description{
A data frame with 20 rows and 4 columns.
}

\usage{data(genes)}

\format{
Column names:
  \describe{
    \item{\code{Name}}{Name of gene}
    \item{\code{Start}}{Start position of gene on the genome}
    \item{\code{End}}{End position of gene on the genome}
    \item{\code{Chr}}{Chromosome of gene}
  }
}

\seealso{
\code{\link{SNPs}}
}

\examples{
data(genes)
print(genes[1:5,])
}

\keyword{datasets}
