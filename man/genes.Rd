\name{genes}
\alias{genes}
\docType{data}
\title{
A data-frame of 20 genes with their chromosomes and positions on the genome
}

\description{
A data-frame of 20 genes with their chromosomes and positions on the genome
}

\usage{data(genes)
}

\format{
  A data-frame with 20 observations for the following 4 variables:
  \describe{
    \item{\code{Name}}{Name of gene}
    \item{\code{Start}}{Start position of gene on the genome}
    \item{\code{End}}{End position of gene on the genome}
    \item{\code{Chr}}{Chromosome of gene}
  }
}

\seealso{
\code{\link{SNPs}}, \code{\link{pathways}}
}

\examples{
data(genes)
print(genes[1:5,])
}

\keyword{datasets}
