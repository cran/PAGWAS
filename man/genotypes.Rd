\name{genotypes}
\alias{genotypes}
\docType{data}
\title{
Genotypes for 100 SNPs
}

\description{
A data-frame of 100 individuals (rows) and 100 SNPs (columns). The entries of the genotype matrix are 0, 1 and 2 for the 100 SNPs There are no missing values
}

\usage{data(genotypes)}

\seealso{\code{\link{SNPs}}, \code{\link{genes}}
}

\examples{
data(genotypes)
dim(genotypes)
genotypes[1:5,1:5]
table(genotypes[,1])
}

\keyword{datasets}
