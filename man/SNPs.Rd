\name{SNPs}
\alias{SNPs}
\docType{data}
\title{
A data frame of 100 artificial SNPs with their chromosomes and positions on the genome
}

\description{ 
A data frame with 100 rows and 3 columns.
}

\usage{data(SNPs)}

\format{
	Column names:
  \describe{
    \item{\code{Name}}{SNP name}
    \item{\code{Position}}{Position of SNP on the genome}
    \item{\code{Chr}}{Chromosome of the SNP}
  }
}


\seealso{
\code{\link{genes}}, \code{\link{genotypes}}
}


\examples{
data(SNPs)
print(SNPs[1:5,])
}

\keyword{datasets}
