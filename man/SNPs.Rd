\name{SNPs}
\alias{SNPs}
\docType{data}
\title{
A data-frame of 100 SNPs with their chromosomes and positions on the genome
}

\description{ 
A data-frame of 100 SNPs with their positions on the genome
}

\usage{data(SNPs)}

\format{
  A data-frame with 100 observations for the following 3 variables:
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
