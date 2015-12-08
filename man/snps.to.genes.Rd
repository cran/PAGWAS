\name{snps.to.genes}
\alias{snps.to.genes}
\title{
Assigns SNPs to genes
}
\description{
Assigns SNPs to genes based on their physical distance.
}

\usage{
snps.to.genes(snp.info, gene.info, distance)
}

\arguments{
  \item{snp.info}{A data frame with 3 columns with names:  Name, Position and Chr that correspond to the SNP name, its position on the genome and its chromosome, respectively
 }
  \item{gene.info}{A data frame with 4 columns with names: Name, Start, End and Chr that correspond to the gene name, start and end positions on the genome and its chromosome, respectively
}
  \item{distance}{A number that corresponds to the distance below and above the Start and End positions of the gene that all SNPs in that region should be assigned to the gene
}

}

\value{
A list of the same size as the number of genes of the gene.info data frame. The names of the SNPs assigned to each gene are returned
}


\seealso{
\code{\link{SNPs}}, \code{\link{genes}}, \code{\link{snps.to.pathways}}
}

\examples{
data(SNPs)
data(genes)
snps.to.genes(snp.info=SNPs,gene.info=genes,distance=50)
}