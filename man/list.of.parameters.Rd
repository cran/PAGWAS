\name{list.of.parameters}
\alias{list.of.parameters}
\docType{data}
\title{
A list with possible hyper-parameters for NBF
}

\description{
A list with 8 different combinations of the hyper-parameters of NBF. Possible values for the hyper-parameters a, b, s2_0 and nu_0 are given in each entry of the list
}

\usage{data(list.of.parameters)
}


\examples{
data(list.of.parameters)
a=list.of.parameters[[1]][1]
b=list.of.parameters[[1]][2]
s2_0=list.of.parameters[[1]][3]
nu_0=list.of.parameters[[1]][4]
}

\keyword{datasets}
