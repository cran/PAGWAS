\name{roc.convex}
\alias{roc.convex}
\title{Computes the area under a ROC curve
}

\description{
Computes the area under a ROC curve using the convex hull of the curve
}

\usage{
roc.convex(sens, spec)
}

\arguments{
  \item{sens}{Vector with the values of the recorded sensitivity (true positive rate)
}
  \item{spec}{Vector with the values of the recorded specificity (1-false positive rate)
}
}

\value{Returns the computed area under the ROC curve
}

\references{
Fawcett, T. (2006). An introduction to ROC analysis. Pattern Recognition Letters, 27

Evangelou, M., Dudbridge, F., Wernisch, L. (2013). Two novel pathway analysis methods based on a hierarchical model. Bioinformatics (to appear)
}

\author{
Marina Evangelou, Lorenz Wernisch}

\examples{
roc.convex(sens=c(0.1,0.5,1),spec=c(0.1,0.3,0.7))
}
