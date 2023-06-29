\name{expressionInAnnotations}
\alias{expressionInAnnotations}
\title{Calculate the prototypes from the mixture matrix}
\description{
	Given a matrix of gene expression per spots, and the 
}
\usage{
NNMFproto(X, mix, mccores=0)
}
\arguments{
	\item{X}{Matrix of expression (genes as columns, spots as rows), in non-log space (counts). }
	\item{mix}{The contribution of the prototypes. It's a matrix with as many rows as spot and as many columns as prototypes.}
	\item{mccores}{If 0, then use lapply, otherwise use mclapply with the requested number of core}
}
\value{
    A list containing the following components:
    \item{proto}{The prototypes.}
    \item{sigma}{The extra dispersion from the negative binomial.}
}

\seealso{\link{STstuff}, \link{clustersInOtherDS}, \link{NNMFfromKM}, \link{expressionInAnnotations},
  \link{deconvoluteClusters}
}