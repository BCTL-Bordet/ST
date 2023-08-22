\name{expressionInAnnotations}
\alias{expressionInAnnotations}
\title{Calculate the prototypes from the mixture matrix}
\description{
	Given a matrix of gene expression per spots, and the 
}
\usage{
expressionInAnnotations(x, f, sigma, maxIter = 5, mccores = TRUE, quiet = FALSE)
}
\arguments{
	\item{x}{Matrix of expression (genes as columns, spots as rows), in non-log space (counts). }
	\item{f}{The contribution of the prototypes. It's a matrix with as many rows as spot and as many columns as prototypes.}
	\item{sigma}{Starting value for the extra variance by gene.}
	\item{maxIter}{Maximum number of iterations (between sigma / counts estimates)}
	\item{mccores}{If 0, then use lapply, otherwise use mclapply with the requested number of core}
	\item{quiet}{If FALSE, print some information while running}
}
\value{
    A list containing the following components:
    \item{W}{The prototypes.}
    \item{Wp}{The prototypes at each iteration.}
    \item{sigma}{The extra dispersion from the negative binomial.}
    \item{sc}{The fit quality for each spot.}
    \item{qual}{The total fit quality.}
    \item{sp}{sigma at each iteration.}
    \item{scp}{sc at each iteration.}
}

\seealso{\link{STstuff}, \link{clustersInOtherDS}, \link{expressionInAnnotations},
  \link{deconvoluteClusters}
}