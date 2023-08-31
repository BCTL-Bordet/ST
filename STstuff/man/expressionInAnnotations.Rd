\name{expressionInAnnotations}
\alias{expressionInAnnotations}
\title{Calculate the gene expression related to a morphological annotation}
\description{
	Given a matrix of gene expression per spots, and the proportion of each
	annotation per spot, will compute the gene expression representative of each annotation.
}
\usage{
expressionInAnnotations(x, f, sigma, maxIter = 5, mccores = TRUE, quiet = FALSE)
}
\arguments{
	\item{x}{Matrix of expression (genes as columns, spots as rows), in non-log space (counts). }
	\item{f}{The proportion of each annotation in each spot.
	    It's a matrix with as many rows as spots and as many columns as annotations.}
	\item{sigma}{Starting value for the extra variance by gene.}
	\item{maxIter}{Maximum number of iterations (between \code{sigma} / counts estimates)}
	\item{mccores}{If 0, then use \code{lapply}, otherwise use \code{mclapply} with the requested number of cores.}
	\item{quiet}{If \code{FALSE}, print some information while running.}
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

\seealso{\link{STstuff}, \code{\link{clustersInOtherDS}},
  \code{\link{deconvoluteClusters}}
}