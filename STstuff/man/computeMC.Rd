\name{computeMC}
\alias{computeMC}
\title{Compute the MC from gene expression in another TNBC dataset.}
\description{
	Allows to estimate the contribution on other datasets of megaclusterings from ST.
}
\usage{
computeMC(x, maxIter=3, rescale=TRUE, mccores=options("mc.cores"), quiet=FALSE)
}
\arguments{
	\item{x}{The bulk data matrix (genes as rows, samples as columns). Must be in count-like space (not log). }
	\item{maxIter}{Maximum number of iteration (rounds of optimizing sigma and scaling factors).}
	\item{rescale}{If TRUE, tries to rescale each gene so that it fits the model as well as possible.}
	\item{mccores}{If 0, then use lapply, otherwise use mclapply with the requested number of core}
	\item{quiet}{If TRUE, runs silently}
}
\value{
    A list containing the following
    components:
    \item{W}{The contribution of each megacluster in each sample.}
    \item{Wp}{Same as W from previous iteration.}
    \item{sigma}{The inferred extra variance by gene.}
    \item{sc}{The scaling factor by genes (if rescale=TRUE).}
    \item{qualb}{The evolution of the model log-likelihood during the iterations.}
}
\seealso{\link{STstuff}, \link{clustersInOtherDS}, \link{deconvoluteClusters},
  \link{computeET}
}