\name{clustersInOtherDS}
\alias{clustersInOtherDS}
\title{Deconvolve the clusterings from other datasets}
\description{
	Allows to estimate the contribution on other datasets of megaclusterings from ST.
}
\usage{
getClustFrom(x, prot, sigma, maxIter=3, rescale=TRUE, mccores=options("mc.cores"), quiet=FALSE)
}
\arguments{
	\item{x}{The bulk data matrix (genes as rows, samples as columns). Must be in count-like space (not log). }
	\item{prot}{The prototypes of the megaclusters obtained from ST. Must be count-like. For instance can use a form of
	  harmonic mean between the clusters. }
	\item{sigma}{Matrix of extra variance (as obtained from NNMFfromKM) by gene. Genes are in rows, clusters in columns.
	  Just used as an initialisation.}
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
\seealso{\link{STstuff}, \link{clustersInOtherDS}, \link{NNMFfromKM}, \link{NNMFproto}, \link{deconvoluteClusters}
}