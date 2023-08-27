\name{clustersInOtherDS}
\alias{clustersInOtherDS}
\title{Deconvolve the clusterings from bulk datasets}
\description{
	Allows to estimate the proportion of each megacluster (MC) in another bulk gene expression dataset.
}
\usage{
clustersInOtherDS(x, prot, sigma, maxIter=3, rescale=TRUE,
  mccores=options("mc.cores"), quiet=FALSE)
}
\arguments{
	\item{x}{The bulk data matrix (genes as rows, samples as columns). Must be in count-like space (not log). }
	\item{prot}{The prototypes of the megaclusters obtained from ST. Must be count-like. For instance can use a form of
	  harmonic mean of the cluster prototypes. }
	\item{sigma}{Matrix of extra variance (as obtained from \code{\link{deconvoluteClusters}}) by gene.
	  Genes are in rows, clusters in columns.
	  Just used as an initialisation.}
	\item{maxIter}{Maximum number of iteration (rounds of alternatively optimizing \code{sigma} and scaling factors).}
	\item{rescale}{If \code{TRUE}, tries to rescale each gene so that it fits the model as well as possible.}
	\item{mccores}{If 0, then use \code{lapply}, otherwise use \code{mclapply} with the requested number of core}
	\item{quiet}{If \code{TRUE}, runs silently}
}
\value{
    A list containing the following
    components:
    \item{W}{The contribution of each megacluster in each sample.}
    \item{Wp}{Same as \code{W} from previous iteration.}
    \item{sigma}{The inferred extra variance by gene.}
    \item{sc}{The scaling factor by genes (if \code{rescale==TRUE}).}
    \item{qualb}{The evolution of the model log-likelihood during the iterations.}
}
\seealso{\link{STstuff}, \code{\link{computeMC}}, \code{\link{deconvoluteClusters}}
}