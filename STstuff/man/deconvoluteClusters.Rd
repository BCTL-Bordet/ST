\name{deconvoluteClusters}
\alias{deconvoluteClusters}
\title{Get the protoypes from a clustering}
\description{
	Get the protoypes from a clustering. Remove prototypes that are a mix of two other protoypes.
}
\usage{
deconvoluteClusters(X, k, rmCutOff=.5, clean=TRUE, quiet=FALSE, killK = c(), Niter=10, mccores=options("mc.cores"))
}
\arguments{
	\item{X}{The data matrix (genes are columns, spots as rows).}
	\item{k}{The clustering, as from a kmeans - so a numerical vector giving the cluster ID for each spot. }
	\item{rmCutOff}{The cutoff chosen to decide which cluster to remove. Higher values will remove more clusters.}
	\item{clean}{If TRUE, removes low expressed. Otherwise keeps all genes. Will still keep only the top 500 genes that fit best the clustering.}
	\item{quiet}{If TRUE, runs silently}
	\item{killK}{List of clusters that should not be used for the prototypes.}
	\item{Niter}{Number of iterations used by NNMFfromKM}
	\item{mccores}{If 0, then use lapply, otherwise use mclapply with the requested number of core}
}
\value{
    A list containing the following
    components:
    \item{proto}{The prototypes. Itself a list, with compenents proto (the mean) and sigma (the extra variance). }
    \item{prevProto}{The prototypes obtained during the algorithm, with less clusters removed.}
    \item{fit}{The mixture matrix, giving the contribution of each prototype in each spot.}
}
\seealso{\link{STstuff}, \link{getClustFrom}, \link{NNMFfromKM}, \link{NNMFproto}, \link{deconvoluteClusters}
}

