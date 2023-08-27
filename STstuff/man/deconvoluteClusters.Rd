\name{deconvoluteClusters}
\alias{deconvoluteClusters}
\title{Get the protoypes from a clustering}
\description{
	Deconvolve the prototypes from a k-means clustering, so that they are less correlated
	while still fitting the original clustering. The method is similar in principle to the 
	c-means clustering in that elements to cluster can belong to multiple clusters.
	Optionally, remove prototypes that are a mix of two other protoypes.
}
\usage{
deconvoluteClusters(X, k, rmCutOff=0, clean=TRUE, quiet=FALSE,
  killK = c(), Niter=10, mccores=options("mc.cores"))
}
\arguments{
	\item{X}{The data matrix (genes are columns, spots as rows).}
	\item{k}{The clustering, as from the \code{kmeans} function:
	  a numerical vector giving the cluster ID for each spot. }
	\item{rmCutOff}{The cutoff chosen to decide which cluster to remove (between 0 and 1).
	  Higher values will remove more clusters, default will not remove any cluster.}
	\item{clean}{If \code{TRUE}, removes low expressed features. Otherwise keeps all features.
	  Analysis will be done on the top 500 features that fit best the clustering among the surviving features.}
	\item{quiet}{If \code{TRUE}, runs silently}
	\item{killK}{List of clusters that should not be used for the prototypes.}
	\item{Niter}{Number of iterations for the optimization}
	\item{mccores}{If 0, then use \code{lapply}, otherwise use \code{mclapply} with the requested number of cores}
}
\value{
    A list containing the following
    components:
    \item{proto}{The prototypes. Itself a list, with compenents \code{proto} (the mean)
      and \code{sigma} (the extra variance). }
    \item{prevProto}{The prototypes obtained during the algorithm, with less clusters removed.}
    \item{fit}{The mixture matrix, giving the contribution of each prototype in each spot.}
}
\seealso{\link{STstuff}
}

