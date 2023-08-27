\name{computeET}
\alias{computeET}
\title{Compute the ecotypes in another TNBC dataset.}
\description{
	Assign ecotypes (ET) to each sample in another TNBC bulk gene expression dataset.
	It starts from the estimate of the megaclusters with \code{\link{computeMC}}.
}
\usage{
computeET(MC)
}
\arguments{
	\item{MC}{The estimate of the presence of each megacluster, as computed by \code{\link{computeMC}}.}
}
\value{
    A vector with the ecotype of the samples (from 1 to 9).
}
\seealso{\link{STstuff}, \link{clustersInOtherDS}, \link{deconvoluteClusters},
  \link{computeMC}
}