\name{computeET}
\alias{computeET}
\title{Compute the ecotypes in another TNBC dataset.}
\description{
	Assign ET to each sample in other TNBC datasets, from the estimate of the MC.
}
\usage{
computeET(MC)
}
\arguments{
	\item{MC}{The estimate of the presence of each megacluster, as computed by computeMC.}
}
\value{
    A vector with the ecotype of the samples (from 1 to 9).
}
\seealso{\link{STstuff}, \link{clustersInOtherDS}, \link{NNMFfromKM}, \link{NNMFproto}, \link{deconvoluteClusters},
  \link{computeMC}
}