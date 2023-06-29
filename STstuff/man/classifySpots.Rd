\name{classifySpots}
\alias{classifySpots}
\title{Classify each spot of a spatial slide}
\description{
	Assign to each spot of a spatial slide the proportion of various annotations.
	
}
\usage{
classifySpots(x, type="TNBC")
}
\arguments{
	\item{x}{The count matrix, not normalized. Genes are in rows, samples in columns. }
	\item{type}{The classifiers to use (at this stage only TNBC is allowed). }
}
\value{
    A matrix giving the results of the various classifiers / regressors at the spot level.
}
\seealso{\link{STstuff},
  \link{computeMC}
}