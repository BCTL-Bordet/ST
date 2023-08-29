\name{classifySpots}
\alias{classifySpots}
\title{Classify each spot of a spatial slide}
\description{
	Assigns to each spot of a spatial slide the proportion of various annotations.
	
}
\usage{
classifySpots(x, type="TNBC")
}
\arguments{
	\item{x}{The count matrix, not normalized. Genes are in rows, samples in columns. }
	\item{type}{The classifiers to use (at this stage only \code{"TNBC"} is allowed, which
	  relates to classifiers obtained on the triple negative breast cancer study). }
}
\value{
    A matrix giving the results of the various classifiers / regressors at the spot level.
}
\seealso{\link{STstuff},
  \code{\link{computeMC}}
}