\name{TNBCclassif}
\alias{TNBCclassif}
\title{Classifies TNBC samples in subtypes}
\description{
	Using gene expression data of TNBC samples, classifies those samples in subtypes based on
	various schemes.
	The methods are based on our best understanding of the published papers.
}
\usage{
TNBCclassif(x, version=c("lehmann", "bareche", "burstein", "burstein2", "jiang"),
  shortName=FALSE, coef=FALSE, sig=NULL, rescale=TRUE)
}
\arguments{
	\item{x}{Matrix of gene expression, with genes in rows and samples in columns.
	  The genes IDs are the matrix row names. }
	\item{version}{The subtyping method used.}
	\item{shortName}{Should an abbreviated version of the group names be reported?}
	\item{coef}{Should the function report the subtypes or the coefficients used to predict those subtypes?}
	\item{sig}{Signatures to use instead of the signatures by default for the different methods (experts only).}
	\item{rescale}{In many methods, a first step is to rescale the data. This implies a subtype
	  distribution similar to the global distribution. In the case where it is known not to be the case,
	  it is possible to rescale \code{x} by hand and set \code{rescale==FALSE} to avoid further rescaling.}
}
\value{
A vector with the subtype of each sample, or (if \code{coef==TRUE}) a matrix with the coefficients
for the subtypes for each sample.
}
\seealso{\link{nonSTstuff}, \code{\link{entrez2Symb}}, \code{\link{symb2Entrez}}, \code{\link{calcSig}}
}