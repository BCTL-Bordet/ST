\name{STstuff}
\alias{STstuff}
\title{The STstuff package}
\description{
	Package with functions used for the analysis of spatial transcriptomics data,
	in particular related to metaclustering, deconvolution of gene expression between annotations
	and specifically triple negative breast cancer (TNBC) data.
}
\details{
	\code{\link{deconvoluteClusters}} is used to find prototypes from a k-means that are less
	  correlated.
	  
	\code{\link{expressionInAnnotations}} is used to find gene expression profile that correspond
	  to morphological annotations in an ST experiment. For instance, having the proportion of
	  tumor cells in each spot, the function will output the gene expression related to the tumor.
	  
	\code{\link{clustersInOtherDS}} is used to find megaclusters (MCs, clusters of clusters) in
	another bulk gene expression dataset.
	
	The following functions are specific to the ST TNBC study.
	\code{\link{computeMC}} computes the MCs in another bulk gene expression TNBC dataset.
	\code{\link{computeET}} computes the ecotypes (groups of samples with similar MCs)
	in another bulk gene expression TNBC dataset.
	\code{\link{classifySpots}} estimates the proportion of each annotation in each spot of
	an ST TNBC dataset, as defined in the ST TNBC study (use on other tissue types at your own
	risk).
}
\seealso{\link{STstuff}, \code{\link{clustersInOtherDS}}, \code{\link{expressionInAnnotations}},
  \code{\link{deconvoluteClusters}},  \code{\link{computeET}},  \code{\link{computeMC}},
  \code{\link{classifySpots}}
}