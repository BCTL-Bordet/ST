\name{calcSig}
\alias{calcSig}
\title{Calculate the value of signatures in gene expression data}
\description{
	Calculates the value of signatures in gene expression data as weighted means
	of the gene expressions. The function will find the genes in common between the signatures and the
	expression data, and compute the signature using those genes.
}
\usage{
calcSig(d, sig, dropEmpty=TRUE, balanced=FALSE, isCount=FALSE)
}
\arguments{
	\item{d}{Matrix of gene expression, with genes in rows and samples in columns.
	  The genes IDs are the matrix row names. }
	\item{sig}{List of signatures, or alternatively a single signature. Each signature is a
	  dataframe with 3 columns: \code{"name"} that contains the gene names, \code{"entrez"}
	  that contains the gene Entrez ID and \code{"coefficient"} that contains the gene coefficient
	  (i.e.\ weight, potentially negative) in the signature.}
	\item{dropEmpty}{If \code{TRUE}, signatures that have no genes in common with the gene expression
	  matrix are dropped.}
	\item{balanced}{If \code{TRUE}, positive and negative signatures are calculated separately as means, and
	  their difference is reported. This ensures a similar importance for the positive and negative part
	  of the signature.}
	\item{isCount}{If \code{TRUE}, uses a modified version for when \code{d} is a count matrix, and not
	  in log space.}
}
\value{
A vector or matrix with the estimate of the signature(s) expression in each sample.
}
\seealso{\link{nonSTstuff}, \code{\link{entrez2Symb}}, \code{\link{symb2Entrez}}, \code{\link{TNBCclassif}}
}